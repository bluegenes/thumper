#! /usr/bin/env python
"""
At whole-MAG or contig-level:
1. do search w/containment
2. aggregate matched hashes to taxonomic rank (optional)
3. save matched annotation info to CSV, signatures to matches sigfiles

modified from charcoal contigs_search.py c7fd192
"""
import sys
import argparse
import os.path
import json
import csv

import screed

import sourmash
from sourmash.lca.command_index import load_taxonomy_assignments
from sourmash.lca import LCA_Database, lca_utils

from .lineage_db import LineageDB
from .version import version
from thumper.charcoal_utils import (gather_at_rank, get_ident, ContigGatherInfo)
from .search_utils import (gather_guess_tax_at_each_rank, search_containment_at_rank, SearchFiles)

def main(args):
    "Main entry point for scripting. Use cmdline for command line entry."
    genomebase = os.path.basename(args.genome)
    match_rank = 'genus'

    # load taxonomy CSV
    tax_assign, _ = load_taxonomy_assignments(args.lineages_csv,
                                              start_column=2)
    print(f'loaded {len(tax_assign)} tax assignments.')

    # load the genome signature
    genome_sig = sourmash.load_one_signature(args.genome_sig, select_moltype=args.alphabet, ksize=args.ksize)

    # load all of the matches from search --containment in the database
    with open(args.matches_sig, 'rt') as fp:
        try:
            siglist = list(sourmash.load_signatures(fp, do_raise=True,
                                                    quiet=False))
        except sourmash.exceptions.SourmashError:
            siglist = []
    print(f"loaded {len(siglist)} matches from '{args.matches_sig}'")

    # Hack for examining members of our search database: remove exact matches.
    new_siglist = []
    for ss in siglist:
        if genome_sig.similarity(ss) == 1.0:
            print(f'removing an identical match: {ss.name()}')
        else:
            new_siglist.append(ss)
    siglist = new_siglist

    # init search and gather tax dicts
    search_tax, gather_tax,genome_gather_tax = {},{},{}
    if not siglist:
        # write empty files so snakemake workflows don't complain; exit.
        print('no non-identical matches for this genome, exiting.')
        if not args.no_search_contigs:
            sf = SearchFiles(args.output_prefix, not args.no_search, args.gather, contigs=True)
            sf.write_taxonomy_json(search_tax, result_type="search")
            sf.write_taxonomy_json(gather_tax, result_type="gather")
            sf.close()
        if args.search_genome:
            gf = SearchFiles(args.output_prefix, not args.no_search, args.gather, contigs=False)
            gf.write_taxonomy_json(genome_gather_tax, result_type="gather")
            gf.close()
        return 0

    # construct a template minhash object that we can use to create new 'uns
    empty_mh = siglist[0].minhash.copy_and_clear()
    ksize = empty_mh.ksize
    scaled = empty_mh.scaled
    moltype = empty_mh.moltype

    # create empty LCA database to populate...
    lca_db = LCA_Database(ksize=ksize, scaled=scaled, moltype=moltype)
    lin_db = LineageDB()

    # ...with specific matches.
    for ss in siglist:
        ident = get_ident(ss)
        lineage = tax_assign[ident]

        lca_db.insert(ss, ident=ident)
        lin_db.insert(ident, lineage)

    print(f'loaded {len(siglist)} signatures & created LCA Database')
    print('')
    print(f'reading contigs from {genomebase}')

    screed_iter = screed.open(args.genome)
    genome_len = 0

    search_tax = {}
    gather_tax = {}

    if not args.no_search_contigs:
        sf = SearchFiles(args.output_prefix, not args.no_search, args.gather, contigs=True)

        for n, record in enumerate(screed_iter):
            # look at each contig individually
            mh = empty_mh.copy_and_clear()
            mh.add_sequence(record.sequence, force=True)
            # search, optionally aggregate matched hashes to get containment at rank

            seq_len = len(record.sequence)
            genome_len+=seq_len
            num_hashes = len(mh.hashes)

            if not args.no_search:
                search_results, search_rank_results = search_containment_at_rank(mh, lca_db, lin_db, match_rank)

                if not search_results:
                    # write to unclassified
                    sf.unmatched.write(">" + record.name + "\n" + record.sequence + "\n")
                    continue # if no search results, don't bother with gather
                else:
                    # first, print normal search --containment results
                    for sr in search_results:
                        sf.write_result(sr, record.name, seq_len, result_type="search")
                    # now, print containment at rank results
                    for sr in search_rank_results:
                        sf.write_result(sr, record.name, seq_len, result_type="ranksearch")

                    #info = ContigRankSearchInfo(len(record.sequence), len(mh), results)
                    #search_tax[record.name] = info

            if args.gather:
                # first, gather at match rank (default genus)
                gather_results = list(gather_at_rank(mh, lca_db, lin_db, match_rank))
                # store dict of gather_results
                info = ContigGatherInfo(len(record.sequence), len(mh), gather_results)
                gather_tax[record.name] = info

                if not gather_results:
                    # write to unclassified. should only get here if no search OR gather results
                    sf.unmatched.write(">" + record.name + "\n" + record.sequence + "\n")
                else:
                    # next, summarize at higher ranks
                    gather_taxonomy_per_rank = gather_guess_tax_at_each_rank(gather_results, num_hashes, \
                                                                             minimum_matches=args.gather_min_matches, \
                                                                             lowest_rank=match_rank, \
                                                                             taxlist=lca_utils.taxlist(include_strain=False))
                    # write taxonomy out
                    for gr in gather_taxonomy_per_rank:
                        sf.write_result(gr, record.name, seq_len, result_type="rankgather")

        print(f"Processed {n+1} contigs.")
        # close contig files
        sf.write_taxonomy_json(gather_tax, result_type ="gather")
        sf.close()

    if args.search_genome:
        gf = SearchFiles(args.output_prefix, not args.no_search, args.gather, contigs=False)
        # MAG workflow
        entire_mh = genome_sig.minhash
        genome_name = genome_sig.name()
        num_hashes = len(entire_mh.hashes)
        if not genome_len:
            for record in screed_iter:
                genome_len+=len(record.sequence)
        if not args.no_search:
            #results are guaranteed, otherwise would have exited before searching
            search_results, search_rank_results = search_containment_at_rank(entire_mh, lca_db, lin_db, match_rank)
            for sr in search_results:
                gf.write_result(sr, genome_name, genome_len, result_type="search")
            for sr in search_rank_results:
                gf.write_result(sr, genome_name, genome_len, result_type="ranksearch")
        if args.gather:
            gather_results = list(gather_at_rank(entire_mh, lca_db, lin_db, match_rank))
            genome_gather_info = ContigGatherInfo(genome_len, len(entire_mh), gather_results)
            genome_gather_tax[genome_name] = genome_gather_info
            # next, summarize at higher ranks
            gather_taxonomy_per_rank = gather_guess_tax_at_each_rank(gather_results, num_hashes, \
                                                                     minimum_matches=args.gather_min_matches, \
                                                                     lowest_rank=match_rank, \
                                                                     taxlist=lca_utils.taxlist(include_strain=False))
            for gather_res in gather_taxonomy_per_rank:
                gf.write_result(gather_res, genome_name, genome_len, result_type="rankgather")
        # close genome files
        gf.write_taxonomy_json(genome_gather_tax, result_type="gather")
        gf.close()

    return 0


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser(sys_args)
    p.add_argument('--genome', help='genome file', required=True)
    p.add_argument('--genome-sig', help='genome sig', required=True)
    p.add_argument('--matches-sig', help='all relevant matches', required=True)
    p.add_argument('--lineages-csv', help='lineage spreadsheet', required=True)
    p.add_argument('--alphabet', help='alphabet', required=True)
    p.add_argument('--ksize', help='ksize', required=True)

    p.add_argument('--force', help='continue past survivable errors',
                   action='store_true')

    # search type switches
    p.add_argument('--no-search', help='do not run search with containment',
                   action='store_true')
    p.add_argument('--gather', help='run sourmash gather',
                   action='store_true')
    p.add_argument('--gather-min-matches', type=int, default=3)
    p.add_argument('--no-search-contigs', action="store_true")
    p.add_argument('--search-genome', action="store_true")

    # output options:
    # build outputs based on the options above.
    p.add_argument('--output-prefix',
                    help='prefix for search or gather results',
                    required=True)
    args = p.parse_args()

    return main(args)


# execute this, when run with `python -m`.
if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
