#! /usr/bin/env python
"""
At whole-MAG or contig-level:
1. do search w/containment
2. aggregate matched hashes to taxonomic rank (optional)
3. save matched hashes to JSON, annotation info to CSV

started from charcoal contigs_search.py c7fd192
"""
import sys
import argparse
import os.path
import json
import csv

import screed

import sourmash
from sourmash.lca.command_index import load_taxonomy_assignments
from sourmash.lca import LCA_Database

from .lineage_db import LineageDB
from .version import version
from thumper.charcoal_utils import (gather_at_rank, get_ident)
from .search_utils import (gather_guess_tax_at_each_rank, search_containment_at_rank)


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

    # if, after removing exact match(es), there is nothing left, quit.
    # (but write an empty JSON file so that snakemake workflows don't
    # complain.)
    if not siglist:
        print('no non-identical matches for this genome, exiting.')
        contigs_tax = {}
        # TO DO - need to modify this
        with open(args.json_out, 'wt') as fp:
            fp.write(json.dumps(contigs_tax))
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

    out_prefix = args.output_prefix
    # contigs with no matches
    unmatchedF = f"{out_prefix}.unmatched.fq"
    unmatched = open(unmatchedF, "w")


    if not args.no_search:
        # set output filenames
        search_csvF = f"{out_prefix}.search.csv"
        search_matchesF = f"{out_prefix}.search.matches.sig"
        ranksearch_csvF = f"{out_prefix}.search.rank.csv"
        ranksearch_matchesF = f"{out_prefix}.search.rank.matches.sig"

        # open csvs, write headers
        search_csv = open(search_csvF, "w")
        search_fieldnames = ['contig_name', 'length', 'similarity', 'name', 'md5', 'lineage']
        #search_fieldnames = ['similarity', 'name', 'filename', 'md5', 'lineage']
        search_w =csv.DictWriter(search_csv, fieldnames=search_fieldnames)
        search_w.writeheader()

        rank_csv = open(ranksearch_csvF, "w")
        rank_fieldnames = ['contig_name', 'contig_length', 'match_rank', 'lineage', 'contained_at_rank', 'contained_bp']
        rank_w = csv.DictWriter(rank_csv, fieldnames=rank_fieldnames)
        rank_w.writeheader()
        #'similarity, match, md5, filename, name, lineage')
        #'lineage, containment, intersect_bp, match_sig')
        # for rank results: contig_name,contig_length,match_rank,lineage,contained_at_rank,bp_contained #also num_hashes?

       # open matches sig files
        search_matches = open(search_matchesF, "w")
        rank_matches = open(ranksearch_matchesF, "w")

    if args.gather:
        #set output filenames
        contigs_tax = {} # dict for gather results
        gather_csvF = f"{out_prefix}.gather.csv"
        rankgather_csvF = f"{out_prefix}.gather.rank.csv"
        # open gather files
        gather_csv = open(gather_csvF, "w")
        rankgather_csv = open(rankgather_csvF, "w")
        # use csv dictwriter?
        gather_rank_fieldnames = ['contig_name', 'contig_length', 'match_rank', 'lineage', 'f_ident', 'f_major']
        gather_rank_w = csv.DictWriter(rank_csv, fieldnames=gather_rank_fieldnames)

    for n, record in enumerate(screed_iter):
        # look at each contig individually
        mh = empty_mh.copy_and_clear()
        mh.add_sequence(record.sequence, force=True)
        # search, optionally aggregate matched hashes to get containment at rank

        contig_len = len(record.sequence)
        num_hashes = len(mh.hashes)

        if not args.no_search:
            search_results, search_rank_results = search_containment_at_rank(mh, lca_db, lin_db, match_rank)

            if not search_results:
                # write to unclassified
                unmatched.write(record.name + "\n" + record.sequence + "\n")
                continue # if no search results, don't bother with gather
            else:
                # first, print normal search --containment results --> functionize this!
                for sr in search_results:
                    d = dict(sr._asdict())
                    # save match to output matches
                    sourmash.signature.save_signatures([d['match']], search_matches)
                    del d['match']
                    # better way to do this?
                    d["contig_name"] = record.name
                    d["length"] = contig_len
                    search_w.writerow(d)

                # now, print containment at rank results
                for sr in search_rank_results:
                    d = dict(sr._asdict())
                    # save match to output matches
                    sig.save_signatures([d['match_sig']], rank_matches)
                    del d['match_sig']
                    d["contig_name"] = contig_name
                    d["length"] = contig_len
                    d["rank"] = sr.lineage[-1]
                    rank_w.writerow(d)
                    #report_csv.write(res)

        if args.gather:
            # first, gather at match rank (default genus)
            gather_results = list(gather_at_rank(mh, lca_db, lin_db, match_rank))
            # write gather_results??

            if not gather_results:
                # write to unclassified. should only get here if no search OR gather results
                unmatched.write(record.name + "\n" + record.sequence + "\n")
            else:
                # next, summarize at higher ranks
                gather_taxonomy_per_rank = gather_guess_tax_at_each_rank(gather_results, num_hashes, \
                                                                         minimum_matches=args.gather_min_matches, \
                                                                         lowest_rank=match_rank, \
                                                                         taxlist=lca_utils.taxlist(include_strain=False))
                #results = list of RankSumGatherResult = namedtuple('RankSumGatherResult', 'lineage, f_ident, f_major')

                # write taxonomy out
                for gr in gather_taxonomy_per_rank:
                    d = dict(gr._asdict())
                    # save match to output matches
                    d["contig_name"] = contig_name
                    d["length"] = contig_len
                    d["rank"] = gr.lineage[-1]
                    d["major_bp"] = get_match_bp(float(gr.f_major))
                    gather_rank_w.writerow(d)
                    #report_csv.write(res)

    print(f"Processed {n} contigs.")

    # close files
    unmatched.close()
    if not args.no_search:
        search_csv.close()
        rank_csv.close()
        search_matches.close()
        rank_matches.close()
    if args.gather:
        gather_csv.close()
        rankgather_csv.close()

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
    # switch search types
    p.add_argument('--no-search', help='do not run search with containment',
                   action='store_true', default=False)
    p.add_argument('--gather', help='run sourmash gather',
                   action='store_true', default=False)
    p.add_argument('--gather-min-matches', type=int, default=3)

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
