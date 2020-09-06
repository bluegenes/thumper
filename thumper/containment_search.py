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

import screed

import sourmash
from sourmash.lca.command_index import load_taxonomy_assignments
from sourmash.lca import LCA_Database

from .lineage_db import LineageDB
from .version import version
from thumper.charcoal_utils import (gather_at_rank, get_ident, ContigSearchInfo, search_containment_at_rank)


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
    #contigs_tax = {}

    # need five files
    unmatched = open(args.no_matches, "w")
    search_csv = open(args.search_csv, "w")
    #search_fieldnames = ['similarity', 'name', 'filename', 'md5', 'lineage']
    rank_fieldnames = ['contig_name', 'length', 'similarity', 'name', 'md5', 'lineage']
    #'similarity, match, md5, filename, name, lineage')
    #'lineage, containment, intersect_bp, match_sig')
    search_w =csv.DictWriter(search_csv, fieldnames=search_fieldnames)
    search_w.writeheader()

    rank_csv = open(args.ranksum_csv, "w")
    rank_fieldnames = ['contig_name', 'contig_length', 'match_rank', 'lineage', 'contained_at_rank', 'contained_bp']
    # for rank results: contig_name,contig_length,match_rank,lineage,contained_at_rank,bp_contained #also num_hashes?
    rank_w = csv.DictWriter(rank_csv, fieldnames=rank_fieldnames)
    rank_w.writeheader()

    search_matches = open(args.search_matches, "w")
    rank_matches = open(args.ranksum_matches, "w")

    for n, record in enumerate(screed_iter):
        # look at each contig individually
        mh = empty_mh.copy_and_clear()
        mh.add_sequence(record.sequence, force=True)
        # search, optionally aggregate matched hashes to get containment at rank

        contig_name = record.name
        contig_len = len(record.sequence)
        num_hashes = len(mh.hashes)

        search_results, search_rank_results = search_containment_at_rank(mh, lca_db, lin_db, match_rank)

        # first, print normal search --containment results

        if not search_results:
            unmatched.write(record.name + "\n" + record.sequence + "\n")
            # write to unclassified
        else:
            for sr in search_results:
                d = dict(sr._asdict())
                # save match to output matches
                sourmash.signature.save_signatures([d['match']], search_matches)
                del d['match']
                # better way to do this?
                d["contig_name"] = contig_name
                d["length"] = contig_len
                w.writerow(d)

            # now, print containment at rank results
            for sr in search_rank_results:
                d = dict(sr._asdict())
                # save match to output matches
                sig.save_signatures([d['match_sig']], rank_matches)
                del d['match_sig']
                d["contig_name"] = contig_name
                d["length"] = contig_len
                d["rank"] = sr.lineage[-1]
                w.writerow(d)
                report_csv.write(res)

    print(f"Processed {n} contigs.")

    # close five files
    unmatched.close()
    search_csv.close()
    rank_csv.close()
    search_matches.close()
    rank_matches.close()

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

    p.add_argument('--json-out',
                   help='JSON-format output file of all tax results',
                   required=True)
    args = p.parse_args()

    return main(args)


# execute this, when run with `python -m`.
if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
