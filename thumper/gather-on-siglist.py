#! /usr/bin/env python
"""
gather-on-siglist.py:: Gather from a list of signatures. Optionally summarize taxonomy.
Authors:
  - N. Tessa Pierce Ward, github.com/bluegenes, ntpierce@gmail.com

"""
import os
import sys
import argparse
import csv
import json
from collections import defaultdict, namedtuple
import pandas as pd

import screed
import sourmash
from sourmash.logging import notify
from sourmash.lca import LCA_Database, lca_utils
from sourmash.lca.command_index import load_taxonomy_assignments

from .lineage_db import LineageDB
from .version import version
from thumper.charcoal_utils import (gather_at_rank, get_ident, ContigGatherInfo)
from .search_utils import (gather_guess_tax_at_each_rank, search_containment_at_rank, SearchFiles)

rareInfo = namedtuple('RarefactionInfo','num_founders, num_members')

def load_sigs_from_list(siglistfiles, moltype, ksize):
    # input lists of signatures instead
    sigs = []
    for sl in siglistfiles:
        notify(f'loading from {sl}')
        sigfiles = sourmash.sourmash_args.load_file_list_of_signatures(sl)
        new_sigs = load_sigs(sigfiles, moltype, ksize, source_type= "input sigfile list")
        notify(f'...got {len(new_sigs)} signatures from {sl} siglist file.')
        sigs+=new_sigs
    return sigs

def load_sigs(sig_sources, moltype, ksize, source_type="input sigfiles"):
    siglist=[]
    for filename in sig_sources:
        if source_type != "input sigfile list":
            notify(f'loading from {filename}')
        m = 0
        for sig in sourmash.sourmash_args.load_file_as_signatures(filename,
                                           select_moltype=moltype,
                                           ksize=ksize):
            m += 1
            siglist.append((filename, sig))
        if source_type != "input sigfile list":
            notify(f'...got {m} signatures from {source_type}.')
    return siglist

def remove_identical_match(query_sig, siglist):
    # Hack for examining members of our search database: remove exact matches.
    new_siglist = []
    for ss in siglist:
        if query_sig.similarity(ss) == 1.0:
            print(f'removing an identical match: {str(ss)}')
        else:
            new_siglist.append(ss)
    return new_siglist


def main(args):

    genomebase = os.path.basename(args.genome)
    match_rank = 'genus'

    # load taxonomy CSV
    tax_assign, _ = load_taxonomy_assignments(args.lineages_csv,
                                              start_column=2)
    print(f'loaded {len(tax_assign)} tax assignments.')

    # load the genome signature
    genome_sig = sourmash.load_file_as_signatures(args.genome_sig, select_moltype=args.alphabet,ksize=args.ksize)

    # load siglist database signatures (if using prefetch, this will be matches only!)
    siglist=load_sigs_from_list(args.siglist_db, args.alphabet, args.ksize)
    notify(f'loaded {len(siglist)} database signatures.')

    # rm identical match
    siglist = remove_identical_match(genome_sig, siglist)
    gather_tax = {}
    if not siglist:
        # write empty files so snakemake workflows don't complain; exit.
        print('no non-identical matches for this genome, exiting.')
        # todo: write gather json and csv!!

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
    gather_tax = {}

    # do I need genome length here? Can I use a csv of lengths instead?
    if not genome_len:
        for record in screed_iter:
            genome_len+=len(record.sequence)

    # query minhash
    entire_mh = genome_sig.minhash
    scaled = entire_mh.scaled
    genome_name = str(genome_sig)
    num_hashes = len(entire_mh.hashes)

    # run gather (ignore 100% identical matches!)
    gather_results = list(gather_at_rank(entire_mh, lca_db, lin_db, match_rank))

    import pdb;pdb.set_trace()

    genome_gather_info = ContigGatherInfo(genome_len, len(entire_mh), gather_results)
    genome_gather_tax[genome_name] = genome_gather_info
    # next, summarize at higher ranks
    gather_taxonomy_per_rank = gather_guess_tax_at_each_rank(gather_results, num_hashes, scaled, minimum_matches=args.gather_min_matches,lowest_rank=match_rank) #taxlist=lca_utils.taxlist(include_strain=False))


    # output: gather json, gather csv, rankgather csv
    prefix = args.prefix
    with open(f'{prefix}.founders.siglist.txt', 'wt') as fp:
        for (founder_from, founder) in founders:
            fp.write(founder_from + "\n")
    with open(f'{prefix}.founders.siglist.csv', 'wt') as fp:
        for (founder_from, founder) in founders:
            fp.write(f"{str(founder)},{founder_from}\n")
    with open(f'{prefix}.members.siglist.txt', 'wt') as fp:
        for (member_from, member) in members:
            fp.write(member_from + "\n")
    with open(f'{prefix}.members.siglist.csv', 'wt') as fp:
        for (member_from, member) in members:
            fp.write(f"{str(member)},{member_from}\n")



def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser(sys_args)
    p.add_argument("--siglist-db", action="append", help="provide siglist database")
    p.add_argument('--genome', help='genome file', required=True)
    p.add_argument('--genome-sig', help='genome sig', required=True)
    p.add_argument('--lineages-csv', help='lineage spreadsheet', required=True)
    p.add_argument('--alphabet', help='alphabet',  default='protein')
    p.add_argument('-k', '--ksize', type=int, default=10)
    p.add_argument('--gather-min-matches', type=int, default=3)

    p.add_argument('--force', help='continue past survivable errors',
                   action='store_true')
    # output options:
    p.add_argument('--output-prefix',
                    help='prefix for search or gather results',
                    required=True)
    args = p.parse_args()

    return main(args)


# execute this, when run with `python -m`.
if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
