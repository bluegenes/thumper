#! /usr/bin/env python
"""
Create a csv with taxonomic annotation results per contig.
inputs: json gather results against single db (same alpha-ksize-scaled)
outputs: contig annotation csv
script modified from version in charcoal c7fd192
"""

import sys
import argparse
import csv
import os.path
from collections import Counter
import json

import sourmash
from sourmash.lca.command_index import load_taxonomy_assignments
from sourmash.lca import LCA_Database, LineagePair, lca_utils

from thumper import charcoal_utils as utils
from thumper.search_utils import gather_guess_tax_at_rank #gather_guess_tax_at_each_rank
from thumper.lineage_db import LineageDB
from thumper.charcoal_utils import (gather_at_rank, get_ident, summarize_at_rank,
                                    pretty_print_lineage, load_contigs_gather_json,
                                    is_contig_contaminated, is_contig_clean, kb)


# defaults from charcoal
GATHER_MIN_MATCHES=3
F_IDENT_THRESHOLD=0.1
F_MAJOR_THRESHOLD=0.2


def calculate_contam(genome_lin, contigs_d, rank, filter_names=None):
    "Calculate not-bad bp at each rank. Be conservative."
    good_names = dict()
    bad_names = dict()

    for contig_name, gather_info in contigs_d.items():
        contig_taxlist = gather_info.gather_tax
        if filter_names and contig_name in filter_names:
            continue

        if is_contig_contaminated(genome_lin, contig_taxlist, rank, GATHER_MIN_MATCHES):
            bad_names[contig_name] = gather_info
        else:
            good_names[contig_name] = gather_info

    return (good_names, bad_names)


# clean calculates the num bp that match to genome lineage. Useful as a confidence metric?
def calculate_clean(genome_lin, contigs_d, rank):
    "Calculate definitely-clean bp, as opposed to not-bad bp."
    good_names = dict()
    bad_names = dict()

    for contig_name, gather_info in contigs_d.items():
        contig_taxlist = gather_info.gather_tax

        if not is_contig_contaminated(genome_lin, contig_taxlist, rank, GATHER_MIN_MATCHES):
            good_names[contig_name] = gather_info
        else:
            bad_names[contig_name] = gather_info

    return (good_names, bad_names)

def guess_tax_by_gather(gather_results, num_hashes, match_rank, report_fp, minimum_matches=GATHER_MIN_MATCHES):
    "Guess likely taxonomy using gather."
    sum_ident = 0
    first_lin = ()
    first_count = 0
    comment=""
    if num_hashes < minimum_matches:
        comment= "insufficient hashes for classification"
    # if match_rank is genus, this should be same as using gather_results
    rank_gather = summarize_at_rank(gather_results, match_rank)
    for lin, count in rank_gather:
        if count >= minimum_matches:
            # record the first lineage we come across as likely lineage.
            if not first_lin:
                first_lin = lin
                first_count = count
        else:
            comment= "insufficient matched hashes"

        sum_ident += count

    if not first_lin:
        return "", 0.0, 0.0, ""

    f_ident = sum_ident / num_hashes
    f_major = first_count / sum_ident

    return first_lin, f_ident, f_major, comment

def get_genome_taxonomy(genome_name, genome_gather_json_filename, match_rank, min_f_ident, min_f_major):

    guessed_genome_lineage, f_major, f_ident, comment = "", 0.0, 0.0, ""
    # did we get gather results?
    genome_info = utils.load_contigs_gather_json(genome_gather_json_filename)

    if genome_info:
        gather_results = genome_info[genome_name].gather_tax
        genome_len = genome_info[genome_name].length
        genome_hashes = genome_info[genome_name].num_hashes

        # calculate lineage from majority vote on LCA
        guessed_genome_lineage, f_major, f_ident, comment = guess_tax_by_gather(gather_results, genome_hashes, match_rank, sys.stdout)

        print(f'Gather classification on this genome yields: {pretty_print_lineage(guessed_genome_lineage)}')

        if f_major == 1.0 and f_ident == 1.0:
            comment = "All genome hashes belong to one lineage! Nothing to do."
            print(comment)

    return guessed_genome_lineage, f_major, f_ident, comment


# gather results from each alpha-ksize?
def main(args):
    "Main entry point for scripting. Use cmdline for command line entry."

    # load genome gather matches
    jsoninfo = [line.strip().split(',') for line in open(args.jsoninfo_file, 'r')]

    json_dir = os.path.dirname(os.path.dirname(args.jsoninfo_file))

    num_genomes = len(jsoninfo)
    summary_d = {}
    detected_contam = {}

    for n, [genome_name, database, genome_json, contigs_json] in enumerate(jsoninfo[1:]):
        genome_jsonfile = os.path.join(json_dir, genome_json)
        contigs_jsonfile = os.path.join(json_dir, contigs_json)
        alpha,ksize,scaled = database.split(".", 1)[1].split("-")
        scaled = int(scaled.split("scaled")[1])
        ksize = int(ksize.split("k")[1])
        print(f'examining {genome_name} for classification accuracy ({n+1} of {num_genomes})')

        # assign genome lineage -- CURRENTLY NOT USING MIN_F_IDENT, MIN_F_MAJOR?
        genome_lineage,f_ident,f_major, comment = get_genome_taxonomy(genome_name, genome_jsonfile, args.match_rank, args.min_f_ident, args.min_f_major)

        vals = {}
        vals['genome'] = genome_name
        vals['f_ident'] = f_ident
        vals['f_major'] = f_major
        vals['comment'] = comment
        vals['lineage'] = sourmash.lca.display_lineage(genome_lineage)
        #vals['needs_lineage_flag'] = 1 if needs_lineage else 0
        #vals['mismatch_at'] = filter_at

        # calculate summary stats for contigs
        nohash_bp = 0
        nohash_count = 0
        noident_bp = 0
        noident_count = 0
        contigs_n = 0
        contigs_bp = 0

        # load contig gather matches
        contigs_d = load_contigs_gather_json(contigs_jsonfile)

        # guess taxonomy per contig
        for contig_name, gather_info in contigs_d.items():
            contigs_n += 1
            contigs_bp += gather_info.length

            if not gather_info.num_hashes:
                nohash_bp += gather_info.length
                nohash_count += 1
            elif not gather_info.gather_tax or not genome_lineage:
                noident_bp += gather_info.length
                noident_count += 1

        vals['ignored_contigs_n'] = nohash_count
        vals['ignored_contigs_bp'] = nohash_bp
        vals['noident_contigs_n'] = noident_count
        vals['noident_contigs_bp'] = noident_bp
        vals['total_contigs_n'] = contigs_n
        vals['total_contigs_bp'] = contigs_bp

        # track contigs that have been eliminated at various ranks
        for rank in sourmash.lca.taxlist():
            (good_names, bad_names) = calculate_contam(genome_lineage,
                                                       contigs_d,
                                                       rank)

            bad_n = len(bad_names)
            bad_bp = sum([ x.length for x in bad_names.values() ])

            print(f'   {rank}: {len(bad_names)} contigs w/ {kb(bad_bp)}kb')
            vals[f'bad_{rank}_bp'] = bad_bp
            vals[f'bad_{rank}_n'] = bad_n

            (good_names, bad_names) = calculate_clean(genome_lineage,
                                                      contigs_d,
                                                      rank)
            good_n = len(good_names)
            good_bp = sum([ x.length for x in good_names.values() ])
            vals[f'good_{rank}_bp'] = good_bp
            vals[f'good_{rank}_n'] = good_n

            assert bad_bp + good_bp == contigs_bp

            # track contamination between source (genome) / target (contig)
            x = []
            for contig_name in bad_names:
                contig_taxlist = contigs_d[contig_name].gather_tax
                for hit, count in contig_taxlist:
                    if utils.is_lineage_match(genome_lineage, hit, rank):
                        continue

                    # contam!
                    source_lin = utils.pop_to_rank(genome_lineage, rank)
                    target_lin = utils.pop_to_rank(hit, rank)

                    x.append((source_lin, target_lin, count))

#                    target = detected_contam.get(source_lin, Counter())
#                    target[target_lin] += count
#                    detected_contam[source_lin] = target
            detected_contam[genome_name] = x


            if rank == args.match_rank:
                break


        vals['total_bad_bp'] = vals['bad_genus_bp']

        print(f"   (total): {vals['bad_genus_n']} contigs w/ {kb(vals['bad_genus_bp'])}kb")

        summary_d[genome_name] = vals

        ###

    # output a sorted hit list CSV
    fp = open(args.output_csv, 'wt')
    hitlist_w = csv.writer(fp)

    hitlist_w.writerow(['genome', 'total_bad_bp', 'superkingdom_bad_bp', 'phylum_bad_bp',
                        'class_bad_bp', 'order_bad_bp', 'family_bad_bp', 'genus_bad_bp',
                        'f_ident', 'f_major', 'lineage', 'comment'])

    summary_items = list(summary_d.items())
    summary_items.sort(key=lambda x: -x[1]["total_bad_bp"])

    for genome_name, vals in summary_items:
        hitlist_w.writerow([genome_name,
                            #vals['filter_at'], '',
                            vals["total_bad_bp"],
                            vals['bad_superkingdom_bp'],
                            vals['bad_phylum_bp'],
                            vals['bad_class_bp'],
                            vals['bad_order_bp'],
                            vals['bad_family_bp'],
                            vals['bad_genus_bp'],
                            f'{vals["f_ident"]:.03}',
                            f'{vals["f_major"]:.03}',
                            vals["lineage"],
                            vals["comment"]])

    fp.close()

    # output a sorted summary CSV with a lot more information!
    fp = open(args.contig_details_summary, 'wt')

    # build column list; put genome first
    vals = summary_items[0][1]
    all_columns = set(vals.keys())
    all_columns.remove('genome')
    all_columns = list(sorted(all_columns))
    all_columns = ['genome'] + all_columns

    w = csv.DictWriter(fp, fieldnames=all_columns)
    w.writeheader()

    for genome, vals in summary_items:
        w.writerow(vals)

    if args.lineages_for_charcoal:
        with open(args.lineages_for_charcoal, "w") as fp:
            header = ["genome", "lineage"] # is header allowed?
            fp.write(",".join(header) + "\n")
            for genome, vals in summary_items:
                fp.write(genome + "," + vals["lineage"] + "\n")

    fp.close()

    ####

    print(f"processed {num_genomes} genomes.")

    print(f"saving contamination summary to {args.contam_summary_json}")
    with open(args.contam_summary_json, 'wt') as fp:
        utils.save_contamination_summary(detected_contam, fp)

    return 0


# from file: write a csv file of
# genome_name,genome_json,contig_json


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser(sys_args)
    p.add_argument('--jsoninfo-file', required=True) # file with list of json files to aggregate
    p.add_argument('--output-csv', required=True)
    p.add_argument('--lineages-for-charcoal')
    p.add_argument('--contig-details-summary', required=True)
    p.add_argument('--contam-summary-json', required=True)
    p.add_argument('--match-rank', default="genus")
    p.add_argument('--gather_min_matches', type=float, default=GATHER_MIN_MATCHES)
    p.add_argument('--min_f_ident', type=float, default=F_IDENT_THRESHOLD)
    p.add_argument('--min_f_major', type=float, default=F_MAJOR_THRESHOLD)
    args = p.parse_args()

    acceptable_ranks =  ["superkingdom", "phylum", "class", "order", "family", "genus"]
    assert args.match_rank in acceptable_ranks, f"Error, match_rank must be one of {','.join(acceptable_ranks)}"

    return main(args)


# execute this, when run with `python -m`.
if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)

