import os.path
from . import pytest_utils as utils
import csv
from sourmash import signature as sig
from thumper import compare_taxonomy

# to do: use proteome subsets instead for smaller testdata. Then add testdata to repo.

@utils.in_tempdir
def test_compare_taxonomy(location):
    # test an empty set of matches (once self is removed)
    args = utils.Args()
    #args.jsoninfo_file = utils.get_testfile("output.test-prot/classify/test_prot_gtdb-nine.protein-k11-scaled10.gather.txt")
    args.jsoninfo_file = utils.get_testfile("test-data/intermediate/classify/test_prot_gtdb-nine.protein-k11-scaled10.gather.txt")
    args.gather_min_matches = 3
    args.min_f_ident=0.1
    args.min_f_major=0.2
    args.match_rank="genus"
    args.output_csv = os.path.join(location, "GB_GCA_000384615.1_protein.100contigs.faa.gz.x.gtdb-nine.protein-k11-scaled10.gather.classify.csv")
    args.contam_summary_json = os.path.join(location, "GB_GCA_000384615.1_protein.100contigs.faa.gz.x.gtdb-nine.protein-k11-scaled10.contam-summary.json")
    args.contig_details_summary = os.path.join(location, "GB_GCA_000384615.1_protein.100contigs.faa.gz.x.gtdb-nine.protein-k11-scaled10.contig-details.csv")

    status = compare_taxonomy.main(args)

    assert os.path.exists(args.output_csv)

    with open(args.output_csv) as fp:
        out = csv.reader(fp)
        firstline = out.__next__()
        header = ['genome', 'total_mismatched_bp', 'superkingdom_mismatched_bp', 'phylum_mismatched_bp', 'class_mismatched_bp', 'order_mismatched_bp', 'family_mismatched_bp', 'genus_mismatched_bp', 'f_ident', 'f_major', 'lineage', 'comment']
        assert firstline == header

