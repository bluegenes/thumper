import os.path
from . import pytest_utils as utils
import csv
from sourmash import signature as sig
from thumper import search_or_gather

# to do: use proteome subsets instead for smaller testdata. Then add testdata to repo.

@utils.in_tempdir
def test_empty_1(location):
    # test an empty set of matches (once self is removed)
    # TODO: change testdata .. this is not empty
    args = utils.Args()
    args.genome = utils.test_file("test-data/proteomes/GB_GCA_000384615.1_protein.faa.gz")
    args.genome_sig = utils.test_file("test-data/intermediate/GB_GCA_000384615.1_protein.faa.gz.sig")
    args.matches_sig = utils.test_file("test-data/intermediate/GB_GCA_000384615.1_protein.faa.gz.x.gtdb-nine.protein-k11.matches.sig")
    args.lineages_csv = utils.test_file("test-data/databases/gtdb-nine.lineages.csv")
    args.alphabet = "protein"
    args.ksize = 33
    args.output_prefix = "GB_GCA_000384615.1_protein.faa.gz.x.gtdb-nine.protein-k11"
    args.no_search=False
    args.gather=False

    # search_outfiles --> can we use these in all tests?
    search_csv = os.path.join(location, f"{args.output_prefix}.search.csv")
    ranksearch_csv = os.path.join(location, f"{args.output_prefix}.ranksearch.csv")
    search_matches = os.path.join(location, f"{args.output_prefix}.search.matches.sig")
    ranksearch_matches = os.path.join(location, f"{args.output_prefix}.ranksearch.matches.sig")

    outfiles = [search_csv, ranksearch_csv, search_matches, ranksearch_matches]
    status = search_or_gather.main(args)

    for outF in outfiles:
        assert os.path.exists(outF)

    with open(search_csv) as fp:
        sr = csv.reader(fp)
        firstline = sr.__next__()
        assert firstline == ['name', 'length', 'similarity', 'name', 'filename', 'md5', 'lineage']
    with open(ranksearch_csv) as rp:
        sr = csv.reader(rp)
        firstline = sr.__next__()
        assert firstline == ['name', 'length', 'match_rank', 'lineage', 'contained_at_rank', 'contained_bp']


def get_csv_set(f):
    return set(map(tuple, csv.reader(f)))

# test contig search w rank
@utils.in_tempdir
def test_contig_search_1(location):
    # test for same results
    args = utils.Args()
    args.genome = utils.test_file("test-data/proteomes/GB_GCA_002691795.1_protein.faa.gz")
    args.genome_sig = utils.test_file("test-data/intermediate/GB_GCA_002691795.1_protein.faa.gz.sig")
    args.matches_sig = utils.test_file("test-data/intermediate/GB_GCA_002691795.1_protein.faa.gz.x.gtdb-nine.protein-k11.matches.sig")
    args.lineages_csv = utils.test_file("test-data/databases/gtdb-nine.lineages.csv")
    args.alphabet = "protein"
    args.ksize = 33
    args.output_prefix = "GB_GCA_002691795.1_protein.faa.gz.x.gtdb-nine.protein-k11"
    args.no_search=False
    args.gather=False

    search_csv = os.path.join(location, f"{args.output_prefix}.search.csv")
    ranksearch_csv = os.path.join(location, f"{args.output_prefix}.ranksearch.csv")
    search_matches = os.path.join(location, f"{args.output_prefix}.search.matches.sig")
    ranksearch_matches = os.path.join(location, f"{args.output_prefix}.ranksearch.matches.sig")
    outfiles = [search_csv, ranksearch_csv, search_matches, ranksearch_matches]
    status = search_or_gather.main(args)
    assert status == 0

    for outF in outfiles:
        assert os.path.exists(outF)

    saved_search_csv = \
    utils.test_file("test-data/intermediate/GB_GCA_002691795.1_protein.faa.gz.x.gtdb-nine.protein-k11.search.csv")
    with open(saved_search_csv) as fp:
        saved_search_csvset = get_csv_set(fp)
    with open(search_csv) as fp:
        this_search_csvset = get_csv_set(fp)
    assert saved_search_csvset == this_search_csvset

    saved_ranksearch_csv = \
    utils.test_file("test-data/intermediate/GB_GCA_002691795.1_protein.faa.gz.x.gtdb-nine.protein-k11.ranksearch.csv")
    with open(saved_ranksearch_csv) as fp:
        saved_ranksearch_csvset = get_csv_set(fp)
    with open(ranksearch_csv) as fp:
        this_ranksearch_csvset = get_csv_set(fp)
    assert saved_ranksearch_csvset == this_ranksearch_csvset

    saved_search_matches = \
    utils.test_file("test-data/intermediate/GB_GCA_002691795.1_protein.faa.gz.x.gtdb-nine.protein-k11.search.matches.sig")
    with open(saved_search_matches) as sm:
        saved_search_sigs = set(sig.load_signatures(sm))
    with open(search_matches) as sm:
        this_search_sigs = set(sig.load_signatures(sm))
    assert saved_search_sigs == this_search_sigs

    saved_ranksearch_matches = \
    utils.test_file("test-data/intermediate/GB_GCA_002691795.1_protein.faa.gz.x.gtdb-nine.protein-k11.ranksearch.matches.sig")
    with open(saved_ranksearch_matches) as rm:
        saved_ranksearch_sigs = set(sig.load_signatures(rm))

    with open(ranksearch_matches) as rm:
        this_ranksearch_sigs = set(sig.load_signatures(rm))
    assert saved_ranksearch_sigs == this_ranksearch_sigs


# test contig search w rank + gather

# test contig gather only (no search)

# test MAG search

# test MAG search + gather

# test MAG gather only





