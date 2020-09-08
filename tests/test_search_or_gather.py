import os.path
from . import pytest_utils as utils
import json

from thumper import search_or_gather


@utils.in_tempdir
def test_1(location):
    # test an empty set of matches (once self is removed)
    args = utils.Args()
    args.genome = utils.test_file("test-data/proteomes/GB_GCA_002691795.1_protein.faa.gz")
    args.genome_sig = utils.test_file("test-data/GB_GCA_002691795.1_protein.faa.gz.sig")
    args.matches_sig = utils.test_file("test-data/GB_GCA_002691795.1_protein.faa.gz.x.gtdb-nine.protein-k11.matches.sig")
    args.lineages_csv = utils.test_file("test-data/databases/gtdb-nine.lineages.csv")
    args.alphabet = "protein"
    args.ksize = 33
    args.output_prefix = "GB_GCA_002691795.1_protein.faa.gz.x.gtdb-nine.protein-k11"
    args.no_search=False
    args.gather=False

    #args.json_out = os.path.join(location, 'tax.json')

    status = search_or_gather.main(args)
    import pdb;pdb.set_trace()
    assert status == 0
    assert os.path.exists(args.json_out)

    with open(args.json_out, 'rt') as fp:
        results = json.load(fp)
        assert results == {}

