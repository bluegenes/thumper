"""
test thumper/search_utils.py
"""
import json
from collections import defaultdict, Counter, namedtuple
from operator import itemgetter
import csv

import sourmash
from sourmash.lca import lca_utils, LineagePair, LCA_Database, taxlist
from thumper.lineage_db import LineageDB
from thumper.charcoal_utils import get_ident, pop_to_rank
# to do: specific imports!
from thumper.search_utils import *


def make_mh(hashvals, ksize=3, scaled=1):
    mh = sourmash.MinHash(n=0, scaled=1, ksize=3)
    mh.add_many(hashvals)
    return mh


def make_sig_and_lin(hashvals, ident, lin, ksize=3, scaled=1):
    mh = make_mh(hashvals)
    sig = sourmash.SourmashSignature(mh, name=ident)
    lineage = lca_utils.make_lineage(lin)
    return mh, sig, lineage


def test_gen_mh():
    mh = make_mh([12345678])
    return mh.copy_and_clear()


def test_add_hashes_at_ranks_1():
    lin1 = lca_utils.make_lineage('a')
    hashval1  = 12345678
    lineage_hashD=defaultdict(test_gen_mh)
    # manually add hashval1 to a
    lineage_hashD[lin1].add_many([hashval1])

    # test that hashval2 gets added to lin1
    lin2 = lca_utils.make_lineage('a;b')
    hashval2 = 87654321
    match_rank="genus"
    lineage_hashD = add_hashes_at_ranks(lineage_hashD, [hashval2], lin2, match_rank)

    assert set(lineage_hashD[lin1].hashes) == set([hashval1, hashval2])
    assert set(lineage_hashD[lin2].hashes) == set([hashval2])


def test_add_hashes_at_ranks_2():
    lin1 = lca_utils.make_lineage('a;b')
    lin2 = lca_utils.make_lineage('a')
    hashval1  = 12345678
    lineage_hashD=defaultdict(test_gen_mh)
    # manually add hashval1 to a;, a;b
    lineage_hashD[lin1].add_many([hashval1])
    lineage_hashD[lin2].add_many([hashval1])

    # test that hashval2 gets added to *both* lin1, lin2
    lin3 = lca_utils.make_lineage('a;b;c')
    hashval2 = 87654321
    match_rank="genus"
    lineage_hashD = add_hashes_at_ranks(lineage_hashD, [hashval2], lin3, match_rank)

    assert set(lineage_hashD[lin1].hashes) == set([hashval1, hashval2])
    assert set(lineage_hashD[lin2].hashes) == set([hashval1, hashval2])
    assert set(lineage_hashD[lin3].hashes) == set([hashval2])


def test_add_hashes_at_ranks_3():
    lin1 = lca_utils.make_lineage('a;b')
    lin2 = lca_utils.make_lineage('a')
    hashval1  = 12345678
    lineage_hashD=defaultdict(test_gen_mh)
    # manually add hashval1 to a;, a;b
    lineage_hashD[lin1].add_many([hashval1])
    lineage_hashD[lin2].add_many([hashval1])

    # test that hashval2 is added separately
    lin3 = lca_utils.make_lineage('d;e;f')
    lin4 = lca_utils.make_lineage('d;e')
    lin5 = lca_utils.make_lineage('d')
    hashval2 = 87654321
    match_rank="genus"
    lineage_hashD = add_hashes_at_ranks(lineage_hashD, [hashval2], lin3, match_rank)

    assert set(lineage_hashD[lin1].hashes) == set([hashval1])
    assert set(lineage_hashD[lin2].hashes) == set([hashval1])
    assert set(lineage_hashD[lin3].hashes) == set([hashval2])
    assert set(lineage_hashD[lin4].hashes) == set([hashval2])
    assert set(lineage_hashD[lin5].hashes) == set([hashval2])


def test_add_hashes_at_ranks_4():
    lin1 = lca_utils.make_lineage('a;b')
    lin2 = lca_utils.make_lineage('a')
    hashval1  = 12345678
    lineage_hashD=defaultdict(test_gen_mh)
    # manually add hashval1 to a;, a;b
    lineage_hashD[lin1].add_many([hashval1])
    lineage_hashD[lin2].add_many([hashval1])

    # test that hashval2 is added appropriately
    lin3 = lca_utils.make_lineage('a;d')
    hashval2 = 87654321
    match_rank="genus"
    lineage_hashD = add_hashes_at_ranks(lineage_hashD, [hashval2], lin3, match_rank)

    assert set(lineage_hashD[lin1].hashes) == set([hashval1])
    assert set(lineage_hashD[lin2].hashes) == set([hashval1, hashval2])
    assert set(lineage_hashD[lin3].hashes) == set([hashval2])


def test_get_lineage_at_match_rank():
    hashval  = 12345678
    ident1 = 'first'
    mh1, sig1, lin1 = make_sig_and_lin([hashval], ident1, 'a;b;c')
    # make lin_db
    lin_db = LineageDB()
    lin_db.insert(ident1, lin1)
    superkingdom = lineage = lca_utils.make_lineage('a')
    phylum = lineage = lca_utils.make_lineage('a;b')

    assert get_lineage_at_match_rank(lin_db, sig1, "phylum") == phylum
    assert get_lineage_at_match_rank(lin_db, sig1, "superkingdom") == superkingdom
    assert get_lineage_at_match_rank(lin_db, sig1, "class") == lin1


def test_calculate_containment_at_rank_1():
    # one lineage, calculate query containment at rank
    hashval1  = 12345678
    ident = 'uniq'
    mh1, sig1, lin1 = make_sig_and_lin([hashval1], ident, 'a;b')
    match_rank="genus"
    # make lineage hashD
    lineage_hashD=defaultdict(test_gen_mh)
    lineage_hashD = add_hashes_at_ranks(lineage_hashD, [hashval1], lin1, match_rank)

    # calculate containment
    containmentD = calculate_containment_at_rank(lineage_hashD, sig1, match_rank)
    # superkingdom lineage that should have 100% containment
    assert containmentD["superkingdom"][0][1] == 1.0
    assert containmentD["superkingdom"][0][2] == 3
    assert containmentD["phylum"][0][1] == 1.0
    assert containmentD["phylum"][0][2] == 3


def test_calculate_containment_at_rank_2():
    # two lineages, match at phylum level
    hashval1  = 12345678
    ident = 'uniq'
    mh1, sig1, lin1 = make_sig_and_lin([hashval1], ident, 'a;b;c')
    lin2 = lca_utils.make_lineage('a;d')
    hashval2 = 87654321
    match_rank="genus"
    # make lineage hashD
    lineage_hashD=defaultdict(test_gen_mh)
    lineage_hashD = add_hashes_at_ranks(lineage_hashD, [hashval1], lin1, match_rank)
    lineage_hashD = add_hashes_at_ranks(lineage_hashD, [hashval2], lin2, match_rank)

    # calculate containment
    containmentD = calculate_containment_at_rank(lineage_hashD, sig1, match_rank)

    # superkingdom lineage that should have 100% containment
    lin3 = lca_utils.make_lineage('a')
    assert containmentD["superkingdom"][0][1] == 1.0
    assert containmentD["class"][0][1] == 1.0
    phylum_containment = set([containmentD["phylum"][0][1], containmentD["phylum"][1][1]])
    assert set([0.0, 1.0]) == phylum_containment


def test_calculate_containment_at_rank_3():
    # two lineages with overlapping hashes (50% containment)
    hashval1  = 12345678
    ident = 'uniq'
    mh1, sig1, lin1 = make_sig_and_lin([hashval1], ident, 'a;b;c')
    lin2 = lca_utils.make_lineage('a;d')
    hashval2 = 87654321
    match_rank="genus"
    # make lineage hashD
    lineage_hashD=defaultdict(test_gen_mh)
    lineage_hashD = add_hashes_at_ranks(lineage_hashD, [hashval1], lin1, match_rank)
    lineage_hashD = add_hashes_at_ranks(lineage_hashD, [hashval2], lin2, match_rank)

    # make query sig
    mh = make_mh([hashval1,hashval2])
    query_sig = sourmash.SourmashSignature(mh, name='query')

    # calculate containment
    containmentD = calculate_containment_at_rank(lineage_hashD, query_sig, match_rank)

    # superkingdom lineage that should have 100% containment
    lin3 = lca_utils.make_lineage('a')
    assert containmentD["superkingdom"][0][1] == 1.0
    # class should have 50% containment
    assert containmentD["class"][0][1] == 0.5
    assert containmentD["class"][0][2] == 3
    phylum_containment = [containmentD["phylum"][0][1], containmentD["phylum"][1][1]]
    assert [0.5, 0.5] == phylum_containment


def test_calculate_containment_at_rank_4():
    # add two (nonmatching) hashvals to query
    hashval1  = 12345678
    ident = 'uniq'
    mh1, sig1, lin1 = make_sig_and_lin([hashval1], ident, 'a;b;c')
    lin2 = lca_utils.make_lineage('a;d')
    hashval2 = 87654321
    match_rank="genus"
    # make lineage hashD
    lineage_hashD=defaultdict(test_gen_mh)
    lineage_hashD = add_hashes_at_ranks(lineage_hashD, [hashval1], lin1, match_rank)
    lineage_hashD = add_hashes_at_ranks(lineage_hashD, [hashval2], lin2, match_rank)

    # make query sig
    mh = make_mh([hashval1,hashval2, 33333333, 44444444])
    query_sig = sourmash.SourmashSignature(mh, name='query')

    # calculate containment
    containmentD = calculate_containment_at_rank(lineage_hashD, query_sig, match_rank)

    # superkingdom lineage that should have 50% containment
    lin3 = lca_utils.make_lineage('a')
    assert containmentD["superkingdom"][0][1] == 0.5
    assert containmentD["superkingdom"][0][2] == 6
    # each class should have 25% containment
    assert containmentD["class"][0][1] == 0.25
    assert containmentD["class"][0][2] == 3
    assert [containmentD["phylum"][0][1], containmentD["phylum"][1][1]] == [0.25, 0.25]


def test_sort_by_rank_and_containment_1():
    # 1. three results, check that they sort by rank, containment
    hashval1  = 12345678
    ident = 'uniq'
    mh1, sig1, lin1 = make_sig_and_lin([hashval1], ident, 'a;b;c')
    lin2 = lca_utils.make_lineage('a;d')
    hashval2 = 87654321
    match_rank="genus"
    # make lineage hashD
    lineage_hashD=defaultdict(test_gen_mh)
    lineage_hashD = add_hashes_at_ranks(lineage_hashD, [hashval1], lin1, match_rank)
    lineage_hashD = add_hashes_at_ranks(lineage_hashD, [hashval2], lin2, match_rank)
    # make query sig
    mh = make_mh([hashval1,hashval2, 33333333, 44444444])
    query_sig = sourmash.SourmashSignature(mh, name='query')
    superK_lin = lca_utils.make_lineage('a')
    phylum_match_lin = lca_utils.make_lineage('a;b')
    # calculate containment
    containmentD = calculate_containment_at_rank(lineage_hashD, query_sig, match_rank)
    sorted_results = sort_by_rank_and_containment(containmentD, match_rank)
    assert sorted_results[0].lineage == superK_lin
    assert sorted_results[0].contained_at_rank == 0.5
    # phylum results
    assert sorted_results[1].lineage[-1].rank == "phylum"
    assert sorted_results[2].lineage[-1].rank == "phylum"
    assert set([sorted_results[1].lineage, sorted_results[2].lineage]) == set([lin2,phylum_match_lin])
    assert sorted_results[1].contained_at_rank == 0.25
    assert sorted_results[2].contained_at_rank == 0.25
    # class results
    assert sorted_results[3].lineage[-1].rank == "class"
    assert sorted_results[3].contained_at_rank == 0.25


def test_sort_by_rank_and_containment_2():
    # 1. three results, check that they sort by rank, containment
    hashval1  = 12345678
    ident = 'uniq'
    mh1, sig1, lin1 = make_sig_and_lin([hashval1], ident, 'a;b;c')
    lin2 = lca_utils.make_lineage('a;d')
    hashval2 = 87654321
    hashval3 = 33333333
    match_rank="genus"
    # make lineage hashD
    lineage_hashD=defaultdict(test_gen_mh)
    lineage_hashD = add_hashes_at_ranks(lineage_hashD, [hashval1, hashval3], lin1, match_rank)
    lineage_hashD = add_hashes_at_ranks(lineage_hashD, [hashval2], lin2, match_rank)
    # make query sig
    mh = make_mh([hashval1,hashval2, hashval3, 44444444])
    query_sig = sourmash.SourmashSignature(mh, name='query')
    superK_lin = lca_utils.make_lineage('a')
    phylum_match_lin = lca_utils.make_lineage('a;b')
    # calculate containment
    containmentD = calculate_containment_at_rank(lineage_hashD, query_sig, match_rank)
    sorted_results = sort_by_rank_and_containment(containmentD, match_rank)
    assert sorted_results[0].lineage == superK_lin
    assert sorted_results[0].contained_at_rank == 0.75
    # phylum results should also be sorted by containment
    assert sorted_results[1].lineage[-1].rank == "phylum"
    assert sorted_results[1].contained_at_rank == 0.5
    assert sorted_results[2].lineage[-1].rank == "phylum"
    assert sorted_results[2].contained_at_rank == 0.25
    # class results
    assert sorted_results[3].lineage[-1].rank == "class"
    assert sorted_results[3].contained_at_rank == 0.5


def test_contain_at_rank_1():
    # one minhash, one set of ranks
    hashval  = 12345678
    ident = 'uniq'
    mh1, sig1, lin1 = make_sig_and_lin([hashval], ident, 'a;b;c')

    lca_db = LCA_Database(scaled=1, ksize=3)
    lca_db.insert(sig1, ident=ident)

    lin_db = LineageDB()
    lin_db.insert(ident, lin1)

    results, rank_results=search_containment_at_rank(mh1, lca_db, lin_db, "class")
    assert len(results) == 1
    assert results[0].lineage == lin1
    assert results[0].name == ident
    assert results[0].similarity == 1.0

    superk_lin = lca_utils.make_lineage('a')
    phylum_match_lin = lca_utils.make_lineage('a;b')
    assert len(rank_results) == 3
    assert rank_results[0].lineage == superk_lin
    assert rank_results[0].contained_at_rank == 1.0
    assert rank_results[1].lineage == phylum_match_lin
    assert rank_results[1].contained_at_rank == 1.0
    assert rank_results[1].contained_bp == 3
    assert rank_results[2].lineage == lin1
    assert rank_results[2].contained_at_rank == 1.0
    assert rank_results[2].contained_bp == 3

def test_contain_at_rank_2():
    #two minhashes, fully shared ranks

    # first sig
    hashval  = 12345678
    ident1 = 'first'
    mh1, sig1, lin1 = make_sig_and_lin([hashval], ident1, 'a;b;c')

    # second sig
    hashval2 = 87654321
    ident2 = 'second'
    mh2, sig2, lin2 = make_sig_and_lin([hashval2], ident2, 'a;b;c')

    # create lca_db w sigs
    lca_db = LCA_Database(scaled=1, ksize=3)
    lca_db.insert(sig1, ident=ident1)
    lca_db.insert(sig2, ident=ident2)

    # make lin_db
    lin_db = LineageDB()
    lin_db.insert(ident1, lin1)
    lin_db.insert(ident2, lin2)

    # search with combined hashvals
    search_mh = make_mh([hashval, hashval2])
    results, rank_results=search_containment_at_rank(search_mh, lca_db, lin_db, "class")

    assert len(results) == 2
    assert set([results[0].lineage, results[1].lineage])  == set([lin1])
    assert set([results[0].similarity, results[1].similarity])  == set([0.5])
    assert set([results[0].name, results[1].name]) == set([ident1, ident2])

    superk_lin = lca_utils.make_lineage('a')
    phylum_match_lin = lca_utils.make_lineage('a;b')
    assert len(rank_results) == 3
    assert rank_results[0].lineage == superk_lin
    assert rank_results[0].contained_at_rank == 1.0
    assert rank_results[1].lineage == phylum_match_lin
    assert rank_results[1].contained_at_rank == 1.0
    assert rank_results[1].contained_bp == 6
    assert rank_results[2].lineage == lin1
    assert rank_results[2].contained_at_rank == 1.0
    assert rank_results[2].contained_bp == 6

def test_contain_at_rank_3():
    # two minhashes, totally distinct ranks
    # first sig
    hashval1  = 12345678
    ident1 = 'first'
    mh1, sig1, lin1 = make_sig_and_lin([hashval1], ident1, 'a;b;c')

    # second sig
    hashval2 = 87654321
    ident2 = 'second'
    mh2, sig2, lin2 = make_sig_and_lin([hashval2], ident2, 'd;e;f')

    # create lca_db w sig1
    lca_db = LCA_Database(scaled=1, ksize=3)
    lca_db.insert(sig1, ident=ident1)
    lca_db.insert(sig2, ident=ident2)

    # next, make lin_db
    lin_db = LineageDB()
    lin_db.insert(ident1, lin1)
    lin_db.insert(ident2, lin2)

    # search with combined hashvals
    search_mh = make_mh([hashval1, hashval2])
    results, rank_results=search_containment_at_rank(search_mh, lca_db, lin_db, "class")

    assert len(results) == 2
    assert set([results[0].lineage, results[1].lineage])  == set([lin1, lin2])
    assert set([results[0].similarity, results[1].similarity])  == set([0.5])
    assert set([results[0].name, results[1].name]) == set([ident1, ident2])

    superk_lin1 = lca_utils.make_lineage('a')
    superk_lin2 = lca_utils.make_lineage('d')
    phylum_lin1 = lca_utils.make_lineage('a;b')
    phylum_lin2 = lca_utils.make_lineage('d;e')
    assert len(rank_results) == 6

    assert set([rank_results[0].lineage, rank_results[1].lineage])== set([superk_lin1, superk_lin2])
    assert set([rank_results[0].contained_at_rank, rank_results[1].contained_at_rank])== set([0.5])
    assert set([rank_results[0].contained_bp, rank_results[1].contained_bp])== set([3])

    assert set([rank_results[2].lineage, rank_results[3].lineage])== set([phylum_lin1, phylum_lin2])
    assert set([rank_results[2].contained_at_rank, rank_results[3].contained_at_rank])== set([0.5])
    assert set([rank_results[2].contained_bp, rank_results[3].contained_bp])== set([3])

    assert set([rank_results[4].lineage, rank_results[5].lineage])== set([lin1, lin2])
    assert set([rank_results[4].contained_at_rank, rank_results[5].contained_at_rank])== set([0.5])
    assert set([rank_results[4].contained_bp, rank_results[5].contained_bp])== set([3])

def test_contain_at_rank_4():
    # two minhashes, share ranks at phylum level
    hashval  = 12345678
    ident1 = 'first'
    mh1, sig1, lin1 = make_sig_and_lin([hashval], ident1, 'a;b;c')
    hashval2 = 87654321
    ident2 = 'second'
    mh2, sig2, lin2 = make_sig_and_lin([hashval2], ident2, 'a;b;f')
    # create lca_db w sigs
    lca_db = LCA_Database(scaled=1, ksize=3)
    lca_db.insert(sig1, ident=ident1)
    lca_db.insert(sig2, ident=ident2)
    # make lin_db
    lin_db = LineageDB()
    lin_db.insert(ident1, lin1)
    lin_db.insert(ident2, lin2)

    # search with combined hashvals
    search_mh = make_mh([hashval, hashval2])
    results, rank_results=search_containment_at_rank(search_mh, lca_db, lin_db, "class")

    assert len(results) == 2
    assert set([results[0].lineage, results[1].lineage])  == set([lin1, lin2])
    assert set([results[0].similarity, results[1].similarity])  == set([0.5])
    assert set([results[0].name, results[1].name]) == set([ident1, ident2])

    superk_lin = lca_utils.make_lineage('a')
    phylum_match_lin = lca_utils.make_lineage('a;b')

    assert len(rank_results) == 4
    # superk and phylum aggregate
    assert rank_results[0].lineage == superk_lin
    assert rank_results[0].contained_at_rank == 1.0
    assert rank_results[1].lineage == phylum_match_lin
    assert rank_results[1].contained_at_rank == 1.0
    assert rank_results[1].contained_bp == 6
    # different results at class
    assert set([rank_results[2].lineage, rank_results[3].lineage])== set([lin1, lin2])
    assert set([rank_results[2].contained_at_rank, rank_results[3].contained_at_rank])== set([0.5])
    assert set([rank_results[2].contained_bp, rank_results[3].contained_bp])== set([3])

