"""
test thumper/search_utils.py
"""
import os
import json
from collections import defaultdict, Counter, namedtuple
from operator import itemgetter
import csv

import sourmash
from sourmash.lca import lca_utils, LineagePair, LCA_Database, taxlist
from thumper.lineage_db import LineageDB
from thumper.charcoal_utils import get_ident, pop_to_rank, gather_at_rank
# to do: specific imports!
from thumper.search_utils import *
from . import pytest_utils as utils


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
    assert containmentD["phylum"][0][1] == 1.0


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
    # each class should have 25% containment
    assert containmentD["class"][0][1] == 0.25
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
    assert rank_results[2].lineage == lin1
    assert rank_results[2].contained_at_rank == 1.0

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
    assert rank_results[2].lineage == lin1
    assert rank_results[2].contained_at_rank == 1.0

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

    assert set([rank_results[2].lineage, rank_results[3].lineage])== set([phylum_lin1, phylum_lin2])
    assert set([rank_results[2].contained_at_rank, rank_results[3].contained_at_rank])== set([0.5])

    assert set([rank_results[4].lineage, rank_results[5].lineage])== set([lin1, lin2])
    assert set([rank_results[4].contained_at_rank, rank_results[5].contained_at_rank])== set([0.5])

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
    # different results at class
    assert set([rank_results[2].lineage, rank_results[3].lineage])== set([lin1, lin2])
    assert set([rank_results[2].contained_at_rank, rank_results[3].contained_at_rank])== set([0.5])

def test_gather_at_rank_1():
    # one minhash, one set of ranks
    hashval  = 12345678
    ident = 'uniq'
    mh1, sig1, lin1 = make_sig_and_lin([hashval], ident, 'a;b;c')

    lca_db = LCA_Database(scaled=1, ksize=3)
    lca_db.insert(sig1, ident=ident)

    lin_db = LineageDB()
    lin_db.insert(ident, lin1)

    gather_results=list(gather_at_rank(mh1, lca_db, lin_db, "class"))
    assert len(gather_results) == 1
    assert gather_results[0][0] == lin1
    assert gather_results[0][1] == 1


def test_gather_at_rank_2():
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
    gather_results=list(gather_at_rank(search_mh, lca_db, lin_db, "class"))
    assert len(gather_results) == 1
    assert gather_results[0][0] == lin1
    assert gather_results[0][1] == 2


def test_gather_at_rank_3():
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
    gather_results=list(gather_at_rank(search_mh, lca_db, lin_db, "class"))

    assert len(gather_results) == 2
    assert set([gather_results[0][0],gather_results[1][0]]) == set([lin1, lin2])
    assert set([gather_results[0][1],gather_results[1][1]]) == set([1])


def test_gather_guess_tax_at_rank_1():
    # one minhash, one set of ranks
    hashval  = 12345678
    ident = 'uniq'
    mh1, sig1, lin1 = make_sig_and_lin([hashval], ident, 'a;b;c')

    lca_db = LCA_Database(scaled=1, ksize=3)
    lca_db.insert(sig1, ident=ident)

    lin_db = LineageDB()
    lin_db.insert(ident, lin1)

    num_hashes = 1
    phylum_match_lin = lca_utils.make_lineage('a;b')

    gather_results=list(gather_at_rank(mh1, lca_db, lin_db, "class"))
    phylum_results = gather_guess_tax_at_rank(gather_results, num_hashes, "phylum", minimum_matches=1)

    assert len(phylum_results) == 3

    assert phylum_results[0] == phylum_match_lin
    assert phylum_results[1] == 1.0

def test_gather_guess_tax_at_rank_2():
    # one minhash, one set of ranks
    hashval  = 12345678
    ident = 'uniq'
    mh1, sig1, lin1 = make_sig_and_lin([hashval], ident, 'a;b;c')

    lca_db = LCA_Database(scaled=1, ksize=3)
    lca_db.insert(sig1, ident=ident)

    lin_db = LineageDB()
    lin_db.insert(ident, lin1)

    num_hashes = 1
    phylum_match_lin = lca_utils.make_lineage('a;b')

    gather_results=list(gather_at_rank(mh1, lca_db, lin_db, "class"))
    phylum_results = gather_guess_tax_at_rank(gather_results, num_hashes, "phylum", minimum_matches=3)

    assert len(phylum_results) == 3
    assert set(phylum_results) == set([""])


def test_gather_guess_tax_at_each_rank_1():
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

    num_hashes = 2
    superk_lin = lca_utils.make_lineage('a')
    phylum_lin = lca_utils.make_lineage('a;b')

    # search with combined hashvals
    search_mh = make_mh([hashval, hashval2])
    gather_results=list(gather_at_rank(search_mh, lca_db, lin_db, "class"))
    rank_results=gather_guess_tax_at_each_rank(gather_results, num_hashes, minimum_matches=1, \
                                               lowest_rank="class",
                                               taxlist=lca_utils.taxlist(include_strain=False))

    assert len(rank_results) == 3

    assert rank_results[0] == RankSumGatherResult(lineage=superk_lin,f_ident=1.0, f_major=1.0)
    assert rank_results[1] == RankSumGatherResult(lineage=phylum_lin,f_ident=1.0, f_major=1.0)
    assert rank_results[2] == RankSumGatherResult(lineage=lin1,f_ident=1.0, f_major=1.0)


def test_gather_guess_tax_at_each_rank_2():
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

    # some useful bits
    num_hashes = 2
    superk_lin = lca_utils.make_lineage('a')
    phylum_lin = lca_utils.make_lineage('a;b')

    # search with combined hashvals
    search_mh = make_mh([hashval, hashval2])
    gather_results=list(gather_at_rank(search_mh, lca_db, lin_db, "class"))
    rank_results=gather_guess_tax_at_each_rank(gather_results, num_hashes, minimum_matches=1, \
                                               lowest_rank="class", \
                                               taxlist=lca_utils.taxlist(include_strain=False))

    assert len(rank_results) == 3

    assert rank_results[0] == RankSumGatherResult(lineage=superk_lin,f_ident=1.0, f_major=1.0)
    assert rank_results[1] == RankSumGatherResult(lineage=phylum_lin,f_ident=1.0, f_major=1.0)

    assert rank_results[2].lineage in [lin1,lin2]
    assert rank_results[2].f_ident == 1.0
    assert rank_results[2].f_major == 0.5


def test_gather_guess_tax_at_each_rank_3():
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

    num_hashes = 2
    #winner seems to be def lineage.. will this remain true always?
    superk_lin = lca_utils.make_lineage('d')
    phylum_lin = lca_utils.make_lineage('d;e')

    # search with combined hashvals
    search_mh = make_mh([hashval1, hashval2])
    gather_results=list(gather_at_rank(search_mh, lca_db, lin_db, "class"))
    rank_results=gather_guess_tax_at_each_rank(gather_results, num_hashes, minimum_matches=1, \
                                               lowest_rank="class",
                                               taxlist=lca_utils.taxlist(include_strain=False))
    assert len(rank_results) == 3

    assert rank_results[0] == RankSumGatherResult(lineage=superk_lin,f_ident=1.0, f_major=0.5)
    assert rank_results[1] == RankSumGatherResult(lineage=phylum_lin,f_ident=1.0, f_major=0.5)
    assert rank_results[2] == RankSumGatherResult(lineage=lin2,f_ident=1.0, f_major=0.5)

def test_get_match_bp_1():
    num_matched_hashes=2
    scaled=1
    ksize=3
    match_bp = get_match_bp(scaled, ksize, num_matched_hashes=num_matched_hashes)
    assert match_bp == 2.0

def test_get_match_bp_2():
    total_num_hashes=4
    f_major = 0.5
    scaled=1
    ksize=3
    match_bp = get_match_bp(scaled, ksize, match_percent=f_major, total_num_hashes=total_num_hashes)
    assert match_bp == 2.0

def test_get_match_bp_3():
    total_num_hashes=4
    f_major = 0.5
    scaled=1
    ksize=3
    match_bp = get_match_bp(scaled, ksize, match_percent=f_major)
    assert match_bp == "NA"

@utils.in_tempdir
def test_searchfiles_contigs_1(location):
    prefix = os.path.join(location, "pref")
    filelist = [f"{prefix}.contigs.ranksearch.csv",
                f"{prefix}.contigs.ranksearch.matches.sig",
                f"{prefix}.contigs.search.csv",
                f"{prefix}.contigs.search.matches.sig",
                f"{prefix}.contigs.unmatched.fq"]

    sf = SearchFiles(prefix)
    sf.close()
    for f in filelist:
        assert os.path.exists(f)


@utils.in_tempdir
def test_searchfiles_contigs_2(location):
    prefix = os.path.join(location, "pref")
    filelist = [f"{prefix}.contigs.rankgather.csv",
                f"{prefix}.contigs.unmatched.fq"]

    sf = SearchFiles(prefix, search=False, gather=True)
    sf.close()
    for f in filelist:
        assert os.path.exists(f)


@utils.in_tempdir
def test_searchfiles_contigs_3(location):
    prefix = os.path.join(location, "prefix")
    filelist = [f"{prefix}.contigs.ranksearch.csv",
                f"{prefix}.contigs.ranksearch.matches.sig",
                f"{prefix}.contigs.search.csv",
                f"{prefix}.contigs.search.matches.sig",
                f"{prefix}.contigs.rankgather.csv",
                f"{prefix}.contigs.unmatched.fq"]

    sf = SearchFiles(prefix, search=True, gather=True)
    sf.close()
    for f in filelist:
        assert os.path.exists(f)


def get_csv_set(f):
    return set(map(tuple, csv.reader(f)))


@utils.in_tempdir
def test_searchfiles_contigs_just_search(location):
    prefix = os.path.join(location, "pref")
    filelist = [f"{prefix}.contigs.ranksearch.csv",
                f"{prefix}.contigs.ranksearch.matches.sig",
                f"{prefix}.contigs.search.csv",
                f"{prefix}.contigs.search.matches.sig",
                f"{prefix}.contigs.unmatched.fq"]

    sf = SearchFiles(prefix, search=True, gather=True)

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
    num_hashes = 2
    # search with combined hashvals
    search_mh = make_mh([hashval, hashval2])
    results, rank_results=search_containment_at_rank(search_mh, lca_db, lin_db, "class")
    gather_results=list(gather_at_rank(search_mh, lca_db, lin_db, "class"))


    #write search results
    name = 'name'
    seq_len = 6
    for res in results:
        sf.write_result(res, name, seq_len, result_type="search")
    for res in rank_results:
        sf.write_result(res, name, seq_len, result_type="ranksearch")

    sf.close()

    # check results are in files
    for f in filelist:
        assert os.path.exists(f)

    with open(f"{prefix}.contigs.search.csv", "r") as searchres:
         this_search_csvset = get_csv_set(searchres)
    with open(utils.test_file("test-data/test.contigs.search.csv"), "r") as searchres:
         saved_search_csvset = get_csv_set(searchres)
    assert saved_search_csvset == this_search_csvset

    with open(f"{prefix}.contigs.ranksearch.csv", "r") as searchres:
         this_ranksearch_csvset = get_csv_set(searchres)
    with open(utils.test_file("test-data/test.contigs.ranksearch.csv"), "r") as searchres:
         saved_ranksearch_csvset = get_csv_set(searchres)
    assert saved_ranksearch_csvset == this_ranksearch_csvset

@utils.in_tempdir
def test_searchfiles_contigs_just_gather(location):
    prefix = os.path.join(location, "pref")
    filelist = [f"{prefix}.contigs.rankgather.csv",
                f"{prefix}.contigs.unmatched.fq"]

    sf = SearchFiles(prefix, search=True, gather=True)

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
    num_hashes = 2
    # search with combined hashvals
    search_mh = make_mh([hashval, hashval2])
    gather_results=list(gather_at_rank(search_mh, lca_db, lin_db, "class"))

    gather_rank_results=gather_guess_tax_at_each_rank(gather_results, num_hashes, minimum_matches=1, \
                                                      lowest_rank="class", \
                                                      taxlist=lca_utils.taxlist(include_strain=False))

    #write search results
    name = 'name'
    seq_len = 6
    for gres in gather_rank_results:
        sf.write_result(gres, name, seq_len, result_type="rankgather")

    sf.close()

    # check results are in files
    for f in filelist:
        assert os.path.exists(f)

    with open(f"{prefix}.contigs.rankgather.csv", "r") as gatherres:
         this_gather_csvset = get_csv_set(gatherres)
    with open(utils.test_file("test-data/test.contigs.rankgather.csv"), "r") as searchres:
         saved_gather_csvset = get_csv_set(searchres)
    assert saved_gather_csvset == this_gather_csvset


@utils.in_tempdir
def test_searchfiles_contigs_search_and_gather(location):
    prefix = os.path.join(location, "pref")
    filelist = [f"{prefix}.contigs.ranksearch.csv",
                f"{prefix}.contigs.ranksearch.matches.sig",
                f"{prefix}.contigs.search.csv",
                f"{prefix}.contigs.search.matches.sig",
                f"{prefix}.contigs.rankgather.csv",
                f"{prefix}.contigs.unmatched.fq"]

    sf = SearchFiles(prefix, search=True, gather=True)

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
    num_hashes = 2
    # search with combined hashvals
    search_mh = make_mh([hashval, hashval2])
    results, rank_results=search_containment_at_rank(search_mh, lca_db, lin_db, "class")
    gather_results=list(gather_at_rank(search_mh, lca_db, lin_db, "class"))

    gather_rank_results=gather_guess_tax_at_each_rank(gather_results, num_hashes, minimum_matches=1, \
                                                      lowest_rank="class", \
                                                      taxlist=lca_utils.taxlist(include_strain=False))

    #write search results
    name = 'name'
    seq_len = 6
    for res in results:
        sf.write_result(res, name, seq_len, result_type="search")
    for res in rank_results:
        sf.write_result(res, name, seq_len, result_type="ranksearch")
    for gres in gather_rank_results:
        sf.write_result(gres, name, seq_len, result_type="rankgather")

    sf.close()

    # check results are in files
    for f in filelist:
        assert os.path.exists(f)

    with open(f"{prefix}.contigs.search.csv", "r") as searchres:
         this_search_csvset = get_csv_set(searchres)
    with open(utils.test_file("test-data/test.contigs.search.csv"), "r") as searchres:
         saved_search_csvset = get_csv_set(searchres)
    assert saved_search_csvset == this_search_csvset

    with open(f"{prefix}.contigs.ranksearch.csv", "r") as searchres:
         this_ranksearch_csvset = get_csv_set(searchres)
    with open(utils.test_file("test-data/test.contigs.ranksearch.csv"), "r") as searchres:
         saved_ranksearch_csvset = get_csv_set(searchres)
    assert saved_ranksearch_csvset == this_ranksearch_csvset

    with open(f"{prefix}.contigs.rankgather.csv", "r") as gatherres:
         this_gather_csvset = get_csv_set(gatherres)
    with open(utils.test_file("test-data/test.contigs.rankgather.csv"), "r") as searchres:
         saved_gather_csvset = get_csv_set(searchres)
    assert saved_gather_csvset == this_gather_csvset


@utils.in_tempdir
def test_searchfiles_contigs_and_genome_1(location):
    prefix = os.path.join(location, "cg1")
    filelist = [f"{prefix}.ranksearch.csv",
                f"{prefix}.ranksearch.matches.sig",
                f"{prefix}.search.csv",
                f"{prefix}.search.matches.sig",
                f"{prefix}.contigs.ranksearch.csv",
                f"{prefix}.contigs.ranksearch.matches.sig",
                f"{prefix}.contigs.search.csv",
                f"{prefix}.contigs.search.matches.sig",
                f"{prefix}.contigs.unmatched.fq"]

    sf = SearchFiles(prefix)
    sf.close()
    gf = SearchFiles(prefix, contigs=False)
    gf.close()

    for f in filelist:
        assert os.path.exists(f)


@utils.in_tempdir
def test_searchfiles_contigs_and_genome_2(location):
    prefix = os.path.join(location, "cg2")
    filelist = [f"{prefix}.rankgather.csv",
                f"{prefix}.contigs.rankgather.csv",
                f"{prefix}.contigs.unmatched.fq"]

    sf = SearchFiles(prefix, search=False, gather=True)
    sf.close()
    gf = SearchFiles(prefix, search=False, gather=True, contigs=False)
    gf.close()
    for f in filelist:
        assert os.path.exists(f)


@utils.in_tempdir
def test_searchfiles_contigs_and_genome_3(location):
    prefix = os.path.join(location, "cg3")
    filelist = [f"{prefix}.ranksearch.csv",
                f"{prefix}.ranksearch.matches.sig",
                f"{prefix}.search.csv",
                f"{prefix}.search.matches.sig",
                f"{prefix}.contigs.ranksearch.csv",
                f"{prefix}.contigs.ranksearch.matches.sig",
                f"{prefix}.contigs.search.csv",
                f"{prefix}.contigs.search.matches.sig",
                f"{prefix}.contigs.rankgather.csv",
                f"{prefix}.contigs.unmatched.fq"]

    sf = SearchFiles(prefix, search=True, gather=True)
    sf.close()
    gf = SearchFiles(prefix, search=True, gather=True, contigs=False)
    gf.close()
    for f in filelist:
        assert os.path.exists(f)
