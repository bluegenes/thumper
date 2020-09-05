"""
utility functions for charcoal.
"""
import json
from collections import defaultdict, Counter, namedtuple
from operator import itemgetter
import csv

import sourmash
from sourmash.lca import lca_utils, LineagePair, taxlist
from thumper.charcoal_utils import *


# generic SearchResult WITH lineage
SearchResult = namedtuple('SearchResult',
                          'similarity, match, md5, filename, name, lineage')

# generic RankSearchResult = summarized containment at rank
RankSumSearchResult = namedtuple('RankSumSearchResult',
                                 'lineage, containment, match_sig')

ContigSearchInfo = namedtuple('ContigSearchInfo',
                              ['length', 'num_hashes', 'search_containment', 'contained_at_rank'])


def add_hashes_at_ranks(lineage_hashD, hashes_to_add, lineage, match_rank):
    # first add full lineage
    lineage_hashD[lineage].add_many(hashes_to_add)
    for rank in lca_utils.taxlist(include_strain=False):
        # TODO: add check to pop ONLY if needed (no need to pop at genus if lineage only has superk, phyl)
        lin_at_rank = pop_to_rank(lineage, rank)
        lineage_hashD[lin_at_rank].add_many(hashes_to_add)
        if rank == match_rank:
            break
    return lineage_hashD


def get_lineage_at_match_rank(linDB, sig, match_rank):
    match_ident = get_ident(sig)
    match_lineage = linDB.ident_to_lineage[match_ident]
    match_lineage = pop_to_rank(match_lineage, match_rank)
    return match_lineage


def calculate_containment_at_rank(lineage_hashD, query_sig, match_rank):
    # calculate containment for each lineage match at each rank
    summarized_results = defaultdict(list)
    for lin, matched_hashes in lineage_hashD.items():
        rank = lin[-1].rank
        print(rank)
        # maybe also calculate matched_bp?
        #intersect_bp = scaled * len(intersected_hashes)
        linmatch_sig = sourmash.SourmashSignature(matched_hashes)
        containment = query_sig.contained_by(linmatch_sig)
        summarized_results[rank].append((lin, containment, linmatch_sig)) # optionally don't keep track of sig here
    return summarized_results


def sort_by_rank_and_containment(summarized_results, match_rank):
    sorted_results = []
    # iterate superkingdom --> match_rank
    for rank in lca_utils.taxlist(include_strain=False):
        rank_res = summarized_results[rank]
        rank_res.sort(key=itemgetter(1), reverse=True)  # sort by containment
        #print(rank_res)
        for (lin, containment, match_sig) in rank_res:
            sorted_results.append(RankSumSearchResult(lineage=lin, containment=containment, match_sig=match_sig))
        if rank == match_rank:
            break
    return sorted_results

# test_sort_by_rank_and_containment
# 1. three results, check that they sort by rank, containment


def sort_and_store_search_results(res):
    # sort normal search --containment results on similarity (reverse)
    sorted_rs = []
    res.sort(key=itemgetter(0), reverse=True)
    for (similarity, match, filename, match_lineage) in res:
       sorted_rs.append(SearchResult(similarity=similarity,
                                     match=match,
                                     md5=match.md5sum(),
                                     filename=filename,
                                     name=match.name(),
                                     lineage=match_lineage))
    return sorted_rs

# test_sort_and_store_search_results
# 1. three results, check that they sort by containment

def search_containment_at_rank(mh, lca_db, lin_db, match_rank, ignore_abundance=False, summarize_at_ranks=True):
    "Run search --containment, and aggregate at given rank and above."

    results=[]
    found_md5=set()
    def gen_mh():
        return mh.copy_and_clear()
    lin_hashes=defaultdict(gen_mh) #defaultdict requires function that defines an empty minhash
    query_hashes = set(mh.hashes)
    query_sig = sourmash.SourmashSignature(mh)

    # search
    search_iter = lca_db.search(query_sig, threshold=0, do_containment=True, \
                  ignore_abundance=ignore_abundance, best_only=False, unload_data=False)

    # iterate through matches
    for (similarity, match_sig, filename) in search_iter:
        md5 = match_sig.md5sum()
        if md5 not in found_md5:
            found_md5.add(md5)
            match_lineage = get_lineage_at_match_rank(lin_db, match_sig, match_rank)
            results.append((similarity, match_sig, filename, match_lineage)) #store search results + lineage

            if summarize_at_ranks:
                # Keep track of matched hashes at higher taxonomic ranks
                intersected_hashes = query_hashes.intersection(set(match_sig.minhash.hashes))
                lin_hashes = add_hashes_at_ranks(lin_hashes, intersected_hashes, match_lineage, match_rank)

    # sort and store results
    search_results = sort_and_store_search_results(results)
    search_results_at_rank = []
    if summarize_at_ranks:
        rank_containment = calculate_containment_at_rank(lin_hashes, query_sig, match_rank)
        search_results_at_rank = sort_by_rank_and_containment(rank_containment, match_rank)

    return search_results, search_results_at_rank





## tests (to do - move to separate testing file)
## TO DO:
  # - move to separate testing file


from sourmash.lca import LCA_Database
from thumper.lineage_db import LineageDB

def make_mh(hashvals, ksize=3, scaled=1):
    mh = sourmash.MinHash(n=0, scaled=1, ksize=3)
    mh.add_many(hashvals)
    return mh

def make_sig_and_lin(hashvals, ident, lin, ksize=3, scaled=1):
    mh = make_mh(hashvals)
    sig = sourmash.SourmashSignature(mh, name=ident)
    lineage = lca_utils.make_lineage('a;b;c')
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


def test_contain_at_rank_1():
    # one minhash, one set of ranks

    # create mh, sig w/hashval
    hashval  = 12345678
    ident = 'uniq'
    mh1, sig1, lin1 = make_sig_and_lin([hashval], ident, 'a;b;c')

    # create lca_db w sig1
    lca_db = LCA_Database(scaled=1, ksize=3)
    lca_db.insert(sig1, ident=ident)

    # make lin_db
    lin_db = LineageDB()
    lin_db.insert(ident, lin1)
    # lineage db created properly?
    # assert 'uniq' in lin_db.lineage_to_idents[lin1]

    results, rank_results=search_containment_at_rank(mh1, lca_db, lin_db, "class")
    #print(results)
    #print(rank_results)

    #assert 'uniq' in ldb.lineage_to_idents[lineage]
    #assert ldb.ident_to_lineage['uniq'] == lineage
    #assert assignments[hashval] == set([ lin ])

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

    #print(results)
    #print(rank_results)

def test_contain_at_rank_3():
    # two minhashes, totally distinct ranks
    # first sig
    hashval  = 12345678
    ident1 = 'first'
    mh1, sig1, lin1 = make_sig_and_lin([hashval], ident1, 'a;b;c')

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
    search_mh = make_mh([hashval, hashval2])
    results, rank_results=search_containment_at_rank(search_mh, lca_db, lin_db, "class")

    # check results
    #print(results)
    #print(rank_results)

def test_contain_at_rank_4():
    # two minhashes, share ranks at phylum level

    # first sig
    hashval  = 12345678
    ident1 = 'first'
    mh1, sig1, lin1 = make_sig_and_lin([hashval], ident1, 'a;b;c')

    # second sig
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

    #print(results)
    #print(rank_results)

