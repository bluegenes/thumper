"""
Perform search --containment and optionally aggregate at rank
"""
import json
from collections import defaultdict, Counter, namedtuple
from operator import itemgetter
import csv

import sourmash
from sourmash.lca import lca_utils, LineagePair, taxlist
from thumper.charcoal_utils import get_ident, pop_to_rank


# generic SearchResult WITH lineage
SearchResultLin = namedtuple('SearchResult',
                          'similarity, match, md5, filename, name, lineage')

# generic RankSearchResult = summarized containment at rank
RankSumSearchResult = namedtuple('RankSumSearchResult',
                                 #'lineage, containment, intersect_bp, match_sig')
                                 'lineage, contained_at_rank, contained_bp, match_sig')

# use csv instead of json for now
#ContigSearchInfo = namedtuple('ContigSearchInfo',
#                              ['length', 'num_hashes', 'contained_at_rank', 'bp_contained'])

# for rank results: contig_name,length,num_hashes,match_rank,lineage,contained_at_rank,bp_contained


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
    scaled_val = int(query_sig.minhash.scaled)
    ksize = int(query_sig.minhash.ksize)
    for lin, matched_hashes in lineage_hashD.items():
        rank = lin[-1].rank
        # TODO; check this. just scaled_val, or scaled * ksize * num matched hashes?
        intersect_bp = scaled_val * len(matched_hashes) * ksize
        linmatch_sig = sourmash.SourmashSignature(matched_hashes)
        containment = query_sig.contained_by(linmatch_sig)
        summarized_results[rank].append((lin, containment, intersect_bp, linmatch_sig)) # optionally don't keep track of sig here
    return summarized_results


def sort_by_rank_and_containment(summarized_results, match_rank):
    sorted_results = []
    # iterate superkingdom --> match_rank
    for rank in lca_utils.taxlist(include_strain=False):
        rank_res = summarized_results[rank]
        rank_res.sort(key=itemgetter(1), reverse=True)  # sort by containment
        for (lin, containment, intersect_bp, match_sig) in rank_res:
            sorted_results.append(RankSumSearchResult(lineage=lin, contained_at_rank=containment, contained_bp=intersect_bp, match_sig=match_sig))
        if rank == match_rank:
            break
    return sorted_results


def sort_and_store_search_results(res): # same as in charcoal_utils
    # sort normal search --containment results on similarity (reverse)
    sorted_rs = []
    res.sort(key=itemgetter(0), reverse=True)
    for (similarity, match, filename, match_lineage) in res:
       sorted_rs.append(SearchResultLin(similarity=similarity,
                                        match=match,
                                        md5=match.md5sum(),
                                        filename=filename,
                                        name=match.name(),
                                        lineage=match_lineage))
    return sorted_rs


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
