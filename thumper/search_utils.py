"""
utility functions for charcoal.
"""
import json
from collections import defaultdict, Counter, namedtuple
from operator import itemgetter
import csv

import sourmash
from sourmash.lca import lca_utils, LineagePair, taxlist
#from thumper.charcoal_utils import get_ident


## NOTE: see here for intersection mins
# https://github.com/dib-lab/sourmash/blob/7df142e4fe4de2bda4891936ea51cbc237168bcf/sourmash/search.py
#found_mins = set(_filter_max_hash(found_mins, new_max_hash))
#intersect_mins = query_mins.intersection(found_mins)
#intersect_orig_query_mins = orig_query_mins.intersection(found_mins)
#intersect_bp = cmp_scaled * len(intersect_orig_query_mins)


# generic SearchResult WITH lineage
SearchResult = namedtuple('SearchResult',
                          'similarity, match, md5, filename, name, lineage')

# generic RankSearchResult = summarized containment at rank
RankSumSearchResult = namedtuple('RankSumSearchResult',
                          'lineage, similarity')

ContigSearchInfo = namedtuple('ContigSearchInfo',
                              ['length', 'num_hashes', 'search_containment', 'contained_at_rank'])

def search_containment_at_rank(mh, lca_db, lin_db, match_rank, ignore_abundance=False, summarize_at_ranks=True):
    "Run search --containment, and aggregate at given rank and above."
    # do we need to copy mh if not modifying it?
    import copy
    minhash = copy.copy(mh)
    query_sig = sourmash.SourmashSignature(minhash)
    results=[]
    lin_hashes={}
    found_md5=set()
    search_iter = lca_db.search(query_sig, threshold=0, do_containment=True, ignore_abundance=ignore_abundance, best_only=False, unload_data=False)
    for (similarity, match_sig, filename) in search_iter:
        md5 = match_sig.md5sum()
        if md5 not in found_md5:
            found_md5.add(md5)

            # get lineage
            match_ident = get_ident(match_sig)
            match_lineage = lin_db.ident_to_lineage[match_ident]
            match_lineage = pop_to_rank(match_lineage, match_rank)

            results.append((similarity, match_sig, filename, match_lineage))

            if summarize_at_ranks: #and len(search_iter) >1: # no need to summarize if we just have one hit
                # add the match_sig hashes so we can calculate containment of contig by this genome
                for rank in lca_utils.taxlist(include_strain=False):
                    lin_at_rank = pop_to_rank(match_lineage, rank)
                    if lin_at_rank not in lin_hashes.keys():
                        # would be neat to have a query_sig.return_common(match_sig) to get *just* the hashes in common
                        fresh_mh = mh.copy_and_clear()
                        fresh_mh.add_many(match_sig.minhash.hashes)
                        lin_hashes[lin_at_rank] = (fresh_mh, rank)
                        #lin_hashes[lin_at_rank] = (match_sig.minhash.hashes, similarity)
                    else:
                        current_mh = lin_hashes[lin_at_rank][0]
                        current_mh.add_many(match_sig.minhash.hashes)
                        lin_hashes[lin_at_rank] = (current_mh, rank)
                        # could also calculate + store containment here. But would recalculated every time we add more hashes.. too many containment operations?
                        #similarity = query_sig.contained_by(current_hashes)
                        #lin_hashes[lin_at_rank] = (match_sig.minhash.hashes, similarity)
                    if rank == match_rank:
                        break

    # sort normal search --containment results on similarity (reverse)
    results.sort(key=itemgetter(0), reverse=True)

    x = []
    for (similarity, match, filename, match_lineage) in results:

        x.append(SearchResult(similarity=similarity,
                              match=[],# match
                              md5=match.md5sum(),
                              filename=filename,
                              name=match.name(),
                              lineage=match_lineage))


    # now, calculate containment for each lineage match at each rank
    summarized_results = defaultdict(list)
    y = []
    if summarize_at_ranks:
        for lin, (match_hashes, rank) in lin_hashes.items():
            linmatch_sig = sourmash.SourmashSignature(match_hashes)
            containment = query_sig.contained_by(linmatch_sig)
            summarized_results[rank].append((rank, lin, containment)) #, linmatch_sig))

        # sort and store results

        # iterate superkingdom --> match_rank
        for rank in lca_utils.taxlist(include_strain=False):
            rank_res = summarized_results[rank]
            # sort by containment
            #rank_res.sort(key=lambda x: -x[2])
            rank_res.sort(key=itemgetter(2), reverse=True)
            for (rank, lin, containment) in rank_res:
                y.append(RankSumSearchResult(lineage=lin, similarity=containment)) #, match=linmatch_sig))
            if rank == match_rank:
                break
    return x,y


## make test fixtures? lca db, lindb
#class FakeLCA_Database(object):
#    def __init__(self):
#        self._assignments = {}
#
#    def _set_lineage_assignment(self, hashval, assignment):
#        self._assignments[hashval] = assignment
#
#    def get_lineage_assignments(self, hashval):
#        if hashval in self._assignments:
            #return self._assignments[hashval]
        #else:
        #    return None

#@pytest.fixture
#def lin_db():
# build lindb
#    return lindb
    # make a mh

    #mh = sourmash.MinHash(n=1, ksize=20, track_abundance=True)
    #mh.add_kmer("AT" * 10)
    #sig1 = SourmashSignature(mh, name='foo')

    #lineage = ((LineagePair('rank1', 'name1'),
    #            LineagePair('rank2', 'name2')))

#from sourmash.lca.command_index import load_taxonomy_assignments

# from charcoal_utils:
#def get_ident(sig):
#    "Hack and slash identifiers."
#    ident = sig.name()
#    ident = ident.split()[0]
#    ident = ident.split('.')[0]
#    return ident

#def get_idents_for_hashval(lca_db, hashval):
#    "Get the identifiers associated with this hashval."
#    idx_list = lca_db.hashval_to_idx.get(hashval, [])
#    for idx in idx_list:
#        ident = lca_db.idx_to_ident[idx]
#        yield ident


def test_contain_at_rank_1():
    from sourmash.lca import LCA_Database
    from lineage_db import LineageDB
    # like, one minhash, one set of ranks

    # create mh, sig w/hashval
    hashval  = 12345678
    mh = sourmash.MinHash(n=0, scaled=1, ksize=3)
    mh.add_hash(hashval) # make another test using add_hash_with_abundance?
    sig1 = sourmash.SourmashSignature(mh)

    # create lca_db w sig1
    ident = 'uniq'
    lca_db = LCA_Database(scaled=1, ksize=3)
    lca_db.insert(sig1, ident=ident)

    # next, make lin_db
    lin = lca_utils.make_lineage('a;b;c')
    lin_db = LineageDB()
    lin_db.insert(ident, lin)
    # lineage db created properly?
    assert 'uniq' in lin_db.lineage_to_idents[lin]

    results, rank_results =search_containment_at_rank(mh, lca_db, lin_db, "order")
    print(results)

    print(rank_results)


    #assert 'uniq' in ldb.lineage_to_idents[lineage]
    #assert ldb.ident_to_lineage['uniq'] == lineage
    #assert assignments[hashval] == set([ lin ])


#def test_contain_at_rank_2():
     #two minhashes, fully shared ranks


#def test_contain_at_rank_3():
    # two minhashes, totally distinct ranks

#def test_contain_at_rank_4():
    # two minhashes, totally distinct ranks
