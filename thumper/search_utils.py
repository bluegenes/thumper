"""
Perform search --containment and optionally aggregate at rank
"""
import json
from collections import defaultdict, Counter, namedtuple
from operator import itemgetter
import csv

import sourmash
from sourmash.lca import lca_utils, LineagePair, taxlist
from thumper.charcoal_utils import get_ident, pop_to_rank, gather_at_rank, summarize_at_rank

#GATHER_MIN_MATCHES=3

# generic SearchResult WITH lineage
SearchResultLin = namedtuple('SearchResult',
                          'similarity, match, md5, filename, name, lineage')

# generic RankSearchResult = summarized containment at rank
RankSumSearchResult = namedtuple('RankSumSearchResult',
                                 #'lineage, containment, intersect_bp, match_sig')
                                 'lineage, contained_at_rank, contained_bp, match')

RankSumGatherResult = namedtuple('RankSumGatherResult', 'lineage, f_ident, f_major')

#RankSumGatherResult = namedtuple('RankSumGatherResult',
#                                 'lineage, num_hashes, matched_bp')

# use csv instead of json for now
#ContigSearchInfo = namedtuple('ContigSearchInfo',
#                              ['length', 'num_hashes', 'contained_at_rank', 'bp_contained'])

# for rank results: name,length,num_hashes,match_rank,lineage,contained_at_rank,bp_contained


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
        #intersect_bp = scaled_val * len(matched_hashes) * ksize
        intersect_bp = get_match_bp(scaled_val, ksize, num_matched_hashes=len(matched_hashes))
        linmatch_sig = sourmash.SourmashSignature(matched_hashes) #ADD MORE INFO (e.g. name/ident?) HERE IF KEEPING SIG?
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
            sorted_results.append(RankSumSearchResult(lineage=lin, contained_at_rank=containment, contained_bp=intersect_bp, match=match_sig))
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


def get_match_bp(scaled, ksize, num_matched_hashes=None, match_percent=None, total_num_hashes=None):
    # TO DO: Check this.
    if match_percent and total_num_hashes:
        return (float(match_percent)*int(total_num_hashes) * int(scaled))
    elif num_matched_hashes:
        return (float(num_matched_hashes) * int(scaled))
    else:
        # to be safe, should probably return something useful if we don't have these...
        print("Can't calculate matched bp. Please make sure you've provided all the right info.")
        return "NA"


def gather_guess_tax_at_rank(gather_results, num_hashes, rank, minimum_matches=3):
    # modified from charcoal_utils
    "Guess likely taxonomy using gather."
    sum_ident = 0
    first_lin = ()
    first_count = 0
    # summarize to rank
    rank_gather = summarize_at_rank(gather_results, rank)
    for lin, count in rank_gather:
        if count >= minimum_matches:
            # record the first lineage we come across as likely lineage.
            if not first_lin:
                first_lin = lin
                first_count = count

        sum_ident += count

    if not first_lin:
        return "","",""

    f_ident = sum_ident / num_hashes
    f_major = first_count / sum_ident

    return first_lin, f_ident, f_major


def gather_guess_tax_at_each_rank(gather_results, num_hashes, taxlist=lca_utils.taxlist(include_strain=False), minimum_matches=3, lowest_rank="genus"):
    rank_results = []
    prev_lineage=""
    top_lineage=""
    for rank in taxlist:
        top_lineage, f_ident, f_major = gather_guess_tax_at_rank(gather_results, num_hashes, rank, minimum_matches=minimum_matches)

        # summarizing at a lower rank than exists will yield same result as prev. break!
        if not top_lineage or top_lineage == prev_lineage:
            break
        rank_results.append(RankSumGatherResult(lineage=top_lineage, f_ident=f_ident, f_major=f_major))
        prev_lineage = top_lineage
        if rank == lowest_rank:
            break

    return rank_results


class SearchFiles:
    """
    Class to handle all the files created during search or gather
    """
    def __init__(self, out_prefix, search=True, gather=False, contigs=True):

        self.gather = gather
        self.search = search
        self.contigs = contigs

        if self.contigs:
            out_prefix = out_prefix + ".contigs"

            self.unmatched = open(f"{out_prefix}.unmatched.fq", "w")

        if self.search:
            self.search_csv = open(f"{out_prefix}.search.csv", "w")
            self.search_matches = open(f"{out_prefix}.search.matches.sig", "w")
            self.ranksearch_csv = open(f"{out_prefix}.ranksearch.csv", "w")
            self.ranksearch_matches = open(f"{out_prefix}.ranksearch.matches.sig", "w")
            self.search_sigs = []
            self.ranksearch_sigs = []

            search_fieldnames = ['name', 'length', 'similarity', 'name', 'filename', 'md5', 'lineage']
            #search_fieldnames = ['similarity', 'name', 'filename', 'md5', 'lineage']
            self.search_w = csv.DictWriter(self.search_csv, fieldnames=search_fieldnames)
            self.search_w.writeheader()

            rank_fieldnames = ['name', 'length', 'match_rank', 'lineage', 'contained_at_rank', 'contained_bp']
            self.rank_w = csv.DictWriter(self.ranksearch_csv, fieldnames=rank_fieldnames)
            self.rank_w.writeheader()

        if self.gather:
            self.rankgather_csv = open(f"{out_prefix}.rankgather.csv", "w")
            gather_rank_fieldnames = ['name', 'length', 'match_rank', 'lineage', 'f_ident', 'f_major']
            self.gather_rank_w = csv.DictWriter(self.rankgather_csv, fieldnames=gather_rank_fieldnames)
            self.gather_rank_w.writeheader()

    def write_result(self, result, name, length, result_type="search"):
        # write single result
        d = dict(result._asdict())
        d["name"] = name
        d["length"] = length
        d["lineage"] = lca_utils.display_lineage(result.lineage)

        if self.search and result_type == "search":
            self.search_sigs.append(d['match'])
            del d['match']
            self.search_w.writerow(d)
        elif self.search and result_type == "ranksearch":
            self.ranksearch_sigs.append(d['match'])
            del d['match']
            d["match_rank"] = result.lineage[-1].rank
            self.rank_w.writerow(d)
        elif self.gather and result_type == "rankgather":
            d["match_rank"] = result.lineage[-1].rank
            #d["major_bp"] = get_match_bp(float(gr.f_major))
            self.gather_rank_w.writerow(d)


    def close(self):
        # close files
        if self.contigs:
            self.unmatched.close()
        if self.search:
            self.search_csv.close()
            self.ranksearch_csv.close()
            sourmash.signature.save_signatures(self.search_sigs, fp=self.search_matches)
            sourmash.signature.save_signatures(self.ranksearch_sigs, fp=self.ranksearch_matches)
            self.search_matches.close()
            self.ranksearch_matches.close()
        if self.gather:
            self.rankgather_csv.close()
