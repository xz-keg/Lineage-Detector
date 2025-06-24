import collections,json,copy
import select
import numpy as np
import re
from numba import njit, float64, int64,types
from numba.typed import Dict
from concurrent.futures import ThreadPoolExecutor

# Score constants
breakpoint_score = 1.2
mismatch_score = 1.0

# Define a robust parser and mappings for nucleotide data
NUC_TO_INT = {'A': 0, 'T': 1, 'C': 2, 'G': 3,  '-': 4,'REF': 5} 
REF_INT = 5# REF is the reference/wild-type state
INT_TO_NUC = {v: k for k, v in NUC_TO_INT.items()}
mutation_parser = re.compile(r"([ATCG])(\d+)([ATCG])")

def lineage_reader(path):
    # read from lineage_notes.txt and extract all parent-children relationships of designations.
    f=open(path,'r')
    lines=f.readlines()
    f.close()
    parent_relation={}
    prev_variant=''
    def comp(ida,idb):
        if ida==idb:
            return False
        if 'X' in ida:
            return False
        if 'X' in idb:
            return True
        if len(ida)>len(idb):
            return True
        if len(ida)<len(idb):
            return False
        if ida>idb:
            return True
        return False
    for line in lines:
        line=line.strip()
        if len(line)>5 and line[0]!='*':
            linsp=line.split('\t')[0]
            this_id=linsp.split('.')[0]
            prevsp=prev_variant.split('.')
            prev_id=prevsp[0]
            
            if len(prevsp)>=3 and comp(this_id,prev_id):
                parent_relation[this_id]=prev_variant
            prev_variant=linsp
        
    return parent_relation




def parse_mutation(mut_str):
    #print(mut_str)
    match = mutation_parser.match(mut_str)
    if not match: return None, None, None
    ref_nuc, pos, alt_nuc = match.groups()
    if int(pos)<50 or int(pos)>29700:
        return None,None,None
    
    return int(pos), ref_nuc, alt_nuc

@njit(float64(types.DictType(int64, int64), types.DictType(int64, int64), types.DictType(int64, int64)), nogil=True)
def _find_score_for_pair_jit(r_dict, p1_dict, p2_dict):
    """
    Calculates min score for a pair. This version correctly uses only Numba-supported
    operations to find the union of positions.
    """
    # 1. Get keys from each dict into a NumPy array (the Numba-safe way)
    r_keys = np.empty(len(r_dict), dtype=np.int64)
    i = 0
    for k in r_dict.keys():
        r_keys[i] = k
        i += 1

    p1_keys = np.empty(len(p1_dict), dtype=np.int64)
    i = 0
    for k in p1_dict.keys():
        p1_keys[i] = k
        i += 1

    p2_keys = np.empty(len(p2_dict), dtype=np.int64)
    i = 0
    for k in p2_dict.keys():
        p2_keys[i] = k
        i += 1

    # 2. Concatenate keys and find the sorted, unique positions using np.unique
    # This is fast and fully supported by Numba.
    combined_keys = np.concatenate((r_keys, p1_keys, p2_keys))
    if combined_keys.size == 0:
        return 0.0
    all_pos = np.unique(combined_keys)
    
    seq_len = len(all_pos)
    dp = np.zeros((seq_len, 2), dtype=float64)

    # 3. Run the DP, fetching from Numba dicts inside the compiled loop
    for i in range(seq_len):
        pos = all_pos[i]
        r_nuc = r_dict.get(pos, REF_INT)
        p1_nuc = p1_dict.get(pos, REF_INT)
        p2_nuc = p2_dict.get(pos, REF_INT)
        
        cost1 = mismatch_score if r_nuc != p1_nuc else 0.0
        cost2 = mismatch_score if r_nuc != p2_nuc else 0.0
        
        if i == 0:
            dp[0, 0], dp[0, 1] = cost1, cost2
            continue

        stay_p1 = dp[i-1, 0]
        switch_to_p1 = dp[i-1, 1] + breakpoint_score
        dp[i, 0] = cost1 + min(stay_p1, switch_to_p1)
        
        stay_p2 = dp[i-1, 1]
        switch_to_p2 = dp[i-1, 0] + breakpoint_score
        dp[i, 1] = cost2 + min(stay_p2, switch_to_p2)

    return min(dp[seq_len-1, 0], dp[seq_len-1, 1])
parent_relation=lineage_reader('lineage_notes.txt')
def min_recombination_score(v_set, r):
    # --- 1. GLOBAL SETUP ---
    global parent_relation
    r_dict_numba = Dict.empty(key_type=types.int64, value_type=types.int64)
    ref_nuc_map = {}
    for mut in r:
        pos, ref, alt = parse_mutation(mut)
        if pos is not None:
            r_dict_numba[pos] = NUC_TO_INT[alt]
            if pos not in ref_nuc_map: ref_nuc_map[pos] = ref

    variants_parsed = []
    for name, mut_list in v_set.items():
        v_dict_numba = Dict.empty(key_type=types.int64, value_type=types.int64)
        for mut in mut_list:
            pos, ref, alt = parse_mutation(mut)
            if pos is not None:
                v_dict_numba[pos] = NUC_TO_INT[alt]
                if pos not in ref_nuc_map: ref_nuc_map[pos] = ref
        variants_parsed.append({'name': name, 'data_numba': v_dict_numba})
    #print(v_set,variants_parsed)
    n_variants = len(variants_parsed)
    if n_variants == 0: return "N/A", "N/A", [], r, len(r) * mismatch_score

    # --- 2. PARALLEL PROCESSING ---
    def process_pair(pair_indices):
        i, j = pair_indices
        p1_dict = variants_parsed[i]['data_numba']
        p2_dict = variants_parsed[j]['data_numba']
        score = _find_score_for_pair_jit(r_dict_numba, p1_dict, p2_dict)
        return score, i, j

    tasks = [(i, j) for i in range(n_variants) for j in range(i+1, n_variants)]
    min_overall_score, best_pair_indices = float('inf'), []

    with ThreadPoolExecutor() as executor:
        results = executor.map(process_pair, tasks)
        for score, i, j in results:
            if score ==min_overall_score:
                best_pair_indices.append((i,j))
            if score < min_overall_score:
                min_overall_score = score
                best_pair_indices = [(i, j)]
            
            
    
    if best_pair_indices == []: return "N/A", "N/A", [], [], float('inf')
    def check_parent(v1,v2):
        #check if v1 is parent of v2
        #print(variants_parsed)
        v1=variants_parsed[v1]['name']
        v2=variants_parsed[v2]['name']
        #print(v1,v2)
        if v1 in v2:
            return True
        v2_prefix=v2.split('.')[0]
        while v2_prefix in parent_relation:
            v2=v2.replace(v2_prefix,parent_relation[v2_prefix])
            v2_prefix=v2.split('.')[0]
            if v1 in v2:
                return True
        return False
    def merge_best_pair(base_pair_indices):
        # merge recombinant candidates if two pairs can be reduced(variants in pair 1 are same/parent of variants pair 2 ), using parent_relation for variants with different prefixs. 
        returned_indices=[]
        #do brute force approach, if one option is dominated by another then remove it
        for i in range(len(best_pair_indices)):
            removed=False
            option=best_pair_indices[i]
            for j in range(len(best_pair_indices)):
                if i==j:
                    continue
                option2=best_pair_indices[j]
                #print(option,option2)
                if check_parent(option2[0],option[0]) and check_parent(option2[1],option[1]):
                    removed=True
                    break
                if check_parent(option2[1],option[0]) and check_parent(option2[0],option[1]):
                    removed=True
                    break 
            if not removed:
                returned_indices.append(option)
        return returned_indices


    returned_indices=merge_best_pair(best_pair_indices)

    # --- 3. RECONSTRUCTION (in Python, called only once) ---
    returned_tuples=[]
    for p1_idx, p2_idx in returned_indices:
    #p1_idx, p2_idx = best_pair_indices
        v1 = {k: v for k, v in variants_parsed[p1_idx]['data_numba'].items()}
        v2 = {k: v for k, v in variants_parsed[p2_idx]['data_numba'].items()}
        r_dict_py = {k: v for k, v in r_dict_numba.items()}
        
        all_pos = sorted(list(r_dict_py.keys() | v1.keys() | v2.keys()))
        if not all_pos:
            return variants_parsed[p1_idx]['name'], variants_parsed[p2_idx]['name'], [], [], 0.0
        
        seq_len = len(all_pos)
        dp, path = np.zeros((seq_len, 2)), np.zeros((seq_len, 2), dtype=int)

        for i, pos in enumerate(all_pos):
            r_nuc, p1_nuc, p2_nuc = r_dict_py.get(pos, REF_INT), v1.get(pos, REF_INT), v2.get(pos, REF_INT)
            cost1, cost2 = (mismatch_score if r_nuc != p1_nuc else 0.0), (mismatch_score if r_nuc != p2_nuc else 0.0)
            
            if i == 0:
                dp[0, 0], dp[0, 1] = cost1, cost2
                continue

            stay1, switch1 = dp[i-1, 0], dp[i-1, 1] + breakpoint_score
            dp[i, 0], path[i, 0] = (cost1 + stay1, 0) if stay1 <= switch1 else (cost1 + switch1, 1)

            stay2, switch2 = dp[i-1, 1], dp[i-1, 0] + breakpoint_score
            dp[i, 1], path[i, 1] = (cost2 + stay2, 1) if stay2 <= switch2 else (cost2 + switch2, 0)

        breakpoints, add_mutations = [], []
        current_parent = 0 if dp[seq_len-1, 0] <= dp[seq_len-1, 1] else 1
        parent_map = {0: v1, 1: v2}
        
        for i in range(seq_len - 1, -1, -1):
            pos = all_pos[i]
            if r_dict_py.get(pos, REF_INT) != parent_map[current_parent].get(pos, REF_INT):
                ref, r_nuc_int = ref_nuc_map.get(pos, 'N'), r_dict_py.get(pos, REF_INT)
                prev_nuc_int=parent_map[current_parent].get(pos, REF_INT)
                alt = INT_TO_NUC.get(r_nuc_int, '?') if r_nuc_int != REF_INT else ref
                prev = INT_TO_NUC.get(prev_nuc_int, '?') if prev_nuc_int != REF_INT else ref
                if r_nuc_int != REF_INT: add_mutations.append(f"{prev}{pos}{alt}")
                else: add_mutations.append(f"{prev}{pos}{alt}(r)")

            if i > 0:
                prev_parent = path[i, current_parent]
                if current_parent != prev_parent: breakpoints.append(all_pos[i])
                current_parent = prev_parent
        
        c1_name = variants_parsed[p1_idx]['name'] if current_parent == 0 else variants_parsed[p2_idx]['name']
        c2_name = variants_parsed[p2_idx]['name'] if current_parent == 0 else variants_parsed[p1_idx]['name']
        returned_tuples.append((c1_name, c2_name, sorted(breakpoints), list(set(add_mutations)), min_overall_score))
    return returned_tuples

def read_ref():
    ref=open("reference_seq.txt",'r')
    q=ref.readlines()
    seq=''
    for l in q:
        for w in l:
            if w in ['a','t','c','g']:
                seq=seq+w
    return seq.upper()



def designation_browser(current_node,current_mut):
    mut=[]
    if 'branch_attrs' in current_node:
        #print(current_node['branch_attrs'])
        if 'nuc' in current_node['branch_attrs']['mutations']:
            mut=current_node['branch_attrs']['mutations']['nuc']
    all_mutations=copy.deepcopy(current_mut)
    for item in mut:
        nuc_pos=int(item[1:-1])
        already=False
        for item2 in current_mut:
            if int(item2[1:-1])==nuc_pos:
                already=True
                correct_mut=item2[0]+item[1:]
                all_mutations.remove(item2)
        if not(already):
            all_mutations.append(item)
        else:
            if correct_mut[0]!=correct_mut[-1]:
                all_mutations.append(correct_mut)
    name=current_node['name']
    if not("NODE" in name):
        variant_mutation_dic[name]=all_mutations
    if 'children' in current_node:
        for child in current_node['children']:
            w=designation_browser(child,all_mutations)

    return 0

variant_mutation_dic={}
def read_designation():
    w=open("des.json",'r')
    q=json.load(w)
    w.close()
    sp=designation_browser(q['tree'],[])
    return 0
read_designation()

ref=read_ref()


seq=copy.deepcopy(variant_mutation_dic['XFG'])
selected_variants={}
for key in variant_mutation_dic:
    if not('internal' in key):    
        if not('XFG' in key):
            selected_variants[key]=copy.deepcopy(variant_mutation_dic[key])

print(min_recombination_score(selected_variants,seq))

