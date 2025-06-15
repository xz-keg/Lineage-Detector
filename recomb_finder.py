from collections import namedtuple
import copy,json
#from nt import read

# --- Configuration ---
# A Mutation is defined by its genomic position (pos) and its specific change (alt)


REVERSION_COST = 1
NEW_MUTATION_COST = 1
BREAKPOINT_COST = 1.5
convergent_pos_list=[670,22032,22995,27810]

def extract_position(mutation):
    """
    Extracts the position from a mutation string.
    """
    return int(mutation[1:-1])

def extract_final_allele(mutation):
    """
    Extracts the final allele from a mutation string.
    """
    return mutation[-1]

def extract_first_allele(mutation):
    """
    Extracts the first allele from a mutation string.
    """
    return mutation[0]
def find_minimum_recombination_cost(v1, v2, r,max_allowed_score):
    """
    Computes the minimum of "breakpoints plus additional mutations" to form a
    recombinant from two parent variants.
    format of mutation is 'A23598G', 1st letter is the reference allele, followed by position, followed by alternative allele
    Args:
        v1 (list[Mutation]): A list of mutations for the first parent variant.
        v2 (list[Mutation]): A list of mutations for the second parent variant.
        r (list[Mutation]): A list of mutations for the recombinant variant.

    Returns:
        float: The minimum possible cost to derive the recombinant.
        ending_variant (int): 1 if v1 is the ending variant, 2 if v2 is the ending variant.
        b (list[Breakpoint]): A list of breakpoints for the minimum recombination cost.
        (Also returns the DP table for debugging and path reconstruction if needed).
        alist (list[Mutation]): A list of additional mutations for the additional mutations.
    """
    # 1. Preprocessing
    position_set=set()
    for mutation in v1+v2+r:
        position_set.add(extract_position(mutation))

    mut_dic_v1={}
    mut_dic_v2={}
    mut_dic_r={}
    mut_dic_original={}
    for mutation in v1:
        pos=extract_position(mutation)
        if pos>50:
            mut_dic_v1[extract_position(mutation)]=extract_final_allele(mutation)
            mut_dic_original[extract_position(mutation)]=extract_first_allele(mutation)
    for mutation in v2:
        pos=extract_position(mutation)
        if pos>50:
            mut_dic_v2[extract_position(mutation)]=extract_final_allele(mutation)
            mut_dic_original[extract_position(mutation)]=extract_first_allele(mutation)
    for mutation in r:
        pos=extract_position(mutation)
        if pos>50:
            mut_dic_r[extract_position(mutation)]=extract_final_allele(mutation)
            mut_dic_original[extract_position(mutation)]=extract_first_allele(mutation)

    # Collect all unique mutation positions from all three variants and sort them.
    # These positions define the boundaries of our segments.
    all_positions = sorted(list(position_set))

    dp=[0,0]  # dp[0]=currently at v1 dp[1]=currently at v2
    best_bp=[[],[]]
    best_additional_list=[[],[]]
    for i in range(1,len(all_positions)):
        #consider whether there to add a breakpoint before the current position
        pos=all_positions[i]
        if pos<=50:
            continue
        if pos in mut_dic_r:
            allele_r=mut_dic_r[pos]
        else:
            allele_r=mut_dic_original[pos]
        if pos in mut_dic_v1:
            allele_v1=mut_dic_v1[pos]
        else:
            allele_v1=mut_dic_original[pos]
        if pos in mut_dic_v2:
            allele_v2=mut_dic_v2[pos]
        else:
            allele_v2=mut_dic_original[pos]
        all_original=mut_dic_original[pos]
        v1_cost=0
        v2_cost=0
        if allele_r!=allele_v1:
            v1_cost=NEW_MUTATION_COST
            if allele_r==all_original and not(pos in convergent_pos_list):
                v1_cost=REVERSION_COST
        if allele_r!=allele_v2:
            v2_cost=NEW_MUTATION_COST
            if allele_r==all_original and not(pos in convergent_pos_list):
                v2_cost=REVERSION_COST
        dp0_new=dp[0]+v1_cost
        best_bp0_new=copy.deepcopy(best_bp[0])
        best_additional_list0=copy.deepcopy(best_additional_list[0])
        

        if dp0_new>dp[1]+BREAKPOINT_COST+v1_cost:
            dp0_new=dp[1]+BREAKPOINT_COST+v1_cost
            best_bp0_new=copy.deepcopy(best_bp[1])+[i]
            best_additional_list0=copy.deepcopy(best_additional_list[1])
        if v1_cost>0:
            addmut=allele_v1+str(pos)+allele_r
            best_additional_list0.append(addmut)
        dp1_new=dp[1]+v2_cost
        best_bp1_new=copy.deepcopy(best_bp[1])
        best_additional_list1=copy.deepcopy(best_additional_list[1])
        if dp1_new>dp[0]+BREAKPOINT_COST+v2_cost:
            dp1_new=dp[0]+BREAKPOINT_COST+v2_cost
            best_bp1_new=copy.deepcopy(best_bp[0])+[i]
            best_additional_list1=copy.deepcopy(best_additional_list[0])
        if v2_cost>0:
            addmut=allele_v2+str(pos)+allele_r
            best_additional_list1.append(addmut)
        best_additional_list[0]=best_additional_list0
        best_additional_list[1]=best_additional_list1
        dp[0]=dp0_new
        dp[1]=dp1_new
        best_bp[0]=best_bp0_new
        best_bp[1]=best_bp1_new
        if dp[0]>max_allowed_score and dp[1]>max_allowed_score:
            return max_allowed_score+1,None,None,None
    bp0_list=[]
    for item in best_bp[0]:
        bp0_list.append((all_positions[item-1],all_positions[item]))
    bp1_list=[]
    for item in best_bp[1]:
        bp1_list.append((all_positions[item-1],all_positions[item]))
    if dp[0]<dp[1]:
        return dp[0],1,bp0_list,best_additional_list[0]
    return dp[1],2,bp1_list,best_additional_list[1]


def find_best(all_candidates,r):
    min_recomb_cost=1000
    bp,additional_list=[],[]
    endv=None
    donors=[]
    namelist=[]
    for item in all_candidates:
        namelist.append(item)
    checked=0
    print(find_minimum_recombination_cost(all_candidates['LS.2.1.1'],all_candidates['LF.7'],r,1000))
    for i in range(len(namelist)):
        for j in range(i+1,len(namelist)):
            checked+=1
            if checked%1000==0:
                print(checked)
            v1=all_candidates[namelist[i]]
            v2=all_candidates[namelist[j]]
            cost,ending_variant,bp_temp,additional_list_temp=find_minimum_recombination_cost(v1,v2,r,max_allowed_score=min_recomb_cost)
            if cost<min_recomb_cost:
                min_recomb_cost=cost
                bp=bp_temp
                additional_list=additional_list_temp
                endv=ending_variant
                donors=[namelist[i],namelist[j]]
                print(min_recomb_cost,endv,bp,additional_list,donors)

    return min_recomb_cost,endv,bp,additional_list,donors

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


def fasta_reader(fasta_file):
    """
    Reads a FASTA file and returns a dictionary of sequence IDs and their corresponding sequences.
    """
    sequences = {}
    with open(fasta_file, 'r') as file:
        current_id = None
        lines=file.readlines()
        for line in lines:
            line = line.strip()
            if line.startswith('>'):
                current_id = line[1:]
                sequences[current_id] = ''
            else:
                sequences[current_id] += line
    return sequences

#judge recomb cost for a sequence
def recomb_cost(mutations):
    print("successfully decided mutations",mutations)
    cost,endv,bp,additional_list,donors=find_best(selected_variants,mutations)
    return cost,endv,bp,additional_list,donors
seq=copy.deepcopy(variant_mutation_dic['XFP'])
selected_variants={}
for key in variant_mutation_dic:
    if len(variant_mutation_dic[key])>180:
        selected_variants[key]=copy.deepcopy(variant_mutation_dic[key])
print(len(variant_mutation_dic),len(selected_variants))
selected_variants.pop('XFP')

cost,endv,bp,additional_list,donors=recomb_cost(seq)
print(key,cost,endv,bp,additional_list,donors)


