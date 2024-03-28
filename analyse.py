import json

import argparse

parser = argparse.ArgumentParser(description='Demo of argparse')
parser.add_argument('--filename', type=str,default='usher')
parser.add_argument('--important-threshold', type=int,default=3)
parser.add_argument('--important-ct-threshold', type=int,default=2)
parser.add_argument('--output-highlighted', type=bool,default=False)

args = parser.parse_args()
important_threshold=args.important_threshold
important_country=args.important_ct_threshold
filename=args.filename
opt=args.output_highlighted

def read_ref():
    ref=open("reference_seq.txt",'r')
    q=ref.readlines()
    seq=''
    for l in q:
        for w in l:
            if w in ['a','t','c','g']:
                seq=seq+w
    return seq.upper()




def read_table():
    fs=open("table.txt",'r')
    flines=fs.readlines()
    trans_table={}
    for l in flines:
        linsp=l.split()
        encoding=linsp[0]
        translated=linsp[2]
        trans_table[encoding]=translated
    return trans_table  

def split_ref(ref,trans_table):
    already=[]
    
    for i in range(len(ref)-3):
        if trans_table[ref[i:i+3]]=='M':
            if not(i in already):
                j=i
                #print(j,i)
                stop=False
                translated=''
                while not(stop):
                    if trans_table[ref[j:j+3]]!='O':
                        translated=translated+trans_table[ref[j:j+3]]
                    else:
                        stop=True
                    if trans_table[ref[j:j+3]]=='M':
                        already.append(j)
                        #print(j)
                    j=j+3
                    if j+3>len(ref):
                        stop=True
                print(i,j,len(translated),translated)
                
                
                
           
    return 0



ref=read_ref()
import re
def comp(a,b):
    diff=0
    for i in range(len(a)):
        if a[i]!=b[i]:
            diff+=1
    return diff
'''
for i in range(len(ref)):
    if ref[i:i+8]=='AAACGAAC':
        print(i)
    if comp(ref[i:i+8],'AAACGAAC')==1:
        print(i,ref[i:i+8])
'''


table=read_table()
#m=split_ref(ref,table)
import copy
def highlight_browser(node):
    node['node_attrs']['userOrOld']['value']='highlighted sample'
    if 'children' in node:
        for child in node['children']:
            highlight_browser(child)
    return 0
def node_browser(node,current_lineage,current_seq,mutation_from_last,backcount):
    lineage=''
    if 'pango_lineage_usher' in node['node_attrs']:
        lineage=node['node_attrs']['pango_lineage_usher']['value']
    else:
        lineage=current_lineage
    is_terminal_lineage=True
    if lineage=='':
        is_terminal_lineage=False
    count=0
   # print("node attr:",node['node_attrs'])
    #print("branch attr:",node['branch_attrs'])
    mut=[]
    if 'branch_attrs' in node:
        mut=node['branch_attrs']['mutations']['nuc']
    this_seq=current_seq
    important_mut=False
   
    if lineage!=current_lineage:
        current_mut=[]
        back_mutation_count=0
    else:
        current_mut=copy.deepcopy(mutation_from_last)
        back_mutation_count=backcount
    
    designated_mutations=[]
    
    if lineage in variant_mutation_dic:
        designated_mutations=variant_mutation_dic[lineage]
    lineage_ref=ref
    for item in designated_mutations:
        ids=int(item[1:-1])
        lineage_ref=lineage_ref[:ids-1]+item[-1]+lineage_ref[ids:]
    i=0
    while (i<len(mut)):
        item=mut[i]
        ids=int(item[1:-1])
        this_seq=this_seq[:ids-1]+item[-1]+this_seq[ids:]
        #if (lineage=='JN.1.7'):
        #   print(item,ids)
        if (item[-1]!=ref[ids-1]):
            # ignore reversions
            # check if S or Orf9b
           
            if (item[-1]!=lineage_ref[ids-1]) and (lineage_ref[ids-1]!='-'):
                for annoitem in anno:
                    if 'start' in anno[annoitem]:
                        
                        if ids>=anno[annoitem]['start'] and ids<=anno[annoitem]['end']:
                            start=anno[annoitem]['start']
                            nuc_st=ids-(ids-start)%3
                            aa=table[current_seq[nuc_st-1:nuc_st+2]]
                            old_aa=table[this_seq[nuc_st-1:nuc_st+2]]
                            if aa!=old_aa:
                                # credible mutations: mutations on S, Orf9b, Orf9c, normal to O, or M to others
                                if (annoitem in ['S','ORF9b']) or (aa=='O') or (nuc_st==start):
                                    important_mut=True
            else:
                
                mut.remove(item)
                i-=1
        else:
            # if the position is not already reverted in designated
            if lineage_ref[ids-1]!=ref[ids-1] and lineage_ref[ids-1]!='-':
                back_mutation_count+=1
            else:
                mut.remove(item)
                i-=1
        i+=1
        
    #if lineage=='JN.1.7':
    #    print("show off:", current_mut, mut)
    current_mut.extend(mut)

    country_list=[]
    date_list=[]
    if 'hCoV' in node['name']:
        count+=1
        #print(node['name'],count)
        if back_mutation_count>5:
            print(lineage,node['name'],current_mut)
        ct=node['name'].split('/')[1]
        if not(ct in country_list):
            country_list.append(ct)
        if 'GBW' in node['name']:
            country_list.append('GBW')
        date=node['name'].split('|')[-1]
        date_list.append(date)
        
    else:
        if 'children' in node:
            for child in node['children']:
                ret_info=node_browser(child,lineage,this_seq,current_mut,back_mutation_count)
                if ret_info[0]!=lineage:
                    is_terminal_lineage=False
                    
                count+=ret_info[1]
                for ct in ret_info[2]:
                    if not(ct in country_list):
                        country_list.append(ct)
                if len(date_list)<=1:
                    for ct in ret_info[3]:
                        if not(ct in date_list):
                            date_list.append(ct)
            #print(node['name'],important_mut)
    if not(is_terminal_lineage):
        lineage='not terminal'
    if important_mut:
        if is_terminal_lineage:
            if (count>=important_threshold or (len(country_list)>=important_country and count>=2 and len(date_list)>=2) ):
                imp_mut=[]
                for item in current_mut:
                    ids=int(item[1:-1])
                    for annoitem in anno:
                        if ('start' in anno[annoitem]) and (annoitem!='nuc') and (this_seq[ids-1]==item[-1]):
                            if ids>=anno[annoitem]['start'] and ids<=anno[annoitem]['end']:
                                #print(item,annoitem,this_seq[ids-1])
                                start=anno[annoitem]['start']
                                nuc_st=ids-(ids-start)%3
                                old_aa=table[ref[nuc_st-1:nuc_st+2]]
                                aa=table[this_seq[nuc_st-1:nuc_st+2]]
                                
                                #print(old_aa,aa)
                                if aa!=old_aa:
                                    # credible mutations: mutations on S, Orf9b, Orf9c, normal to O, or M to others
                                    if (annoitem in ['S','ORF9b']) or (aa=='O') or (nuc_st==start):
                                        muta=annoitem+':'+old_aa+str(int((nuc_st-start)/3)+1)+aa
                                        imp_mut.append(muta)
                if len(imp_mut)>0:
                    print(node['name'],lineage,','.join(current_mut),count,len(country_list),imp_mut)
                w=highlight_browser(node)
    return [lineage,count,copy.deepcopy(country_list),copy.deepcopy(date_list)]

variant_mutation_dic={}

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

def read_designation():
    w=open("des.json",'r')
    q=json.load(w)
    w.close()
    sp=designation_browser(q['tree'],[])
    return 0
read_designation() 
#print(len(variant_mutation_dic))

f=open(filename+".json",'r')
js=json.load(f)
f.close()
anno=js['meta']['genome_annotations']
js['meta']['colorings'][0]['scale'].append(['highlighted sample','#CCCC00'])

anno['ORF1a']=anno['ORF1ab']['segments'][0]
anno['ORF1b']=anno['ORF1ab']['segments'][1]
anno['ORF9b']={'start':28284,'end':28577}
anno['ORF9c']={'start':28734,'end':28955}
anno['ORF3c']={'start':25457,'end':25582}

tree=js['tree']


root=tree['children'][0]

l,count,countries,dates=node_browser(root,'',ref,[],0)
outjs=open(filename+".json",'w')
json.dump(js,outjs)

'''
while 'children' in node_main:
    max_length=0
    for child in node_main['children']:
        lenth=len(str(child))
        if lenth>max_length:
            best_child=child
            max_length=lenth
    node_main=best_child
    print(node_main['name'])
    print(node_main['branch_attrs'])
    print(node_main['node_attrs'])
'''
