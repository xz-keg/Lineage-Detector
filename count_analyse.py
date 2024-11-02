import json
import os
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
deweighted_countries=["england",'scotland','wales','canada']

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
    if 'userOrOld' in node['node_attrs']:
        node['node_attrs']['userOrOld']['value']='highlighted sample'
    if 'children' in node:
        for child in node['children']:
            highlight_browser(child)
    return 0

def node_browser(node):
    highlight_count=0
    total_count=0
    if 'node_attrs' in node:
        if 'userOrOld' in node['node_attrs']:
            val=node['node_attrs']['userOrOld']['value']
            if val=='highlighted sample':
                highlight_count+=1
            if val in ['highlighted sample','uploaded sample']:
                total_count+=1

    if 'children' in node:
        for child in node['children']:
            ret_info=node_browser(child)
            highlight_count+=ret_info[1]
            total_count+=ret_info[0]
    
    return [total_count,highlight_count]

all_files=os.listdir()
keys=['2024-7','2024-8','2024-9','2024-10']
sum_total={}
sum_ht={}

for item in all_files:
    if '.json' in item:
        for key in keys:
            if key in item:

                f=open(item,'r')
                js=json.load(f)
                f.close()
                tree=js['tree']
                root=tree['children'][0]

                total_count,highlight_count=node_browser(root)
                print(item,total_count,highlight_count)
                if not(key in sum_total):
                    sum_total[key]=0
                    sum_ht[key]=0
                sum_total[key]+=total_count
                sum_ht[key]+=highlight_count
print(sum_total, sum_ht)