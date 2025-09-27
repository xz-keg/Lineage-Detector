
import json
import copy
import numpy as np
def read_ref():
    ref=open("reference_seq.txt",'r')
    q=ref.readlines()
    seq=''
    for l in q:
        for w in l:
            if w in ['a','t','c','g']:
                seq=seq+w
    return seq.upper()

all_count={}
rcount={}
ycount={}
pcount={}


def designation_browser(current_node,parent_mut,anno,current_ref):
    mut=[]
    all_mutations=[]
    if "designation_date" in current_node['node_attrs']:
        date=current_node['node_attrs']['designation_date']
        yy=date['value'].split('-')[0]
        if not(yy in ycount):
            ycount[yy]=1
        else:
            ycount[yy]+=1
        name=current_node['name']
        all_mutations=copy.deepcopy(current_node['branch_attrs']['mutations'])
        mut=all_mutations
        if len(mut)==0:
            #print(name,parent_mut,mut)
            mut=parent_mut
            
        if 'X' in name and not('.' in name):
            if not(yy in rcount):
                rcount[yy]=1
            else:
                rcount[yy]+=1
        else:
           
            for item in mut:
                if not(item in pcount):
                    pcount[item]={}
                if not(yy in pcount[item]):
                    pcount[item][yy]=1
                else:
                    pcount[item][yy]+=1
            if not('nuc' in mut):
                print(name)
           
    all_mutations=current_node['branch_attrs']['mutations']

  
    
    
    if 'children' in current_node:
        for child in current_node['children']:
            w=designation_browser(child,copy.deepcopy(all_mutations),anno,current_ref)

    return 0

def read_designation():
    w=open("des.json",'r')
    js=json.load(w)
    w.close()
    
    anno=js['meta']['genome_annotations']
    print(anno)
    js['meta']['colorings'][0]['scale'].append(['highlighted sample','#CCCC00'])

    anno['ORF9c']={'start':28734,'end':28955}
    anno['ORF3c']={'start':25457,'end':25582}
    anno['ORF3b']={'start':25814,'end':25882}
    anno['ORF3d']={'start':25524,'end':25697}
    anno['ORF3d-2']={'start':25968,'end':26069}
    anno['ORF0']={'start':107,'end':136}
    ref=read_ref()
    sp=designation_browser(js['tree'],[],anno,ref)

    return 0
read_designation() 

print(ycount,rcount,pcount)