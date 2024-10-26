ycount={}
import json
import copy
def read_ref():
    ref=open("reference_seq.txt",'r')
    q=ref.readlines()
    seq=''
    for l in q:
        for w in l:
            if w in ['a','t','c','g']:
                seq=seq+w
    return seq.upper()

def designation_browser(current_node,current_mut,anno,current_ref):
    mut=[]
    if 'branch_attrs' in current_node:
        #print(current_node['branch_attrs'])
        if 'nuc' in current_node['branch_attrs']['mutations']:
            mut=current_node['branch_attrs']['mutations']['nuc']
    all_mutations=copy.deepcopy(current_mut)
    for item in mut:
        nuc_pos=int(item[1:-1])
        current_ref[:nuc_pos]=
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
    
    if "designation_date" in current_node['node_attrs']:
        #print(current_node['node_attrs']['designation_date'])
        date=current_node['node_attrs']['designation_date']
        yy=date['value'].split('-')[0]
        if not(yy in ycount):
            ycount[yy]=1
        else:
            ycount[yy]+=1
    if 'children' in current_node:
        for child in current_node['children']:
            w=designation_browser(child,all_mutations,anno,current_ref)

    return 0

def read_designation():
    w=open("des.json",'r')
    js=json.load(w)
    w.close()
    anno=js['meta']['genome_annotations']
    js['meta']['colorings'][0]['scale'].append(['highlighted sample','#CCCC00'])

    anno['ORF1a']=anno['ORF1ab']['segments'][0]
    anno['ORF1b']=anno['ORF1ab']['segments'][1]
    anno['ORF9b']={'start':28284,'end':28577}
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
print(ycount)