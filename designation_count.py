ycount={}
import json
import copy
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
            w=designation_browser(child,all_mutations)

    return 0

def read_designation():
    w=open("des.json",'r')
    q=json.load(w)
    w.close()
    sp=designation_browser(q['tree'],[])
    return 0
read_designation() 
print(ycount)