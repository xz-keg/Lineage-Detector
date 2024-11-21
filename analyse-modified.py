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
deweighted_countries=['united_kingdom','canada']

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
def read_tsv():
	ans=[]
	for l in tsvlines[1:-1]:
		linsp=l.split("\t")
		exmut=[]
		name=linsp[0]
		insertion=linsp[10].split(", ")
		deletion=linsp[12].split(", ")
		masked=linsp[13].split(",")
		# for m in linsp[12].split(", "):
		# 	for n in range(int(m.split(":")[0].split("-")[1])-int(m.split(":")[0].split("-")[0])+1):
		# 		exmut.append(ref[int(m.split(":")[0].split("-")[0])+n+1]+str(int(m.split(":")[0].split("-")[0])+n)+"-")
		# for m in linsp[13].split(", "):
		# 	if len(m)>0:
		# 		exmut.append(m)
		ans.append([name,insertion,deletion,masked])
	return ans 
	
def search_tsv(name):
	ans=[]
	if len(tsvlines)>0:
		for l in tsvlines[1:]:
			linsp=l.split("\t")
			exmut=[]
			if name ==linsp[0]:			
				insertion=linsp[10].split(", ")
				deletion=linsp[12].split(", ")
				masked=linsp[13].split(",")
				return [insertion,deletion,masked]
	else:
		return[[],[],[]]
	 

import copy
def highlight_browser(node):
	if 'userOrOld' in node['node_attrs']:
		node['node_attrs']['userOrOld']['value']='highlighted sample'
	if 'children' in node:
		for child in node['children']:
			highlight_browser(child)
	return 0
def node_browser(node,current_lineage,current_seq,mutation_from_last,backcount):
	lineage=''
	if 'pango_lineage_usher' in node['node_attrs']:
		lineage=node['node_attrs']['pango_lineage_usher']['value'].split("_")[0]
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
		mut=copy.deepcopy(node['branch_attrs']['mutations']['nuc'])
	this_seq=current_seq
	important_mut=False
   
	if lineage!=current_lineage:
		current_mut=[]
		back_mutation_count=0
	else:
		current_mut=copy.deepcopy(mutation_from_last)
		back_mutation_count=backcount
	
	designated_mutations=[]
		
	while lineage not in variant_mutation_dic:
		if "." in lineage:
			lineage=lineage.rsplit(".",1)[0]
		elif lineage in list(alias):
			lineage = alias[lineage]
		else:
			break
	else:
		designated_mutations=variant_mutation_dic[lineage]
	lineage_ref=ref
	for item in designated_mutations:
		if item[0]!= "i":
			ids=int(item[1:-1])
			lineage_ref=lineage_ref[:ids-1]+item[-1]+lineage_ref[ids:]

	country_list=[]
	date_list=[]
	exmut_dict=dict()
	i=0
	n_glycan=False
	while (i<len(mut) and len(mut)>0):
		item=mut[i]
		ids=int(item[1:-1])
		this_seq=this_seq[:ids-1]+item[-1]+this_seq[ids:]
		#if (lineage=='JN.1.7'):
		#   print(item,ids)
   
		if (item[-1]!=ref[ids-1]):
			# ignore reversions
			# check if S or Orf9b

			if (item[-1]!=lineage_ref[ids-1]) and (lineage_ref[ids-1]!='-'):				
				ins_artefact = False
				for itemi in designated_mutations:
					if itemi[0]=="i":
						idi=int(itemi[3:].split(":")[0])
						inserted=itemi[3:].split(":")[1]
						if ids > idi and ids <= idi + len(inserted):
							if item[-1]==inserted[ids-idi-1]:
								ins_artefact = True
								mut.remove(item)
								i-=1
								break 
				if ins_artefact:
					continue
				for annoitem in anno:
					if 'start' in anno[annoitem]:
						
						if ids>=anno[annoitem]['start'] and ids<=anno[annoitem]['end']:
							start=anno[annoitem]['start']
							nuc_st=ids-(ids-start)%3
							old_aa=table[current_seq[nuc_st-1:nuc_st+2]]
							aa=table[this_seq[nuc_st-1:nuc_st+2]]
							if aa!=old_aa:
								# credible mutations: mutations on S, Orf9b, Orf9c, normal to O, or M to others
								if (annoitem in ['S','ORF9b']) or (aa=='O') or (nuc_st==start):
									important_mut=True
									if annoitem=='S' and not(n_glycan):
										pot=40
										local_ref=lineage_ref[nuc_st-41:nuc_st+21]
										ii=0
										while ii<len(local_ref):
											if local_ref[ii]!='-':
												ii+=1
											else:
												local_ref=local_ref[:ii]+local_ref[ii+1:]
												if ii<pot:
													pot-=1				
										aa_1=table[local_ref[pot-3:pot]]								
										aa_2=table[local_ref[pot-6:pot-3]]
										aa_n1=table[local_ref[pot+3:pot+6]]
										aa_n2=table[local_ref[pot+6:pot+9]]
										if aa_2=='N' and aa_1!='P' and (aa in ['S','T'] and not(old_aa in ['S','T'])):
											n_glycan=True
										if aa_1=='N' and old_aa=='P' and aa_n1 in ['S','T']:
											n_glycan=True
										if aa=='N' and aa_n1!='P' and aa_n2 in ['S','T']:
											n_glycan=True
													


			elif  (lineage_ref[ids-1]!='-'):
				# if 'hCoV-19/' in node['name'] and (lineage_ref[ids-1]=='-'):
				# 	deletion=search_tsv(node['name'])[1]
				# 	for itemd in deletion:
				# 		if ids >= int(itemd.split(":")[0].split("-")[0])-3 or ids <= int(itemd.split(":")[0].split("-")[1])+3:
				# 			if item in list(exmut_dict):
				# 				exmut_dict[item]+=1
				# 			else:
				# 				exmut_dict.update({item:1})	
				mut.remove(item)
				i-=1
		else:
			# if the position is not already reverted in designated
			#if lineage_ref[ids-1]!=ref[ids-1] and lineage_ref[ids-1]!='-':
			if lineage_ref[ids-1]!=ref[ids-1]:
				back_mutation_count+=1
			else:
				mut.remove(item)
				i-=1
		i+=1
		
	#if lineage=='JN.1.7':
	#	print("show off:", current_mut, mut)	
	
	current_mut.extend(mut)

	if 'hCoV-19/' == node['name'][:8]:
		count+=1
		exmut=[]
		#print(node['name'])		
		insertion=search_tsv(node['name'])[0]
		deletion=search_tsv(node['name'])[1]
		masked=search_tsv(node['name'])[2]
		# 	for n in range(int(m.split(":")[0].split("-")[1])-int(m.split(":")[0].split("-")[0])+1):
		# 		exmut.append(ref[int(m.split(":")[0].split("-")[0])+n+1]+str(int(m.split(":")[0].split("-")[0])+n)+"-")
		#this_seq=lineage_ref
		i=0
		while (i<len(insertion)):
			item=insertion[i]
			if len(item)>0:
				ids=int(item.split(":")[0].split("-")[0])
				exists_ref = False
				for itemi in designated_mutations:
					if itemi[0]=="i":
						idi=int(itemi[3:].split(":")[0])
						inserted=itemi[3:].split(":")[1]
						if ids > idi - len(inserted) and ids <= idi + len(inserted):
							exists_ref = True
							break 
				if exists_ref:
					if len(item.split(":")[1].replace("A","").replace("C","").replace("G","").replace("T",""))==0:
						if ids < idi:
							local_ref=lineage_ref[ids-3:idi]+inserted+lineage_ref[idi:idi+3]
							thisseq=lineage_ref[ids-3:ids]+item.split(":")[1]+lineage_ref[ids:idi+3]
						else:
							local_ref=lineage_ref[idi-3:idi]+inserted+lineage_ref[idi:ids+3]
							thisseq=lineage_ref[idi-3:ids]+item.split(":")[1]+lineage_ref[ids:ids+3]	
						if local_ref==thisseq:
							insertion.remove(item)
							i-=1
						else:
							if "ins"+item.split(":")[0].split("-")[0]+":"+item.split(":")[1] not in exmut:
								exmut.append("ins"+item.split(":")[0].split("-")[0]+":"+item.split(":")[1])
				else:
					if "ins"+item.split(":")[0].split("-")[0]+":"+item.split(":")[1] not in exmut:
						exmut.append("ins"+item.split(":")[0].split("-")[0]+":"+item.split(":")[1])
			i+=1
		i=0
		while (i<len(deletion)):
			item=deletion[i]
			if len(item)>0:
				local_ref=lineage_ref[int(item.split(":")[0].split("-")[0])-4: int(item.split(":")[0].split("-")[1])+3]				
				for ids in range(int(item.split(":")[0].split("-")[0])-4, int(item.split(":")[0].split("-")[1])+3):
					if lineage_ref[ids] == "-":
						this_seq=this_seq[:ids]+ref[ids]+this_seq[ids+1:]
				for ids in range(int(item.split(":")[0].split("-")[0]), int(item.split(":")[0].split("-")[1])+1):
					this_seq=this_seq[:ids-1]+"-"+this_seq[ids:]
				thisseq=this_seq[int(item.split(":")[0].split("-")[0])-4: int(item.split(":")[0].split("-")[1])+3]
				if local_ref.replace("-","")==thisseq.replace("-",""):
					deletion.remove(item)
					i-=1
				else:
					if "del"+item.split(":")[0] not in exmut:
						exmut.append("del"+item.split(":")[0])
			i+=1
		i=0
		while (i<len(masked)):
			item=masked[i]
			if len(item)>0:
				ids=int(item[1:-1])
				if lineage_ref[ids-1]==item[-1]:
					masked.remove(item)
					i-=1
				else:
					if item not in exmut:
						exmut.append(item)
			i+=1
		i=0
		while (i<len(mut)):
			item=mut[i]
			ids=int(item[1:-1])
			if lineage_ref[ids-1]=='-':
				if this_seq[ids-1]!='-':
					if item not in exmut:
						exmut.append(item)
				mut.remove(item)
				i-=1
			i+=1
		i=0
		while i<len(current_mut):
			item=current_mut[i]
			ids=int(item[1:-1])
			if lineage_ref[ids-1]=='-':
				if this_seq[ids-1]!='-':
					exmut.append(item)
				current_mut.remove(item)
				i-=1
			i+=1
		#print(node['name'],count)
		if back_mutation_count>=5:
			temp_lin,temp_mut,bkct = recomb_helper(copy.deepcopy(lineage),copy.deepcopy(current_mut))
			for i in range(len(temp_mut)):
				ids = int(temp_mut[i][1:-1])
				if ref[ids-1]==temp_mut[i][-1]:
					temp_mut[i]+="r"
			print(lineage,node['name'],temp_lin,",".join(temp_mut),bkct,",".join(exmut))
			w=highlight_browser(node)
		ct=node['name'].split('/')[1]
		if ct in ["env", "ENV", "deer", "bat"]:
			ct=node['name'].split('/')[2]
		if ct in ["England","Scotland","Northern_Ireland","Wales"]:
			ct = "United_Kingdom"
		if ct in ["CHN","Beijing","Tianjin","Tianjn","Hebei","Henan","Shanxi","Shandong","Heilongjiang","Jilin","Liaoning","Inner_Mongolia","Neimenggu","Shaanxi","Ningxia","Gansu","Xinjiang","Tibet","Xizang","Qinghai","Shanghai","Jiangsu","Zhejiang","Anhui","Fujian","Jiangxi","Hunan","Hubei","Guangdong","Hainan","Guangxi","Yunnan","Guizhou","Sichuan","Chongqing"]:
			ct = "China"
		if ct == "Macao":
			ct = "Macau"
		if ct == "Crimea":
			ct = "Ukraine"
		if ct == "Kosovo":
			ct == "Serbia"
		if not(ct in country_list):
			country_list.append(ct)
		if 'GBW' in node['name']:
			country_list.append('GBW')
		date=node['name'].split('|')[-1]
		date_list.append(date)
		for m in exmut:
			if m in list(exmut_dict):
				exmut_dict[m]+=1
			else:
				exmut_dict.update({m:1})				
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
				for ct in list(ret_info[4]):
					if ct in list(exmut_dict):
						exmut_dict[ct]+=ret_info[4][ct]
					else:
						exmut_dict.update({ct:ret_info[4][ct]})
			#print(node['name'],important_mut)
	i=0
	while (i<len(mut)):
		item=mut[i]
		ids=int(item[1:-1])
		if lineage_ref[ids-1]=='-':
			mut.remove(item)
			i-=1
		i+=1
	if not(is_terminal_lineage):
		lineage='not terminal'
	for m in list(exmut_dict):
		if exmut_dict[m]>=2 and exmut_dict[m] >= count/2:
			if m[0]=="d":
				for annoitem in anno:
					if ('start' in anno[annoitem]) and (annoitem!='nuc'):
						if int(m[3:].split("-")[1])>=anno[annoitem]['start'] and int(m[3:].split("-")[0])<=anno[annoitem]['start']:
								important_mut = True
						if int(m[3:].split("-")[1])>=anno[annoitem]['end'] and int(m[3:].split("-")[0])<=anno[annoitem]['end']:
								important_mut = True
				if int(m[3:].split("-")[1])>=anno['S']['start'] and int(m[3:].split("-")[0])<=anno['S']['end']:
							important_mut = True
				if int(m[3:].split("-")[1])>=anno['ORF9b']['start'] and int(m[3:].split("-")[0])<=anno['ORF9b']['end']:
							important_mut = True
			if m[0]=="i":
				if int(m[3:].split(":")[0])>=anno['S']['start'] and int(m[3:].split(":")[0])<=anno['S']['end']:
							important_mut = True
				if int(m[3:].split(":")[0])>=anno['ORF9b']['start'] and int(m[3:].split(":")[0])<=anno['ORF9b']['end']:
							important_mut = True
	if important_mut:
		if is_terminal_lineage:
			if len(country_list)==1:
				if country_list[0].lower().strip() in deweighted_countries:	 #ignore groups from uk and canada.
					count=count*0.5
					for item in list(exmut_dict):
						exmut_dict[item]=exmut_dict[item]*0.5
			if (count>=important_threshold or (len(country_list)>=important_country and count>=2 and len(date_list)>=2) or (count>=2 and n_glycan)):
				i=0
				while (i<len(current_mut)):
					item=current_mut[i]
					ids=int(item[1:-1])
					if lineage_ref[ids-1]=='-':
						current_mut.remove(item)
						i-=1
					i+=1
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
								ref_aa=""
								if not "-" in lineage_ref[nuc_st-1:nuc_st+2]:
									ref_aa=table[lineage_ref[nuc_st-1:nuc_st+2]]
								elif lineage_ref[nuc_st-1:nuc_st+2]=="---":
									ref_aa="-"
								elif lineage_ref[nuc_st-1]=="-":
									i=0
									while lineage_ref[nuc_st-1-3*i]=="-":
										i+=1
									temp=lineage_ref[nuc_st-1-3*i:nuc_st+2].replace("-","")
									if len(temp)==3:
										ref_aa=table[temp]
								elif lineage_ref[nuc_st+1]=="-":
									i=0
									while lineage_ref[nuc_st+1+3*i]=="-":
										i+=1
									temp=lineage_ref[nuc_st-1:nuc_st+2+3*i].replace("-","")
									if len(temp)==3:
										ref_aa=table[temp]
								aa=table[this_seq[nuc_st-1:nuc_st+2]]								
								#print(old_aa,aa)
								if aa!=old_aa and aa!=ref_aa:
									# credible mutations: mutations on S, Orf9b, Orf9c, normal to O, or M to others
									#if (annoitem in ['S','ORF9b']) or (aa=='O') or (nuc_st==start):
										muta=annoitem+':'+old_aa+str(int((nuc_st-start)/3)+1)+aa
										if muta not in imp_mut:
											imp_mut.append(muta)
				for item in list(exmut_dict):
					if exmut_dict[item] >= count/2:						
						if item[0]=="d":
							ids=int(item[3:].split("-")[0])
							ide=int(item[3:].split("-")[1])
							if (ide-ids)%3 == 2:
								for annoitem in anno:
									if ('start' in anno[annoitem]) and (annoitem!='nuc') :		
										if ids <= anno[annoitem]['end'] and ide >= anno[annoitem]['start']:				
											start=anno[annoitem]['start']
											#print(item,annoitem,anno[annoitem]['start'],anno[annoitem]['end'],ids,ide)
											nuc_st=ids-(ids-start)%3-3
											nuc_en=ide-(ide-start)%3+3
											thisseq=this_seq
											for ids in range(len(lineage_ref)):
												if lineage_ref[ids] == "-":
													thisseq=thisseq[:ids]+ref[ids]+thisseq[ids+1:]
											for ids in range(int(item[3:].split(":")[0].split("-")[0]), int(item[3:].split(":")[0].split("-")[1])+1):
												thisseq=thisseq[:ids-1]+"-"+thisseq[ids:]
											for item2 in list(exmut_dict):
												if item2[0]!="d" and item2[0]!="i" and exmut_dict[item2]>=exmut_dict[item]/2:
													id=int(item2[1:-1])
													if id in range(nuc_st-1,nuc_en+2) and thisseq[id-1]!="-":
														thisseq=thisseq[:id-1]+item2[-1]+thisseq[id:]
											#print(nuc_st,nuc_en,lineage,ref[nuc_st-1:nuc_en+2],lineage_ref[nuc_st-1:nuc_en+2],thisseq[nuc_st-1:nuc_en+2])
											old_aa=""
											ref_aa=""
											aa=""
											temp=ref[nuc_st-1:nuc_en+2]
											for i in range(0,len(temp),3):
												old_aa+=table[temp[i:i+3]]
											temp=lineage_ref[nuc_st-1:nuc_en+2]
											temp=temp.replace("-","")
											if len(temp)%3==0:
												for i in range(0,len(temp),3):
													ref_aa+=table[temp[i:i+3]]
											else:
												print("Error:", lineage, annoitem, nuc_st, nuc_en+2, "lineage_ref", lineage_ref[nuc_st-1:nuc_en+2])
											temp=thisseq[nuc_st-1:nuc_en+2]
											temp=temp.replace("-","")
											if len(temp)%3==0:
												for i in range(0,len(temp),3):
													aa+=table[temp[i:i+3]]
											else:
												print("Error:", lineage, annoitem, nuc_st, nuc_en+2, "thisseq", thisseq[nuc_st-1:nuc_en+2])
											#print(old_aa,ref_aa,aa)
											if aa!=old_aa and aa!= ref_aa:
												# credible mutations: mutations on S, Orf9b, Orf9c, normal to O, or M to others
												if len(aa)==len(ref_aa):
													for i in range(int((nuc_en-nuc_st)/3+1)):
														if lineage_ref[nuc_st-1+3*i]=="-":
															ref_aa=ref_aa[:i]+"-"+ref_aa[i:]
															aa=aa[:i]+"-"+aa[i:]
													for i in range(len(old_aa)):
														if (old_aa[i]!=aa[i]) and (ref_aa[i]!=aa[i]):
															muta=annoitem+':'+old_aa[i]+str(int((nuc_st-start)/3)+1+i)+aa[i]
															if muta not in imp_mut:
																imp_mut.append(muta)
												else:
													while len(aa)>0:
														if (old_aa[-1]==aa[-1]) or (ref_aa[-1]==aa[-1]):
															old_aa=old_aa[:-1]
															ref_aa=ref_aa[:-1]
															aa=aa[:-1]
															nuc_en-=3
														elif (old_aa[0]==aa[0]) or (ref_aa[0]==aa[0]):
															old_aa=old_aa[1:]
															ref_aa=ref_aa[1:]
															aa=aa[1:]
															nuc_st+=3
														else:
															break
													aa=aa.replace("-","")
													if len(aa)==0:
														aa="-"
													#if (annoitem in ['S','ORF9b']) or (aa=='O'):
													if nuc_st==nuc_en:
														muta=annoitem+':'+old_aa+str(int((nuc_st-start)/3)+1)+aa
													else:
														muta=annoitem+':'+old_aa+str(int((nuc_st-start)/3)+1)+"-"+str(int((nuc_en-start)/3)+1)+aa
													if muta not in imp_mut:
														imp_mut.append(muta)
												# if nuc_st<=start and nuc_en>=start:
												# 	muta = annoitem+':M1-'
												# 	if muta not in imp_mut:
												# 		imp_mut.append(muta)
						elif item[0]=="i":
							ids=int(item[3:].split(":")[0])
							inserted=item[3:].split(":")[1]
							if len(inserted)%3 == 0:
								for annoitem in anno:
									if ('start' in anno[annoitem]) and (annoitem!='nuc'):
										if ids>=anno[annoitem]['start'] and ids<=anno[annoitem]['end']:
											thisseq=this_seq
											start=anno[annoitem]['start']
											nuc_st=ids-(ids-start)%3
											nuc_en=nuc_st
											for item2 in list(exmut_dict):
												if item2[0]!="d" and item2[0]!="i" and exmut_dict[item2]>=exmut_dict[item]/2:
													id=int(item2[1:-1])
													if id in range(nuc_st-1,nuc_st+2) and id != ids and thisseq[id-1]!="-":
														thisseq=thisseq[:id-1]+item2[-1]+thisseq[id:]
											ref_aa=""
											exists_ref=False
											for itemi in designated_mutations:
												if itemi[0]=="i":
													idi=int(itemi[3:].split(":")[0])
													ref_ins=itemi[3:].split(":")[1]
													if ids > idi - len(ref_ins) and ids <= idi + len(ref_ins):
														exists_ref = True
														break 
											if exists_ref:
												if ids < idi:
													nuc_en=idi-(idi-start)%3													
												else:
													nuc_st=idi-(idi-start)%3													
											temp=lineage_ref[nuc_st-1:idi]+ref_ins+lineage_ref[idi:nuc_en+2]
											if temp[0]=="-":
												i=0
												while lineage_ref[nuc_st-1-3*i]=="-":
													i+=1
												temp=lineage_ref[nuc_st-1-3*i:nuc_st-1]+temp
											if temp[-1]=="-":
												i=0
												while lineage_ref[nuc_en+1+3*i]=="-":
													i+=1
												temp=temp+lineage_ref[nuc_en+2:nuc_en+2+3*i]
											if len(temp)%3==0:
												for i in range(0,len(temp),3):
													ref_aa+=table[temp[i:i+3]]
											old_aa=""
											temp=ref[nuc_st-1:nuc_en+2]
											for i in range(0,len(temp),3):
												old_aa+=table[temp[i:i+3]]
											aa=""
											temp=thisseq[nuc_st-1:ids]+inserted+thisseq[ids:nuc_en+2]
											if temp[0]=="-":
												i=0
												while thisseq[nuc_st-1-3*i]=="-":
													i+=1
												temp=thisseq[nuc_st-1-3*i:nuc_st-1]+temp
											if temp[-1]=="-":
												i=0
												while thisseq[nuc_en+1+3*i]=="-":
													i+=1
												temp=temp+thisseq[nuc_en+2:nuc_en+2+3*i]
											if len(temp)%3==0:
												for i in range(0,len(temp),3):
													if temp[i:i+3] in table:
														aa+=table[temp[i:i+3]]
													else:
														aa+="X"
											#print(old_aa, ref_aa, aa)
											if aa!= old_aa and aa!=ref_aa:
												if len(old_aa)>1:
													if ids < idi:
														while nuc_st < nuc_en:
															if aa[0]!=old_aa[0] and aa[0]!=ref_aa[0]:
																muta=annoitem+':'+old_aa[0]+str(int((nuc_st-start)/3)+1)+aa[0]
																if muta not in imp_mut:
																	imp_mut.append(muta)
															old_aa=old_aa[1:]
															ref_aa=ref_aa[1:]
															aa=aa[1:]
															nuc_st+=3
													else:
														while nuc_st < nuc_en:
															if aa[-1]!=old_aa[-1] and aa[-1]!=ref_aa[-1]:
																muta=annoitem+':'+old_aa[-1]+str(int((nuc_en-start)/3)+1)+aa[-1]
																if muta not in imp_mut:
																	imp_mut.append(muta)
															old_aa=old_aa[:-1]
															ref_aa=ref_aa[:-1]
															aa=aa[:-1]
															nuc_en-=3
												if ref_aa[0]==old_aa or aa[0]==old_aa:
													if aa[0]!=old_aa:
														muta=annoitem+':'+old_aa[0]+str(int((nuc_st-start)/3)+1)+aa[0]
														if muta not in imp_mut:
															imp_mut.append(muta)
													if aa[1:]!=ref_aa[1:]:
														muta=annoitem+':ins'+str(int((nuc_st-start)/3)+1)+aa[1:]
														if muta not in imp_mut:
															imp_mut.append(muta)
												elif ref_aa[-1]==old_aa or aa[-1]==old_aa:
													if aa[-1]!=old_aa:
														muta=annoitem+':'+old_aa[-1]+str(int((nuc_st-start)/3)+1)+aa[-1]
														if muta not in imp_mut:
															imp_mut.append(muta)
													if aa[:-1]!=ref_aa[:-1]:
														muta=annoitem+':ins'+str(int((nuc_st-start)/3))+aa[:-1]
														if muta not in imp_mut:
															imp_mut.append(muta)
												else:
													if aa!=ref_aa:
														muta=annoitem+':'+old_aa+str(int((nuc_st-start)/3))+aa
														if muta not in imp_mut:
															imp_mut.append(muta)
						else:
							ids=int(item[1:-1])
							for annoitem in anno:
								if ('start' in anno[annoitem]) and (annoitem!='nuc'):
									if ids>=anno[annoitem]['start'] and ids<=anno[annoitem]['end']:
										thisseq=this_seq
										thisseq=thisseq[:ids-1]+item[-1]+thisseq[ids:]
										#print(item,annoitem,this_seq[ids-1])
										start=anno[annoitem]['start']
										nuc_st=ids-(ids-start)%3
										for item2 in list(exmut_dict):
											if item2[0]!="d" and item2[0]!="i" and exmut_dict[item2]>=exmut_dict[item]/2:
												id=int(item2[1:-1])
												if id in range(nuc_st-1,nuc_st+2) and id != ids and thisseq[id-1]!="-":
													thisseq=thisseq[:id-1]+item2[-1]+thisseq[id:]
										old_aa=table[ref[nuc_st-1:nuc_st+2]]
										ref_aa=""
										if not "-" in lineage_ref[ids-4:ids+3]:
											ref_aa=table[lineage_ref[nuc_st-1:nuc_st+2]]
										else:
											ref_aa="-"
										# elif lineage_ref[nuc_st+1]=="-":
										# 	ref_aa="="
										# elif lineage_ref[nuc_st-1]=="-":
										# 	i=0
										# 	while lineage_ref[nuc_st+1+3*i]=="-":
										# 		i+=1
										# 	temp=lineage_ref[nuc_st-1:nuc_st+2+3*i].replace("-","")
										# 	if len(temp)==3:
										# 		ref_aa=table[temp]
										aa=table[thisseq[nuc_st-1:nuc_st+2]]
										#print(old_aa,aa)
										if aa!=old_aa and aa!=ref_aa and ref_aa!="-":
											# credible mutations: mutations on S, Orf9b, Orf9c, normal to O, or M to others
											#if (annoitem in ['S','ORF9b']) or (aa=='O') or (nuc_st==start):
											muta=annoitem+':'+old_aa+str(int((nuc_st-start)/3)+1)+aa
											if muta not in imp_mut:
												imp_mut.append(muta)
				if len(imp_mut)>0:
					important_mut = False
					imp_mut_dict={}
					imp_mut_list=[]
					for item in imp_mut:
						if item.split(":")[0] in list(imp_mut_dict):
							imp_mut_dict[item.split(":")[0]].append(item.split(":")[1])
						else:
							imp_mut_dict.update({item.split(":")[0]:[item.split(":")[1],]})
					for item in list(imp_mut_dict):
						if item in ["S", "ORF9b"]:
							important_mut=True
						mutlist=imp_mut_dict[item]
						idlist=[]
						for item2 in mutlist:
							if "-" in item2[:-1]:
								while item2[-1]in"ABCDEFGHIJKLMNOPQRSTUVWXYZ-":
									item2 = item2[:-1]
								while item2[0]in"ABCDEFGHIJKLMNOPQRSTUVWXYZ-":
									item2 = item2[1:]
								idlist.append([int(item2.split("-")[0]),int(item2.split("-")[1])])
							else:
								if item2[0]=="i":
									while item2[-1]in"ABCDEFGHIJKLMNOPQRSTUVWXYZ-ins":
										item2 = item2[:-1]
									while item2[0]in"ABCDEFGHIJKLMNOPQRSTUVWXYZ-ins":
										item2 = item2[1:]
									idlist.append([int(item2)+0.5,int(item2)+0.5])
								else:
									while item2[-1]in"ABCDEFGHIJKLMNOPQRSTUVWXYZ-":
										item2 = item2[:-1]
									while item2[0]in"ABCDEFGHIJKLMNOPQRSTUVWXYZ-":
										item2 = item2[1:]
									idlist.append([int(item2),int(item2)])
						i=0
						while i < len(idlist):
							j=i+1
							while j < len(idlist):
								#print(idlist[i][0],anno[item]['start']+3*(idlist[i][0]-1)-1, lineage_ref[anno[item]['start']+3*(idlist[i][0]-1)-1],idlist[j][0],anno[item]['start']+3*(idlist[j][0]-1)-1, lineage_ref[anno[item]['start']+3*(idlist[j][0]-1)-1])
								if idlist[j][0]>=idlist[i][0] and idlist[j][1]<=idlist[i][1]:
									idlist.remove(idlist[j])
									mutlist.remove(mutlist[j])
									i-=1
									j-=1
								elif idlist[i][0]>=idlist[j][0] and idlist[i][1]<=idlist[j][1]:
									idlist.remove(idlist[i])
									mutlist.remove(mutlist[i])
									i-=1
									j-=1
									break								
								elif lineage_ref[anno[item]['start']+3*(int(idlist[j][0])-1)-1]=="-":
									temp=idlist[j][0]
									while lineage_ref[anno[item]['start']+3*(int(temp)-1)-1]=="-":
										temp-=1
									#print(idlist[j][0],temp,idlist[i][0],int(temp) == int(idlist[i][0]))
									if temp == idlist[i][0]:
										idlist.remove(idlist[j])
										mutlist.remove(mutlist[j])
										i-=1
										j-=1
									# temp=idlist[j][0]
									# while lineage_ref[anno[item]['start']+3*(temp-1)-1]=="-":
									# 	temp+=1
									# if temp == idlist[i][0]:
									# 	idlist.remove(idlist[j])
									# 	mutlist.remove(mutlist[j])
									# 	i-=1
									# 	j-=1
								elif lineage_ref[anno[item]['start']+3*(int(idlist[i][0])-1)-1]=="-":
									temp=idlist[i][0]
									while lineage_ref[anno[item]['start']+3*(int(temp)-1)-1]=="-":
										temp-=1
									#print(idlist[j][0],temp,idlist[i][0],int(temp) == int(idlist[j][0]))
									if temp == idlist[j][0]:
										idlist.remove(idlist[i])
										mutlist.remove(mutlist[i])
										i-=1
										j-=1
										break
									# temp=idlist[i][0]
									# while lineage_ref[anno[item]['start']+3*(temp-1)-1]=="-":
									# 	temp+=1
									# if temp == idlist[j][0]:
									# 	idlist.remove(idlist[i])
									# 	mutlist.remove(mutlist[i])
									# 	i-=1
									# 	j-=1
									# 	break
								j+=1
							i+=1						
						for i in range(len(idlist)):
							idlist[i]=idlist[i][0]
						if 1 in idlist:
							important_mut=True
						newidlist=sorted(idlist)
						newmutlist=[]						
						for i in range(len(idlist)):
							for j in range(len(idlist)):
								if newidlist[i]==idlist[j]:
									if "O" in mutlist[j]:
										important_mut==True
									newmutlist.append(mutlist[j].replace("O","*"))													
						imp_mut_dict[item]=",".join(newmutlist)
						imp_mut_list.append(item+":"+",".join(newmutlist))
					exmut_list=[]
					for item in list(exmut_dict):
						#if exmut_dict[item]>=count/2:
							exmut_list.append(item.replace(":","")+":"+str(exmut_dict[item]))
					if important_mut:
						print(node['name'],lineage,','.join(current_mut),count,len(country_list),max(date_list),';'.join(imp_mut_list),','.join(exmut_list),n_glycan)
						w=highlight_browser(node)
	return [lineage,count,copy.deepcopy(country_list),copy.deepcopy(date_list),copy.deepcopy(exmut_dict)]

variant_mutation_dic={}
def read_alias():
	w=open("alias_key.json",'r')
	q=json.load(w)
	w.close()
	rev=dict()
	for item in list(q):
		if item !="A" and item !="B" and item[0] != "X":
			rev.update({q[item]:item})
	alias=dict()	
	for item in list(q):
		if item !="A" and item !="B" and item[0] != "X":
			temp=q[item].rsplit(".",3)
			if temp[0] !="A" and temp[0] !="B" and temp[0][0] != "X":
				alias.update({item:rev[temp[0]]+"."+".".join(temp[1:])})
			else:
				alias.update({item:q[item]})
		elif item[0] == "X":
			temp1=q[item][0].removesuffix("*")+"."
			temp2=q[item][1].removesuffix("*")+"."
			while temp1.split(".",1)[0] != "A" and temp1.split(".",1)[0] != "B":
				temp1=alias[temp1.split(".",1)[0]]+"."+temp1.split(".",1)[1]
			while temp2.split(".",1)[0] != "A" and temp2.split(".",1)[0] != "B":
				temp2=alias[temp2.split(".",1)[0]]+"."+temp2.split(".",1)[1]
			temp=""
			i=0
			while temp1.split(".")[i]==temp2.split(".")[i]:
				temp+=("."+temp1.split(".")[i])
				i+=1
				if i == len(temp1) or i == len(temp2):
					break
			temp=temp[1:]
			while len(temp.split("."))%3!=0:
				temp+=("."+item)
			temp+=("."+item)
			rev.update({temp:item})
			if len(temp.split("."))>4:
				temp=rev[temp.rsplit(".",3)[0]]+"."+".".join(temp.rsplit(".",3)[1:])
			alias.update({item:temp})
	for item in list(alias):
		temp=alias[item].split(".")
		while len(temp)>1 and temp[-1][0]=="X":
			temp=temp[:-1]
		if len(temp) ==1:
			if temp[0] != "A" and temp[0] != "B":
				temp=alias[temp[0]]
			else:
				temp=temp[0]
		else:
			temp=".".join(temp)
		alias.update({item:temp})
	return alias

alias=read_alias()

def read_insertion():
	insref=open("insertion.txt",'r')
	flines=insref.readlines()
	ins_dic={}
	for l in flines:
		linsp=l.split()
		lineage=linsp[0]
		insertions=[]
		for i in range(1, len(linsp),2):
			site=linsp[i]
			inserted=linsp[i+1]
			insertions.append("ins"+site+":"+inserted)
		ins_dic[lineage]=",".join(insertions)
	return ins_dic
	
ins_dic=read_insertion()
#print(ins_dic)

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
	if name.split(".",1)[0] in alias:
		lineage=name
		while lineage !="A" and lineage !="B":
			if lineage in ins_dic:
				for item in ins_dic[lineage].split(","):
					if item not in all_mutations:
						all_mutations.append(item)
			if "." in lineage:
				lineage=lineage.rsplit(".",1)[0]
			elif lineage in list(alias):
				lineage = alias[lineage]
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

#print(variant_mutation_dic)
#print(len(variant_mutation_dic))

def recomb_helper(lineage,current_mut):
	temp=[]
	designated_mutations=variant_mutation_dic[lineage]
	lineage_ref=ref
	for item in designated_mutations:
		if item[0]!="i":
			ids=int(item[1:-1])
			lineage_ref=lineage_ref[:ids-1]+item[-1]+lineage_ref[ids:]
	this_seq=lineage_ref
	for item in current_mut:
		ids=int(item[1:-1])
		this_seq=this_seq[:ids-1]+item[-1]+this_seq[ids:]
	for lin in list(variant_mutation_dic):
		temp_lin=lin
		while lineage != temp_lin:
			if "." in temp_lin:
				temp_lin=temp_lin.rsplit(".",1)[0]
			elif temp_lin in list(alias):
				temp_lin = alias[temp_lin]
			else:
				break
		else:
			designated_mutations=variant_mutation_dic[lin]
			lineage_ref=ref
			for item in designated_mutations:
				if item[0]!="i":
					ids=int(item[1:-1])
					lineage_ref=lineage_ref[:ids-1]+item[-1]+lineage_ref[ids:]
			current_mut = []
			for ids in range(len(this_seq)):
				if this_seq[ids] != lineage_ref[ids] and this_seq[ids]!= "-" and lineage_ref[ids]!= "-":
					current_mut.append(lineage_ref[ids]+str(ids+1)+this_seq[ids])
			backcount=0
			for item in current_mut:
				ids = int(item[1:-1])
				if ref[ids-1]==item[-1]:
					backcount+=1
			temp.append([lin,current_mut,backcount])		
	while True:
		current_mut = []
		for ids in range(len(this_seq)):
			if this_seq[ids] != lineage_ref[ids] and this_seq[ids]!= "-":
				current_mut.append(lineage_ref[ids]+str(ids+1)+this_seq[ids])
		backcount=0
		for item in current_mut:
			ids = int(item[1:-1])
			if ref[ids-1]==item[-1]:
				backcount+=1
		temp.append([lineage,current_mut,backcount])
		if lineage == "A" or lineage == "B" or backcount == 0:
			break
		if "." in lineage:
				lineage=lineage.rsplit(".",1)[0]
		elif lineage[0] == "X":
			lineage = alias[lineage]
		else:
			lineage = alias[lineage].rsplit(".",1)[0]
		while lineage not in variant_mutation_dic:
			if "." in lineage:
				lineage=lineage.rsplit(".",1)[0]
			elif lineage in list(alias):
				lineage = alias[lineage]
			else:
				break
		else:
			designated_mutations=variant_mutation_dic[lineage]
		lineage_ref=ref
		for item in designated_mutations:
			if item[0]!="i":
				ids=int(item[1:-1])
				lineage_ref=lineage_ref[:ids-1]+item[-1]+lineage_ref[ids:]
	mutcount=len(temp[0][1])+0*temp[0][2]
	for item in temp:
		if len(item[1])+0*item[2]<mutcount:
			mutcount = len(item[1])+0*item[2]
	for item in temp:
		if len(item[1])+0*item[2]==mutcount:
			return item
f=open(filename+".json",'r')
js=json.load(f)
f.close()
if os.path.exists(filename+".tsv"):
	tsv=open(filename+".tsv",'r')
	tsvlines=tsv.readlines()
	tsv.close()
else:
	tsvlines=""
#print(read_tsv())
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

tree=js['tree']


root=tree['children'][0]

l,count,countries,dates,exmuts=node_browser(root,'',ref,[],0)
outjs=open(filename+"-out.json",'w')
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
