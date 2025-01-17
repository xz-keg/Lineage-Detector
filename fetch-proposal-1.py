import json
import os
import argparse
import re
import urllib.request
import urllib3
import markdown
import html2text
from bs4 import BeautifulSoup

# parser = argparse.ArgumentParser(description='Demo of argparse')
# parser.add_argument('--issue', type=str,default='usher')
# args = parser.parse_args()
# issue=args.issue
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
table=read_table()
variant_mutation_dic={}

w=open("alias_key.json",'r')
alias=list(json.load(w))
w.close()

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


import copy
ins_dic=read_insertion()
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
	if (name.split(".",1)[0] in (["A", "B"]+list(alias))) and ("_" not in name) :
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

anno=dict()
anno['ORF0']={'start':107,'end':136}
anno['ORF1a']={"start": 266, "end": 13468}
anno['ORF1b']={"start": 13468, "end": 21555}
anno['S']={"start": 21563, "end": 25384}
anno['ORF3a']={"start": 25393, "end": 26220}
anno['ORF3b']={'start':25814,'end':25882}
anno['ORF3c']={'start':25457,'end':25582}
anno['ORF3d']={'start':25524,'end':25697}
anno['ORF3d-2']={'start':25968,'end':26069}
anno['E']={"start": 26245, "end": 26472}
anno['M']={"start": 26523, "end": 27191}
anno['ORF6']={"start": 27202, "end": 27387}
anno['ORF7a']={"start": 27394, "end": 27759}
anno['ORF7b']={"start": 27756, "end": 27887}
anno['ORF8']={"start": 27894, "end": 28259}
anno['N']={"start": 28274, "end": 29533}
anno['ORF9b']={'start':28284,'end':28577}
anno['ORF9c']={'start':28734,'end':28955}
anno['ORF10']={"start": 29558, "end": 29674}


sc2="https://api.github.com/repos/sars-cov-2-variants/lineage-proposals/issues"
pango="https://api.github.com/repos/cov-lineages/pango-designation/issues"

query_re=re.compile(r"((?:(?:\-\ *)?(?:ins[0-9]*[ACGTacgt]*|del[0-9_]*|[ACGTacgt]?[0-9]+[ACGTacgt]?)(?:\ *,\ *))+(?:(?:\-\ *)?(?:ins[0-9]*[ACGTacgt]*|del[0-9_]*|[ACGTacgt]?[0-9]+[ACGTacgt]?)))(?![0-9]* seq)")
#path_re=re.compile(r"((?:[A-HJ-NP-WYZ]+(?:\.[0-9]+){1,3}|X[A-HJ-NP-WYZ]+(?:\.[0-9]+){0,3})(?:\ *(?:[>,+]|\>\ ?\>)\ *(?:ins[0-9]+[ACGTacgt]+(?:rev)?|del[0-9_-]+|[ACGTacgt][0-9]+[ACGTacgt])(?:\ *\((?:(?:[,=]\ *)?(?:(?:S|N|M|E|[Oo][Rr][Ff][0-9]+[a-z](?:\-[0-9]+)?|NSP[0-9]+)[:_]\ *)?(?:fs|frameshift|(?:[A-Z*]\>)?[A-Z*][0-9]+(?:[A-Z*-]|del|fs|frameshift)(?:rev)?|ins[0-9]+[A-Z*]+|del[0-9]+)+\ *)+\))?+)+)")
path_re=re.compile(r"((?:[A-HJ-NP-WYZa-hj-np-wyz]+(?:\.[0-9]+){1,3}|[Xx][A-HJ-NP-WYZa-hj-np-wyz]+(?:\.[0-9]+){0,3})(?:\[.*\]|\(.*\))?(?:\ *(?:(?:[>,+→]|\-\>|\>\ ?\>)\ *)?(?:(?:ins[0-9]+[ACGTacgt]+(?:\ *(?:r|rev))?|(?:del|Δ)[0-9_-]+|[ACGTacgt]?[0-9]+[ACGTacgt])(?:\ *\((?:(?:[,=]\ *)?(?:(?:S|N|M|E|[Oo][Rr][Ff][0-9]+[a-z]?(?:\-[0-9]+)?|NSP[0-9]+|PLPRO)[:_]\ *)?(?:fs|frameshift|(?:[A-Z*]\>)?[A-Z*]?[0-9]+(?:[A-Z*-]|del|Δ|fs|frameshift)(?:\ *(?:r|rev))?|ins[0-9]+[A-Z*]+|(?:del|Δ)[0-9]+)+\ *)+\))?+|(?:(?:S|N|M|E|[Oo][Rr][Ff][0-9]+[a-z]?(?:\-[0-9]+)?|NSP[0-9]+|PLPRO)[:_]\ *)?(?:fs|frameshift|(?:[A-Z*]\>)?[A-Z*]?[0-9]+(?:[A-Z*-]|del|Δ|fs|frameshift)(?:\ *(?:r|rev))?|ins[0-9]+[A-Z*]+|(?:del|Δ)[0-9]+)(?:\ *\((?:(?:[,=]\ *)?(?:(?:r|rev)\ *)?(?:ins[0-9]+[ACGTacgt]+(?:\ *(?:r|rev))?|(?:del|Δ)[0-9_-]+|[ACGTacgt]?[0-9]+[ACGTacgt])\ *)+\))?+))+)")
seqid_re=re.compile(r"((?:EPI_ISL_|(?:C_)?[A-Z]{2})[0-9]+(?:\.[0-9]+)?)")
date_re=re.compile(r"((?:[0-9]{2,4}\-)?[0-9]{1,2}\-[0-9]{1,2}|[0-9]{4}\-[0-9]{1,2})")


f=open("token.txt","r")
token=f.read()
f.close
http=urllib3.PoolManager(timeout=urllib3.Timeout(connect=10,read=10))
MyHeader={"User-Agent":"Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/131.0.0.0 Safari/537.36 Edg/131.0.0.0","Authorization":"Bearer "+token}

def markdown_to_text(markdown_text):
    # 将 Markdown 转换为 HTML
    html = markdown.markdown(markdown_text)#.replace("\n","")
    soup=BeautifulSoup(html, "html.parser")
    plain_text = "".join(soup.find_all(text=True))
    return plain_text.strip()

def query_normalize(q):
	mut=q.split(",")
	ids=[0 for i in range(len(mut))]
	for i in range(len(mut)):
		mut[i]=mut[i].strip().rstrip().replace("a","A").replace("c","C").replace("g","G").replace("t","T")
		temp=re.search(r"[0-9]+",mut[i])
		if temp == None:
			return q
		ids[i]=int(temp[0])
		if mut[i][0]=="-":
			ids[i]+=100000
	tempdict=dict(zip(ids,mut))
	sortids=sorted(ids)
	sortmut=["" for i in range(len(mut))]
	for i in range(len(mut)):
		sortmut[i]=tempdict[sortids[i]]
	return ", ".join(sortmut)

def path_normalize(p):
	lineage=re.findall(r"([A-HJ-NP-WYZa-hj-np-wyz]+(?:\.[0-9]+){1,3}|[Xx][A-HJ-NP-WYZa-hj-np-wyz]+(?:\.[0-9]+){0,3})",p)[0].upper()	
	lineage_ref=ref
	if lineage in variant_mutation_dic:
		designated_mutations=variant_mutation_dic[lineage]
		for item in designated_mutations:
			if item[0]!= "i":
				ids=int(item[1:-1])
				lineage_ref=lineage_ref[:ids-1]+item[-1]+lineage_ref[ids:]
	else:
		print("ERROR:Wrong Lineage")	
	spl=re.split(r"\ *(?:[>+→]|\-\>|\>\ ?\>)\ *",p)
	parts=[]
	for part in spl:
		muts=[]
		ids=[]
		for mut in re.findall(r"(?![:]\ *)([ACGTacgt][0-9]+[ACGTacgt]|ins[0-9]+[ACGTacgt]+|(?:del|Δ)[0-9_-]+)",part):
			if ("del" in mut) or ("Δ" in mut):
				id=re.findall(r"[0-9]+",mut)
				if len(id)>1:
					if len(id[1])<len(id[0]):
						id[1]=id[0][:len(id[0])-len(id[1])]+id[1]
					muts.append("del"+id[0]+"-"+id[1])
				else:
					muts.append("del"+id[0])
				ids.append(float(id[0]))
			elif "ins" in mut:
				id=re.findall(r"[0-9]+",mut)[0]
				insertion=re.findall(r"[ACGTacgt]+",mut)[0].upper()
				ids.append(float(id)+0.5)
				muts.append("ins"+id+insertion)
			else:
				if mut[0] in "ACGTacgt":
					id=mut[1:-1]
				else:
					id=mut[:-1]
					mut=lineage_ref[int(id)-1]+mut
				muts.append(mut.upper())
				ids.append(float(id))
		if len(muts)>0:
			tempdict=dict(zip(ids,muts))
			sortids=sorted(ids)
			sortmut=["" for i in range(len(muts))]
			for i in range(len(muts)):
				sortmut[i]=tempdict[sortids[i]]
			parts.append(", ".join(sortmut))
	if len(parts)>0:
		return lineage+"> "+"> ".join(parts)

def date_normalize(d,proposingtime):
	pt=proposingtime.split("T")[0].split("-")
	spl=d.split("-")
	if len(spl)==3:
		if len(spl[0])==2:
			spl[0]="20"+spl[0]
		spl[1]=spl[1].zfill(2)
		spl[2]=spl[2].zfill(2)
	elif len(spl[0])==4 or int(spl[0])>=12:
		if len(spl[0])==2:
			spl[0]="20"+spl[0]
		spl[1]=spl[1].zfill(2)
		spl.append["00"]
	else:
		spl.append(spl[1])
		spl[1]=spl[0]
		if int(spl[1])*32+int(spl[2])<int(pt[1])*32+int(pt[2]):
			spl[0]=pt[0]
		else:
			spl[0]=str(int(pt[0])-1)
		spl[1]=spl[1].zfill(2)
		spl[2]=spl[2].zfill(2)
	return "-".join(spl)


def findloc(l):
	spl=re.split(r"(?:EPI_ISL_|(?:C_)?[A-Z]{2})[0-9]+(?:\.[0-9]+)?|(?:[0-9]{2,4}\-)?[0-9]{1,2}\-[0-9]{1,2}|[0-9]{4}\-[0-9]{1,2}",l)
	loc=""
	for x in spl:
		if len(x)>0:
			if not re.fullmatch(r".+\:\ *",x):
				loc+=x.replace("|","").strip(", ()").rstrip(", ()")
	if not (("from" in loc) or ("import" in loc) or re.search("[A-Z ]via", loc)):
		if "Puerto Rico" in loc:
			return "PR"
		elif "Guam" in loc:
			return "GU"
		elif "American Samoa" in loc:
			return "AS"
		elif ("United States" in loc) or ("US" in loc) or (not ("th America" in loc or "Americas" in loc) and "America" in loc):
			return "US"
		elif ("Gerogia") in loc:
			return "GE"
		if ("HK" in loc) or ("Hong Kong" in loc):
			return "HK"
		elif ("Macao" in loc) or ("Macau" in loc):
			return "MO"
		elif "Taiwan" in loc:
			return "TW"
		elif ("Namsai" in loc) or ("Changlang" in loc) or ("Tirap" in loc) or ("Longding" in loc):
			return "IN"
		elif ("CN" in loc) or ("China" in loc) or ("PRC" in loc) or ("Arunachal" in loc) or ("Anhui" in loc) or ("Beijing" in loc) or ("Peking" in loc) or ("Chongqing" in loc) or ("Chungching" in loc) or ("Fujian" in loc) or ("Guangdong" in loc) or ("Canton" in loc) or ("Gansu" in loc) or ("Guangxi" in loc) or ("Guizhou" in loc) or ("Henan" in loc) or ("Hubei" in loc) or ("Hainan" in loc) or ("Heilongjiang" in loc) or ("Hunan" in loc) or ("Jilin" in loc) or ("Jiangsu" in loc) or ("Jiangxi" in loc) or ("Liaoning" in loc) or ("Neimenggu" in loc) or ("Inner Mongolia" in loc) or ("Ningxia" in loc) or ("Qinghai" in loc) or ("Sichuan" in loc) or ("Szechuan" in loc) or ("Shandong" in loc) or ("Shanghai" in loc) or ("Shaanxi" in loc) or ("Shanxi" in loc) or ("Tianjin" in loc) or ("Tientsin" in loc) or ("Xinjiang" in loc) or ("Xizang" in loc) or ("Tibet" in loc) or ("Yunnan" in loc) or ("Zhejiang" in loc):
			return "CN"
		if ("North Korea" in loc) or ("DPRK" in loc):
			return "KP"
		elif "Korea" in loc:
			return "KR"
		if ("Ukraine" in loc) or ("Ukraina" in loc) or ("Crimea" in loc) or ("Krym" in loc) or ("Sevastopol" in loc):
			return "UA"
		elif "Russia" in loc:
			return "RU"
		if "Ascension" in loc:
			return "AC"
		elif "Anguilla" in loc:
			return "AI"
		elif "Bermuda" in loc:
			return "BM"
		elif ("Britain" in loc) or ("GB" in loc) or ("UK" in loc) or ("United Kingdom" in loc) or ("England" in loc) or ("Scotland" in loc) or ("Wales" in loc) or ("Northern Ireland" in loc):
			return "GB"
		elif "Ireland" in loc:
			return "IE"
		if ("Saint Barthélemy" in loc) or ("St. Bart" in loc):
			return "BL"
		elif "France" in loc:
			return "FR"
		if "Åland" in loc:
			return "AX"
		elif ("Finland" in loc) or ("Suomi" in loc):
			return "FI"
		if "Aruba" in loc:
			return "AW"
		elif ("Netherlands" in loc) or ("Nederland" in loc) or ("Holland" in loc):
			return "NL"
		if ("Ceuta" in loc) or ("Melilla" in loc):
			return "EA"
		elif ("Spain" in loc) or ("España" in loc):
			return "ES"
		elif ("Morocco" in loc):
			return "MA"
		if "Andorra" in loc:
			return "AD"
		if ("UAE" in loc) or ("United Arab Emirates" in loc):
			return "AE"
		if "Afghanistan" in loc:
			return "AF"
		if "Antigua and Barbuda" in loc:
			return "AG"
		if ("Albania" in loc) or ("Shqipëri" in loc):
			return "AL"
		if ("Armenia" in loc) or ("Hayastan" in loc):
			return "AM"
		if "Angola" in loc:
			return "AO"
		if "Argentina" in loc:
			return "AR"
		if ("Austria" in loc) or ("Österreich" in loc):
			return "AT"
		if ("Austalia" in loc) or ("AU" in loc):
			return "AU"
		elif ("SA" in loc):
			return "ZA"
		if "Azerbaijan" in loc:
			return "AZ"
		if ("Bosnia" in loc) or ("Herzegovina" in loc) or ("Srpska" in loc):
			return "BA"
		if "Barbados" in loc:
			return "BB"
		if "Bangladesh" in loc:
			return "BD"
		if "Belg" in loc:
			return "BE"
		if "Burkina Faso" in loc:
			return "BF"
		if "Bulgaria" in loc:
			return "BG"
		if "Bahrain" in loc:
			return "BH"
		if "Burundi" in loc:
			return "BI"
		if "Benin" in loc:
			return "BJ"
		if "Brunei" in loc:
			return "BN"
		if "Bolivia" in loc:
			return "BO"
		if ("Brazil" in loc) or ("BR" in loc):
			return "BR"
		if "Bahamas" in loc:
			return "BS"
		if "Bhutan" in loc:
			return "BT"
		if "BV" in loc:
			return "BO"
	return loc

def myJoin(sep,stringlist):
	if len(stringlist)>0:
		return sep.join(stringlist)
	else:
		return ""

def fetch(url):
	data=json.loads(http.request("GET",url, headers=MyHeader).data)
	shortaddr=url.removeprefix("https://api.github.com/repos/").replace("/issues/","#")
	print(shortaddr)
	title=data["title"]
	print(title)
	proposer=data["user"]["login"]	
	print(proposer)	
	proposingtime=data["created_at"]	
	print(proposingtime)
	statue=data["state"]
	if data["state_reason"]:
		statue+=": "+data["state_reason"]
	print(statue)
	labels=[]
	for label in data["labels"]:
		labels.append(label["name"])
	print(myJoin("; ",labels))
	queries=[]
	paths=[]
	source=""
	firstseqsid=[]
	firstseqstime=[]
	firstseqsloc=[]	
	if not "pull_request" in data:		
		rec = "recomb" in title.lower()
		text=markdown_to_text(data["body"])
		#print(markdown.markdown(text))
		#print(text)
		spl=text.splitlines()
		for i in range(len(spl)):
			#print(spl[i])
			z=spl[i]
			#print(z)
			if len(z)>0:
				if (z.split(".")[0] in alias) or (z[0]=="X"):						
					for path in path_re.findall(z):
						print(path)
						if path_normalize(path):
							print(path_normalize(path))
							paths.append(path_normalize(path))
					# spl1=z.split("(")
					# spl2=[]
					# for j in range(len(spl1)):
					# 	if ")" in spl1[j]:						
					# 		spl2.append(spl1[j].rsplit(")",1)[0])
					# 		spl1[j]=spl1[j].rsplit(")",1)[-1]
					# temp=(",".join(spl1+spl2)).replace(" ","")
					# temp=temp.replace(">",",").replace("+",",").split(",")
					# print(temp)
			if ("quer" in z.lower()) or ("gisaid:" in z.lower()):
				#tempq=re.findall(r"((?:(?:\-\ *)?(?:ins|del)?[ACGTacgt0-9_][ACGTacgt0-9 _]*,)(?:\ *\-?(?:ins|del)?[ACGTacgt0-9 _]+,)*(?:\ *\-?(?:ins|del)?[ACGTacgt0-9 _]*[ACGTacgt0-9_]))(?![0-9]* seq)",z)
				tempq=query_re.findall(z.replace("&",",").replace("!","-"))
				for query in tempq:
					if query_normalize(query) not in queries:
						print(query_normalize(query))
						queries.append(query_normalize(query))
			elif query_re.fullmatch(z.replace("&",",").replace("!","-")):
				if query_normalize(z.replace("&",",").replace("!","-")) not in queries:
					print(query_normalize(z.replace("&",",").replace("!","-")))
					queries.append(query_normalize(z.replace("&",",").replace("!","-")))
			#if path_re.fullmatch(z):
				#print(z)
			if ("first" in z.lower()) or ("earliest" in z.lower()):
				if seqid_re.search(z):
					print(z)
					print(seqid_re.findall(z)[0])
					firstseqsid.append(seqid_re.findall(z)[0])
					if date_re.search(z):
						print(date_normalize(date_re.findall(z)[0],proposingtime))
						firstseqstime.append(date_normalize(date_re.findall(z)[0],proposingtime))
					print(findloc(z))
					firstseqsloc.append(findloc(z))
			if ("mutation" in z.lower()):
				if not "cov-spectrum" in z.lower():
					if path_re.search(z):
						for path in path_re.findall(z):								
							print(path)
							if path_normalize(path):
								print(path_normalize(path))
								paths.append(path_normalize(path))
					else:
						print(z)
						if z[-1]==":" and i < len(spl)-1:
							if path_re.search(spl[i+1]):
								for path in path_re.findall(z):								
									print(path)
									if path_normalize(path):
										print(path_normalize(path))
										paths.append(path_normalize(path))
							else:
								spl1=spl[i+1].split("(")
								spl2=[]
								for j in range(len(spl1)):
									if ")" in spl1[j]:						
										spl2.append(spl1[j].rsplit(")",1)[0])
										spl1[j]=spl1[j].rsplit(")",1)[-1]
								temp=(",".join(spl1+spl2)).replace(" ","")
								temp=temp.replace(">",",").replace("+",",").split(",")
								print(temp)
			if (pango in url) and ("sars-cov-2-variants/lineage-proposals" in z):
				source=re.findall(r"(sars-cov\-2\-variants/lineage\-proposals(?:/issues/|\#)[0-9]+)",z)[0].replace("/issues/","#")
				print(source)
			if rec:
				if ("breakpoint" in z.lower()) or ("brpt" in z.lower()):
					print(z)
		if len(queries)==0:
			comments=json.loads(http.request("GET",url+"/comments", headers=MyHeader).data)
			for x in comments:
				text=x["body"]
				spl=text.splitlines()
				for z in spl:
					if ("quer" in z.lower()) or ("gisaid:" in z.lower()):
						tempq=re.findall(r"((?:(?:\-\ *)?(?:ins[0-9]*[ACGTacgt]*|del[0-9_]*|[ACGTacgt]?[0-9]+[ACGTacgt]?)(?:\ *,\ *))+(?:(?:\-\ *)?(?:ins[0-9]*[ACGTacgt]*|del[0-9_]*|[ACGTacgt]?[0-9]+[ACGTacgt]?)))(?![0-9]* seq)",z.replace("&",",").replace("!","-"))
						for query in tempq:
							if query_normalize(query) not in queries:
								print(query_normalize(query))
								queries.append(query_normalize(query))
					elif re.fullmatch(r"((?:(?:\-\ *)?(?:ins[0-9]*[ACGTacgt]*|del[0-9_]*|[ACGTacgt]?[0-9]+[ACGTacgt]?)(?:\ *,\ *))+(?:(?:\-\ *)?(?:ins[0-9]*[ACGTacgt]*|del[0-9_]*|[ACGTacgt]?[0-9]+[ACGTacgt]?)))(?![0-9]* seq)",z.replace("&",",").replace("!","-")):
						if query_normalize(z.replace("&",",").replace("!","-")) not in queries:
							print(query_normalize(z.replace("&",",").replace("!","-")))
							queries.append(query_normalize(z.replace("&",",").replace("!","-")))
	else:
		print("Pull Request")
		source="Pull Request"
	print()
	return "\t".join([shortaddr,title,statue,myJoin("; ",labels),proposer,proposingtime,source,myJoin("; ",queries),myJoin("; ",paths),myJoin("; ",firstseqsid),myJoin("; ",firstseqstime),myJoin("; ",firstseqsloc)])+"\n"	

outputtext=""
import time
page=json.loads(http.request("GET",sc2, headers=MyHeader).data)
maxissue=int(page[0]["number"])
for issue in range(2100,maxissue+1):
	time.sleep(1)
	outputtext+=fetch(sc2+"/"+str(issue))

page=json.loads(http.request("GET",pango, headers=MyHeader).data)
maxissue=int(page[0]["number"])
for issue in range(2790,maxissue+1):	
	time.sleep(1)
	outputtext+=fetch(pango+"/"+str(issue))

outputfile=open("fetch.tsv","w", encoding="utf-8")
outputfile.write(outputtext)
outputfile.close()