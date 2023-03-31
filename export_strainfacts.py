import os
import lz4.frame
import io
import sys
import csv
import gzip
from Bio import SeqIO
import pickle

## Usage:
## export_sf.py [SAMPLE_DIR] [MERGE_DIR] [CANDIDATE_SPECIES] [DB_DIR]

sample_dir=sys.argv[1]
merge_dir=sys.argv[2]
species=sys.argv[3]
db_dir=sys.argv[4]

samples=os.listdir(sample_dir)

## Make global position
## Calculate length of contigs
## load gene annotations fna and save the length
prepend_dic={}
db_species = db_dir+"/gene_annotations/"+species
db_species_list = os.listdir(db_species)
k=0
seqlen=[]
for db in db_species_list:
	fa = db_species+"/"+db+"/"+db+".fna"
	with open(fa, "rt") as handle:
		for record in SeqIO.parse(handle, "fasta"):
			id=record.id.split("|")[2]
			if k==0:
				prepend_dic[id]=0
			else:
				prepend_dic[id]=sum(seqlen)
			seqlen.append(len(record.seq)+1)
			k=k+1

## Load merge dir and preserve bi-allelic positions across samples
merge_path=merge_dir+"/snps/"+species+"/"+species+".snps_info.tsv.lz4"
file = lz4.frame.open(merge_path).read().decode('utf-8')
df=csv.DictReader(io.StringIO(file), delimiter='\t')
dic={}
reference="major"
for row in df:
	if row["snp_type"]=="bi":
		sid=row["site_id"]
		contig=sid.split("|")[2]
		pos=sid.split("|")[3]
		if reference=="ref":
			ref=sid.split("|")[4]
		elif reference=="major":
			ref=row["major_allele"]
		else:
			sys.exit("reference must be ref or major")
		min=row["minor_allele"]
		dic[contig+"|"+pos]=[ref,min]


## Recursively search per-sample SNV for depth
export=[]
export.append("sample\tposition\tallele\tmetagenotype\n")
snum=0
snumdic={}
for sample in samples:
	sample_dic={}
	print(sample)
	snumdic[sample]=snum
	if sample!="merge":
		sample_list=os.listdir(sample_dir+"/"+sample+"/snps")
		species_list=[j.split(".")[0] for j in sample_list]
		if species in species_list:
			path=sample_dir+"/"+sample+"/snps/"+species+".snps.tsv.lz4"
			file = lz4.frame.open(path).read().decode('utf-8')
			df=csv.DictReader(io.StringIO(file), delimiter='\t')

			for row in df:
				tmp_key=row["ref_id"].split("|")[2]+"|"+row["ref_pos"]
				sample_dic[tmp_key]=row

			for key in dic.keys():
				if key in sample_dic.keys():
					ref=dic[key][0]
					alt=dic[key][1]					
					ref_count=float(sample_dic[key]["count_"+ref.lower()])
					alt_count=float(sample_dic[key]["count_"+alt.lower()])
					export.append(str(snum)+"\t"+str(int(sample_dic[key]["ref_pos"])+prepend_dic[sample_dic[key]["ref_id"].split("|")[2]])+"\t"+"alt"+"\t"+str(alt_count)+"\n")
					export.append(str(snum)+"\t"+str(int(sample_dic[key]["ref_pos"])+prepend_dic[sample_dic[key]["ref_id"].split("|")[2]])+"\t"+"ref"+"\t"+str(ref_count)+"\n")
				else:
					export.append(str(snum)+"\t"+str(int(key.split("|")[1])+prepend_dic[key.split("|")[0]])+"\talt\t0.0\n")
					export.append(str(snum)+"\t"+str(int(key.split("|")[1])+prepend_dic[key.split("|")[0]])+"\tref\t0.0\n")
	snum=snum+1

with open(species+"_metagenotype.tsv", "w") as f:
	f.writelines(export)
with open(species+"_sample.pickle","wb") as f:
	pickle.dump(snumdic,f)
