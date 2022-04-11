import pandas as pd
import os
import glob
from sys import exit


######  reading and filtering sample data ######
### read in data, set index to lib
df = pd.read_csv(config["samples"], sep="\t").set_index("lib", drop=False)
## print("####\nraw df")
## print(df)

#reduce sample list on only those listed via '--config select=""'
if "select" in config:
	df = df[df["sample"].isin(config["select"].split(","))]
##	print("###\nselect df")
##	print(df)
	if df.empty:
		print("Error: Specified id(s) not in "+config["samples"], file=sys.stderr)
		sys.exit()

### create non redundant list of sample IDS
all_samples = list(set(df.set_index(['sample']).index.values.tolist()))
## print (all_samples)

#####   partition data by input data type   #####
### create dataframe containing just the samples that have data in the 'f_fastq' column
df_fastq = df[df['f_fastq'].notna()]
## print(df_fastq["lib"].tolist())

### create dataframe containing just the samples that have data in the 'bam' column
df_bam = df[df['bam'].notna()]
print("df_bam: \n")
print(df_bam)
## print("|".join(df_bam["lib"].tolist()))

### make dataframe with only the samples that contain Illumina data, i.e. data in columnis f_fastq r_fastq or bam
Illumina_process_df = df[df["lib"].isin(df_fastq["lib"].tolist()+df_bam["lib"].tolist())]
## print("\n####THESE ARE THE SAMPLES")
## print(Illumina_process_df["sample"].tolist())
## print(Illumina_process_df["lib"].tolist())




#fill dictionary
dic = {'sample': [], 'lib': []}
unitdict = {}
for lib in set(df.index.values.tolist()):
    sample = str(df.loc[lib, ["sample"]].values[0])
    dic["sample"].append(str(sample))
    dic["lib"].append(str(lib))
    if not sample in unitdict:
        unitdict[sample] = [str(lib)]
    else:
        unitdict[sample].append(str(lib))

print("UNITDICT: "+str(unitdict))
units = pd.DataFrame(dic).set_index(['sample','lib'], drop=False)
print(str(units))
units.index = units.index.set_levels(
    [i.astype(str) for i in units.index.levels])  # enforce str in index

#def get_sample_for_lib(wildcards):
#	return pd.read_csv(config["samples"], dtype=str, sep="\t").set_index(["lib"], drop=False).loc[(wildcards.lib), ["sample"]].dropna()


#Illumina_process_df = df[df["lib"].isin(df_fastq["lib"].tolist()+df_bam["lib"].tolist())]
#print(Illumina_process_df.index)
#Illumina_dict = {}
#for k in unitdict.keys():
#	for i in range(len(unitdict[k])):
#		if unitdict[k][i] in Illumina_process_df["lib"]:
#			print(k,unitdict[k][i])
#			if not k in Illumina_dict:
#				Illumina_dict[k] = [str(unitdict[k][i])]
#			else:
#				Illumina_dict[k].append(str(unitdict[k][i]))

#print(Illumina_dict)

illumina_units = units[units["lib"].isin(df_fastq["lib"].tolist()+df_bam["lib"].tolist())]
print("Illumina_units:")
print(illumina_units)
print()

#long reads
#df = pd.read_csv(config["samples"], sep="\t").set_index("lib", drop=False)
#only keep rows in dataframe that do not have na in column fast5_dir
df_fast5 = df[df['fast5_dir'].notna()]
fast5_units = units[units["lib"].isin(df_fast5["lib"].tolist())]
#print("\ndf_fast5")
#print(df_fast5["sample"].tolist())
#print("\nunits - sample")
#print(units["sample"])
#print("\nfast5units")
#print(fast5_units.itertuples())

#only keep rows in dataframe that do not have na in column long
df_long = df[df['long'].notna()]

def input_for_trimgalore_f(wildcards):
	print("Input for trimming: "+str(wildcards.lib))
	if wildcards.lib in df_fastq.index:
		return "results/{sample}/Illumina/raw_reads/from_fastq/{lib}/{sample}.{lib}.raw.1.fastq.gz".format(lib=wildcards["lib"], sample=wildcards["sample"])
	elif wildcards.lib in df_bam.index:
		return "results/{sample}/Illumina/raw_reads/from_bam/{lib}/{sample}.{lib}.raw.1.fastq.gz".format(lib=wildcards["lib"], sample=wildcards["sample"])

def input_for_trimgalore_r(wildcards):
	if wildcards.lib in df_fastq.index:
		return "results/{sample}/Illumina/raw_reads/from_fastq/{lib}/{sample}.{lib}.raw.2.fastq.gz".format(lib=wildcards["lib"], sample=wildcards["sample"])
	elif wildcards.lib in df_bam.index:
		return "results/{sample}/Illumina/raw_reads/from_bam/{lib}/{sample}.{lib}.raw.2.fastq.gz".format(lib=wildcards["lib"], sample=wildcards["sample"])
		
def input_for_clean_trim_fp(wildcards):
	lis = []
#	print(wildcards.sample)
#	print(Illumina_process_df.set_index("sample").index)
#	print(Illumina_process_df.set_index("sample").loc[(wildcards.sample), ["lib"]].dropna())
	for l in Illumina_process_df.set_index("sample").loc[(wildcards.sample), ["lib"]].values:
		if type(l) is str:
			lis.append("results/{sample}/trimming/trimgalore/{lib}/{sample}.{lib}.1.fastq.gz".format(sample=wildcards["sample"], lib=l))
		else:
			lis.append("results/{sample}/trimming/trimgalore/{lib}/{sample}.{lib}.1.fastq.gz".format(sample=wildcards["sample"], lib=l[0]))
#	print(lis)
	return lis
def input_for_clean_trim_rp(wildcards):
	lis = []
#	print(wildcards.sample)
#	print(Illumina_process_df.set_index("sample").index)
	for l in Illumina_process_df.set_index("sample").loc[(wildcards.sample), ["lib"]].values:
		if type(l) is str:
			lis.append("results/{sample}/trimming/trimgalore/{lib}/{sample}.{lib}.2.fastq.gz".format(sample=wildcards["sample"], lib=l))
		else:
			lis.append("results/{sample}/trimming/trimgalore/{lib}/{sample}.{lib}.2.fastq.gz".format(sample=wildcards["sample"], lib=l[0]))
#	print(lis)
	return lis
def input_for_clean_trim_fo(wildcards):
	lis = []
#	print(wildcards.sample)
#	print(Illumina_process_df.set_index("sample").index)
	for l in Illumina_process_df.set_index("sample").loc[(wildcards.sample), ["lib"]].values:
		if type(l) is str:
			lis.append("results/{sample}/trimming/trimgalore/{lib}/{sample}.{lib}.unpaired.1.fastq.gz".format(sample=wildcards["sample"], lib=l))
		else:
			lis.append("results/{sample}/trimming/trimgalore/{lib}/{sample}.{lib}.unpaired.1.fastq.gz".format(sample=wildcards["sample"], lib=l[0]))
#	print(lis)
	return lis
def input_for_clean_trim_ro(wildcards):
	lis = []
#	print(wildcards.sample)
#	print(Illumina_process_df.set_index("sample").index)
	for l in Illumina_process_df.set_index("sample").loc[(wildcards.sample), ["lib"]].values:
		if type(l) is str:
			lis.append("results/{sample}/trimming/trimgalore/{lib}/{sample}.{lib}.unpaired.2.fastq.gz".format(sample=wildcards["sample"], lib=l))
		else:
			lis.append("results/{sample}/trimming/trimgalore/{lib}/{sample}.{lib}.unpaired.2.fastq.gz".format(sample=wildcards["sample"], lib=l[0]))
#	print(lis)
	return lis

#get paths to the fastq read files
def get_raw_f_fastqs(wildcards):
	return df_fastq.loc[(wildcards.lib), ["f_fastq"]].dropna()
#	return pd.read_csv(config["samples"], dtype=str, sep="\t").set_index(["lib"], drop=False).loc[(wildcards.lib), ["f_fastq"]].dropna()
def get_raw_r_fastqs(wildcards):
	return df_fastq.loc[(wildcards.lib), ["r_fastq"]].dropna()
#	return pd.read_csv(config["samples"], dtype=str, sep="\t").set_index(["lib"], drop=False).loc[(wildcards.lib), ["r_fastq"]].dropna()
#function that gets the path to the bam files
def get_bam(wildcards):
	print("get_bam: "+str(wildcards.lib))
	if wildcards.lib in df_bam.index:
		return df_bam.loc[(wildcards.lib), ["bam"]].dropna()
#	else:
#		return []
#	return pd.read_csv(config["samples"], dtype=str, sep="\t").set_index(["lib"], drop=False).loc[(wildcards.lib), ["bam"]].dropna()
#function that gets the path to the long reads
def get_long(wildcards):
	return df_long.loc[(wildcards.lib), ["long"]].dropna()
#function that gets the path to the fast5 directory


def get_fast5_dir(wildcards):
	return df_fast5.loc[(wildcards.lib), ["fast5_dir"]]

#def get_ont_called_by_lib(wildcards):
#	print("wildcards for get called reads individual by library input function: "+str(wildcards))

def get_called_by_sample_by_lib(wildcards):
#	print("\nwildcards for 'get_called_by_sample_by_lib' basecalled input function: "+str(wildcards))
	lis = []
	for u in fast5_units.itertuples():
#		print("sample: %s, lib: %s" %(u.sample, u.lib))
		if str(u.sample) == wildcards.sample:
			lib = u.lib
			if str(lib) == wildcards.lib:
				for unit in flappie_unit_list:
					lis.append("results/{sample}/reads/ont/{basecaller}/{lib}/{sample}.{basecaller}.{unit}.fastq.gz".format(sample=wildcards.sample, lib=lib, basecaller=wildcards.basecaller, unit=unit))
#	print(len(lis))
	return lis
	

## this function gathers called reads for direct correction
def get_called_by_sample(wildcards):
#	print("\nwildcards for 'get_called_by_sample' basecalled input function: "+str(wildcards))
	lis = []
	for u in fast5_units.itertuples():
#		print("sample: %s, lib: %s" %(u.sample, u.lib))
		if str(u.sample) == wildcards.sample:
			lib = u.lib
			print("\tinput for %s and %s" %(u.sample, u.lib))
#		for basecaller in config["basecaller"]:
			for unit in flappie_unit_list:
				lis.append("results/{sample}/reads/ont/{basecaller}/{lib}/{sample}.{basecaller}.{unit}.fastq.gz".format(sample=wildcards.sample, lib=lib, basecaller=wildcards.basecaller, unit=unit))
#	print(len(lis))
	return lis


## this function gathers corrected reads from each library
def gather_corrected_by_lib(wildcards):
#	print("\nwildcards for 'gather_corrected' input function: "+str(wildcards))
#	print("longcorrection: "+str(wildcards.longcorrection))
	lis = []
	for u in fast5_units.itertuples():
#		print("sample: %s, lib: %s" %(u.sample, u.lib))
		if str(u.sample) == wildcards.sample:
			lib = u.lib
			if wildcards.longcorrection in ["ratatosk"]:
				lis.append("results/{sample}/errorcorrection/{longcorrection}/{trimmer}-{corrector}-{merger}-{basecaller}/{lib}/{sample}.{lib}.{basecaller}.{longcorrection}.fastq.gz".format(sample=wildcards.sample, lib=lib, trimmer=wildcards.trimmer, corrector=wildcards.corrector, merger=wildcards.merger, longcorrection=wildcards.longcorrection, basecaller=wildcards.basecaller))
#	print(lis)
	return lis

#creat list for paralllel steps
#ec_unit_list = range(1,ec_concurrency+1)
flappie_unit_list = range(1,config["ont_basecalling"]["concurrency"]+1)

def control_illumina_ec(wildcards):
	lis = []
	for ass in config["illumina_correction"]:
		if ass == "spades":
			lis.append("results/{sample}/errorcorrection/spades/{trimmer}/spades-correct.ok".format(sample=wildcards.sample, trimmer=wildcards.trimmer))
		if ass == "bless":
			lis.append("results/{sample}/errorcorrection/bless/{trimmer}/bless-kbest-se.done".format(sample=wildcards.sample, trimmer=wildcards.trimmer))
			lis.append("results/{sample}/errorcorrection/bless/{trimmer}/bless-kbest-pe.done".format(sample=wildcards.sample, trimmer=wildcards.trimmer))
#	print(lis)
	return lis

def to_merge(wildcards):
	lis = []
	if wildcards.corrector == "None":
		if wildcards.trimmer == "trimgalore":
			lis.append("results/{sample}/trimming/{trimmer}/{sample}-full/{sample}.trimgalore.1.fastq.gz".format(sample=wildcards.sample, trimmer=wildcards.trimmer))
			lis.append("results/{sample}/trimming/{trimmer}/{sample}-full/{sample}.trimgalore.2.fastq.gz".format(sample=wildcards.sample, trimmer=wildcards.trimmer))
	else:
		if wildcards.corrector == "bless":
			lis.append("results/{sample}/errorcorrection/bless/{trimmer}/bless-kbest-corrected/bless.corrected.1.fastq.gz".format(sample=wildcards.sample, trimmer=wildcards.trimmer))
			lis.append("results/{sample}/errorcorrection/bless/{trimmer}/bless-kbest-corrected/bless.corrected.2.fastq.gz".format(sample=wildcards.sample, trimmer=wildcards.trimmer))
		elif wildcards.corrector == "spades":
			lis.append("results/{sample}/errorcorrection/spades/{trimmer}/spades-corrected/corrected/spades.corrected.1.fastq.gz".format(sample=wildcards.sample, trimmer=wildcards.trimmer))
			lis.append("results/{sample}/errorcorrection/spades/{trimmer}/spades-corrected/corrected/spades.corrected.2.fastq.gz".format(sample=wildcards.sample, trimmer=wildcards.trimmer))
	return lis

def get_illumina_assembly_input(wildcards):
	lis = []
#	print("wildcards for 'get_illumina_assembly_input' input function: "+str(wildcards))
	if wildcards.trimmer != "None" and wildcards.corrector == "None" and wildcards.merger == "None":
		lis.append("results/{sample}/trimming/{trimmer}/{sample}-full/{sample}.{trimmer}.1.fastq.gz".format(sample=wildcards["sample"], trimmer=wildcards.trimmer))
		lis.append("results/{sample}/trimming/{trimmer}/{sample}-full/{sample}.{trimmer}.2.fastq.gz".format(sample=wildcards["sample"], trimmer=wildcards.trimmer))
		lis.append("results/{sample}/trimming/{trimmer}/{sample}-full/{sample}.{trimmer}.se.fastq.gz".format(sample=wildcards["sample"], trimmer=wildcards.trimmer))
	elif wildcards.trimmer != "None" and wildcards.corrector != "None" and wildcards.merger == "None":
		if wildcards.corrector == "bless":
			lis.append("results/{sample}/errorcorrection/bless/{trimmer}/bless-kbest-corrected/bless.corrected.1.fastq.gz".format(sample=wildcards.sample, trimmer=wildcards.trimmer))
			lis.append("results/{sample}/errorcorrection/bless/{trimmer}/bless-kbest-corrected/bless.corrected.2.fastq.gz".format(sample=wildcards.sample, trimmer=wildcards.trimmer))
			lis.append("results/{sample}/errorcorrection/bless/{trimmer}/bless-kbest-corrected/bless.corrected.se.fastq.gz".format(sample=wildcards.sample, trimmer=wildcards.trimmer))
		elif wildcards.corrector == "spades":
			lis.append("results/{sample}/errorcorrection/spades/{trimmer}/spades-corrected/corrected/spades.corrected.1.fastq.gz".format(sample=wildcards.sample, trimmer=wildcards.trimmer))
			lis.append("results/{sample}/errorcorrection/spades/{trimmer}/spades-corrected/corrected/spades.corrected.2.fastq.gz".format(sample=wildcards.sample, trimmer=wildcards.trimmer))
			lis.append("results/{sample}/errorcorrection/spades/{trimmer}/spades-corrected/corrected/spades.corrected.se.fastq.gz".format(sample=wildcards.sample, trimmer=wildcards.trimmer))
	elif (wildcards.trimmer != "None" and wildcards.corrector != "None" and wildcards.merger != "None") or (wildcards.trimmer != "None" and wildcards.corrector == "None" and wildcards.merger != "None"):
			lis.append("results/{sample}/read-merging/{merger}/{trimmer}-{corrector}/{sample}.notmerged.1.fastq.gz".format(sample=wildcards.sample, trimmer=wildcards.trimmer, corrector=wildcards.corrector, merger=wildcards.merger))
			lis.append("results/{sample}/read-merging/{merger}/{trimmer}-{corrector}/{sample}.notmerged.2.fastq.gz".format(sample=wildcards.sample, trimmer=wildcards.trimmer, corrector=wildcards.corrector, merger=wildcards.merger))
			lis.append("results/{sample}/read-merging/{merger}/{trimmer}-{corrector}/{sample}.merged.fastq.gz".format(sample=wildcards.sample, trimmer=wildcards.trimmer, corrector=wildcards.corrector, merger=wildcards.merger))
	return lis

longcordic = {'methods':[]}
for c in config["long_correction"]:
#	print(c)
	longcordic["methods"].append(str(c))
longcor = pd.DataFrame(longcordic).set_index(['methods'], drop=False)


def find_assemblies(wildcards):
#function finds all assemblies for the assembler specified in the config file
	lis = []
	for ass in config["assemble"]["assembler"]:
		if ass == "spades":
			lis.extend(glob.glob("results/"+wildcards.sample+"/assembly/spades/*/*/"+wildcards.sample+"/scaffolds.fasta"))
		if ass == "abyss":
			lis.extend(glob.glob("results/"+wildcards.sample+"/assembly/abyss/*/bestk/"+wildcards.sample+"-scaffolds.fa"))
		if ass == "minia":
			lis.extend(glob.glob("results/"+wildcards.sample+"/assembly/minia/*/bestk/*[0-9].contigs.fa"))
		if ass == "platanus":
			lis.extend(glob.glob("results/"+wildcards.sample+"/assembly/platanus/*/auto/*_gapClosed.fa"))
	return lis

def find_new_assemblies(wildcards):
#function finds all assemblies for the assembler specified in the config file
	lis = []
	for ass in config["assemble"]["assembler"]:
		if ass == "spades":
			lis.extend(glob.glob("results/"+wildcards.sample+"/assembly/spades/*/*/"+wildcards.sample+"/scaffolds.fasta"))
		if ass == "abyss":
			lis.extend(glob.glob("results/"+wildcards.sample+"/assembly/abyss/*/bestk/"+wildcards.sample+"-scaffolds.fa"))
		if ass == "minia":
			lis.extend(glob.glob("results/"+wildcards.sample+"/assembly/minia/*/bestk/*[0-9].contigs.fa"))
		if ass == "platanus":
			lis.extend(glob.glob("results/"+wildcards.sample+"/assembly/platanus/*/auto/*_gapClosed.fa"))
	print("ALL ASSEMBLIES FOUND: "+str(lis))
	for i in reversed(range(len(lis))):
		prefix="-".join(lis[i].split("/")[3:6])
#		print(prefix+" - results/"+wildcards.sample+"/assembly_evaluation/assemblies/"+prefix+".fasta")
		if glob.glob("results/"+wildcards.sample+"/assembly_evaluation/assemblies/"+prefix+".min"+str(config["evaluate_assemblies"]["minlength"])+".fasta"):
			print("found: 'results/"+wildcards.sample+"/assembly_evaluation/assemblies/"+prefix+".fasta'")
#			print(str(os.path.getctime("results/"+wildcards.sample+"/assembly_evaluation/assemblies/"+prefix+".fasta"))+" results/"+wildcards.sample+"/assembly/evaluation/assemblies/"+prefix+".fasta")
#			print(str(os.path.getctime(lis[i]))+" "+lis[i])
#			print("older file: "+str(min([lis[i], "results/"+wildcards.sample+"/assembly_evaluation/assemblies/"+prefix+".fasta"], key=os.path.getctime)))
			#delete from list if the original assembly file is older than the new one - nothing should be done - othterwise it will be returned by the input function
		
			if os.path.getctime(lis[i]) <= os.path.getctime("results/"+wildcards.sample+"/assembly_evaluation/assemblies/"+prefix+".min"+str(config["evaluate_assemblies"]["minlength"])+".fasta"):
				del(lis[i])
	print("NEW ASSEMBLIES: "+str(lis))
	return lis

def find_samples_with_assemblies(all_samples):
	lis = []
	for s in all_samples:
		nfiles = 0
		for ass in config["assemble"]["assembler"]:
			file = ass+".ok"
			for root, dirs, files in os.walk("results/"+s+"/assembly/"):
				if file in files:
					nfiles+=1
					break
		if nfiles > 0:
			lis.append(s)
#	print("Assemblies found for the following samples: "+str(lis))
	return lis
			
			

#def gather_assemblies(wildcards):
#	lis = []
#	print("wildcards for quast input function: "+str(wildcards))
#	for input in config["illumina_assembly_input"]:
#		if wildcards.sample in df_fastq["sample"].tolist() or wildcards.sample in df_bam["sample"].tolist():
#			if "abyss" in config["assemble"]["assembler"]:
#				lis.append("results/{sample}/assembly/{assinput}/abyss/bestk/abyss.ok".format(sample=wildcards["sample"], assinput=input))
#			if "minia" in config["assemble"]["assembler"]:
#				lis.append("results/{sample}/assembly/{assinput}/minia/bestk/minia.ok".format(sample=wildcards["sample"], assinput=input))
#			if "platanus" in config["assemble"]["assembler"]:
#				lis.append("results/{sample}/assembly/{assinput}/platanus/auto/platanus.ok".format(sample=wildcards["sample"], assinput=input))
#			if "spades" in config["assemble"]["assembler"]:
#				for kmode in ['default', 'bestk']:
#					lis.append("results/{sample}/assembly/{assinput}/spades/{kmode}/spades.ok".format(sample=wildcards["sample"], assinput=input, kmode=kmode))

def gather_assemblies(wildcards):
	lis = []
#	print("wildcards for quast input function: "+str(wildcards))
	if wildcards.sample in df_fastq["sample"].tolist() or wildcards.sample in df_bam["sample"].tolist():
		if "abyss" in config["assemble"]["assembler"]:
			for trimmer in trim_list:
				for corrector in correct_list:
					for merger in merge_list:
						lis.append("results/{sample}/assembly/abyss/{trimmer}-{corrector}-{merger}/bestk/abyss.ok".format(sample=wildcards["sample"], trimmer=trimmer, corrector=corrector, merger=merger))
		if "minia" in config["assemble"]["assembler"]:
			for trimmer in trim_list:
				for corrector in correct_list:
					for merger in merge_list:
						lis.append("results/{sample}/assembly/minia/{trimmer}-{corrector}-{merger}/bestk/minia.ok".format(sample=wildcards["sample"], trimmer=trimmer, corrector=corrector, merger=merger))
		if "platanus" in config["assemble"]["assembler"]:
			for trimmer in trim_list:
				for corrector in correct_list:
					for merger in merge_list:
						lis.append("results/{sample}/assembly/platanus/{trimmer}-{corrector}-{merger}/auto/platanus.ok".format(sample=wildcards["sample"], trimmer=trimmer, corrector=corrector, merger=merger))
		if "spades" in config["assemble"]["assembler"]:
			for kmode in ['default', 'bestk']:
				for trimmer in trim_list:
					for corrector in correct_list:
						for merger in merge_list:
							lis.append("results/{sample}/assembly/spades/{trimmer}-{corrector}-{merger}/{kmode}/spades.ok".format(sample=wildcards["sample"], trimmer=trimmer, corrector=corrector, merger=merger, kmode=kmode))

		#further assemblers that require fast5 data as input AND fastq reads for correction or alike
		if wildcards.sample in df_fast5["sample"].tolist():
			for basecaller in config["ont_basecalling"]["basecaller"]:
				if "haslr" in config["assemble"]["assembler"]:
					for trimmer in trim_list:
						for corrector in correct_list:
							for merger in ["None"]:
								lis.append("results/{sample}/assembly/haslr/{trimmer}-{corrector}-{merger}-{basecaller}/haslr.ok".format(sample=wildcards["sample"], trimmer=trimmer, corrector=corrector, merger=merger, basecaller=basecaller))
				if "abyss_scaffold" in config["assemble"]["assembler"]:
					for trimmer in trim_list:
						for corrector in correct_list:
							for merger in ["None"]:
								for l in config["long_correction"]:
									lis.append("results/{sample}/assembly/abyss_scaffold/{trimmer}-{corrector}-{merger}-{basecaller}-{longcorrection}/abyss.ok".format(sample=wildcards["sample"], trimmer=trimmer, corrector=corrector, merger=merger, basecaller=basecaller, longcorrection=l))
				if "spades_hybrid" in config["assemble"]["assembler"]:
					for trimmer in trim_list:
						for corrector in correct_list:
							for merger in ["None"]:
								if wildcards.sample in Illumina_process_df["sample"].tolist():
									for l in config["long_correction"]:
										lis.append("results/{sample}/assembly/spades_hybrid/{trimmer}-{corrector}-{merger}-{basecaller}-{longcorrection}/spades.ok".format(sample=wildcards["sample"], trimmer=trimmer, corrector=corrector, merger=merger, basecaller=basecaller, longcorrection=l))

	#assemblers that do only with fast5 data
	for basecaller in config["ont_basecalling"]["basecaller"]:
		if "canu" in config["assemble"]["assembler"]:
			if wildcards.sample in df_fast5["sample"].tolist():
				lis.append("results/{sample}/assembly/canu/raw/{basecaller}/canu.ok".format(sample=wildcards["sample"], basecaller=basecaller))
				for l in config["long_correction"]:
					if l in ["consent", "canucorrect"]:
						lis.append("results/{sample}/assembly/canu/None-None-None-{basecaller}-{longcorrection}/canu.ok".format(sample=wildcards["sample"], longcorrection=l, basecaller=basecaller))
					if l in ["ratatosk"]:
						if wildcards.sample in Illumina_process_df["sample"].tolist():
							for trimmer in trim_list:
								for corrector in correct_list:
									for merger in ["None"]:
										lis.append("results/{sample}/assembly/canu/{trimmer}-{corrector}-{merger}-{basecaller}-{longcorrection}/canu.ok".format(sample=wildcards["sample"], trimmer=trimmer, corrector=corrector, merger=merger, longcorrection=l, basecaller=basecaller))
		if "flye" in config["assemble"]["assembler"]:
			if wildcards.sample in df_fast5["sample"].tolist():
				lis.append("results/{sample}/assembly/flye/raw/{basecaller}/flye.ok".format(sample=wildcards["sample"], basecaller=basecaller))
				for l in config["long_correction"]:
					if l in ["consent", "canucorrect"]:
						lis.append("results/{sample}/assembly/flye/None-None-None-{basecaller}-{longcorrection}/flye.ok".format(sample=wildcards["sample"], longcorrection=l, basecaller=basecaller))
					if l in ["ratatosk"]:
						if wildcards.sample in Illumina_process_df["sample"].tolist():
							for trimmer in trim_list:
								for corrector in correct_list:
									for merger in ["None"]:
										lis.append("results/{sample}/assembly/flye/{trimmer}-{corrector}-{merger}-{basecaller}-{longcorrection}/flye.ok".format(sample=wildcards["sample"], trimmer=trimmer, corrector=corrector, merger=merger, longcorrection=l, basecaller=basecaller))
#	print(lis)
	return lis

def get_long_assembly_input(wildcards):
#	print("\nwildcards for long basecalled input function: "+str(wildcards))
	lis = []
	for u in fast5_units.itertuples():
#		print("sample: %s, lib: %s" %(u.sample, u.lib))
		if str(u.sample) == wildcards.sample:
			lib = u.lib
#			print("\tinput for %s and %s" %(u.sample, u.lib))
#			for basecaller in config["ont_basecalling"]["basecaller"]:
			for unit in flappie_unit_list:
				lis.append("results/{sample}/reads/ont/{basecaller}/{lib}/{sample}.{basecaller}.{unit}.fastq.gz".format(sample=wildcards.sample, lib=lib, basecaller=wildcards.basecaller, unit=unit))
#	print(lis)
	return lis

def get_long_corrected_assembly_input(wildcards):
	lis = []
#	print("\nwildcards for 'get_long_corrected_assembly_input' input function: "+str(wildcards))
	if wildcards.longcorrection == "ratatosk":
		lis.append("results/{sample}/errorcorrection/ratatosk/{trimmer}-{corrector}-None-{basecaller}/{sample}.{basecaller}.ratatosk.fastq.gz".format(sample=wildcards.sample, trimmer=wildcards.trimmer, corrector=wildcards.corrector, basecaller=wildcards.basecaller))
	if wildcards.longcorrection == "canucorrect":
		lis.append("results/{sample}/errorcorrection/canucorrect/None-None-None-{basecaller}/{sample}.{basecaller}.canucorrect.fastq.gz".format(sample=wildcards.sample, basecaller=wildcards.basecaller))
	if wildcards.longcorrection == "consent":
		lis.append("results/{sample}/errorcorrection/consent/None-None-None-{basecaller}/{sample}.{basecaller}.consent.fastq.gz".format(sample=wildcards.sample, basecaller=wildcards.basecaller))
#	print(lis)
	return lis
