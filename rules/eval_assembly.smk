rule ass_gather_assemblies:
	input:
		gather_assemblies
	output:
		"results/{sample}/assembly/all_assemblies.done",
	shell:
		"""
		touch {output}
		"""

s_with_ass = find_samples_with_assemblies(all_samples)
print(s_with_ass)

rule eva_prepare_assemblies:
	input:
		find_new_assemblies
	output:
		"results/{sample}/assembly_evaluation/assemblies/prepare_min"+str(config["evaluate_assemblies"]["minlength"])+".done"
	params:
		minlength = config["evaluate_assemblies"]["minlength"]
	singularity:
		"docker://reslp/biopython_plus:1.77"
	shell:
		"""
		for f in {input}
		do
			prefix=$(echo $f | cut -d "/" -f 4-6 | sed 's/\//-/g' | tr '\\n' ',' | sed 's/,$//')
			bin/lengthfilter.py $f {params.minlength} > results/{wildcards.sample}/assembly_evaluation/assemblies/$prefix.min{params.minlength}.fasta
		done
		touch {output}
		"""	

def gather_busco(wildcards):
	lisout=[]
	for l in glob.glob("results/"+wildcards.sample+"/assembly_evaluation/assemblies/*.min"+str(config["evaluate_assemblies"]["minlength"])+".fasta"):
		print(l)
		lisout.append("/".join(l.split("/")[:3])+"/busco/"+config["evaluate_assemblies"]["busco"]["set"]+"/"+".".join(l.split("/")[-1].split(".")[:-2])+".min"+str(config["evaluate_assemblies"]["minlength"])+".busco.done")
	print(wildcards.sample,str(lisout))
	return lisout

def gather_blobtools(wildcards):
	lisout=[]
	for l in glob.glob("results/"+wildcards.sample+"/assembly_evaluation/assemblies/*.min"+str(config["evaluate_assemblies"]["minlength"])+".fasta"):
		print(l)
		lisout.append("/".join(l.split("/")[:3])+"/blobtools/"+config["evaluate_assemblies"]["busco"]["set"]+"/"+".".join(l.split("/")[-1].split(".")[:-2])+".min"+str(config["evaluate_assemblies"]["minlength"])+".blobtools.done")
	print(wildcards.sample,str(lisout))
	return lisout

def gather_assemblies(wildcards):
	lis = glob.glob("results/"+wildcards.sample+"/assembly_evaluation/assemblies/*.min"+str(config["evaluate_assemblies"]["minlength"])+".fasta")
	print("\n\nList of file to be processed by quast:\n"+str(lis))
	return lis

rule eva_quast:
	input:
		assemblies = gather_assemblies
#		assemblies = find_assemblies,
#		assemblies = expand("results/{sample}/assembly_evaluation/assemblies/{combination}.fasta", sample=s_with_ass, combination=combinations),
#		assemblies = lambda wildcards: expand("results/{{sample}}/assembly_evaluation/assemblies/{combination}.fasta", zip, sample=samps, combination=cs),
#		assemblies = expand("results/{sample}/assembly_evaluation/assemblies/{combination}.fasta", sample=s_with_ass, combination=get_combs_per_sample(s_with_ass)),
	output:
		ok = "results/{sample}/assembly_evaluation/quast/quast.done",
	log:
		log = "results/{sample}/logs/quast.log.txt",
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		dir = "results/{sample}/assembly_evaluation/quast",
		minlen = config["evaluate_assemblies"]["minlength"]
	threads: config["threads"]["quast"]
	singularity: "docker://reslp/quast:5.0.2"
	shell:
		"""
		labels=$(echo "{input.assemblies}" | tr ' ' '\\n' | cut -d "/" -f 2,5 | sed 's/\//-/g' | sed 's/\.min[0-9]*\.fasta$//' | tr '\n' ',' | sed 's/,$//')
		quast -o {params.dir} -m {params.minlen} --labels $labels -t {threads} {input.assemblies} 2>&1 | tee {log.log}
		touch {output.ok}
		"""

rule download_busco_set:
	output:
		busco_set = directory("dbs/busco/busco_set/"+config["evaluate_assemblies"]["busco"]["set"]),
		done = "dbs/busco/busco_set/download_busco_set-"+config["evaluate_assemblies"]["busco"]["set"]+".done"
	params:
		set = config["evaluate_assemblies"]["busco"]["set"],
	log:
		log = "dbs/logs/download_busco_set-"+config["evaluate_assemblies"]["busco"]["set"]+".log.txt"
	shell:
		"""
		echo -e "[$(date)]\\tBUSCO set specified: {params.set}" 2>&1 | tee {log}
		if [ -d {output.busco_set} ]; then rm -rf {output.busco_set}; fi
		mkdir {output.busco_set}

		base_url="https://busco-data.ezlab.org/v5/data/lineages"
		current=$(curl -s $base_url/ | grep "{params.set}" | cut -d ">" -f 2 | sed 's/<.*//')
		echo -e "[$(date)]\\tCurrent version is: $current" 2>&1 | tee -a {log}
		echo -e "[$(date)]\\tDownloading .." 2>&1 | tee -a {log}
		wget -q -c $base_url/$current -O - --no-check-certificate | tar -xz --strip-components 1 -C {output.busco_set}/

		echo -ne "[$(date)]\\tDone!\\n" 2>&1 | tee -a {log}
		touch {output.done}
		"""

rule eva_busco:
	input:
		assembly = "results/{sample}/assembly_evaluation/assemblies/{combination}.fasta",
		busco_set = rules.download_busco_set.output.busco_set,
		script = "bin/tar_folder.sh"
	output:
		done = "results/{sample}/assembly_evaluation/busco/"+config["evaluate_assemblies"]["busco"]["set"]+"/{combination}.busco.done",
		output = "results/{sample}/assembly_evaluation/busco/"+config["evaluate_assemblies"]["busco"]["set"]+"/{combination}/run_busco/software_output.tar.gz",
		logs = "results/{sample}/assembly_evaluation/busco/"+config["evaluate_assemblies"]["busco"]["set"]+"/{combination}/run_busco/logs.tar.gz",
		full_table = "results/{sample}/assembly_evaluation/busco/"+config["evaluate_assemblies"]["busco"]["set"]+"/{combination}/run_busco/full_table_busco.tsv",
		short_summary ="results/{sample}/assembly_evaluation/busco/"+config["evaluate_assemblies"]["busco"]["set"]+"/{combination}/run_busco/short_summary_busco.txt",
		missing_busco_list ="results/{sample}/assembly_evaluation/busco/"+config["evaluate_assemblies"]["busco"]["set"]+"/{combination}/run_busco/missing_busco_list_busco.tsv",
		single_copy_buscos = "results/{sample}/assembly_evaluation/busco/"+config["evaluate_assemblies"]["busco"]["set"]+"/{combination}/run_busco/single_copy_busco_sequences.tar",
		single_copy_buscos_tarlist = "results/{sample}/assembly_evaluation/busco/"+config["evaluate_assemblies"]["busco"]["set"]+"/{combination}/run_busco/single_copy_busco_sequences.txt"

	threads: int(config["threads"]["busco"])
	shadow: "minimal"
	log:
		log = "results/{sample}/logs/busco."+config["evaluate_assemblies"]["busco"]["set"]+".{combination}.log.txt"
	params:
		wd = os.getcwd(),
		sp = config["evaluate_assemblies"]["busco"]["augustus_species"],
		additional_params = config["evaluate_assemblies"]["busco"]["additional_parameters"],
		mode = "genome",
		augustus_config_in_container = "/usr/local/config",
		set = config["evaluate_assemblies"]["busco"]["set"]
	singularity:
		"docker://ezlabgva/busco:v5.2.1_cv1"
	shell:
		"""
		# prepare stripped down version auf augustus config path.
		# this is introduced to lower the number of files.
		mkdir augustus
		cp -R {params.augustus_config_in_container}/cgp augustus
		cp -R {params.augustus_config_in_container}/extrinsic augustus
		cp -R {params.augustus_config_in_container}/model augustus
		cp -R {params.augustus_config_in_container}/profile augustus
		mkdir augustus/species
		cp -R {params.augustus_config_in_container}/species/generic augustus/species/
		
		if [ -d {params.augustus_config_in_container}/species/{params.sp} ]
		then
			cp -R {params.augustus_config_in_container}/species/{params.sp} augustus/species
		fi		

		export AUGUSTUS_CONFIG_PATH=$(pwd)/augustus
	
		threads=$(( {threads} / 2 ))
		if [ "$threads" -eq 0 ]; then threads=1; fi

		augustus_params=""	
		echo "Assembly used for BUSCO is {input.assembly}" 2>&1 | tee {log}
		if [[ -z "{params.sp}" ]] || [[ "{params.sp}" == "None" ]]
		then
			echo "using metaeuk" 2>&1 | tee {log}
		else
			echo "using augustus (species: {params.sp})" 2>&1 | tee {log}
			augustus_params="--augustus --augustus_species {params.sp}"
		fi
		busco -i {input.assembly} -f --out {wildcards.combination} -c $threads $augustus_params --lineage_dataset $(pwd)/{input.busco_set} -m {params.mode} {params.additional_params} 2>&1 | tee -a {log}
		# do some cleanup to save space
		echo -e "\\n[$(date)]\\tCleaning up after BUSCO to save space" 2>&1 | tee -a {log}
		basedir=$(pwd)
		cd {wildcards.combination}/run_{params.set}
		mkdir software_outputs
		mv *_output software_outputs
		$basedir/bin/tar_folder.sh $basedir/{output.output} software_outputs 2>&1 | tee -a $basedir/{log}
		cd ..
		$basedir/bin/tar_folder.sh $basedir/{output.logs} logs 2>&1 | tee -a $basedir/{log}
		cd ..
		tar -pcf {output.single_copy_buscos} -C {wildcards.combination}/run_{params.set}/busco_sequences single_copy_busco_sequences 
		tar -tvf {output.single_copy_buscos} > {output.single_copy_buscos_tarlist} 2>&1 | tee -a $basedir/{log}

		#move output files:
		mv {wildcards.combination}/run_{params.set}/full_table.tsv {output.full_table}
		mv {wildcards.combination}/run_{params.set}/short_summary.txt {output.short_summary}
		mv {wildcards.combination}/run_{params.set}/missing_busco_list.tsv {output.missing_busco_list}
		
		#touch checkpoint
		touch {output.done}
		"""

rule eva_gather_busco:
	input:
		gather_busco
#		lambda wildcards: expand("results/{{sample}}/assembly_evaluation/busco/"+config["evaluate_assemblies"]["busco"]["set"]+"/{combination}.busco.done", zip, sample=samps, combination=cs)
#		expand("results/{sample}/assembly_evaluation/busco/{combination}.busco.done", sample=s_with_ass, combination=combinations)
#		expand("results/{sample}/assembly_evaluation/busco/{combination}.busco.done", sample=s_with_ass, combination=get_combs_per_sample(s_with_ass))
	output:
		"results/{sample}/assembly_evaluation/busco/busco.done"
	shell:
		"""
		touch {output}
		"""

rule eva_blobtools:
	input:
		"results/{sample}/assembly_evaluation/mapping/{combination}/{sample}.indexing_durmvd.done"
	output:
		"results/{sample}/assembly_evaluation/blobtools/{combination}.blobtools.done"
	shell:
		"""
		touch {output}
		"""

rule eva_gather_blobtools:
	input:
		gather_blobtools
	output:
		"results/{sample}/assembly_evaluation/blobtools/blobtools.done"
	shell:
		"""
		touch {output}
		"""

rule eva_x_gather_evaluations:
	input: 
		expand("results/{sample}/assembly_evaluation/{evaluator}/{evaluator}.done", sample=s_with_ass, evaluator=config["evaluate_assemblies"]["evaluator"])
	output:
		"results/{sample}/assembly_evaluation.done"
	shell:
		"""
		touch {output}
		"""

