rule eva_quast:
	input:
		gather_assemblies
	output:
		ok = "results/{sample}/assembly/quast/quast.ok",
	log:
		stdout = "results/{sample}/logs/quast.stdout.txt",
		stderr = "results/{sample}/logs/quast.stderr.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		dir = "results/{sample}/assembly/quast",
	singularity: "docker://reslp/quast:5.0.2"
	shell:
		"""
		echo "{input}" 1> {log.stdout} 2> {log.stderr}
		touch {output.ok}
		"""

rule eva_just_quast:
	input:
		find_assemblies
	output:
		ok = "results/{sample}/assembly/quast/just_quast.ok",
	log:
		log = "results/{sample}/logs/just_quast.log.txt",
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		dir = "results/{sample}/assembly/quast",
		minlen = 500
	singularity: "docker://reslp/quast:5.0.2"
	shell:
		"""
		labels=$(echo "{input}" | tr ' ' '\\n' | cut -d "/" -f 2,4-6 | sed 's/\//-/g' | tr '\\n' ',' | sed 's/,$//')
		quast -o {params.dir} -m {params.minlen} --labels $labels -t {threads} {input} 2>&1 | tee {log.log}
		touch {output.ok}
		"""
#		quast -o {params.dir} -m {params.minlen} -t {threads} --labels {input}

##from $BINFL/../reslp/platyhelminthes_phylo/rules/busco.smk
#rule download_busco_set:
#	output:
#		busco_set = directory("results/busco_set"),
#		checkpoint = "results/checkpoints/download_busco_set.done"
#	params:
#		set = config["busco"]["set"]
#	benchmark:
#		"results/statistics/benchmarks/setup/dowload_busco_set.txt"
#	shell:
#		"""
#		echo {params.set}
#		wget http://busco.ezlab.org/v2/datasets/{params.set}.tar.gz
#		tar xfz {params.set}.tar.gz
#		rm {params.set}.tar.gz
#		if [ -d {output.busco_set} ]; then rm -rf {output.busco_set}; fi
#		mv {params.set} {output.busco_set}
#		touch {output.checkpoint}
#		"""
#rule run_busco:
#	input:
#		assembly = "results/assemblies/{species}.fna",
#		augustus_config_path = rules.prepare_augustus.output.augustus_config_path,
#		busco_set = rules.download_busco_set.output.busco_set,
#	output:
#		checkpoint = "results/checkpoints/busco/busco_{species}.done",
#		augustus_output = "results/busco/{species}/run_busco/augustus_output.tar.gz",
#		blast_output = "results/busco/{species}/run_busco/blast_output.tar.gz",
#		hmmer_output = "results/busco/{species}/run_busco/hmmer_output.tar.gz",
#		full_table = "results/busco/{species}/run_busco/full_table_busco.tsv",
#		short_summary ="results/busco/{species}/run_busco/short_summary_busco.txt",
#		missing_busco_list ="results/busco/{species}/run_busco/missing_busco_list_busco.tsv",
#		single_copy_buscos = directory("results/busco/{species}/run_busco/single_copy_busco_sequences")
#	benchmark: "results/statistics/benchmarks/busco/run_busco_{species}.txt"
#	threads: config["busco"]["threads"]
#	shadow: "shallow"
#	params:
#		wd = os.getcwd(),
#		sp = config["busco"]["augustus_species"],
#		additional_params = config["busco"]["additional_parameters"],
#		species = lambda wildcards: wildcards.species
#	singularity:
#		"docker://reslp/busco:3.0.2"
#	shell:
#		"""
#		mkdir -p log
#		dir=results/busco/{params.species}
#		# prepare stripped down version auf augustus config path.
#		# this is introduced to lower the number of files.
#		mkdir augustus
#		cp -R /opt/conda/config/cgp augustus
#		cp /opt/conda/config/config.ini augustus
#		cp -R /opt/conda/config/extrinsic augustus
#		cp -R /opt/conda/config/model augustus
#		cp -R /opt/conda/config/profile augustus
#		mkdir augustus/species
#		
#		if [ -d /opt/conda/config/species/{params.sp} ]
#		then
#			cp -R /opt/conda/config/species/{params.sp} augustus/species
#		fi		
#
#		export AUGUSTUS_CONFIG_PATH=$(pwd)/augustus
#		echo $AUGUSTUS_CONFIG_PATH
#		run_busco -i {input.assembly} -f --out busco -c {threads} -sp {params.sp} --lineage_path {input.busco_set} -m genome {params.additional_params}
#		# do some cleanup to save space
#		bin/tar_folder.sh {output.blast_output} run_busco/blast_output
#		bin/tar_folder.sh {output.hmmer_output} run_busco/hmmer_output
#		bin/tar_folder.sh {output.augustus_output} run_busco/augustus_output
