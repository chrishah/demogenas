rule cor_bless_by_k:
	input:
		forward = "results/{sample}/trimming/{trimmer}/{sample}-full/{sample}.{trimmer}.1.fastq.gz",
		reverse = "results/{sample}/trimming/{trimmer}/{sample}-full/{sample}.{trimmer}.2.fastq.gz",
		orphans = "results/{sample}/trimming/{trimmer}/{sample}-full/{sample}.{trimmer}.se.fastq.gz",
	output:
		ok = "results/{sample}/errorcorrection/bless/{trimmer}/bless-k{blessk}/bless-k{blessk}.done"
	log:
		"results/{sample}/logs/bless-k{blessk}-{trimmer}.log.txt",
	benchmark: "results/{sample}/benchmarks/bless-k{blessk}-{trimmer}.txt"
	params:
		dir = "results/{sample}/errorcorrection/bless/{trimmer}/bless-k{blessk}",
		k = "{blessk}"
	singularity:
		"docker://chrishah/bless:v1.02"
	resources:
		mem_gb=config["max_mem_in_GB"]["bless"]
	threads: config["threads"]["bless"]
	shell:
		"""
		#set stage
		if [[ ! -d {params.dir} ]]; then mkdir -p {params.dir}; else rm -rf {params.dir}/*; fi
		#create empty file to suck up corrected reads to save space at this stage
		ln -s /dev/null {params.dir}/bless-k{params.k}.corrected.fastq.00000
		#concatenate all reads
		cat {input.forward} {input.reverse} {input.orphans} > {params.dir}/all_reads.fastq.gz

		#run bless
		bless -read {params.dir}/all_reads.fastq.gz -kmerlength {params.k} -max_mem {resources.mem_gb} -notrim -smpthread {threads} -prefix {params.dir}/bless-k{params.k} 2>&1 | tee {log}
		#clean up
		rm {params.dir}/all_reads.fastq.gz

		touch {output.ok}
		"""


#rule cor_find_best_bless:
#	input:
#		expand("results/{{sample}}/errorcorrection/bless/{{trimmer}}/bless-k{blessk}/bless-k{blessk}.done", sample=Illumina_process_df["sample"], blessk=config["bless_k"], trimmer=config["illumina_trimming"])
#	output:
#		ok = "results/{sample}/errorcorrection/bless/{trimmer}/bless-bestk-corrected/bless-bestk-correct.ok"
#	log:
#		stdout = "results/{sample}/logs/bless-bestk-{trimmer}.stdout.txt",
#		stderr = "results/{sample}/logs/bless-bestk-{trimmer}.stderr.txt"
#	benchmark: "results/{sample}/benchmarks/bless-bestk-{trimmer}.txt"
#	params:
#		wd = os.getcwd(),
#		sample = "{sample}",
#	singularity:
#		"docker://chrishah/bless:v1.02"
#	threads: 2
#	shadow: "minimal"
#	shell:
#		"""
#		#read in outputs from the previous bless runs and decide on best k
#		touch {output.ok}
#		"""

rule cor_bless_pe:
	input:
		forward = "results/{sample}/trimming/{trimmer}/{sample}-full/{sample}.{trimmer}.1.fastq.gz",
		reverse = "results/{sample}/trimming/{trimmer}/{sample}-full/{sample}.{trimmer}.2.fastq.gz",
		bless_runs = expand("results/{{sample}}/errorcorrection/bless/{{trimmer}}/bless-k{blessk}/bless-k{blessk}.done", sample=Illumina_process_df["sample"], blessk=config["bless_k"], trimmer=config["illumina_trimming"])
	wildcard_constraints:
		trimmer="trimgalore.*",
	output:
		ok = "results/{sample}/errorcorrection/bless/{trimmer}/bless-kbest-pe.done",
		forward = "results/{sample}/errorcorrection/bless/{trimmer}/bless-kbest-corrected/bless.corrected.1.fastq.gz",
		reverse = "results/{sample}/errorcorrection/bless/{trimmer}/bless-kbest-corrected/bless.corrected.2.fastq.gz",
	log:
		"results/{sample}/logs/bless-kbest-pe-{trimmer}.log.txt",
	benchmark: "results/{sample}/benchmarks/bless-kbest-pe-{trimmer}.txt"
	singularity:
		"docker://chrishah/bless:v1.02"
	params:
		wd = os.getcwd(),
		dir = "results/{sample}/errorcorrection/bless/{trimmer}"
	resources:
		mem_gb=config["max_mem_in_GB"]["bless"]
	threads: config["threads"]["bless"]
	shell:
		"""
		./bin/gather_bless_stats.sh {wildcards.sample} {input.bless_runs} > {params.dir}/bless_stats_pe.tsv

		bestk=$(tail -n 1 {params.dir}/bless_stats_pe.tsv | cut -f 1)
		echo -e "\\n\\n[$(date)]\tCorrecting paired end reads with bless - k: $bestk" 2>&1 | tee {log}

		#set stage
		if [[ ! -d {params.dir}/bless-kbest-corrected ]]; then mkdir -p {params.dir}/bless-kbest-corrected; else rm -rf {params.dir}/bless-kbest-corrected/bless-k$bestk-pe*; fi
		#run bless
		bless -read1 {input.forward} -read2 {input.reverse} -kmerlength $bestk -max_mem {resources.mem_gb} -load {params.dir}/bless-k$bestk/bless-k$bestk -notrim -smpthread {threads} -gzip -prefix {params.dir}/bless-kbest-corrected/bless-k$bestk-pe 2>&1 | tee {log}
		cd {params.dir}/bless-kbest-corrected/
		ln -s bless-k$bestk-pe.1.corrected.fastq.gz $(basename {output.forward})
		ln -s bless-k$bestk-pe.2.corrected.fastq.gz $(basename {output.reverse})
		cd {params.wd}
		touch {output.ok}
		"""
rule cor_bless_se:
	input:
		orphans = "results/{sample}/trimming/{trimmer}/{sample}-full/{sample}.{trimmer}.se.fastq.gz",
		bless_runs = expand("results/{{sample}}/errorcorrection/bless/{{trimmer}}/bless-k{blessk}/bless-k{blessk}.done", sample=Illumina_process_df["sample"], blessk=config["bless_k"], trimmer=config["illumina_trimming"])
	wildcard_constraints:
		trimmer="trimgalore.*",
	output:
		orphans = "results/{sample}/errorcorrection/bless/{trimmer}/bless-kbest-corrected/bless.corrected.se.fastq.gz",
		ok = "results/{sample}/errorcorrection/bless/{trimmer}/bless-kbest-se.done"
	log:
		"results/{sample}/logs/bless-kbest-se-{trimmer}.log.txt",
	benchmark: "results/{sample}/benchmarks/bless-kbest-se-{trimmer}.txt"
	singularity:
		"docker://chrishah/bless:v1.02"
	params:
		wd = os.getcwd(),
		dir = "results/{sample}/errorcorrection/bless/{trimmer}",
	resources:
		mem_gb=config["max_mem_in_GB"]["bless"]
	threads: config["threads"]["bless"]
	shell:
		"""
		./bin/gather_bless_stats.sh {wildcards.sample} {input.bless_runs} > {params.dir}/bless_stats-se.tsv

		bestk=$(tail -n 1 {params.dir}/bless_stats-se.tsv | cut -f 1)
		echo -e "\\n\\n[$(date)]\tCorrecting single end reads with bless - k: $bestk" 2>&1 | tee {log}
		
		#set stage
		if [[ ! -d {params.dir}/bless-kbest-corrected ]]; then mkdir -p {params.dir}/bless-kbest-corrected; else rm -rf {params.dir}/bless-kbest-corrected/bless-k$bestk-se*; fi
		#run bless
		bless -read {input.orphans} -kmerlength $bestk -max_mem {resources.mem_gb} -load {params.dir}/bless-k$bestk/bless-k$bestk -notrim -smpthread {threads} -gzip -prefix {params.dir}/bless-kbest-corrected/bless-k$bestk-se 2>&1 | tee {log}

		ln -s {params.wd}/{params.dir}/bless-kbest-corrected/bless-k$bestk-se.corrected.fastq.gz {output.orphans}
		touch {output.ok}
		"""


rule cor_correct_spades:
	input:
		forward = "results/{sample}/trimming/{trimmer}/{sample}-full/{sample}.{trimmer}.1.fastq.gz",
		reverse = "results/{sample}/trimming/{trimmer}/{sample}-full/{sample}.{trimmer}.2.fastq.gz",
		orphans = "results/{sample}/trimming/{trimmer}/{sample}-full/{sample}.{trimmer}.se.fastq.gz",
	wildcard_constraints:
		trimmer="trimgalore.*",
	output:
		forward = "results/{sample}/errorcorrection/spades/{trimmer}/spades-corrected/corrected/spades.corrected.1.fastq.gz",
		reverse = "results/{sample}/errorcorrection/spades/{trimmer}/spades-corrected/corrected/spades.corrected.2.fastq.gz",
		orphans = "results/{sample}/errorcorrection/spades/{trimmer}/spades-corrected/corrected/spades.corrected.se.fastq.gz",
		ok = "results/{sample}/errorcorrection/spades/{trimmer}/spades-correct.ok",
	log:
		stdout = "results/{sample}/logs/correct-spades-{trimmer}.stdout.txt",
		stderr = "results/{sample}/logs/correct-spades-{trimmer}.stderr.txt"
	benchmark: "results/{sample}/benchmarks/correct-spades-{trimmer}.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		dir = "results/{sample}/errorcorrection/spades/{trimmer}/spades-corrected",
	singularity:
		"docker://reslp/spades:3.15.3"
	threads: config["threads"]["spades_correct"]
	resources:
		mem_gb=150
	shell:
		"""
		echo "Host: $HOSTNAME" 1> {log.stdout} 2> {log.stderr}
		if [[ ! -d {params.dir} ]]; then mkdir -p {params.dir}; else rm -rf {params.dir}/*; fi
		spades.py \
		--only-error-correction \
		-o {params.dir} \
		-1 {input.forward} -2 {input.reverse} -s {input.orphans} \
		-t $(( {threads} - 1 )) \
		-m $(( {resources.mem_gb} - 5 )) 1>> {log.stdout} 2>> {log.stderr}

		cd {params.dir}/corrected
		ln -s {wildcards.sample}.{wildcards.trimmer}.se.fastq.*.cor.fastq.gz spades.corrected.se.fastq.gz
		ln -s {wildcards.sample}.{wildcards.trimmer}.1.fastq.*.cor.fastq.gz spades.corrected.1.fastq.gz
		ln -s {wildcards.sample}.{wildcards.trimmer}.2.fastq.*.cor.fastq.gz spades.corrected.2.fastq.gz
		cd - &> /dev/null
		touch {output.ok}
		"""

rule cor_gather_illumina_corrected:
	input:
		control_illumina_ec
	wildcard_constraints:
		trimmer="trimgalore.*",
	output:
		ok = "results/{sample}/errorcorrection/illumina-correction-{trimmer}.ok",
	threads: 1
	shell:
		"""
		touch {output.ok}
		"""
