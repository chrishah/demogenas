rule bless_by_k:
	input:
		forward = "results/{sample}/trimming/trim_galore/{sample}-full/{sample}.trimgalore.1.fastq.gz",
		reverse = "results/{sample}/trimming/trim_galore/{sample}-full/{sample}.trimgalore.2.fastq.gz",
		orphans = "results/{sample}/trimming/trim_galore/{sample}-full/{sample}.trimgalore.se.fastq.gz",
	output:
		ok = "results/{sample}/errorcorrection/bless/bless-k{blessk}/bless-k{blessk}.done"
	log:
		"results/{sample}/logs/bless-k{blessk}.log.txt",
	benchmark: "results/{sample}/benchmarks/bless-k{blessk}.txt"
	params:
		dir = "results/{sample}/errorcorrection/bless/bless-k{blessk}",
		k = "{blessk}"
	singularity:
		"docker://chrishah/bless:v1.02"
	resources:
		mem_gb=20
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


rule find_best_bless:
	input:
		expand("results/{{sample}}/errorcorrection/bless/bless-k{blessk}/bless-k{blessk}.done", sample=Illumina_process_df["sample"], blessk=config["bless_k"])
	output:
		ok = "results/{sample}/errorcorrection/bless/bless-bestk-corrected/bless-bestk-correct.ok"
	log:
		stdout = "results/{sample}/logs/bless-bestk.stdout.txt",
		stderr = "results/{sample}/logs/bless-bestk.stderr.txt"
	benchmark: "results/{sample}/benchmarks/bless-bestk.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		dir = "results/{sample}/errorcorrection/bless-bestk/",
	singularity:
		"docker://chrishah/bless:v1.02"
	threads: 2
	shadow: "minimal"
	shell:
		"""
		#read in outputs from the previous bless runs and decide on best k
		touch {output.ok}
		"""
rule bless_pe:
	input:
		forward = "results/{sample}/trimming/trim_galore/{sample}-full/{sample}.trimgalore.1.fastq.gz",
		reverse = "results/{sample}/trimming/trim_galore/{sample}-full/{sample}.trimgalore.2.fastq.gz",
		bless_runs = expand("results/{{sample}}/errorcorrection/bless/bless-k{blessk}/bless-k{blessk}.done", sample=Illumina_process_df["sample"], blessk=config["bless_k"])
	output:
		ok = "results/{sample}/errorcorrection/bless/bless-kbest-pe.done",
		forward = "results/{sample}/errorcorrection/bless/bless-kbest-corrected/bless.corrected.1.fastq.gz",
		reverse = "results/{sample}/errorcorrection/bless/bless-kbest-corrected/bless.corrected.2.fastq.gz",
	log:
		"results/{sample}/logs/bless-kbest-pe.log.txt",
	benchmark: "results/{sample}/benchmarks/bless-kbest-pe.txt"
	singularity:
		"docker://chrishah/bless:v1.02"
	params:
		wd = os.getcwd(),
		dir = "results/{sample}/errorcorrection/bless"
	resources:
		mem_gb=20
	threads: config["threads"]["bless"]
	shell:
		"""
		./bin/gather_bless_stats.sh {wildcards.sample} {input.bless_runs} > {params.dir}/bless_stats.tsv

		bestk=$(tail -n 1 {params.dir}/bless_stats.tsv | cut -f 1)
		echo -e "\\n\\n[$(date)]\tCorrecting paired end reads with bless - k: $bestk" 2>&1 | tee {log}

		#set stage
		if [[ ! -d {params.dir}/bless-kbest-corrected ]]; then mkdir -p {params.dir}/bless-kbest-corrected; else rm -rf {params.dir}/bless-kbest-corrected/bless-k$bestk-pe*; fi
		#run bless
		bless -read1 {input.forward} -read2 {input.reverse} -kmerlength $bestk -max_mem {resources.mem_gb} -load {params.dir}/bless-k$bestk/bless-k$bestk -notrim -smpthread {threads} -gzip -prefix {params.dir}/bless-kbest-corrected/bless-k$bestk-pe 2>&1 | tee {log}

		ln -s {params.wd}/{params.dir}/bless-kbest-corrected/bless-k$bestk-pe.1.corrected.fastq.gz {output.forward}
		ln -s {params.wd}/{params.dir}/bless-kbest-corrected/bless-k$bestk-pe.2.corrected.fastq.gz {output.reverse}
		touch {output.ok}
		"""
rule bless_se:
	input:
		orphans = "results/{sample}/trimming/trim_galore/{sample}-full/{sample}.trimgalore.se.fastq.gz",
		bless_runs = expand("results/{{sample}}/errorcorrection/bless/bless-k{blessk}/bless-k{blessk}.done", sample=Illumina_process_df["sample"], blessk=config["bless_k"])
	output:
		orphans = "results/{sample}/errorcorrection/bless/bless-kbest-corrected/bless.corrected.se.fastq.gz",
		ok = "results/{sample}/errorcorrection/bless/bless-kbest-se.done"
	log:
		"results/{sample}/logs/bless-kbest-se.log.txt",
	benchmark: "results/{sample}/benchmarks/bless-kbest-se.txt"
	singularity:
		"docker://chrishah/bless:v1.02"
	params:
		wd = os.getcwd(),
		dir = "results/{sample}/errorcorrection/bless"
	resources:
		mem_gb=20
	threads: config["threads"]["bless"]
	shell:
		"""
		./bin/gather_bless_stats.sh {input.bless_runs} > bless_stats.tsv

		bestk=$(tail -n 1 bless_stats.tsv | cut -f 1)
		echo -e "\\n\\n[$(date)]\tCorrecting single end reads with bless - k: $bestk" 2>&1 | tee {log}
		
		#set stage
		if [[ ! -d {params.dir}/bless-kbest-corrected ]]; then mkdir -p {params.dir}/bless-kbest-corrected; else rm -rf {params.dir}/bless-kbest-corrected/bless-k$bestk-pe*; fi
		#run bless
		bless -read {input.orphans} -kmerlength $bestk -max_mem {resources.mem_gb} -load {params.dir}/bless-k$bestk/bless-k$bestk -notrim -smpthread {threads} -gzip -prefix {params.dir}/bless-kbest-corrected/bless-k$bestk-se 2>&1 | tee {log}

		ln -s {params.wd}/{params.dir}/bless-kbest-corrected/bless-k$bestk-se.corrected.fastq.gz {output.orphans}
		touch {output.ok}
		"""


rule correct_spades:
	input:
		forward = "results/{sample}/trimming/trim_galore/{sample}-full/{sample}.trimgalore.1.fastq.gz",
		reverse = "results/{sample}/trimming/trim_galore/{sample}-full/{sample}.trimgalore.2.fastq.gz",
		orphans = "results/{sample}/trimming/trim_galore/{sample}-full/{sample}.trimgalore.se.fastq.gz",
	output:
		ok = "results/{sample}/errorcorrection/spades/spades-correct.ok",
		results = directory("results/{sample}/errorcorrection/spades/spades-corrected")
	log:
		stdout = "results/{sample}/logs/correct-spades.stdout.txt",
		stderr = "results/{sample}/logs/correct-spades.stderr.txt"
	benchmark: "results/{sample}/benchmarks/correct-spades.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		dir = "results/{sample}/errorcorrection/spades",
	singularity:
		"docker://chrishah/spades:v3.14.0"
	shadow: "minimal"
	threads: config["threads"]["spades_correct"]
	resources:
		mem_gb=90
	shell:
		"""
		echo "Host: $HOSTNAME" 1> {log.stdout} 2> {log.stderr}

		spades.py \
		--only-error-correction \
		-o {output.results} \
		-1 {input.forward} -2 {input.reverse} -s {input.orphans} \
		-t $(( {threads} - 1 )) \
		-m $(( {resources.mem_gb} - 5 )) 1>> {log.stdout} 2>> {log.stderr}

		ln -s {params.wd}/{output.results}/corrected/{wildcards.sample}.trimgalore.se.fastq.*.cor.fastq.gz {params.wd}/{output.results}/corrected/spades.corrected.se.fastq.gz
		ln -s {params.wd}/{output.results}/corrected/{wildcards.sample}.trimgalore.1.fastq.*.cor.fastq.gz {params.wd}/{output.results}/corrected/spades.corrected.1.fastq.gz
		ln -s {params.wd}/{output.results}/corrected/{wildcards.sample}.trimgalore.2.fastq.*.cor.fastq.gz {params.wd}/{output.results}/corrected/spades.corrected.2.fastq.gz
		touch {output.ok}
		"""

rule gather_illumina_corrected:
	input:
		control_illumina_ec
	output:
		ok = "results/{sample}/errorcorrection/illumina-correction.ok",
	threads: 1
	shell:
		"""
		touch {output.ok}
		"""
