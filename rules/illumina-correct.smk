#rule bless_by_k:
#	input:
#		reads = get_illumina_assembly_input,
#	output:
#		ok = "results/{sample}/errorcorrection/bless/bless.k-{blessk}.ok"
#	log:
#		stdout = "results/{sample}/logs/bless.k-{blessk}.stdout.txt",
#		stderr = "results/{sample}/logs/bless.k-{blessk}.stderr.txt"
#	benchmark: "results/{sample}/benchmarks/bless.k-{blessk}.txt"
#	params:
#		wd = os.getcwd(),
#		sample = "{sample}",
#		dir = "results/{sample}/errorcorrection/bless/",
#		k = "{blessk}"
#	threads: 90
#	shadow: "minimal"
#	shell:
#		"""
#		"""

rule correct_spades:
	input:
		forward = "results/{sample}/trimming/trim_galore/{sample}-full/{sample}.trimgalore.1.fastq.gz",
		reverse = "results/{sample}/trimming/trim_galore/{sample}-full/{sample}.trimgalore.2.fastq.gz",
		orphans = "results/{sample}/trimming/trim_galore/{sample}-full/{sample}.trimgalore.se.fastq.gz",
	output:
		ok = "results/{sample}/errorcorrection/spades/spades-correct.ok",
		results = directory("results/{sample}/errorcorrection/spades/{sample}-spades-corrected")
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
	threads: 90
	resources:
		mem_gb=750
	shell:
		"""
		echo "Host: $HOSTNAME" 1> {log.stdout} 2> {log.stderr}

		spades.py \
		--only-error-correction \
		-o {output.results} \
		-1 {input.forward} -2 {input.reverse} -s {input.orphans} \
		-t $(( {threads} - 1 )) \
		-m $(( {resources.mem_gb} - 5 )) 1>> {log.stdout} 2>> {log.stderr}

		touch {output.ok}
		"""
