#import pandas as pd
#import sys

#genome_data = pd.read_table(config["samples"]).set_index("assembly_name", drop=False)
#print(genome_data.loc["CAE_full", ["assembly_path"]].to_list())
#print(genome_data.index.tolist())
#def get_assembly_path(wildcards):
#	print ("wildcards:", wildcards)
#        return genome_data.loc[wildcards.assembly_name, ["assembly_path"]].to_list()

#sample_data = pd.read_table(config["samples"]).set_index("sample", drop=False)

#def get_FR(wildcards):
#	return sample_data.loc[wildcards.sample, ["FR"]].to_list()

#def get_RR(wildcards):
#	return sample_data.loc[wildcards.sample, ["RR"]].to_list()

#rule all:
#	input:
##		expand("results/{samp}/mapping/{samp}_sam2bam.done", samp=sample_data.index.tolist()),
##		expand("results/{samp}/mapping/{samp}.duprmvd.done", samp=sample_data.index.tolist())
##		expand("results/genome_index/genome_index.{assembly_name}.done", assembly_name=genome_data.index.tolist()),
#		expand("results/{samp}/mapping/{assembly_name}/{samp}_calculate_depth.done", zip, samp=sample_data.index.tolist(), assembly_name=genome_data.index.tolist()),
#		expand("results/{samp}/mapping/{assembly_name}/{samp}_indexing_durmvd.done", zip, samp=sample_data.index.tolist(), assembly_name=genome_data.index.tolist())

#rule allclean:
#	input: 
#		checksam2bam = expand("results/{samp}/mapping/{samp}_sam2bam.done", samp=sample_data.index.tolist()),
#		checkdepth = expand("results/{samp}/mapping/{samp}_calculate_depth.done", samp=sample_data.index.tolist()),
#		checkindexing = expand("results/{samp}/mapping/{samp}_indexing_durmvd.done", samp=sample_data.index.tolist()),
#		checkrmbam = expand("results/{samp}/mapping/{samp}.rmbam.done", samp=sample_data.index.tolist())
#rule rmbam:
#	input:
#		expand("results/{samp}/mapping/{samp}.duprmvd.done", samp=sample_data.index.tolist())
#	output:
#		expand("results/{samp}/mapping/{samp}.rmbam.done", samp=sample_data.index.tolist())
#	params:
#		bam = expand("results/{samp}/mapping/{samp}.bam", samp=sample_data.index.tolist())
#		
#	threads: 1
#	shell:
#		"""
#		rm {params.bam}
#		touch {output}
#		"""

rule eva_bowtie2_index:
	input:
		assembly = "results/{sample}/assembly_evaluation/assemblies/{combination}.fasta",
#		assembly = get_assembly_path,
	singularity:
		"docker://reslp/bowtie2:2.3.5" 
	output: 
		"results/{sample}/assembly_evaluation/genome_index/genome_index.{combination}.done"
	params:
		wd = os.getcwd(),
		assemblyID = "{combination}",
		outdir = "results/{sample}/assembly_evaluation/genome_index"
	shell:
		"""
		bowtie2-build {input.assembly} {params.outdir}/{params.assemblyID}.index -q
		touch {output}
		"""

rule eva_bowtie2_mapping:
	input:
		forward_reads = "results/{sample}/trimming/trimgalore/{sample}-full/{sample}.trimgalore.1.fastq.gz",
		reverse_reads = "results/{sample}/trimming/trimgalore/{sample}-full/{sample}.trimgalore.2.fastq.gz",
		check = rules.eva_bowtie2_index.output
	output:
		"results/{sample}/assembly_evaluation/mapping/{combination}/{sample}.{combination}.mapping.done"
	params:
		sample = "{sample}",
		assemblyID = "{combination}",
		indir = "results/{sample}/assembly_evaluation/genome_index",
		outdir = "results/{sample}/assembly_evaluation/mapping/{combination}",
		additional_params = config["evaluate_assemblies"]["bowtie2"]["additional_parameters"]
	singularity:
		"docker://reslp/bowtie2:2.3.5"
	threads: config["threads"]["bowtie2"]
	shell:
		"""
		bowtie2 -1 {input.forward_reads} -2 {input.reverse_reads} -p {threads} -q -x {params.indir}/{params.assemblyID}.index -S {params.outdir}/{params.sample}.sam {params.additional_params} 
		touch {output}
		"""

rule eva_conversion_samtools:
	input:
		rules.eva_bowtie2_mapping.output
	output:
		"results/{sample}/assembly_evaluation/mapping/{combination}/{sample}.sam2bam.done"
	params:
		sample = "{sample}",
		outdir = "results/{sample}/assembly_evaluation/mapping/{combination}"
	singularity:
		"docker://reslp/samtools:1.9"
	threads: config["threads"]["samtools"]
	shell:
		"""
		samtools view -bS {params.outdir}/{params.sample}.sam -o {params.outdir}/{params.sample}.bam -@ {threads}
		samtools sort -o {params.outdir}/{params.sample}.bam {params.outdir}/{params.sample}.bam -@ {threads}
		samtools index {params.outdir}/{params.sample}.bam -@ {threads}
		rm {params.outdir}/{params.sample}.sam 
		touch {output} 
		"""

rule eva_remove_duplicates:
	input:
		rules.eva_conversion_samtools.output
	output:
		bam = "results/{sample}/assembly_evaluation/mapping/{combination}/{sample}.sorted.duprmvd.bam",
		done = "results/{sample}/assembly_evaluation/mapping/{combination}/{sample}.duprmvd.done"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		outdir = "results/{sample}/assembly_evaluation/mapping/{combination}"
	singularity:
		"docker://broadinstitute/picard:2.20.6"
	log:
		stdout = "results/{sample}/logs/remove_duplicates.{combination}.stdout.txt",
		stderr = "results/{sample}/logs/remove_duplicates.{combination}.stderr.txt"
	threads: 1
	shell:
		"""
		cd {params.outdir}
		java -jar /usr/picard/picard.jar MarkDuplicates \
		INPUT={params.sample}.bam OUTPUT={params.sample}.sorted.duprmvd.bam \
		METRICS_FILE={params.sample}.sorted.duprmvd.metrics REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 1> {params.wd}/{log.stdout} 2> {params.wd}/{log.stderr}
		rm {params.sample}.bam
		touch {params.wd}/{output.done}
		"""

#rule eva_calculate_depth:
#	input:
#		rules.eva_remove_duplicates.output
#	output:
#		"results/{sample}/assembly_evaluation/mapping/{combination}/{sample}.calculate_depth.done"
#	params:
#		wd = os.getcwd(),
#		sample = "{sample}",
#		outdir = "results/{sample}/assembly_evaluation/mapping/{combination}"
#	singularity:
#		"docker://reslp/samtools:1.9"
#	threads: 1
#	shell:
#		"""
#		samtools depth {params.outdir}/{params.sample}.sorted.duprmvd.bam | cut -f 3 | sort -n | uniq -c > {params.outdir}/{params.sample}.depth.hist.txt
#		touch {output}
#		"""

rule eva_index_duprmvd:
	input:
		rules.eva_remove_duplicates.output
	output:
		bai = "results/{sample}/assembly_evaluation/mapping/{combination}/{sample}.sorted.duprmvd.bam.bai",
		done = "results/{sample}/assembly_evaluation/mapping/{combination}/{sample}.indexing_durmvd.done"
	params:
		sample = "{sample}",
		outdir = "results/{sample}/assembly_evaluation/mapping/{combination}"
	singularity:
		"docker://reslp/samtools:1.9"
	threads: config["threads"]["samtools_index"]
	shell:
		"""
		samtools index {params.outdir}/{params.sample}.sorted.duprmvd.bam -@ {threads}
		touch {output}
		"""
