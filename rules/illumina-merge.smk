#if merging is not to be done at 'corrected' stage, reset the list of correctors to None
if not "corrected" in config["illumina_merge"]["merge_at_stage"]:
	correct_list_for_merging = ["None"]
else:
	correct_list_for_merging = list(config["illumina_correction"])

#if merging should be done at 'trimmed' stage, add None to whatever correctors are in the list at this stage
if "trimmed" in config["illumina_merge"]["merge_at_stage"]:
	trimmed_list_for_merging = list(config["illumina_trimming"])
	if not "None" in config["illumina_correction"]:
		correct_list_for_merging.append(None)

merge_list_for_merging = list(config["illumina_merge"]["merger"])
#print("FINAL - Illumina_trimming for merging: "+str(trimmed_list_for_merging))
#print("FINAL - Illumina_correction for merging: "+str(correct_list_for_merging))
#print("FINAL - Illumina_merging for merging: "+str(merge_list_for_merging))

rule setup_usearch:
	output:
		"bin/usearch"
	shadow: "minimal"
	shell:
		"""
		wget https://www.drive5.com/downloads/usearch11.0.667_i86linux32.gz
		gunzip -v $(find ./ -name "*gz")
		chmod a+x $(find ./ -name "*linux32")
		mv $(find ./ -name "*linux32") {output}
		"""

rule usearch:
	input:
		reads = to_merge,
		usearch = rules.setup_usearch.output,
		script = "bin/usearch_mergepairs.sh"
	output:
		merged = "results/{sample}/read-merging/usearch/{trimmer}-{corrector}/{sample}.merged.fastq.gz",
		nm1 = "results/{sample}/read-merging/usearch/{trimmer}-{corrector}/{sample}.notmerged.1.fastq.gz",
		nm2 = "results/{sample}/read-merging/usearch/{trimmer}-{corrector}/{sample}.notmerged.2.fastq.gz",
		ok = "results/{sample}/read-merging/usearch/{trimmer}-{corrector}/{trimmer}-{corrector}-usearch.ok"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		batchsize = "4000000"
	threads: config["threads"]["mergepairs_usearch"]
	singularity:
		"docker://chrishah/usearch-docker-onbuild:v012020"
	log:
		stdout = "results/{sample}/logs/usearch.{sample}.{trimmer}.{corrector}.stdout.txt",
		stderr = "results/{sample}/logs/usearch.{sample}.{trimmer}.{corrector}.stderr.txt"
	shadow: "minimal"
	shell:
		"""
		mkdir tmp
		export TMPDIR=$(pwd)/tmp
		export PATH=$PATH:$(pwd)/bin
		export HOME=$(pwd)

		{input.script} {input.reads[0]} {input.reads[1]} {params.sample} {threads} {params.batchsize} 1> {log.stdout} 2> {log.stderr}
		cp {wildcards.sample}_1.nm.fastq.gz {output.nm1}
		cp {wildcards.sample}_2.nm.fastq.gz {output.nm2}
		cp {wildcards.sample}.merged.fastq.gz {output.merged}
		
		echo -e "###\\n$(date)\\tLogs from individual usearch runs:\\n" >> {log.stdout}
		cat merging.log >> {log.stdout}
		touch {output.ok}
		"""
rule flash:
	input:
		to_merge
	output:
		"results/{sample}/read-merging/flash/{trimmer}-{corrector}/{trimmer}-{corrector}-flash.ok"
	shell:
		"""
		echo {input}
		touch {output}
		"""
rule merger3:
	output:
		"results/{sample}/read-merging/merger3/{trimmer}-{corrector}/{trimmer}-{corrector}-merger3.ok"
	shell:
		"""
		echo {input}
		touch {output}
		"""
	
rule gather_illumina_merged:
	input:
#		control_illumina_me
		expand("results/{{sample}}/read-merging/{merger}/{trimmer}-{corrector}/{trimmer}-{corrector}-{merger}.ok",  
			trimmer=trimmed_list_for_merging,
			corrector=correct_list_for_merging,
			merger=merge_list_for_merging)
	output:
		ok = "results/{sample}/read-merging/illumina-merging.ok",
	threads: 1
	shell:
		"""
		touch {output.ok}
		"""
