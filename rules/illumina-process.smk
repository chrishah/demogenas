rule pre_prepare_fastq:
	input:
		forward = get_raw_f_fastqs,
		reverse = get_raw_r_fastqs
	output:
		forward = "results/{sample}/Illumina/raw_reads/from_fastq/{lib}/{sample}.{lib}.raw.1.fastq.gz",
		reverse = "results/{sample}/Illumina/raw_reads/from_fastq/{lib}/{sample}.{lib}.raw.2.fastq.gz"
	shell:
		"""
		ln -s $(pwd)/{input.forward} $(pwd)/{output.forward}
		ln -s $(pwd)/{input.reverse} $(pwd)/{output.reverse}
		"""
rule pre_sort_bam:
	input:
		get_bam
	output:
		ok = "results/{sample}/Illumina/raw_reads/from_bam/{lib}/{sample}.{lib}.sorting.ok",
	log:
		stdout = "results/{sample}/logs/sort_bam.{lib}.stdout.txt",
		stderr = "results/{sample}/logs/sort_bam.{lib}.stderr.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		lib = "{lib}",
		targetdir = "results/{sample}/Illumina/raw_reads/from_bam/{lib}/"
	threads: config["threads"]["samtools"]
	singularity: "docker://reslp/samtools:1.11"
	resources:
		mem_gb = config["threads"]["samtools"]
	shadow: "minimal"
	shell:
		"""
		samtools sort -@ $(( {threads} - 1 )) -m 1G -n {input} -o {params.sample}.{params.lib}.sorted.bam 1> {log.stdout} 2> {log.stderr}
		mv *.bam {params.wd}/{params.targetdir}
		touch {output.ok}
		"""
rule pre_bam2fastq:
	input:
		rules.pre_sort_bam.output
	output:
		forward = "results/{sample}/Illumina/raw_reads/from_bam/{lib}/{sample}.{lib}.raw.1.fastq.gz",
		reverse = "results/{sample}/Illumina/raw_reads/from_bam/{lib}/{sample}.{lib}.raw.2.fastq.gz"
	log:
		stdout = "results/{sample}/logs/bam2fastq.{lib}.stdout.txt",
		stderr = "results/{sample}/logs/bam2fastq.{lib}.stderr.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		lib = "{lib}",
		targetdir = "results/{sample}/Illumina/raw_reads/from_bam/{lib}/"
	threads: 2
	singularity: "docker://reslp/bedtools:2.29.2"
	shadow: "minimal"
	shell:
		"""
		bedtools bamtofastq -i {params.wd}/{params.targetdir}/{params.sample}.{params.lib}.sorted.bam -fq {params.targetdir}/{params.sample}.{params.lib}.raw.1.fastq -fq2 {params.targetdir}/{params.sample}.{params.lib}.raw.2.fastq 1> {log.stdout} 2> {log.stderr}
		gzip -v {params.targetdir}*.fastq 1>> {log.stdout} 2>> {log.stderr}
		rm {params.wd}/{params.targetdir}/*.bam
		"""		


rule tri_trimgalore:
	input:
		forward = input_for_trimgalore_f,
		reverse = input_for_trimgalore_r,
	params:
		wd = os.getcwd(),
		lib = "{lib}",
		sample = "{sample}",
#		adapters = "--illumina"
	singularity:
		"docker://chrishah/trim_galore:0.6.0"
	log:
		stdout = "results/{sample}/logs/trimgalore.{lib}.stdout.txt",
		stderr = "results/{sample}/logs/trimgalore.{lib}.stderr.txt"
	output:
		ok = "results/{sample}/trimming/trimgalore/{lib}/{sample}.{lib}.status.ok",
		f_trimmed = "results/{sample}/trimming/trimgalore/{lib}/{sample}.{lib}.1.fastq.gz",
		r_trimmed = "results/{sample}/trimming/trimgalore/{lib}/{sample}.{lib}.2.fastq.gz",
		f_orphans = "results/{sample}/trimming/trimgalore/{lib}/{sample}.{lib}.unpaired.1.fastq.gz",
		r_orphans = "results/{sample}/trimming/trimgalore/{lib}/{sample}.{lib}.unpaired.2.fastq.gz"
	shadow: "minimal"
	threads: 4
	shell:
		"""
		trim_galore \
		--paired --length 70 -r1 71 -r2 71 --retain_unpaired --stringency 2 --quality 30 \
		{input.forward} {input.reverse} 1> {log.stdout} 2> {log.stderr}

		mv $(find ./ -name "*_val_1.fq.gz") {output.f_trimmed}
		mv $(find ./ -name "*_val_2.fq.gz") {output.r_trimmed}
	
		if [[ -f $(find ./ -name "*_unpaired_1.fq.gz") ]]; then mv $(find ./ -name "*_unpaired_1.fq.gz") {output.f_orphans}; else touch {output.f_orphans}; fi
		if [[ -f $(find ./ -name "*_unpaired_2.fq.gz") ]]; then mv $(find ./ -name "*_unpaired_2.fq.gz") {output.r_orphans}; else touch {output.r_orphans}; fi

		touch {output.ok}

		"""
rule tri_gather_short_trimmed_by_lib:
	input:
		forward = input_for_clean_trim_fp,
		reverse = input_for_clean_trim_rp,
		forward_orphans = input_for_clean_trim_fo,
		reverse_orphans = input_for_clean_trim_ro
	params:
		wd = os.getcwd(),
		sample = "{sample}",
	singularity:
		"docker://chrishah/trim_galore:0.6.0"
	output:
		ok = "results/{sample}/trimming/trimgalore/{sample}-full/{sample}.cat.status.ok",
		f_trimmed = "results/{sample}/trimming/trimgalore/{sample}-full/{sample}.trimgalore.1.fastq.gz",
		r_trimmed = "results/{sample}/trimming/trimgalore/{sample}-full/{sample}.trimgalore.2.fastq.gz",
		orphans = "results/{sample}/trimming/trimgalore/{sample}-full/{sample}.trimgalore.se.fastq.gz",
	shadow: "minimal"
	threads: 2
	shell:
		"""
		if [ $(echo {input.forward} | wc -w) -gt 1 ]
		then
			cat {input.forward} > {output.f_trimmed}
			cat {input.reverse} > {output.r_trimmed}
		else
			ln -s ../../../../../{input.forward} {output.f_trimmed}
			ln -s ../../../../../{input.reverse} {output.r_trimmed}
		fi
		cat {input.forward_orphans} {input.reverse_orphans} > {output.orphans}
		touch {output.ok}
		"""

rule tri_gather_short_direct_by_lib:
	input:
		forward = input_for_clean_direct_fp,
		reverse = input_for_clean_direct_rp,
	params:
		wd = os.getcwd(),
		sample = "{sample}",
	output:
		ok = "results/{sample}/trimming/None/{sample}-full/{sample}.cat.status.ok",
		f_direct = "results/{sample}/trimming/None/{sample}-full/{sample}.None.1.fastq.gz",
		r_direct = "results/{sample}/trimming/None/{sample}-full/{sample}.None.2.fastq.gz",
		o_direct = "results/{sample}/trimming/None/{sample}-full/{sample}.None.se.fastq.gz",
	shadow: "minimal"
	threads: 2
	shell:
		"""
		if [ $(echo {input.forward} | wc -w) -gt 1 ]
		then
			cat {input.forward} > {output.f_direct}
			cat {input.reverse} > {output.r_direct}
		else
			ln -s ../../../../../{input.forward} {output.f_direct}
			ln -s ../../../../../{input.reverse} {output.r_direct}
		fi
		touch {output.o_direct}
		touch {output.ok}
		"""

if config["skip_trimming"] == "yes":
	rule eva_fastqc_direct:
		input:
			f_paired = input_for_trimgalore_f,
			r_paired = input_for_trimgalore_r,
		params:
			wd = os.getcwd(),
			lib = "{lib}",
			sample = "{sample}",
			dir = "results/{sample}/trimming/direct/{lib}"
		singularity:
			"docker://chrishah/trim_galore:0.6.0"
		log:
			stdout = "results/{sample}/logs/fastqc_direct.{sample}.{lib}.stdout.txt",
			stderr = "results/{sample}/logs/fastqc_direct.{sample}.{lib}.stderr.txt"
		output:
			ok = "results/{sample}/trimming/direct/{lib}/{sample}.{lib}.fastqc.status.ok",
		shadow: "minimal"
		threads: 2
		shell:
			"""
			for f in $(echo "{input}"); do if [[ -s $f ]]; then fastqc -t {threads} -o ./ $f; fi; done 1> {log.stdout} 2> {log.stderr}
			mv *.zip *.html {params.wd}/{params.dir}/
			touch {output}
			
			"""
	rule eva_stats_direct:
		input:
			lambda wildcards: expand(rules.eva_fastqc_direct.output.ok, sample=wildcards.sample, lib=list(set(unitdict[wildcards.sample]) & set(Illumina_process_df["lib"].tolist())))
		singularity:
			"docker://chrishah/r-docker:latest"
		log:
			stdout = "results/{sample}/logs/readstats.{sample}.stdout.txt",
			stderr = "results/{sample}/logs/readstats.{sample}.stderr.txt"
		output: 
			stats = "results/{sample}/trimming/direct/{sample}.readstats.txt"
		threads: 2
		shadow: "shallow"
		shell:
			"""
			echo -e "$(date)\tGetting read stats" 1> {log.stdout}
			bin/get_readstats_from_fastqc.sh $(echo {input} | sed 's/ /,/g' | sed 's/.fastqc.status.ok//g') 1> {output.stats} 2> {log.stderr}
			stats=$(cat {output.stats})
			echo -e "Cummulative length: $(cat {output.stats} | cut -d " " -f 1)" 1>> {log.stdout} 2>> {log.stderr}
			echo -e "Average read length: $(cat {output.stats} | cut -d " " -f 2)\\n" 1>> {log.stdout} 2>> {log.stderr}
			echo -e "$(date)\tDone!" 1>> {log.stdout}
			"""
	rule eva_kmc_direct:
		input:
			f_paired = lambda wildcards: expand("results/{{sample}}/Illumina/raw_reads/from_fastq/{lib}/{{sample}}.{lib}.raw.1.fastq.gz", sample=wildcards.sample, lib=list(set(unitdict[wildcards.sample]) & set(Illumina_process_df["lib"].tolist()))),
			r_paired = lambda wildcards: expand("results/{{sample}}/Illumina/raw_reads/from_fastq/{lib}/{{sample}}.{lib}.raw.2.fastq.gz", sample=wildcards.sample, lib=list(set(unitdict[wildcards.sample]) & set(Illumina_process_df["lib"].tolist()))),
		params:
			sample = "{sample}",
			k = "{k}",
			max_mem_in_GB = config["kmc"]["max_mem_in_GB"],
			mincount = config["kmc"]["mincount"],
			maxcount = config["kmc"]["maxcount"],
			maxcounter = config["kmc"]["maxcounter"],
			nbin = 64,
		threads: config["threads"]["kmc"]
		singularity:
			"docker://chrishah/kmc3-docker:v3.0"
		log:
			stdout = "results/{sample}/logs/kmc.{sample}.k{k}.stdout.txt",
			stderr = "results/{sample}/logs/kmc.{sample}.k{k}.stderr.txt"
		output: 
	#		pre = "results/{sample}/kmc/{sample}.k{k}.kmc_pre",
	#		suf = "results/{sample}/kmc/{sample}.k{k}.kmc_suf",
			hist = "results/{sample}/kmc/{sample}.k{k}.histogram.txt",
			genomescope = "results/{sample}/kmc/{sample}.k{k}.histogram.genomescope.txt"
		shadow: "minimal"
		shell:
			"""
			echo -e "$(date)\tStarting kmc"
			## add all files that are not empty to file of filenames
			for f in $(echo "{input}"); do if [[ -s $f ]]; then echo "$f"; fi; done > fastqs.txt
			mkdir {params.sample}.db
			kmc -k{params.k} -m$(( {params.max_mem_in_GB} - 2 )) -v -sm -ci{params.mincount} -cx{params.maxcount} -cs{params.maxcounter} -n{params.nbin} -t$(( {threads} - 1 )) @fastqs.txt {params.sample} {params.sample}.db 1>> {log.stdout} 2>> {log.stderr}
			#kmc_tools histogram {params.sample} -ci{params.mincount} -cx{params.maxcount} {output.hist}
			kmc_tools histogram {params.sample} -ci{params.mincount} {output.hist} 1> /dev/null 2>> {log.stderr}
			cat {output.hist} | awk '{{print $1" "$2}}' > {output.genomescope}
			"""
	rule eva_plot_k_hist_direct:
		input:
			hist = rules.eva_kmc_direct.output.hist,
			stats = rules.eva_stats_direct.output.stats
		output:
			full = "results/{sample}/plots/{sample}-k{k}-distribution-full.pdf",		
		params:
			sample = "{sample}",
			k = "{k}",
			script = "bin/plot.freq.in.R"
		singularity:
			"docker://chrishah/r-docker:latest"
		log:
			stdout = "results/{sample}/logs/plotkmerhist.{sample}.k{k}.stdout.txt",
			stderr = "results/{sample}/logs/plotkmerhist.{sample}.k{k}.stderr.txt"
		shadow: "shallow"
		shell:
			"""
			stats=$(cat {input.stats})
			echo -e "Cummulative length: $(echo -e "$stats" | cut -d " " -f 1)"
			echo -e "Average read length: $(echo -e "$stats" | cut -d " " -f 2)"
			Rscript {params.script} {input.hist} {params.sample} {params.k} $stats 1> {log.stdout} 2> {log.stderr}
			cp {params.sample}-k{params.k}-distribution* results/{params.sample}/plots/
			"""

else:
	rule eva_fastqc_raw:
		input:
			forward = input_for_trimgalore_f,
			reverse = input_for_trimgalore_r,
		params:
			wd = os.getcwd(),
			lib = "{lib}",
			sample = "{sample}",
		singularity:
			"docker://chrishah/trim_galore:0.6.0"
		log:
			stdout = "results/{sample}/logs/fastqc_raw.{sample}.{lib}.stdout.txt",
			stderr = "results/{sample}/logs/fastqc_raw.{sample}.{lib}.stderr.txt"
		output:
			ok = "results/{sample}/read_qc/fastqc_raw/{lib}/{sample}.{lib}.status.ok",
		shadow: "minimal"
		threads: 2
		shell:
			"""
			fastqc -o ./ {input} 1> {log.stdout} 2> {log.stderr}
			mv *.zip *.html {params.wd}/results/{params.sample}/read_qc/fastqc_raw/{params.lib}/
			touch {output}
			"""

	rule eva_fastqc_trimmed:
		input:
			f_paired = "results/{sample}/trimming/trimgalore/{lib}/{sample}.{lib}.1.fastq.gz",
			r_paired = "results/{sample}/trimming/trimgalore/{lib}/{sample}.{lib}.2.fastq.gz",
			f_unpaired = "results/{sample}/trimming/trimgalore/{lib}/{sample}.{lib}.unpaired.1.fastq.gz",
			r_unpaired = "results/{sample}/trimming/trimgalore/{lib}/{sample}.{lib}.unpaired.2.fastq.gz"
		params:
			wd = os.getcwd(),
			lib = "{lib}",
			sample = "{sample}",
		singularity:
			"docker://chrishah/trim_galore:0.6.0"
		log:
			stdout = "results/{sample}/logs/fastqc_trimgalore.{sample}.{lib}.stdout.txt",
			stderr = "results/{sample}/logs/fastqc_trimgalore.{sample}.{lib}.stderr.txt"
		output:
			ok = "results/{sample}/trimming/trimgalore/{lib}/{sample}.{lib}.fastqc.status.ok",
		shadow: "minimal"
		threads: 2
		shell:
			"""
			for f in $(echo "{input}"); do if [[ -s $f ]]; then fastqc -t {threads} -o ./ $f; fi; done 1> {log.stdout} 2> {log.stderr}
			mv *.zip *.html {params.wd}/results/{params.sample}/trimming/trimgalore/{params.lib}/
			touch {output}
			
			"""
	rule eva_stats:
		input:
			lambda wildcards: expand(rules.eva_fastqc_trimmed.output.ok, sample=wildcards.sample, lib=list(set(unitdict[wildcards.sample]) & set(Illumina_process_df["lib"].tolist())))
		singularity:
			"docker://chrishah/r-docker:latest"
		log:
			stdout = "results/{sample}/logs/readstats.{sample}.stdout.txt",
			stderr = "results/{sample}/logs/readstats.{sample}.stderr.txt"
		output: 
			stats = "results/{sample}/trimming/trimgalore/{sample}.readstats.txt"
		threads: 2
		shadow: "shallow"
		shell:
			"""
			echo -e "$(date)\tGetting read stats" 1> {log.stdout}
			bin/get_readstats_from_fastqc.sh $(echo {input} | sed 's/ /,/g' | sed 's/.fastqc.status.ok//g') 1> {output.stats} 2> {log.stderr}
			stats=$(cat {output.stats})
			echo -e "Cummulative length: $(cat {output.stats} | cut -d " " -f 1)" 1>> {log.stdout} 2>> {log.stderr}
			echo -e "Average read length: $(cat {output.stats} | cut -d " " -f 2)\\n" 1>> {log.stdout} 2>> {log.stderr}
			echo -e "$(date)\tDone!" 1>> {log.stdout}
			"""
	rule eva_kmc:
		input:
			f_paired = lambda wildcards: expand("results/{{sample}}/trimming/trimgalore/{lib}/{{sample}}.{lib}.1.fastq.gz", sample=wildcards.sample, lib=list(set(unitdict[wildcards.sample]) & set(Illumina_process_df["lib"].tolist()))),
			r_paired = lambda wildcards: expand("results/{{sample}}/trimming/trimgalore/{lib}/{{sample}}.{lib}.2.fastq.gz", sample=wildcards.sample, lib=list(set(unitdict[wildcards.sample]) & set(Illumina_process_df["lib"].tolist()))),
			f_unpaired = lambda wildcards: expand("results/{{sample}}/trimming/trimgalore/{lib}/{{sample}}.{lib}.unpaired.1.fastq.gz", sample=wildcards.sample, lib=list(set(unitdict[wildcards.sample]) & set(Illumina_process_df["lib"].tolist()))),
			r_unpaired = lambda wildcards: expand("results/{{sample}}/trimming/trimgalore/{lib}/{{sample}}.{lib}.unpaired.2.fastq.gz", sample=wildcards.sample, lib=list(set(unitdict[wildcards.sample]) & set(Illumina_process_df["lib"].tolist()))),
		params:
			sample = "{sample}",
			k = "{k}",
			max_mem_in_GB = config["kmc"]["max_mem_in_GB"],
			mincount = config["kmc"]["mincount"],
			maxcount = config["kmc"]["maxcount"],
			maxcounter = config["kmc"]["maxcounter"],
			nbin = 64,
		threads: config["threads"]["kmc"]
		singularity:
			"docker://chrishah/kmc3-docker:v3.0"
		log:
			stdout = "results/{sample}/logs/kmc.{sample}.k{k}.stdout.txt",
			stderr = "results/{sample}/logs/kmc.{sample}.k{k}.stderr.txt"
		output: 
	#		pre = "results/{sample}/kmc/{sample}.k{k}.kmc_pre",
	#		suf = "results/{sample}/kmc/{sample}.k{k}.kmc_suf",
			hist = "results/{sample}/kmc/{sample}.k{k}.histogram.txt",
			genomescope = "results/{sample}/kmc/{sample}.k{k}.histogram.genomescope.txt"
		shadow: "minimal"
		shell:
			"""
			echo -e "$(date)\tStarting kmc"
			## add all files that are not empty to file of filenames
			for f in $(echo "{input}"); do if [[ -s $f ]]; then echo "$f"; fi; done > fastqs.txt
			mkdir {params.sample}.db
			kmc -k{params.k} -m$(( {params.max_mem_in_GB} - 2 )) -v -sm -ci{params.mincount} -cx{params.maxcount} -cs{params.maxcounter} -n{params.nbin} -t$(( {threads} - 1 )) @fastqs.txt {params.sample} {params.sample}.db 1>> {log.stdout} 2>> {log.stderr}
			#kmc_tools histogram {params.sample} -ci{params.mincount} -cx{params.maxcount} {output.hist}
			kmc_tools histogram {params.sample} -ci{params.mincount} {output.hist} 1> /dev/null 2>> {log.stderr}
			cat {output.hist} | awk '{{print $1" "$2}}' > {output.genomescope}
			"""
	rule eva_plot_k_hist:
		input:
			hist = rules.eva_kmc.output.hist,
			stats = rules.eva_stats.output.stats
		output:
			full = "results/{sample}/plots/{sample}-k{k}-distribution-full.pdf",		
		params:
			sample = "{sample}",
			k = "{k}",
			script = "bin/plot.freq.in.R"
		singularity:
			"docker://chrishah/r-docker:latest"
		log:
			stdout = "results/{sample}/logs/plotkmerhist.{sample}.k{k}.stdout.txt",
			stderr = "results/{sample}/logs/plotkmerhist.{sample}.k{k}.stderr.txt"
		shadow: "shallow"
		shell:
			"""
			stats=$(cat {input.stats})
			echo -e "Cummulative length: $(echo -e "$stats" | cut -d " " -f 1)"
			echo -e "Average read length: $(echo -e "$stats" | cut -d " " -f 2)"
			Rscript {params.script} {input.hist} {params.sample} {params.k} $stats 1> {log.stdout} 2> {log.stderr}
			cp {params.sample}-k{params.k}-distribution* results/{params.sample}/plots/
			"""

# these parts only work when data is trimmed as part of the pipeline so far

rule fil_kmc_create_db:
	input:
		f_paired = lambda wildcards: expand("results/{{sample}}/trimming/trimgalore/{lib}/{{sample}}.{lib}.1.fastq.gz", sample=wildcards.sample, lib=list(set(unitdict[wildcards.sample]) & set(Illumina_process_df["lib"].tolist()))),
		r_paired = lambda wildcards: expand("results/{{sample}}/trimming/trimgalore/{lib}/{{sample}}.{lib}.2.fastq.gz", sample=wildcards.sample, lib=list(set(unitdict[wildcards.sample]) & set(Illumina_process_df["lib"].tolist()))),
		f_unpaired = lambda wildcards: expand("results/{{sample}}/trimming/trimgalore/{lib}/{{sample}}.{lib}.unpaired.1.fastq.gz", sample=wildcards.sample, lib=list(set(unitdict[wildcards.sample]) & set(Illumina_process_df["lib"].tolist()))),
		r_unpaired = lambda wildcards: expand("results/{{sample}}/trimming/trimgalore/{lib}/{{sample}}.{lib}.unpaired.2.fastq.gz", sample=wildcards.sample, lib=list(set(unitdict[wildcards.sample]) & set(Illumina_process_df["lib"].tolist()))),
	params:
		sample = "{sample}",
		k = "{filter_k}",
		max_mem_in_GB = config["kmc"]["max_mem_in_GB"],
		mincount = config["kmc"]["mincount"],
		maxcount = config["kmc"]["maxcount"],
		maxcounter = config["kmc"]["maxcounter"],
		nbin = 64,
	threads: config["threads"]["kmc"]
	singularity:
		"docker://chrishah/kmc3-docker:v3.0"
	log:
		stdout = "results/{sample}/logs/kmc_db.{sample}.k{filter_k}.stdout.txt",
		stderr = "results/{sample}/logs/kmc_db.{sample}.k{filter_k}.stderr.txt"
	output: 
		db = directory("results/{sample}/kmc/filtered/{sample}.k{filter_k}.kmers"),
	shadow: "minimal"
	shell:
		"""
		echo -e "$(date)\tStarting kmc"
		echo "{input}" | sed 's/ /\\n/g' > fastqs.txt
		mkdir -p {output}/db
		kmc -k{params.k} -m$(( {params.max_mem_in_GB} - 2 )) -v -sm -ci{params.mincount} -cx{params.maxcount} -cs{params.maxcounter} -n{params.nbin} -t$(( {threads} - 1 )) @fastqs.txt {output}/{params.sample} {output}/db 1>> {log.stdout} 2>> {log.stderr}
		"""

rule fil_kmc_filter_forw:
	input:
		db = rules.fil_kmc_create_db.output,
		f_paired = lambda wildcards: expand("results/{{sample}}/trimming/trimgalore/{lib}/{{sample}}.{lib}.1.fastq.gz", sample=wildcards.sample, lib=list(set(unitdict[wildcards.sample]) & set(Illumina_process_df["lib"].tolist()))),
#		r_paired = lambda wildcards: expand("results/{{sample}}/trimming/trimgalore/{lib}/{{sample}}.{lib}.2.fastq.gz", sample=wildcards.sample, lib=unitdict[wildcards.sample]),
#		f_unpaired = lambda wildcards: expand("results/{{sample}}/trimming/trimgalore/{lib}/{{sample}}.{lib}.unpaired.1.fastq.gz", sample=wildcards.sample, lib=unitdict[wildcards.sample]),
#		r_unpaired = lambda wildcards: expand("results/{{sample}}/trimming/trimgalore/{lib}/{{sample}}.{lib}.unpaired.2.fastq.gz", sample=wildcards.sample, lib=unitdict[wildcards.sample]),
	params:
		sample = "{sample}",
		k = "{filter_k}",
		max_mem_in_GB = config["kmc"]["max_mem_in_GB"],
		mincov = "{mincov}",
		minprop = "{minprop}",
		maxcount = config["kmc"]["maxcount"],
		maxcounter = config["kmc"]["maxcounter"],
		nbin = 64,
	threads: config["threads"]["kmc"]
	singularity:
		"docker://chrishah/kmc3-docker:v3.0"
	log:
		stdout = "results/{sample}/logs/kmc_filter.{sample}.k{filter_k}-{mincov}-{minprop}.stdout.txt",
		stderr = "results/{sample}/logs/kmc_filter.{sample}.k{filter_k}-{mincov}-{minprop}.stderr.txt"
	output: 
		filtered = "results/{sample}/kmc/filtered/{filter_k}-{mincov}-{minprop}/{sample}.k{filter_k}.{mincov}.{minprop}.filtered.forw.fastq",
	shadow: "minimal"
	shell:
		"""
		echo -e "$(date)\tStarting kmc filtering"
		echo "{input}" | sed 's/ /\\n/g' | tail -n +2 | sort > fastqs.txt
		kmc_tools filter {input.db}/{wildcards.sample} -ci{params.mincov} @fastqs.txt -ci{params.minprop} {output.filtered} 1>> {log.stdout} 2>> {log.stderr}
		"""
rule fil_kmc_filter_reve:
	input:
		db = rules.fil_kmc_create_db.output,
		r_paired = lambda wildcards: expand("results/{{sample}}/trimming/trimgalore/{lib}/{{sample}}.{lib}.2.fastq.gz", sample=wildcards.sample, lib=list(set(unitdict[wildcards.sample]) & set(Illumina_process_df["lib"].tolist()))),
	params:
		sample = "{sample}",
		k = "{filter_k}",
		max_mem_in_GB = config["kmc"]["max_mem_in_GB"],
		mincov = "{mincov}",
		minprop = "{minprop}",
		maxcount = config["kmc"]["maxcount"],
		maxcounter = config["kmc"]["maxcounter"],
		nbin = 64,
	threads: config["threads"]["kmc"]
	singularity:
		"docker://chrishah/kmc3-docker:v3.0"
	log:
		stdout = "results/{sample}/logs/kmc_filter.{sample}.k{filter_k}-{mincov}-{minprop}.stdout.txt",
		stderr = "results/{sample}/logs/kmc_filter.{sample}.k{filter_k}-{mincov}-{minprop}.stderr.txt"
	output: 
		filtered = "results/{sample}/kmc/filtered/{filter_k}-{mincov}-{minprop}/{sample}.k{filter_k}.{mincov}.{minprop}.filtered.reve.fastq",
	shadow: "minimal"
	shell:
		"""
		echo -e "$(date)\tStarting kmc filtering"
		echo "{input}" | sed 's/ /\\n/g' | tail -n +2 | sort > fastqs.txt
		kmc_tools filter {input.db}/{wildcards.sample} -ci{params.mincov} @fastqs.txt -ci{params.minprop} {output.filtered} 1>> {log.stdout} 2>> {log.stderr}
		"""

rule fil_repair_sort:
	input:
		forw = rules.fil_kmc_filter_forw.output,
		reve = rules.fil_kmc_filter_reve.output
	threads: 1
	singularity:
		"docker://reslp/samtools:1.9"
	log:
		"results/{sample}/logs/kmc_sort.{sample}.k{filter_k}-{mincov}-{minprop}.log.txt",
	output:
		bam = "results/{sample}/kmc/filtered/{filter_k}-{mincov}-{minprop}/{sample}.k{filter_k}.{mincov}.{minprop}.sorted.bam"
	shadow: "minimal"
	shell:
		"""
		echo -e "$(date)\tStarting sort" &> {log}
		cat {input.forw} {input.reve} | perl -ne 'chomp; $h=$_; $s=<>; $p=<>; $q=<>; chomp $s; chomp $q; if (($h =~ / 1/) || ($h =~ /\/1$/)){{$flag=68}}else{{$flag=132}}; ($h,$e)=split(" ", $h); print substr($h, 1)."\\t$flag\\t*\\t0\\t0\\t*\\t*\\t0\\t0\\t$s\\t$q\\tRG:Z:$e\\n"' | samtools view -bS | samtools sort -n -@ {threads} 1> {output.bam} 2>> {log}
		"""

rule fil_repair_extract_pe:
	input:
		bam = rules.fil_repair_sort.output,
		forw = rules.fil_kmc_filter_forw.output,
		reve = rules.fil_kmc_filter_reve.output
	params:
		wd = os.getcwd()
	threads: 1
	singularity:
		"docker://reslp/samtools:1.9"
	log:
		stdout = "results/{sample}/logs/kmc_extract-pe.{sample}.k{filter_k}-{mincov}-{minprop}.stdout.txt",
		stderr = "results/{sample}/logs/kmc_extract-pe.{sample}.k{filter_k}-{mincov}-{minprop}.stderr.txt"
	output:
		forw = "results/{sample}/kmc/filtered/{filter_k}-{mincov}-{minprop}/0000.k{filter_k}.{mincov}.{minprop}.1.fastq",
		reve = "results/{sample}/kmc/filtered/{filter_k}-{mincov}-{minprop}/0000.k{filter_k}.{mincov}.{minprop}.2.fastq",
		singletons = "results/{sample}/kmc/filtered/{filter_k}-{mincov}-{minprop}/{sample}.k{filter_k}.{mincov}.{minprop}.filtered.singletons.txt",
	shadow: "minimal"
	shell:
		"""
		echo -e "$(date)\tStarting extract pe" 1> {log.stdout} 2> {log.stderr}

		samtools view {input.bam} | perl -ne 'chomp; push(@a, $_); if (scalar @a == 2){{@1 = split("\\t", $a[0]); @2 = split("\\t", $a[1]); if ($1[0] eq $2[0]){{$1[11] =~ s/RG:Z://; $2[11] =~ s/RG:Z://; $1[0] = "$1[0] $1[11]"; $2[0] = "$2[0] $2[11]"; print STDOUT "\@$1[0]\\n$1[9]\\n+\\n$1[10]\\n"; print STDOUT "\@$2[0]\\n$2[9]\\n+\\n$2[10]\\n"; @a=()}}else{{$singleton = shift @a; @line = split("\\t", $singleton); $line[11] =~ s/RG:Z://; $line[0] = "$line[0] $line[11]"; print STDERR "$line[0]\\n"}}}}' 1> results/{wildcards.sample}/kmc/filtered/{wildcards.filter_k}-{wildcards.mincov}-{wildcards.minprop}/0000.interleaved.fastq 2> results/{wildcards.sample}/kmc/filtered/{wildcards.filter_k}-{wildcards.mincov}-{wildcards.minprop}/{wildcards.sample}.k{wildcards.filter_k}.{wildcards.mincov}.{wildcards.minprop}.filtered.singletons.txt
		cat results/{wildcards.sample}/kmc/filtered/{wildcards.filter_k}-{wildcards.mincov}-{wildcards.minprop}/0000.interleaved.fastq | perl -ne '$h1=$_; $s1=<>; $p1=<>; $q1=<>; $h2=<>; $s2=<>; $p2=<>; $q2=<>; print STDOUT "$h1$s1$p1$q1"; print STDERR "$h2$s2$p2$q2"' 1> {output.forw} 2> {output.reve}

		rm -rv {params.wd}/{input.forw} {params.wd}/{input.reve} {params.wd}/{input.bam} 1>> {log.stdout} 2>> {log.stderr}
		"""

rule fil_repair_extract_se:
	input:
		singletons = rules.fil_repair_extract_pe.output.singletons,
		f_paired = lambda wildcards: expand("results/{{sample}}/trimming/trimgalore/{lib}/{{sample}}.{lib}.1.fastq.gz", sample=wildcards.sample, lib=list(set(unitdict[wildcards.sample]) & set(Illumina_process_df["lib"].tolist()))),
		r_paired = lambda wildcards: expand("results/{{sample}}/trimming/trimgalore/{lib}/{{sample}}.{lib}.2.fastq.gz", sample=wildcards.sample, lib=list(set(unitdict[wildcards.sample]) & set(Illumina_process_df["lib"].tolist()))),
	params:
		wd = os.getcwd()
	threads: 1
	singularity:
		"docker://chrishah/mira:v4.9.6"
	log:
		stdout = "results/{sample}/logs/kmc_extract-se.{sample}.k{filter_k}-{mincov}-{minprop}.stdout.txt",
		stderr = "results/{sample}/logs/kmc_extract-se.{sample}.k{filter_k}-{mincov}-{minprop}.stderr.txt"
	output:
		forw = "results/{sample}/kmc/filtered/{filter_k}-{mincov}-{minprop}/{sample}.k{filter_k}.{mincov}.{minprop}.filtered.1.fastq.gz",
		reve = "results/{sample}/kmc/filtered/{filter_k}-{mincov}-{minprop}/{sample}.k{filter_k}.{mincov}.{minprop}.filtered.2.fastq.gz",
	resources:
		mem_gb=config["max_mem_in_GB"]["miraconvert"]
#	shadow: "minimal"
	shell:
		"""
		echo -e "$(date)\tStarting extract se" 1> {log.stdout}
		cat {input.singletons} | perl -ne 'chomp; $forw = $_; $reve = $_; if ($_ =~ / 1/){{$reve =~ s/ 1/ 2/;}}else{{$forw =~ s/ 2/ 1/;}} print STDOUT "$forw\\n"; print STDERR "$reve\\n"' 1> results/{wildcards.sample}/kmc/filtered/{wildcards.filter_k}-{wildcards.mincov}-{wildcards.minprop}/{wildcards.sample}.k{wildcards.filter_k}.{wildcards.mincov}.{wildcards.minprop}.1.list.txt 2> results/{wildcards.sample}/kmc/filtered/{wildcards.filter_k}-{wildcards.mincov}-{wildcards.minprop}/{wildcards.sample}.k{wildcards.filter_k}.{wildcards.mincov}.{wildcards.minprop}.2.list.txt
		i=0001
		for f in $(echo "{input.f_paired}" | sort); do cmd="miraconvert -n results/{wildcards.sample}/kmc/filtered/{wildcards.filter_k}-{wildcards.mincov}-{wildcards.minprop}/{wildcards.sample}.k{wildcards.filter_k}.{wildcards.mincov}.{wildcards.minprop}.1.list.txt $f results/{wildcards.sample}/kmc/filtered/{wildcards.filter_k}-{wildcards.mincov}-{wildcards.minprop}/$i.1"; echo $cmd; $cmd; i=$(printf "%04d" $(( i + 1 ))); done 1>> {log.stdout} 2> {log.stderr}
		i=0001
		for f in $(echo "{input.r_paired}" | sort); do cmd="miraconvert -n results/{wildcards.sample}/kmc/filtered/{wildcards.filter_k}-{wildcards.mincov}-{wildcards.minprop}/{wildcards.sample}.k{wildcards.filter_k}.{wildcards.mincov}.{wildcards.minprop}.2.list.txt $f results/{wildcards.sample}/kmc/filtered/{wildcards.filter_k}-{wildcards.mincov}-{wildcards.minprop}/$i.2"; echo $cmd; $cmd; i=$(printf "%04d" $(( i + 1 ))); done 1>> {log.stdout} 2>> {log.stderr}
		cat results/{wildcards.sample}/kmc/filtered/{wildcards.filter_k}-{wildcards.mincov}-{wildcards.minprop}/0*.1.fastq > results/{wildcards.sample}/kmc/filtered/{wildcards.filter_k}-{wildcards.mincov}-{wildcards.minprop}/{wildcards.sample}.k{wildcards.filter_k}.{wildcards.mincov}.{wildcards.minprop}.filtered.1.fastq
		cat results/{wildcards.sample}/kmc/filtered/{wildcards.filter_k}-{wildcards.mincov}-{wildcards.minprop}/0*.2.fastq > results/{wildcards.sample}/kmc/filtered/{wildcards.filter_k}-{wildcards.mincov}-{wildcards.minprop}/{wildcards.sample}.k{wildcards.filter_k}.{wildcards.mincov}.{wildcards.minprop}.filtered.2.fastq
		rm results/{wildcards.sample}/kmc/filtered/{wildcards.filter_k}-{wildcards.mincov}-{wildcards.minprop}/0*.fastq results/{wildcards.sample}/kmc/filtered/{wildcards.filter_k}-{wildcards.mincov}-{wildcards.minprop}/{wildcards.sample}.k{wildcards.filter_k}.{wildcards.mincov}.{wildcards.minprop}.?.list.txt

		gzip -v results/{wildcards.sample}/kmc/filtered/{wildcards.filter_k}-{wildcards.mincov}-{wildcards.minprop}/{wildcards.sample}.k{wildcards.filter_k}.{wildcards.mincov}.{wildcards.minprop}.filtered.?.fastq 1>> {log.stdout} 2>> {log.stderr}
		
		rm -rv {params.wd}/{input.singletons} 1>> {log.stdout} 2>> {log.stderr}
		"""

rule fil_kmc_filtered:
	input:
		reads = rules.fil_repair_extract_se.output
	params:
		sample = "{sample}",
		k = "{filter_k}",
		max_mem_in_GB = config["kmc"]["max_mem_in_GB"],
		mincount = config["kmc"]["mincount"],
		maxcount = config["kmc"]["maxcount"],
		maxcounter = config["kmc"]["maxcounter"],
		nbin = 64,
	threads: config["threads"]["kmc"]
	singularity:
		"docker://chrishah/kmc3-docker:v3.0"
	log:
		stdout = "results/{sample}/logs/kmc_filtered.{sample}.k{filter_k}-{mincov}-{minprop}.stdout.txt",
		stderr = "results/{sample}/logs/kmc_filtered.{sample}.k{filter_k}-{mincov}-{minprop}.stderr.txt"
	output: 
		hist = "results/{sample}/kmc/filtered/{filter_k}-{mincov}-{minprop}/{sample}.k{filter_k}.{mincov}.{minprop}.filtered.histogram.txt",
		genomescope = "results/{sample}/kmc/filtered/{filter_k}-{mincov}-{minprop}/{sample}.k{filter_k}.{mincov}.{minprop}.filtered.histogram.genomescope.txt"
	shadow: "minimal"
	shell:
		"""
		echo -e "$(date)\tStarting kmc"
		echo "{input}" | sed 's/ /\\n/g' > fastqs.txt
		mkdir {params.sample}.db
		kmc -k{params.k} -m$(( {params.max_mem_in_GB} - 2 )) -v -sm -ci{params.mincount} -cx{params.maxcount} -cs{params.maxcounter} -n{params.nbin} -t$(( {threads} - 1 )) @fastqs.txt {params.sample} {params.sample}.db 1>> {log.stdout} 2>> {log.stderr}
		kmc_tools histogram {params.sample} -ci{params.mincount} {output.hist} 1> /dev/null 2>> {log.stderr}
		cat {output.hist} | awk '{{print $1" "$2}}' > {output.genomescope}
		"""

if config["kmer_filtering"]["assemble"] == "yes":
	rule fil_link_kmerfiltered:
		input:
			forw = rules.fil_repair_extract_se.output.forw,
			reve = rules.fil_repair_extract_se.output.reve,
			orphans = rules.tri_gather_short_trimmed_by_lib.output.orphans
		output:
			forw = "results/{sample}/trimming/trimgalore.k{filter_k}.{mincov}.{minprop}/{sample}-full/{sample}.trimgalore.k{filter_k}.{mincov}.{minprop}.1.fastq.gz",
			reve = "results/{sample}/trimming/trimgalore.k{filter_k}.{mincov}.{minprop}/{sample}-full/{sample}.trimgalore.k{filter_k}.{mincov}.{minprop}.2.fastq.gz",
			orphans = "results/{sample}/trimming/trimgalore.k{filter_k}.{mincov}.{minprop}/{sample}-full/{sample}.trimgalore.k{filter_k}.{mincov}.{minprop}.se.fastq.gz",
		shell:
			"""
			ln -s ../../../../../{input.forw} {output.forw}
			ln -s ../../../../../{input.reve} {output.reve}
			head -n 4 <(zcat {input.orphans}) | gzip > {output.orphans}
			"""

rule reformat_read_headers:
	input:
		"results/{sample}/trimming/trimgalore/{lib}/{sample}.{lib}.{pe}.fastq.gz"
	log:
		stdout = "results/{sample}/logs/reformat.{lib}.{pe}.stdout.txt",
		stderr = "results/{sample}/logs/reformat.{lib}.{pe}.stderr.txt"
	params:
		pe = "{pe}"
	threads: 2
	output:
		read = "results/{sample}/errorcorrection/bless/{lib}/{lib}.{pe}.fastq.gz",
		ok = "results/{sample}/errorcorrection/bless/{lib}/headerformat.{lib}.{pe}.ok"
	shell:
		"""
		#check if the read headers contain space - if yes reformat
		if [ "$(zcat {input} | head -n 1 | grep " " | wc -l)" -gt 0 ]
		then
			zcat {input} | sed 's/ {params.pe}.*/\/{params.pe}/' | gzip > {output.read} 
		else
			ln -s $(pwd)/{input} $(pwd)/{output.read}
		fi
		touch {output.ok}
		"""



rule clean_corrected_libs:
	input:
		forward = lambda wildcards: expand("results/{{sample}}/errorcorrection/bless/{lib}/{{sample}}.{lib}.1.corrected.fastq.gz", sample=wildcards.sample, lib=list(set(unitdict[wildcards.sample]) & set(Illumina_process_df["lib"].tolist()))),
		reverse = lambda wildcards: expand("results/{{sample}}/errorcorrection/bless/{lib}/{{sample}}.{lib}.2.corrected.fastq.gz", sample=wildcards.sample, lib=list(set(unitdict[wildcards.sample]) & set(Illumina_process_df["lib"].tolist()))),
		orphans = lambda wildcards: expand("results/{{sample}}/errorcorrection/bless/{lib}/{{sample}}.{lib}.se.corrected.fastq.gz", sample=wildcards.sample, lib=list(set(unitdict[wildcards.sample]) & set(Illumina_process_df["lib"].tolist()))),
	params:
		wd = os.getcwd(),
		sample = "{sample}",
	singularity:
		"docker://chrishah/trim_galore:0.6.0"
	output:
		ok = "results/{sample}/errorcorrection/bless/{sample}-full/{sample}.cat.status.ok",
		f = "results/{sample}/errorcorrection/bless/{sample}-full/{sample}.blesscorrected.1.fastq.gz",
		r = "results/{sample}/errorcorrection/bless/{sample}-full/{sample}.blesscorrected.2.fastq.gz",
		o = "results/{sample}/errorcorrection/bless/{sample}-full/{sample}.blesscorrected.se.fastq.gz"
	shadow: "minimal"
	threads: 2
	shell:
		"""
		if [ $(echo {input.forward} | wc -w) -gt 1 ]
		then
			cat {input.forward} > {output.f}
			cat {input.reverse} > {output.r}
		else
			ln -s ../../../../../{input.forward} {output.f}
			ln -s ../../../../../{input.reverse} {output.r}
		fi
		cat {input.orphans} > {output.o}
		touch {output.ok}
		"""

rule clean_merged_libs:
	input:
		merged = lambda wildcards: expand("results/{{sample}}/readmerging/usearch/{lib}/{{sample}}.{lib}.merged.fastq.gz", sample=wildcards.sample, lib=list(set(unitdict[wildcards.sample]) & set(Illumina_process_df["lib"].tolist()))),
		forward = lambda wildcards: expand("results/{{sample}}/readmerging/usearch/{lib}/{{sample}}.{lib}_1.nm.fastq.gz", sample=wildcards.sample, lib=list(set(unitdict[wildcards.sample]) & set(Illumina_process_df["lib"].tolist()))),
		reverse = lambda wildcards: expand("results/{{sample}}/readmerging/usearch/{lib}/{{sample}}.{lib}_2.nm.fastq.gz", sample=wildcards.sample, lib=list(set(unitdict[wildcards.sample]) & set(Illumina_process_df["lib"].tolist()))),
	params:
		wd = os.getcwd(),
		sample = "{sample}",
	singularity:
		"docker://chrishah/trim_galore:0.6.0"
	output:
		ok = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.cat.status.ok",
		m = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.usearchmerged.fastq.gz",
		f = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.usearchnotmerged.1.fastq.gz",
		r = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.usearchnotmerged.2.fastq.gz",
	shadow: "minimal"
	threads: 2
	shell:
		"""
		if [ $(echo {input.merged} | wc -w) -gt 1 ]
		then
			cat {input.merged} > {output.m}
			cat {input.forward} > {output.f}
			cat {input.reverse} > {output.r}
		else
			ln -s ../../../../../{input.merged} {output.m}
			ln -s ../../../../../{input.forward} {output.f}
			ln -s ../../../../../{input.reverse} {output.r}
		fi
		touch {output.ok}
		"""

#rule merge_libs_merged:
#rule filter_by_kmer_coverage:
#	input:
#	params:
#		sample = "{sample}",
#		k = "{k}",
#		max_mem_in_GB = config["kmc"]["max_mem_in_GB"],
#		mincount = config["kmc"]["mincount"],
#		maxcount = config["kmc"]["maxcount"],
#		maxcounter = config["kmc"]["maxcounter"],
#		nbin = 64,
#	threads: config["threads"]["kmc"]
#	singularity:
#		"docker://chrishah/kmc3-docker:v3.0"
#	log:
#		stdout = "results/{sample}/logs/kmc.{sample}.k{k}.stdout.txt",
#		stderr = "results/{sample}/logs/kmc.{sample}.k{k}.stderr.txt"
#	output: 
##		pre = "results/{sample}/kmc/{sample}.k{k}.kmc_pre",
##		suf = "results/{sample}/kmc/{sample}.k{k}.kmc_suf",
#		hist = "results/{sample}/kmc/{sample}.k{k}.histogram.txt",
#	shadow: "shallow"
#	shell:
#		"""
#		echo "{input}" | sed 's/ /\\n/g' > fastqs.txt
#		kmc_tools filter $db_prefix -ci$min_kmer_cov -cx$max_kmer_cov @fastqs.txt -ci$min_perc_kmers -cx$max_perc_kmers {output.fastq}
#		"""
#rule merge_k_hists:
#pdftk $(ls -1 $s-* | grep "\-full.pdf" | tr '\n' ' ') cat output $s.distribution.full.pdf
#pdftk $(ls -1 $s-* | grep -v "\-full.pdf" | tr '\n' ' ') cat output $s.distribution.pdf
