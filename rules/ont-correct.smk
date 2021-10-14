rule ratatosk:
	input:
		reads = get_illumina_assembly_input,
		long = get_called_by_sample_by_lib
	output:
		fastq = "results/{sample}/errorcorrection/ratatosk/{basecaller}/{lib}/{sample}.{lib}.{basecaller}.ratatosk.fastq.gz"
	log:
		stdout = "results/{sample}/logs/ratatosk.{sample}.{lib}.{basecaller}.stdout.txt",
		stderr = "results/{sample}/logs/ratatosk.{sample}.{lib}.{basecaller}.stderr.txt"
	benchmark: "results/{sample}/benchmarks/ratatosk.{sample}.{lib}.{basecaller}.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		dir = "results/{sample}/errorcorrection/ratatosk/{basecaller}/{lib}",
		isize = 500,
		k1 = 31,
		k2 = 63,
		options = "",
#		options = "--trim-split 20",
	singularity: "docker://chrishah/ratatosk:0.4-kdefault"
	shadow: "minimal"
	threads: 90
	resources:
		mem_gb=750
	shell:
		"""
		#get all raw long reads (even if concurrency setting doesn't reflect the real number)
		echo -e "[$(date)]\\tConcatenating the following files into {wildcards.sample}.{wildcards.lib}.{wildcards.basecaller}.raw.fastq.gz prior to correction" > {log.stdout}
		echo -e "$(ls -1 {params.wd}/results/{wildcards.sample}/reads/ont/{wildcards.basecaller}/{wildcards.lib}/{wildcards.sample}.{wildcards.basecaller}.*.fastq.gz)" >> {log.stdout}
		cat $(echo -e "$(ls -1 {params.wd}/results/{wildcards.sample}/reads/ont/{wildcards.basecaller}/{wildcards.lib}/{wildcards.sample}.{wildcards.basecaller}.*.fastq.gz)") > {wildcards.sample}.{wildcards.lib}.{wildcards.basecaller}.raw.fastq.gz

		Ratatosk \
		-s {input.reads[0]} -s {input.reads[1]} -s {input.reads[2]} \
		-l {wildcards.sample}.{wildcards.lib}.{wildcards.basecaller}.raw.fastq.gz \
		-c {threads} -i {params.isize} -k {params.k1} -K {params.k2} -v {params.options} \
		-o {params.dir}/{wildcards.sample}.{wildcards.lib}.{wildcards.basecaller}.ratatosk 1>> {log.stdout} 2> {log.stderr}
		gzip {params.dir}/{wildcards.sample}.{wildcards.lib}.{wildcards.basecaller}.ratatosk.fastq
		"""

#ruleorder: consent > gather_corrected_by_lib

rule gather_corrected_by_lib:
	input:
		long = gather_corrected_by_lib,
	output:
		fastq = "results/{sample}/errorcorrection/{corrector}/{basecaller}/{sample}.{basecaller}.{corrector}.fastq.gz"
#	wildcard_constraints:
#		corrector = "ratatosk"
	log:
		stdout = "results/{sample}/logs/{corrector}-gather.{basecaller}.stdout.txt",
		stderr = "results/{sample}/logs/{corrector}-gather.{basecaller}.stderr.txt"
	benchmark: "results/{sample}/benchmarks/{corrector}-gather.{basecaller}.txt"
	params:
		wd = os.getcwd(),
	singularity: "docker://reslp/consent:v2.1"
	shadow: "minimal"
	threads: 2
	resources:
		mem_gb=4
	shell:
		"""
		echo -e "[$(date)]\\tStarting on host: $HOSTNAME" 1> {log.stdout} 
		cat {input.long} > {output.fastq} 2> {log.stderr}
		echo -e "[$(date)]\\tDone!" 1> {log.stdout} 
		"""

rule consent:
	input:
		long = get_called_by_sample
	output:
		fastq = "results/{sample}/errorcorrection/consent/{basecaller}/{sample}.{basecaller}.consent.fastq.gz"
	log:
		stdout = "results/{sample}/logs/consent.{basecaller}.stdout.txt",
		stderr = "results/{sample}/logs/consent.{basecaller}.stderr.txt"
	benchmark: "results/{sample}/benchmarks/consent.{basecaller}.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		dir = "results/{sample}/errorcorrection/consent/{basecaller}",
#		tempdir = "/scratch/{sample}-{lib}",
		consent_options = "--minimapIndex 2000M",
	singularity: "docker://reslp/consent:v2.1"
	shadow: "minimal"
	threads: 95
	resources:
		mem_gb=90
	shell:
		"""
#		cat {input.long} > long.fastq.gz
#		CONSENT-correct {params.consent_options} -j {threads} --in long.fastq.gz --out {params.dir}/{wildcards.sample}.consent.fastq --type ONT 1> {log.stdout} 2> {log.stderr}
#		rm long.fastq.gz
		zcat {input.long} | perl -ne 'chomp; $h=$_; $s=<>; $p=<>; $q=<>; $h =~ s/^@/>/; print "$h\\n$s"' > long.fasta
		CONSENT-correct {params.consent_options} -j {threads} --in long.fasta --out {params.dir}/{wildcards.sample}.{wildcards.basecaller}.consent.fasta --type ONT 1> {log.stdout} 2> {log.stderr}
		cat {params.dir}/{wildcards.sample}.{wildcards.basecaller}.consent.fasta | perl -ne 'chomp; $h=$_; $s=<>; chomp $s; $h =~ s/^>/@/; $l = length($s); $qline = print "$h\\n$s\\n+\\n"; print "I"x$l; print "\\n"' | gzip -v > {output.fastq}
		"""

rule canu_correct:
	input:
		long = get_called_by_sample
	output:
		fastq = "results/{sample}/errorcorrection/canu/{basecaller}/{sample}.{basecaller}.canu.fastq.gz"
#		dir = directory("results/{sample}/errorcorrection/canu/{basecaller}/{sample}")
	log:
		stdout = "results/{sample}/logs/canu-correct.{basecaller}.stdout.txt",
		stderr = "results/{sample}/logs/canu-correct.{basecaller}.stderr.txt"
	benchmark:
		"results/{sample}/benchmarks/canu-correct.{basecaller}.benchmark.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		genome_size = "200m",
		dir = "results/{sample}/errorcorrection/canu/{basecaller}",
		options = "" #
	singularity: "docker://chrishah/canu:v2.1.1"
#	shadow: "minimal"
	threads: 94
	resources:
		mem_gb=750
	shell:
		"""
		canu -correct \
		-d {params.dir}/{wildcards.sample} \
		-p {sample} genomeSize={params.genome_size} \
		useGrid=false {params.options} \
		-nanopore {input.long} 1> {log.stdout} 2> {log.stderr}

		touch {output.fastq}
		"""

#rule cdhitest_unitigs:
#	input:
#		unitigs = rules.abyss_tricormer.output.unitigs,
#	output:
#		nr = "results/{sample}/assembly/abyss/bestk/{sample}-unitigs.nr-{similarity}.fa"
#	log:
#		stdout = "results/{sample}/logs/cdhitest-abyss-unitigs-nr-{similarity}.stdout.txt",
#		stderr = "results/{sample}/logs/cdhitest-abyss-unitigs-nr-{similarity}.stderr.txt"
#	params:
#		wd = os.getcwd(),
#		sample = "{sample}",
#		similarity = "{similarity}",
#	singularity: "docker://chrishah/cdhit:v4.8.1"
#	shadow: "minimal"
#	threads: 16
#	resources:
#		mem_gb=8
#	shell:
#		"""
#		cd-hit-est -i {input.unitigs} -o {output.nr} -c {params.similarity} -T {threads} -M $(( ( {resources.mem_gb} - 1 ) * 1000 )) 1> {log.stdout} 2> {log.stderr}
#		"""
#
#rule split:
#	input:
##		unitigs = rules.abyss.output.unitigs,
##		unitigs = rules.cdhitest_unitigs.output.nr,
#		unitigs = lambda wildcards: expand("results/{{sample}}/assembly/abyss/bestk/{{sample}}-unitigs.nr-{{similarity}}.fa", sample="{sample}", similarity=config["ectools"]["unitig_similarity"]),
#		long = expand(rules.flappie.output.fastq, lib="{lib}", sample="{sample}", lib="{lib}", unit=flappie_unit_list)
#	output:
#		ok = "results/{sample}/errorcorrection/ectools/pass1/partition-{similarity}.ok",
#	log:
#		stdout = "results/{sample}/logs/long.split-pass1-{similarity}.stdout.txt",
#		stderr = "results/{sample}/logs/long.split-pass1-{similarity}.stderr.txt"
#	params:
#		wd = os.getcwd(),
#		sample = "{sample}",
#		dir = "results/{sample}/errorcorrection/ectools/pass1",
#		minlength_in = config["ectools"]["minlength_in"],
#		minlength_out = config["ectools"]["minlength_out"],
#		reads_per_file = config["ectools"]["reads_per_file"],
#		files_per_dir = config["ectools"]["files_per_dir"]
#	singularity: "docker://chrishah/ectools-docker:latest"
#	threads: 2
#	resources:
#		mem_gb=4
#	shell:
#		"""
#		cd {params.dir}
#		#make sure there are no directories
#		rm -rf $(find ./ -maxdepth 1 -mindepth 1 -type d)
#		#get correct paths to all files
#		inlong=$(for f in $(echo "{input.long}"); do echo {params.wd}/$f; done | tr '\\n' ' ')
#		cat $inlong > long.fastq.gz
#		#partition (script in container path)
#		partition.py -minlen {params.minlength_in} {params.reads_per_file} {params.files_per_dir} long.fastq.gz 1> {params.wd}/{log.stdout} 2> {params.wd}/{log.stderr}
#		rm long.fastq.gz
#		#filter unitigs by length
#		cat {params.wd}/{input.unitigs} | perl -ne 'chomp; if ($_ =~ />/){{print "\\n$_\\n"}}else{{print "$_"}}' | grep -v "^$" | perl -ne 'chomp; $h=$_; $s=<>; if (length($s) > {params.minlength_in}){{print "$h\\n$s"}}' > unitigs.nr-{wildcards.similarity}-min-{params.minlength_in}.fasta
#		#prepare correction script (original in container)
#		cat /usr/src/ectools/correct.sh | \
#		sed 's?/path/to/ectools?/usr/src/ectools?' | \
#		sed "s?UNITIG_FILE=.*?UNITIG_FILE=$(pwd)/unitigs.nr-{wildcards.similarity}-min-{params.minlength_in}.fasta?" | \
#		sed 's/^source /#source /' | \
#		sed 's/^nucmer /nucmer -t $1 /' | \
#		sed 's/^MIN_READ_LEN=.*/MIN_READ_LEN={params.minlength_out}/' > correct.sh
#		#write checkpoint file
#		touch {params.wd}/{output.ok}
#		"""
#
#rule ectools:
#	input:
##		rules.split.output.ok
#		lambda wildcards: expand("results/{{sample}}/errorcorrection/ectools/pass1/partition-{{similarity}}.ok", sample="{sample}", similarity=config["ectools"]["unitig_similarity"])
#	output:
#		ok = "results/{sample}/errorcorrection/ectools/pass1/ectools.{similarity}.{unit}.ok",
#	log:
#		stdout = "results/{sample}/logs/ectools.pass1.{unit}.{similarity}.stdout.txt",
#		stderr = "results/{sample}/logs/ectools.pass1.{unit}.{similarity}.stderr.txt"
#	params:
#		wd = os.getcwd(),
#		sample = "{sample}",
#		unit = "{unit}",
#		dir = "results/{sample}/errorcorrection/ectools/pass1",
#		nbatches = ec_concurrency,
#	singularity: "docker://chrishah/ectools-docker:latest"
##	shadow: "minimal"
#	threads: 40
#	resources:
#		mem_gb=20
#	shell:
#		"""
#		cd {params.dir}
#		unit={params.unit}
#		step={params.nbatches}
#		unset TMPDIR
#		for c in $(find ./ -mindepth 1 -maxdepth 1 -type d | sed 's?./??' | sort -n )
#		do
#			#remove leading zeros (double curly bracket is just to escape wildcards)
#			#count=$(echo $c | sed 's/^0*//')
#			count=$(printf %0.f $c)
#			times=$(( count / step ))
#			if [ $(( ( count - times * step ) + 1 )) -eq "$unit" ]
#			then
#				echo -e "$(date)\\tprocessing - $count" 
#				cd $(printf "%04d" $count)
#				for suffix in $(ls -1 p* | cut -d "." -f 1 | sort | uniq | sed 's/^p//')
#				do
#					echo -ne "\\t$(date)\\t$count-$suffix - "
#					if [ ! -f "p$suffix.cor.fa" ] && [ ! -f "p$suffix.failed.fa" ]
#					then
#						#make sure to remove any old files
#						if [ $(ls -1 p$suffix* | wc -l) -gt 1 ]; then rm p$suffix.*; fi
#						#create tmpdir
#						if [ -d tmp-$suffix ]; then rm -rf tmp-$suffix; fi
#						mkdir tmp-$suffix
#						export TMPDIR=$(pwd)/tmp-$suffix
#						export SGE_TASK_ID=$(printf %0.f $suffix)
#						bash ../correct.sh $(( {threads} - 1 )) && returncode=$? || returncode=$?
#						if [ $returncode -gt 0 ]
#						then
#							mv p$suffix p$suffix.failed.fa
#							echo -e "$(date)\\tsomething went wrong with p$suffix"
#						else
#							rm $(ls -1 p$suffix* | grep -v "cor.fa")
#							echo -e "$(date)\\tp$suffix.cor.fa done"
#						fi
#						rm -rf tmp-$suffix
#					else
#						echo -e "$(date)\\tp$suffix.cor.fa done previously"
#					fi
#				done
#				cd ..
#			fi
#		done 1> {params.wd}/{log.stdout} 2> {params.wd}/{log.stderr}
#		touch {params.wd}/{output.ok}
#		"""
#
#rule split_pass2:
#	input:
#		pass1 = expand(rules.ectools.output.ok, sample="{sample}", lib="{lib}", unit=flappie_unit_list, similarity=config["ectools"]["unitig_similarity"])
#	output:
#		ok = "results/{sample}/errorcorrection/ectools/pass2/partition-{similarity}.ok",
#	log:
#		stdout = "results/{sample}/logs/long.split-pass2-{similarity}.stdout.txt",
#		stderr = "results/{sample}/logs/long.split-pass2-{similarity}.stderr.txt"
#	params:
#		wd = os.getcwd(),
#		sample = "{sample}",
#		dir = "results/{sample}/errorcorrection/ectools/pass2",
#		minlength_in = config["ectools"]["minlength_in"],
#		minlength_out = config["ectools"]["minlength_out"],
#		reads_per_file = 1,
#		nbatches = ec_concurrency,
#	singularity: "docker://chrishah/ectools-docker:latest"
#	threads: 2
#	resources:
#		mem_gb=4
#	shell:
#		"""
#		cd {params.dir}
#		#make sure there are no directories
#		rm -rf $(find ./ -maxdepth 1 -mindepth 1 -type d)
#
#		#gather all failed reads
#		cat ../pass1/*/*.failed.fa > failed.fa
#
#		#find the number of reads per directory to use the concurrency optimally
#		files_per_dir=$(cat failed.fa | grep ">" | sed -n '1~{params.nbatches}p' | wc -l)
#		#partition (script in container path)
#		partition.py -minlen {params.minlength_in} {params.reads_per_file} $files_per_dir failed.fa 1> {params.wd}/{log.stdout} 2> {params.wd}/{log.stderr}
#		rm failed.fa
#
#		#create symbolic link to unitigs file used for correction
#		ln -s ../pass1/unitigs.nr-{wildcards.similarity}-min-{params.minlength_in}.fasta .
#		
#		#prepare correction script (original in container)
#		cat /usr/src/ectools/correct.sh | \
#		sed 's?/path/to/ectools?/usr/src/ectools?' | \
#		sed "s?UNITIG_FILE=.*?UNITIG_FILE=$(pwd)/unitigs.nr-{wildcards.similarity}-min-{params.minlength_in}.fasta?" | \
#		sed 's/^source /#source /' | \
#		sed 's/^nucmer /nucmer -t $1 /' | \
#		sed 's/^MIN_READ_LEN=.*/MIN_READ_LEN={params.minlength_out}/' > correct.sh
#
#		#write checkpoint file
#		touch {params.wd}/{output.ok}
#		"""
#rule ectools_pass2:
#	input:
#		lambda wildcards: expand("results/{{sample}}/errorcorrection/ectools/pass2/partition-{{similarity}}.ok", sample="{sample}", similarity=config["ectools"]["unitig_similarity"])
#	output:
#		ok = "results/{sample}/errorcorrection/ectools/pass2/ectools.{similarity}.{unit}.ok",
#	log:
#		stdout = "results/{sample}/logs/ectools.pass2.{unit}.{similarity}.stdout.txt",
#		stderr = "results/{sample}/logs/ectools.pass2.{unit}.{similarity}.stderr.txt"
#	params:
#		wd = os.getcwd(),
#		sample = "{sample}",
#		unit = "{unit}",
#		dir = "results/{sample}/errorcorrection/ectools/pass2",
#		nbatches = ec_concurrency,
#	singularity: "docker://chrishah/ectools-docker:latest"
##	shadow: "minimal"
#	threads: 1
#	resources:
#		mem_gb=20
#	shell:
#		"""
#		cd {params.dir}
#		unit={params.unit}
#		step={params.nbatches}
#		unset TMPDIR
#		for c in $(find ./ -mindepth 1 -maxdepth 1 -type d | sed 's?./??' | sort -n )
#		do
#			#remove leading zeros (double curly bracket is just to escape wildcards)
#			#count=$(echo $c | sed 's/^0*//')
#			count=$(printf %0.f $c)
#			times=$(( count / step ))
#			if [ $(( ( count - times * step ) + 1 )) -eq "$unit" ]
#			then
#				echo -e "$(date)\\tprocessing - $count" 
#				cd $(printf "%04d" $count)
#				for suffix in $(ls -1 p* | cut -d "." -f 1 | sort | uniq | sed 's/^p//')
#				do
#					echo -ne "\\t$(date)\\t$count-$suffix - "
#					if [ ! -f "p$suffix.cor.fa" ] && [ ! -f "p$suffix.failed.fa" ]
#					then
#						#make sure to remove any old files
#						if [ $(ls -1 p$suffix* | wc -l) -gt 1 ]; then rm p$suffix.*; fi
#						#create tmpdir
#						if [ -d tmp-$suffix ]; then rm -rf tmp-$suffix; fi
#						mkdir tmp-$suffix
#						export TMPDIR=$(pwd)/tmp-$suffix
#						export SGE_TASK_ID=$(printf %0.f $suffix)
#						bash ../correct.sh $(( {threads} - 1 )) && returncode=$? || returncode=$?
#						if [ $returncode -gt 0 ]
#						then
#							mv p$suffix p$suffix.failed.fa
#							echo -e "$(date)\\tsomething went wrong with p$suffix"
#						else
#							rm $(ls -1 p$suffix* | grep -v "cor.fa")
#							echo -e "$(date)\\tp$suffix.cor.fa done"
#						fi
#						rm -rf tmp-$suffix
#					else
#						echo -e "$(date)\\tp$suffix.cor.fa done previously"
#					fi
#				done
#				cd ..
#			fi
#		done 1> {params.wd}/{log.stdout} 2> {params.wd}/{log.stderr}
#		touch {params.wd}/{output.ok}
#		"""
#rule gather_ectools_corrected:
#	input:
#		lambda wildcards: expand("results/{{sample}}/errorcorrection/ectools/pass2/ectools.{{similarity}}.{unit}.ok", sample=df_fast5["sample"], unit=ec_unit_list, similarity=config["ectools"]["unitig_similarity"])
#	output:
#		"results/{sample}/errorcorrection/ectools/{sample}.ectools.{similarity}.corrected.fa.gz"
#	log:
#		stdout = "results/{sample}/logs/ectools-{similarity}.gather.stdout.txt",
#		stderr = "results/{sample}/logs/ectools-{similarity}.gather.stderr.txt"
#	params:
#		dir = "results/{sample}/errorcorrection/ectools",
#	threads: 2
#	shell:
#		"""
#		cat $(find {params.dir} -name "*.cor.fa") | gzip > {output}
#		"""		

#rule ratatosk:
#	input:
#		reads = get_illumina_assembly_input,
##		f = rules.clean_trimmed_libs.output.f_trimmed,
##		r = rules.clean_trimmed_libs.output.r_trimmed,
##		se = rules.clean_trimmed_libs.output.orphans,
##		merged = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.usearchmerged.fastq.gz",
##		f = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.usearchnotmerged.1.fastq.gz",
##		r = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.usearchnotmerged.2.fastq.gz",
##		se = "results/{sample}/errorcorrection/bless/{sample}-full/{sample}.blesscorrected.se.fastq.gz", 
##		long = get_long_assembly_input
#		long = expand("results/{{sample}}/reads/ont/{basecaller}/{{lib}}/{{sample}}.{basecaller}.{unit}.fastq.gz", basecaller=config["basecaller"], unit=flappie_unit_list)
##		long = expand("results/{{units.sample}}/reads/ont/{basecaller}/{{units.lib}}/{{units.sample}}.{basecaller}.{unit}.fastq.gz", basecaller=config["basecaller"], units=fast5_units.itertuples(), unit=flappie_unit_list),
#	output:
#		fastq = "results/{sample}/errorcorrection/ratatosk/{basecaller}/{lib}/{sample}.{basecaller}.ratatosk.fastq.gz"
#	log:
#		stdout = "results/{sample}/logs/ratatosk.{sample}.{lib}.{basecaller}.stdout.txt",
#		stderr = "results/{sample}/logs/ratatosk.{sample}.{lib}.{basecaller}.stderr.txt"
#	benchmark: "results/{sample}/benchmarks/ratatosk.{sample}.{lib}.{basecaller}.txt"
#	params:
#		wd = os.getcwd(),
#		sample = "{sample}",
#		dir = "results/{sample}/errorcorrection/ratatosk/{basecaller}/{lib}",
#		isize = 500,
#		k1 = 31,
#		k2 = 63,
#		options = "",
##		options = "--trim-split 20",
#	singularity: "docker://chrishah/ratatosk:0.4-kdefault"
##	shadow: "minimal"
#	threads: 90
#	resources:
#		mem_gb=750
#	shell:
#		"""
#		Ratatosk \
#		-s {input.reads[0]} -s {input.reads[1]} -s {input.reads[2]} \
#		$(echo " {input.long}" | sed 's/ / -l /g') \
#		-c {threads} -i {params.isize} -k {params.k1} -K {params.k2} -v {params.options} \
#		-o {params.dir}/{wildcards.sample}.ratatosk 1> {log.stdout} 2> {log.stderr}
#		gzip {params.dir}/{wildcards.sample}.ratatosk.fastq
#		"""

#		$(echo " {input}[:3]" | sed 's/ / -s /g') \
#		-s {input.merged} -s {input.se} -s {input.f} -s {input.r} \

#rule gather_corrected:
#	input:
#		long = gather_corrected_by_lib,
#	output:
#		fastq = "results/{sample}/errorcorrection/ratatosk/{basecaller}/{sample}.{basecaller}.ratatosk.fastq.gz"
#	log:
#		stdout = "results/{sample}/logs/ratatosk-gather.{basecaller}.stdout.txt",
#		stderr = "results/{sample}/logs/ratatosk-gather.{basecaller}.stderr.txt"
#	benchmark: "results/{sample}/benchmarks/ratatosk-gather.{basecaller}.txt"
#	params:
#		wd = os.getcwd(),
#		sample = "{sample}",
#		dir = "results/{sample}/errorcorrection/ratatosk/{basecaller}",
##		tempdir = "/scratch/{sample}-{lib}",
#		consent_options = "--minimapIndex 2000M",
#	singularity: "docker://reslp/consent:v2.1"
#	shadow: "minimal"
#	threads: 95
#	resources:
#		mem_gb=90
#	shell:
#		"""
#		touch {output.fastq}
#		"""
#rule consent:
#	input:
#		long = get_long_assembly_input
#	output:
#		fastq = "results/{sample}/errorcorrection/consent/{basecaller}/{sample}.{basecaller}.consent.fastq.gz"
#	log:
#		stdout = "results/{sample}/logs/consent.{basecaller}.stdout.txt",
#		stderr = "results/{sample}/logs/consent.{basecaller}.stderr.txt"
#	benchmark: "results/{sample}/benchmarks/consent.{basecaller}.txt"
#	params:
#		wd = os.getcwd(),
#		sample = "{sample}",
#		dir = "results/{sample}/errorcorrection/consent/{basecaller}",
##		tempdir = "/scratch/{sample}-{lib}",
#		consent_options = "--minimapIndex 2000M",
#	singularity: "docker://reslp/consent:v2.1"
#	shadow: "minimal"
#	threads: 95
#	resources:
#		mem_gb=90
#	shell:
#		"""
##		cat {input.long} > long.fastq.gz
##		CONSENT-correct {params.consent_options} -j {threads} --in long.fastq.gz --out {params.dir}/{wildcards.sample}.consent.fastq --type ONT 1> {log.stdout} 2> {log.stderr}
##		rm long.fastq.gz
#		zcat {input.long} | perl -ne 'chomp; $h=$_; $s=<>; $p=<>; $q=<>; $h =~ s/^@/>/; print "$h\\n$s"' > long.fasta
#		CONSENT-correct {params.consent_options} -j {threads} --in long.fasta --out {params.dir}/{wildcards.sample}.consent.fasta --type ONT 1> {log.stdout} 2> {log.stderr}
#		cat {params.dir}/{wildcards.sample}.consent.fasta | perl -ne 'chomp; $h=$_; $s=<>; chomp $s; $h =~ s/^>/@/; $l = length($s); $qline = print "$h\\n$s\\n+\\n"; print "I"x$l; print "\\n"' | gzip -v > {output}
#		"""

#rule canu_correct:
#	input:
#		long = get_long_assembly_input
#	output:
#		fastq = "results/{sample}/errorcorrection/canu/{sample}.{basecaller}.canu.fastq.gz",
##		dir = directory("results/{sample}/errorcorrection/canu/{sample}")
#	log:
#		stdout = "results/{sample}/logs/canu-correct.{basecaller}.stdout.txt",
#		stderr = "results/{sample}/logs/canu-correct.{basecaller}.stderr.txt"
#	benchmark:
#		"results/{sample}/benchmarks/canu-correct.{basecaller}.benchmark.txt"
#	params:
#		wd = os.getcwd(),
#		sample = "{sample}",
#		genome_size = "200m",
#		dir = "results/{sample}/errorcorrection/canu/{basecaller}",
#		options = "" #
#	singularity: "docker://chrishah/canu:v2.1.1"
##	shadow: "minimal"
#	threads: 94
#	resources:
#		mem_gb=750
#	shell:
#		"""
#		canu -correct \
#		-d {params.dir}/{wildcards.sample} \
#		-p {sample} genomeSize={params.genome_size} \
#		useGrid=false {params.options} \
#		-nanopore {input.long} 1> {log.stdout} 2> {log.stderr}
#
#		touch {output.fastq}
#		"""
