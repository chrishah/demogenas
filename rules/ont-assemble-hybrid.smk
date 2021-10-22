#
#rule spades_hybrid_ectools:
#	input:
#		merged = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.usearchmerged.fastq.gz",
#		f = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.usearchnotmerged.1.fastq.gz",
#		r = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.usearchnotmerged.2.fastq.gz",
#		se = "results/{sample}/errorcorrection/bless/{sample}-full/{sample}.blesscorrected.se.fastq.gz", 
#		long = "results/{sample}/errorcorrection/ectools/{sample}.ectools.{similarity}.corrected.fa.gz"
#	output:
#		ok = "results/{sample}/assembly/spades/hybrid-ectools-{similarity}/spades.ok"
#	log:
#		stdout = "results/{sample}/logs/spades.hybrid-ectools-{similarity}.stdout.txt",
#		stderr = "results/{sample}/logs/spades.hybrid-ectools-{similarity}.stderr.txt"
#	params:
#		wd = os.getcwd(),
#		sample = "{sample}",
#		dir = "results/{sample}/assembly/spades/hybrid-ectools-{similarity}",
#		mode = "only-assembler" #could be careful, only-error-correction, only-assembler 
#	singularity:
#		"docker://chrishah/spades:v3.14.0"
##	shadow: "minimal"
#	threads: 90
#	resources:
#		mem_gb=750
#	shell:
#		"""
#		echo "Host: $HOSTNAME" 1> {log.stdout} 2> {log.stderr}
#		if [ ! -d {params.dir} ]
#		then
#			mkdir -p {params.dir}
#		else
#			rm -rf {params.dir}/* 
#		fi
#
#		cd {params.dir}
#
#		spades.py \
#		-o ./{params.sample} \
#		-s {params.wd}/{input.se} \
#		--merged {params.wd}/{input.merged} \
#		-1 {params.wd}/{input.f} \
#		-2 {params.wd}/{input.r} \
#		--nanopore {params.wd}/{input.long} \
#		--checkpoints last \
#		--{params.mode} \
#		-t $(( {threads} - 1 )) \
#		-m $(( {resources.mem_gb} - 5 )) 1>> {params.wd}/{log.stdout} 2>> {params.wd}/{log.stderr}
#
#		touch spades.ok
#		"""
#
#rule abyss_rescaffold_ectools:
#	input:
#		merged = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.usearchmerged.fastq.gz",
#		f = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.usearchnotmerged.1.fastq.gz",
#		r = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.usearchnotmerged.2.fastq.gz",
#		se = "results/{sample}/errorcorrection/bless/{sample}-full/{sample}.blesscorrected.se.fastq.gz", 
#		bestk = rules.kmergenie.output.bestk,
#		abyss = "results/{sample}/assembly/abyss/{ass_to_scf}-bestk/abyss.ok", #rules.abyss_tricormer.output.ok,
#		long = rules.gather_ectools_corrected.output
#	output:
#		ok = "results/{sample}/assembly/abyss/ectools-{similarity}-rescaffold-{ass_to_scf}/abyss.ok"
#	log:
#		stdout = "results/{sample}/logs/abyss.ectools-{similarity}-rescaffold-{ass_to_scf}.stdout.txt",
#		stderr = "results/{sample}/logs/abyss.ectools-{similarity}-rescaffold-{ass_to_scf}.stderr.txt"
#	params:
#		wd = os.getcwd(),
#		sample = "{sample}",
#		dir = "results/{sample}/assembly/abyss/ectools-{similarity}-rescaffold-{ass_to_scf}",
#		defaultk = k,
#		ass_to_scaffold = "{ass_to_scf}"
#	singularity: "docker://reslp/abyss:2.2.5"
##	shadow: "minimal"
#	threads: 40
#	resources:
#		mem_gb=40
#	shell:
#		"""
#		export TMPDIR={params.wd}/{params.dir}/tmp
#		echo "Host: $HOSTNAME" 1> {log.stdout} 2> {log.stderr}
#
#		bestk=$(cat {input.bestk})
#		cd {params.dir}
#		#cp all relevant files from Illumina only assembly
#		for f in $(ls -1 ../{params.ass_to_scaffold}/{params.sample}-* | grep -v "\-9\." | grep -v "\-1[0-9]\." | grep -v "\-stats"); do ln -fs $f .; done 
#		if [ $bestk -eq 0 ]
#		then
#			echo -e "Setting k to {params.defaultk}" 1> {params.wd}/{log.stdout}
#			touch k.set.manually.to.{params.defaultk}
#			bestk={params.defaultk}
#		fi
#
#		abyss-pe -C {params.wd}/{params.dir} k=$bestk name={params.sample} np=$(( {threads} - 1 )) in='{params.wd}/{input.f} {params.wd}/{input.r}' se='{params.wd}/{input.merged} {params.wd}/{input.se}' long='longa' longa='{params.wd}/{input.long}' 1> {params.wd}/{log.stdout} 2> {params.wd}/{log.stderr}
#		touch {params.wd}/{output.ok}
#		"""
#
#rule abyss_rescaffold_long:
#	input:
#		merged = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.usearchmerged.fastq.gz",
#		f = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.usearchnotmerged.1.fastq.gz",
#		r = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.usearchnotmerged.2.fastq.gz",
#		se = "results/{sample}/errorcorrection/bless/{sample}-full/{sample}.blesscorrected.se.fastq.gz", 
#		bestk = rules.kmergenie.output.bestk,
#		abyss = rules.abyss_tricormer.output.ok,
#		long = get_long
#	output:
#		ok = "results/{sample}/assembly/abyss/long-rescaffold/abyss.ok"
#	log:
#		stdout = "results/{sample}/logs/abyss.long-rescaffold.stdout.txt",
#		stderr = "results/{sample}/logs/abyss.long-rescaffold.stderr.txt"
#	params:
#		wd = os.getcwd(),
#		sample = "{sample}",
#		dir = "results/{sample}/assembly/abyss/long-rescaffold",
#		defaultk = k
#	singularity: "docker://reslp/abyss:2.2.5"
##	shadow: "minimal"
#	threads: 40
#	resources:
#		mem_gb=40
#	shell:
#		"""
#		export TMPDIR={params.wd}/{params.dir}/tmp
#		echo "Host: $HOSTNAME" 1> {log.stdout} 2> {log.stderr}
#
#		bestk=$(cat {input.bestk})
#		cd {params.dir}
#		#cp all relevant files from Illumina only assembly
#		for f in $(ls -1 ../bestk/{params.sample}-* | grep -v "\-9\." | grep -v "\-1[0-9]\." | grep -v "\-stats"); do ln -fs $f .; done 
#		if [ $bestk -eq 0 ]
#		then
#			echo -e "Setting k to {params.defaultk}" 1> {params.wd}/{log.stdout}
#			touch k.set.manually.to.{params.defaultk}
#			bestk={params.defaultk}
#		fi
#
#		abyss-pe -C {params.wd}/{params.dir} k=$bestk name={params.sample} np=$(( {threads} - 1 )) in='{params.wd}/{input.f} {params.wd}/{input.r}' se='{params.wd}/{input.merged} {params.wd}/{input.se}' long='longa' longa='{params.wd}/{input.long}' 1> {params.wd}/{log.stdout} 2> {params.wd}/{log.stderr}
#		touch {params.wd}/{output.ok}
#		"""
#
#rule abyss_rescaffold_consent:
#	input:
#		merged = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.usearchmerged.fastq.gz",
#		f = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.usearchnotmerged.1.fastq.gz",
#		r = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.usearchnotmerged.2.fastq.gz",
#		se = "results/{sample}/errorcorrection/bless/{sample}-full/{sample}.blesscorrected.se.fastq.gz", 
#		bestk = rules.kmergenie.output.bestk,
#		abyss = rules.abyss_tricormer.output.ok,
#		long = rules.consent.output.fastq
#	output:
#		ok = "results/{sample}/assembly/abyss/consent-rescaffold/abyss.ok"
#	log:
#		stdout = "results/{sample}/logs/abyss.consent-rescaffold.stdout.txt",
#		stderr = "results/{sample}/logs/abyss.consent-rescaffold.stderr.txt"
#	params:
#		wd = os.getcwd(),
#		sample = "{sample}",
#		dir = "results/{sample}/assembly/abyss/consent-rescaffold",
#		defaultk = k
#	singularity: "docker://reslp/abyss:2.2.5"
##	shadow: "minimal"
#	threads: 40
#	resources:
#		mem_gb=40
#	shell:
#		"""
#		export TMPDIR={params.wd}/{params.dir}/tmp
#		echo "Host: $HOSTNAME" 1> {log.stdout} 2> {log.stderr}
#
#		bestk=$(cat {input.bestk})
#		cd {params.dir}
#		#cp all relevant files from Illumina only assembly
#		for f in $(ls -1 ../bestk/{params.sample}-* | grep -v "\-9\." | grep -v "\-1[0-9]\." | grep -v "\-stats"); do ln -fs $f .; done 
#		if [ $bestk -eq 0 ]
#		then
#			echo -e "Setting k to {params.defaultk}" 1> {params.wd}/{log.stdout}
#			touch k.set.manually.to.{params.defaultk}
#			bestk={params.defaultk}
#		fi
#
#		abyss-pe -C {params.wd}/{params.dir} k=$bestk name={params.sample} np=$(( {threads} - 1 )) in='{params.wd}/{input.f} {params.wd}/{input.r}' se='{params.wd}/{input.merged} {params.wd}/{input.se}' long='longa' longa='{params.wd}/{input.long}' 1> {params.wd}/{log.stdout} 2> {params.wd}/{log.stderr}
#		touch {params.wd}/{output.ok}
#		"""
#
rule ahy_abyss_scaffold:
	input:
		reads = get_illumina_assembly_input,
		bestk = rules.ail_kmergenie.output.bestk,
		abyss = rules.ail_abyss.output.ok,
		long = get_long_assembly_input
	output:
		ok = "results/{sample}/assembly/abyss_scaffold/{trimmer}-{corrector}-{merger}-{basecaller}-{longcorrection}/abyss.ok",
	log:
		stdout = "results/{sample}/logs/abyss_scaffold.{trimmer}-{corrector}-{merger}-{basecaller}-{longcorrection}.stdout.txt",
		stderr = "results/{sample}/logs/abyss_scaffold.{trimmer}-{corrector}-{merger}-{basecaller}-{longcorrection}.stderr.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		dir = "results/{sample}/assembly/abyss_scaffold/{trimmer}-{corrector}-{merger}-{basecaller}-{longcorrection}",
		origdir = rules.ail_abyss.output.ok.replace("/abyss.ok",""),
		defaultk = k
	singularity: "docker://reslp/abyss:2.2.5"
#	shadow: "minimal"
	threads: 90
	resources:
		#can problably use much less memory
		mem_gb=700
	shell:
		"""
		export TMPDIR={params.wd}/{params.dir}/tmp
		echo "Host: $HOSTNAME" 1> {log.stdout} 2> {log.stderr}

		bestk=$(cat {input.bestk})

		cd {params.dir}

		#get all raw long reads (even if concurrency setting doesn't reflect the real number)
		echo -e "[$(date)]\\tConcatenating the following files into {wildcards.sample}.{wildcards.basecaller}.fastq.gz prior to assembly" >> {params.wd}/{log.stdout}
		echo -e "$(ls -1 $(echo -e "{input.long}" | tr ' ' '\\n' | cut -d "/" -f 1-6 | sort | uniq | sed 's?^?{params.wd}/?' | sed 's?$?/{wildcards.sample}.{wildcards.basecaller}.*.fastq.gz?g'))" >> {params.wd}/{log.stdout}
		cat $(echo -e "{input.long}" | tr ' ' '\\n' | cut -d "/" -f 1-6 | sort | uniq | sed 's?^?{params.wd}/?' | sed 's?$?/{wildcards.sample}.{wildcards.basecaller}.*.fastq.gz?g') > {wildcards.sample}.{wildcards.basecaller}.fastq.gz

		#cp all relevant files from Illumina only assembly
		for f in $(ls -1 {params.wd}/{params.origdir}/{params.sample}-* | grep -v "\-9\." | grep -v "\-1[0-9]\." | grep -v "\-stats"); do ln -fs $f .; done 
		if [ $bestk -eq 0 ]
		then
			echo -e "Setting k to {params.defaultk}" 1> {params.wd}/{log.stdout}
			touch k.set.manually.to.{params.defaultk}
			bestk={params.defaultk}
		fi

#		inlong=$(for f in $(echo "{input.long}"); do echo {params.wd}/$f; done | tr '\\n' ' ')
		abyss-pe -C {params.wd}/{params.dir} k=$bestk name={params.sample} np=$(( {threads} - 1 )) in='{params.wd}/{input.reads[0]} {params.wd}/{input.reads[1]}' se='{params.wd}/{input.reads[2]}' long='longall' longall="{wildcards.sample}.{wildcards.basecaller}.fastq.gz" 1> {params.wd}/{log.stdout} 2> {params.wd}/{log.stderr}
		
		rm {wildcards.sample}.{wildcards.basecaller}.fastq.gz
		touch {params.wd}/{output.ok}
		"""

#rule abyss_rescaffold_corrected:
#	input:
#		reads = get_illumina_assembly_input,
##		merged = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.usearchmerged.fastq.gz",
##		f = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.usearchnotmerged.1.fastq.gz",
##		r = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.usearchnotmerged.2.fastq.gz",
##		se = "results/{sample}/errorcorrection/bless/{sample}-full/{sample}.blesscorrected.se.fastq.gz", 
#		bestk = rules.kmergenie.output.bestk,
#		abyss = rules.abyss.output.ok,
#		long = "results/{sample}/errorcorrection/{longcor}/{basecaller}/{sample}.{basecaller}.{longcor}.fastq.gz",
##		long = get_long_corrected_assembly_input
##		long = expand("results/{{sample}}/reads/ont/flappie/{lib}/{{sample}}.flappie.{unit}.fastq.gz", sample="{sample}", lib="{lib}", unit=flappie_unit_list)
#	output:
#		ok = "results/{sample}/assembly/{assinput}/abyss/rescaffold-corrected/{basecaller}/{longcor}/abyss.ok",
#	log:
#		stdout = "results/{sample}/logs/abyss.{assinput}.rescaffold-corrected-{basecaller}-{longcor}.stdout.txt",
#		stderr = "results/{sample}/logs/abyss.{assinput}.rescaffold-corrected-{basecaller}-{longcor}.stderr.txt"
#	params:
#		wd = os.getcwd(),
#		sample = "{sample}",
#		dir = "results/{sample}/assembly/{assinput}/abyss/rescaffold-corrected/{basecaller}/{longcor}",
#		origdir = rules.abyss.output.ok.replace("/abyss.ok",""),
#		defaultk = k
#	singularity: "docker://reslp/abyss:2.2.5"
##	shadow: "minimal"
#	threads: 90
#	resources:
#		#can problably use much less memory
#		mem_gb=700
#	shell:
#		"""
#		export TMPDIR={params.wd}/{params.dir}/tmp
#		echo "Host: $HOSTNAME" 1> {log.stdout} 2> {log.stderr}
#
#		bestk=$(cat {input.bestk})
#
#		cd {params.dir}
#		#cp all relevant files from Illumina only assembly
#		for f in $(ls -1 {params.wd}/{params.origdir}/{params.sample}-* | grep -v "\-9\." | grep -v "\-1[0-9]\." | grep -v "\-stats"); do ln -fs $f .; done 
#		if [ $bestk -eq 0 ]
#		then
#			echo -e "Setting k to {params.defaultk}" 1> {params.wd}/{log.stdout}
#			touch k.set.manually.to.{params.defaultk}
#			bestk={params.defaultk}
#		fi
#
#		inlong=$(for f in $(echo "{input.long}"); do echo {params.wd}/$f; done | tr '\\n' ' ')
#		abyss-pe -C {params.wd}/{params.dir} k=$bestk name={params.sample} np=$(( {threads} - 1 )) in='{params.wd}/{input.reads[0]} {params.wd}/{input.reads[1]}' se='{params.wd}/{input.reads[2]}' long='longall' longall="<(cat $inlong)" 1> {params.wd}/{log.stdout} 2> {params.wd}/{log.stderr}
#		
#		touch {params.wd}/{output.ok}
#		"""

rule ahy_spades_hybrid:
	input:
		reads = get_illumina_assembly_input,
		bestk = rules.ail_kmergenie.output.bestk,
		long = get_long_assembly_input
	output:
		ok = "results/{sample}/assembly/spades_hybrid/{trimmer}-{corrector}-{merger}-{basecaller}-{longcorrection}/spades.ok",
	log:
		stdout = "results/{sample}/logs/spades_hybrid.{trimmer}-{corrector}-{merger}-{basecaller}-{longcorrection}.stdout.txt",
		stderr = "results/{sample}/logs/spades_hybrid.{trimmer}-{corrector}-{merger}-{basecaller}-{longcorrection}.stderr.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		dir = "results/{sample}/assembly/spades_hybrid/{trimmer}-{corrector}-{merger}-{basecaller}-{longcorrection}",
		mode = "only-assembler" #could be careful, only-error-correction, only-assembler 
	singularity:
		"docker://chrishah/spades:v3.14.0"
#	shadow: "minimal"
	threads: 90
	resources:
		mem_gb=750
	shell:
		"""
		echo "Host: $HOSTNAME" 1> {log.stdout} 2> {log.stderr}
		if [ ! -d {params.dir} ]
		then
			mkdir -p {params.dir}
		else
			rm -rf {params.dir}/* 
		fi

		cd {params.dir}

		#get all raw long reads (even if concurrency setting doesn't reflect the real number)
		echo -e "[$(date)]\\tConcatenating the following files into {wildcards.sample}.{wildcards.basecaller}.fastq.gz prior to assembly" >> {log.stdout}
		echo -e "$(ls -1 $(echo -e "{input.long}" | tr ' ' '\\n' | cut -d "/" -f 1-6 | sort | uniq | sed 's?^?{params.wd}/?' | sed 's?$?/{wildcards.sample}.{wildcards.basecaller}.*.fastq.gz?g'))" >> {log.stdout}
		cat $(echo -e "{input.long}" | tr ' ' '\\n' | cut -d "/" -f 1-6 | sort | uniq | sed 's?^?{params.wd}/?' | sed 's?$?/{wildcards.sample}.{wildcards.basecaller}.*.fastq.gz?g') > {wildcards.sample}.{wildcards.basecaller}.fastq.gz

		#run spades
		spades.py \
		-o ./{params.sample} \
		-s {params.wd}/{input.reads[2]} \
		-1 {params.wd}/{input.reads[0]} \
		-2 {params.wd}/{input.reads[1]} \
		--nanopore {wildcards.sample}.{wildcards.basecaller}.fastq.gz \
		--checkpoints last \
		--{params.mode} \
		-t $(( {threads} - 1 )) \
		-m $(( {resources.mem_gb} - 5 )) 1>> {params.wd}/{log.stdout} 2>> {params.wd}/{log.stderr}

		rm {wildcards.sample}.{wildcards.basecaller}.fastq.gz
		touch spades.ok
		"""
#		--merged {params.wd}/{input.merged} \
#		$(for f in $(echo "{input.long}"); do echo "--nanopore {params.wd}/$f"; done | tr '\\n' ' ') \

#rule spades_hybrid_corrected:
#	input:
#		reads = get_illumina_assembly_input,
##		merged = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.usearchmerged.fastq.gz",
##		f = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.usearchnotmerged.1.fastq.gz",
##		r = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.usearchnotmerged.2.fastq.gz",
##		se = "results/{sample}/errorcorrection/bless/{sample}-full/{sample}.blesscorrected.se.fastq.gz", 
#		bestk = rules.kmergenie.output.bestk,
#		long = "results/{sample}/errorcorrection/{longcor}/{basecaller}/{sample}.{basecaller}.{longcor}.fastq.gz",
##		long = get_long_corrected_assembly_input
##		long = expand("results/{{sample}}/reads/ont/flappie/{lib}/{{sample}}.flappie.{unit}.fastq.gz", sample="{sample}", lib="{lib}", unit=flappie_unit_list)
#	output:
#		ok = "results/{sample}/assembly/{assinput}/spades-hybrid/corrected/{basecaller}/{longcor}/spades.ok",
#	log:
#		stdout = "results/{sample}/logs/spades.{assinput}.hybrid-corrected-{basecaller}-{longcor}.stdout.txt",
#		stderr = "results/{sample}/logs/spades.{assinput}.hybrid-corrected-{basecaller}-{longcor}.stderr.txt"
#	params:
#		wd = os.getcwd(),
#		sample = "{sample}",
#		dir = "results/{sample}/assembly/{assinput}/spades-hybrid/corrected/{basecaller}/{longcor}",
#		mode = "only-assembler" #could be careful, only-error-correction, only-assembler 
#	singularity:
#		"docker://chrishah/spades:v3.14.0"
##	shadow: "minimal"
#	threads: 90
#	resources:
#		mem_gb=750
#	shell:
#		"""
#		echo "Host: $HOSTNAME" 1> {log.stdout} 2> {log.stderr}
#		if [ ! -d {params.dir} ]
#		then
#			mkdir -p {params.dir}
#		else
#			rm -rf {params.dir}/* 
#		fi
#
#		cd {params.dir}
#
#		spades.py \
#		-o ./{params.sample} \
#		-s {params.wd}/{input.reads[2]} \
#		-1 {params.wd}/{input.reads[0]} \
#		-2 {params.wd}/{input.reads[1]} \
#		$(for f in $(echo "{input.long}"); do echo "--nanopore {params.wd}/$f"; done | tr '\\n' ' ') \
#		--checkpoints last \
#		--{params.mode} \
#		-t $(( {threads} - 1 )) \
#		-m $(( {resources.mem_gb} - 5 )) 1>> {params.wd}/{log.stdout} 2>> {params.wd}/{log.stderr}
#
#		touch spades.ok
#		"""

#rule spades_consent:
#	input:
#		merged = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.usearchmerged.fastq.gz",
#		f = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.usearchnotmerged.1.fastq.gz",
#		r = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.usearchnotmerged.2.fastq.gz",
#		se = "results/{sample}/errorcorrection/bless/{sample}-full/{sample}.blesscorrected.se.fastq.gz", 
#		long = rules.consent.output.fastq
#	output:
#		ok = "results/{sample}/assembly/spades/consent-long/spades.ok"
#	log:
#		stdout = "results/{sample}/logs/spades.consent-long.stdout.txt",
#		stderr = "results/{sample}/logs/spades.consent-long.stderr.txt"
#	params:
#		wd = os.getcwd(),
#		sample = "{sample}",
#		dir = "results/{sample}/assembly/spades/consent-long",
#		mode = "only-assembler" #could be careful, only-error-correction, only-assembler 
#	singularity:
#		"docker://chrishah/spades:v3.14.0"
##	shadow: "minimal"
#	threads: 90
#	resources:
#		mem_gb=750
#	shell:
#		"""
#		echo "Host: $HOSTNAME" 1> {log.stdout} 2> {log.stderr}
#		if [ ! -d {params.dir} ]
#		then
#			mkdir -p {params.dir}
#		else
#			rm -rf {params.dir}/* 
#		fi
#
#		cd {params.dir}
#
#		spades.py \
#		-o ./{params.sample} \
#		-s {params.wd}/{input.se} \
#		--merged {params.wd}/{input.merged} \
#		-1 {params.wd}/{input.f} \
#		-2 {params.wd}/{input.r} \
#		$(for f in $(echo "{input.long}"); do echo "--nanopore {params.wd}/$f"; done | tr '\\n' ' ') \
#		--checkpoints last \
#		--{params.mode} \
#		-t $(( {threads} - 1 )) \
#		-m $(( {resources.mem_gb} - 5 )) 1>> {params.wd}/{log.stdout} 2>> {params.wd}/{log.stderr}
#
#		touch spades.ok
#		"""

rule ahy_haslr:
	input:
		reads = get_illumina_assembly_input,
		bestk = rules.ail_kmergenie.output.bestk,
		long = get_long_assembly_input
	output:
		ok = "results/{sample}/assembly/haslr/{trimmer}-{corrector}-{merger}-{basecaller}/haslr.ok",
		dir = directory("results/{sample}/assembly/haslr/{trimmer}-{corrector}-{merger}-{basecaller}/{sample}")
	log:
		stdout = "results/{sample}/logs/haslr.{trimmer}-{corrector}-{merger}-{basecaller}.stdout.txt",
		stderr = "results/{sample}/logs/haslr.{trimmer}-{corrector}-{merger}-{basecaller}.stderr.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		dir = "results/{sample}/assembly/haslr/{trimmer}-{corrector}-{merger}-{basecaller}",
		genome_size = "350m",
		defaultk = k,
		minia_mode = "unitigs", #could also be contigs
		cov_lr = 50,
	singularity:
		"docker://chrishah/haslr:0.8a1"
	shadow: "minimal"
	threads: 90
	resources:
		mem_gb=750
	shell:
		"""
		bestk=$(cat {input.bestk})
		if [ $bestk -eq 0 ]
		then
			echo -e "[$(date)]\\tSetting k to {params.defaultk} (default)" 1> {log.stdout}
			touch k.set.manually.to.{params.defaultk}
			bestk={params.defaultk}
		else
			echo -e "[$(date)]\\tSetting k to $bestk (determined by kmergenie)" 1> {log.stdout}
		fi

		#get all raw long reads (even if concurrency setting doesn't reflect the real number)
		echo -e "[$(date)]\\tConcatenating the following files into {wildcards.sample}.{wildcards.basecaller}.fastq.gz prior to assembly" >> {log.stdout}
		echo -e "$(ls -1 $(echo -e "{input.long}" | tr ' ' '\\n' | cut -d "/" -f 1-6 | sort | uniq | sed 's?^?{params.wd}/?' | sed 's?$?/{wildcards.sample}.{wildcards.basecaller}.*.fastq.gz?g'))" >> {log.stdout}
		cat $(echo -e "{input.long}" | tr ' ' '\\n' | cut -d "/" -f 1-6 | sort | uniq | sed 's?^?{params.wd}/?' | sed 's?$?/{wildcards.sample}.{wildcards.basecaller}.*.fastq.gz?g') > {wildcards.sample}.{wildcards.basecaller}.fastq.gz

		haslr.py -t {threads} -o {params.dir}/{wildcards.sample} --minia-kmer $bestk --minia-asm {params.minia_mode} --cov-lr {params.cov_lr} -g {params.genome_size} -l {wildcards.sample}.{wildcards.basecaller}.fastq.gz -x nanopore -s {input.reads[0]} {input.reads[1]} {input.reads[2]} 1>> {log.stdout} 2> {log.stderr}

		echo -e "\\nExitcode: $?\\n"

		touch {output.ok}
		"""
