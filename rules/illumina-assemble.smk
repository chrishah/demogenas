rule kmergenie:
	input:
		reads = get_illumina_assembly_input,
##		merged = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.usearchmerged.fastq.gz",
##		f = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.usearchnotmerged.1.fastq.gz",
##		r = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.usearchnotmerged.2.fastq.gz",
##		se = "results/{sample}/errorcorrection/bless/{sample}-full/{sample}.blesscorrected.se.fastq.gz" 
	output:
		report = directory("results/{sample}/assembly/{assinput}/kmergenie/report"),
		bestk = "results/{sample}/assembly/{assinput}/kmergenie/{sample}.bestk",
		bestkcutoff = "results/{sample}/assembly/{assinput}/kmergenie/{sample}.bestk-cutoff",
	log:
		stdout = "results/{sample}/logs/kmergenie.{assinput}.stdout.txt",
		stderr = "results/{sample}/logs/kmergenie.{assinput}.stderr.txt"
	params:
		sample = "{sample}",
		mink = 21,
		stepk = 10,
		maxk = 121
	singularity: "docker://reslp/kmergenie:1.7051"
	threads: 10
	resources:
		mem_gb = 20
	shadow: "minimal"
	shell:
		"""
		echo {input} | tr ' ' '\\n' > fofn.txt
		#capturing the returncode both in case of success or (||) when dying witha n error stops snakemake to exit the shell upon an error. 
		#kmergenie exits with and error if no best kmer could be found - in this case we want the pipeline to go on so I check if this was the cause of the error in the next if
		# and if yes we create the expected output files
		kmergenie fofn.txt -l {params.mink} -k {params.maxk} -s {params.stepk} -o {params.sample} -t $(( {threads} - 1 )) --diploid 1> {log.stdout} 2> {log.stderr} && returncode=$? || returncode=$? 

		if [ $returncode -gt 0 ]
		then
			#report was not written - could be because no best k was found or because of some other error
			#check if best k could be not be found
			if [ "$(cat {log.stderr})" = "No best k found" ]
			then
				echo "\\nBest k could not be found - moving on" 1>> {log.stdout}
				echo "0" > {output.bestk}
				echo "0" > {output.bestkcutoff}
				mkdir -p {output.report}
			else
				exit 1
			fi
		else	
			grep "^best k" {log.stdout} | cut -d " " -f 3 > {output.bestk}
			grep "cut-off for best k" {log.stdout} | cut -d " " -f 7 > {output.bestkcutoff}
			mkdir -p {output.report}
			mv *.html *.pdf {output.report}
		fi
		"""
rule abyss:
	input:
		reads = get_illumina_assembly_input,
#		f = rules.clean_trimmed_libs.output.f_trimmed,
#		r = rules.clean_trimmed_libs.output.r_trimmed,
#		se = rules.clean_trimmed_libs.output.orphans,
		bestk = rules.kmergenie.output.bestk,
		bestkcutoff = rules.kmergenie.output.bestkcutoff
	output:
		ok = "results/{sample}/assembly/{assinput}/abyss/bestk/abyss.ok",
		unitigs = "results/{sample}/assembly/{assinput}/abyss/bestk/{sample}-unitigs.fa",
		scaffolds = "results/{sample}/assembly/{assinput}/abyss/bestk/{sample}-scaffolds.fa",
	log:
		stdout = "results/{sample}/logs/abyss.{assinput}.bestk.stdout.txt",
		stderr = "results/{sample}/logs/abyss.{assinput}.bestk.stderr.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		dir = "results/{sample}/assembly/{assinput}/abyss/bestk",
		defaultk = k
	singularity: "docker://reslp/abyss:2.2.5"
#	shadow: "minimal"
	threads: 90
	resources:
		mem_gb=370
	shell:
		"""
		export TMPDIR={params.wd}/{params.dir}/tmp
		echo "Host: $HOSTNAME" 1> {log.stdout} 2> {log.stderr}

		bestk=$(cat {input.bestk})
		kcutoff=$(cat {input.bestkcutoff})

		cd {params.dir}
		if [ $bestk -eq 0 ]
		then
			echo -e "Setting k to {params.defaultk}" 1> {params.wd}/{log.stdout}
			touch k.set.manually.to.{params.defaultk}
			bestk={params.defaultk}
			kcutoff=2
		fi

		abyss-pe -C {params.wd}/{params.dir} k=$bestk name={params.sample} np=$(( {threads} - 1 )) in='{params.wd}/{input.reads[0]} {params.wd}/{input.reads[1]}' se='{params.wd}/{input.reads[2]}' default 1> {params.wd}/{log.stdout} 2> {params.wd}/{log.stderr}
		touch {params.wd}/{output.ok}
		"""
rule minia:
	input:
		reads = get_illumina_assembly_input,
#		merged = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.usearchmerged.fastq.gz",
#		f = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.usearchnotmerged.1.fastq.gz",
#		r = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.usearchnotmerged.2.fastq.gz",
#		se = "results/{sample}/errorcorrection/bless/{sample}-full/{sample}.blesscorrected.se.fastq.gz",
		bestk = rules.kmergenie.output.bestk,
		bestkcutoff = rules.kmergenie.output.bestkcutoff
	output:
		ok = "results/{sample}/assembly/{assinput}/minia/bestk/minia.ok",
		unitigs = "results/{sample}/assembly/{assinput}/minia/bestk/{sample}_bestk.unitigs.fa",
		contigs = "results/{sample}/assembly/{assinput}/minia/bestk/{sample}_bestk.contigs.fa",
	log:
		stdout = "results/{sample}/logs/minia.{assinput}.bestk.stdout.txt",
		stderr = "results/{sample}/logs/minia.{assinput}.bestk.stderr.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		dir = "results/{sample}/assembly/{assinput}/minia/bestk",
		defaultk = k
	singularity:
		"docker://chrishah/minia:3.2.4"
#	shadow: "minimal"
	threads: 40
	resources:
		mem_gb=20
	shell:
		"""
		bestk=$(cat {input.bestk})
		kcutoff=$(cat {input.bestkcutoff})

		cd {params.dir}
		if [ $bestk -eq 0 ]
		then
			echo -e "Setting k to {params.defaultk}" 1> {params.wd}/{log.stdout}
			touch k.set.manually.to.{params.defaultk}
			bestk={params.defaultk}
			kcutoff=2
		fi

		minia \
		-in {params.wd}/{input.reads[0]} \
		-in {params.wd}/{input.reads[1]} \
		-in {params.wd}/{input.reads[2]} \
		-kmer-size $bestk -abundance-min $kcutoff -max-memory $(( ( {resources.mem_gb} - 5 ) * 1000 )) -out {params.sample}_k$bestk -nb-cores $(( {threads} - 1 )) 1>> {params.wd}/{log.stdout} 2> {params.wd}/{log.stderr}
		
		ln -s {params.sample}_k$bestk.unitigs.fa {params.sample}_bestk.unitigs.fa
		ln -s {params.sample}_k$bestk.contigs.fa {params.sample}_bestk.contigs.fa
		rm *.h5
		touch {params.wd}/{output.ok}
		"""
#		-in {params.wd}/{input.merged} \

rule platanus:
	input:
		reads = get_illumina_assembly_input,
#		merged = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.usearchmerged.fastq.gz",
#		f = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.usearchnotmerged.1.fastq.gz",
#		r = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.usearchnotmerged.2.fastq.gz",
#		se = "results/{sample}/errorcorrection/bless/{sample}-full/{sample}.blesscorrected.se.fastq.gz" 
	output:
		ok = "results/{sample}/assembly/{assinput}/platanus/auto/platanus.ok"
	log:
		stdout = "results/{sample}/logs/platanus.{assinput}.stdout.txt",
		stderr = "results/{sample}/logs/platanus.{assinput}.stderr.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		dir = "results/{sample}/assembly/{assinput}/platanus/auto",
		min = 100
	singularity:
		"docker://chrishah/platanus:v1.2.4"
#	shadow: "minimal"
	threads: 90
	resources:
		mem_gb=365
	shell:
		"""
		echo "Host: $HOSTNAME" 1> {log.stdout} 2> {log.stderr}
		if [ ! -d {params.dir} ]
		then
			mkdir -p {params.dir}
		else
			rm -f {params.dir}/* 
		fi
		
		cd {params.dir}
		echo -e "\n$(date)\tRunning platanus assemble" 1>> {params.wd}/{log.stdout} 2>> {params.wd}/{log.stderr}
		platanus assemble \
		-o {params.sample} \
		-f <(zcat {params.wd}/{input.reads[0]}) <(zcat {params.wd}/{input.reads[1]}) <(zcat {params.wd}/{input.reads[2]}) \
		-t {threads} \
		-m {resources.mem_gb} 1>> {params.wd}/{log.stdout} 2>> {params.wd}/{log.stderr}

		echo -e "\n$(date)\tRunning platanus scaffold" 1>> {params.wd}/{log.stdout} 2>> {params.wd}/{log.stderr}
		platanus scaffold \
		-o {params.sample} \
		-c {params.sample}_contig.fa \
		-b {params.sample}_contigBubble.fa \
		-IP1 <(zcat {params.wd}/{input.reads[0]}) <(zcat {params.wd}/{input.reads[1]}) \
		-t {threads} 1>> {params.wd}/{log.stdout} 2>> {params.wd}/{log.stderr}
#		-IP1 temp.f1.fasta temp.f2.fasta \

		echo -e "\n$(date)\tRunning platanus gap close" 1>> {params.wd}/{log.stdout} 2>> {params.wd}/{log.stderr}
		platanus gap_close \
		-o {params.sample} \
		-c {params.sample}_scaffold.fa \
		-f <(zcat {params.wd}/{input.reads[2]}) \
		-IP1 <(zcat {params.wd}/{input.reads[0]}) <(zcat {params.wd}/{input.reads[1]}) \
		-t {threads} 1>> {params.wd}/{log.stdout} 2>> {params.wd}/{log.stderr}

		touch platanus.ok
#		if [ ! -d {params.dir} ]
#		then
#			mkdir -p {params.dir} && mv {params.sample}.* {params.dir}/
#		else
#			rm -f {params.dir}/* && mv {params.sample}.* {params.dir}/
#		fi
#		touch {params.dir}/platanus.status.ok
		"""
		#filtering paired end reads by length (minimum 100)
#		paste <(zcat {params.wd}/{input.f}) <(zcat {params.wd}/{input.r}) | perl -ne 'chomp; $h=$_; $s=<>; chomp $s; $p=<>; $q=<>; chomp $p; chomp $q; @s=split("\\t",$s); if ((length($s[0]) >= {params.min}) && (length($s[1]) >= {params.min})){{@h=split("\\t",$h); $h[0] =~ s/^@/>/g; $h[1] =~ s/^@/>/g; @q=split("\\t",$q); print STDOUT "$h[0]\\n$s[0]\\n"; print STDERR "$h[1]\\n$s[1]\\n";}}' 1> temp.f1.fasta 2> temp.f2.fasta


rule spades:
	input:
		reads = get_illumina_assembly_input,
#		merged = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.usearchmerged.fastq.gz",
#		f = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.usearchnotmerged.1.fastq.gz",
#		r = "results/{sample}/readmerging/usearch/{sample}-full/{sample}.usearchnotmerged.2.fastq.gz",
#		se = "results/{sample}/errorcorrection/bless/{sample}-full/{sample}.blesscorrected.se.fastq.gz",
		bestk = rules.kmergenie.output.bestk 
	output:
		ok = "results/{sample}/assembly/{assinput}/spades/{kmode}/spades.ok"
	log:
		stdout = "results/{sample}/logs/spades.{assinput}.{kmode}.stdout.txt",
		stderr = "results/{sample}/logs/spades.{assinput}.{kmode}.stderr.txt"
	benchmark: "results/{sample}/benchmarks/spades.{assinput}.{kmode}.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		dir = "results/{sample}/assembly/{assinput}/spades/{kmode}",
		defaultks = "21,33,55,77",
		kmode = "{kmode}",
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
		bestk=$(cat {input.bestk})
		if [ $bestk -eq 0 ]
		then
			echo -e "No optimal k found" 1> {log.stdout} 2> {log.stderr}
			touch {output.ok}
			exit 0
		else
			ks=$(echo -e "{params.defaultks},$(cat {input.bestk})" | tr ',' '\\n' | sort -n | tr '\\n' ',' | sed 's/,$//')
			krange=$(if [[ "{params.kmode}" == "bestk" ]]; then echo "-k $ks"; else echo "-k auto"; fi)
		fi

		if [ ! -d {params.dir} ]
		then
			mkdir -p {params.dir}
		else
			rm -rf {params.dir}/* 
		fi

		cd {params.dir}

		spades.py \
		-o ./{params.sample} \
		-1 {params.wd}/{input.reads[0]} -2 {params.wd}/{input.reads[1]} -s {params.wd}/{input.reads[2]} \
		$krange \
		--checkpoints last \
		--{params.mode} \
		-t $(( {threads} - 1 )) \
		-m $(( {resources.mem_gb} - 5 )) 1>> {params.wd}/{log.stdout} 2>> {params.wd}/{log.stderr}

		touch spades.ok
		"""
#		--merged {params.wd}/{input.merged} \
