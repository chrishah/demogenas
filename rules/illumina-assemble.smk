#the next block determines which combinations of trimmer-corrector-merger are considered for assembly
#it reads the lists of trimmers, correctors and mergers in from the config files
#then adds these if the particular stage (trimmed and/or corrected and/or merged) is specified

if "merged" in config["assemble"]["assemble_at_stage"]:
	merge_list = list(config["illumina_merge"]["merger"])
else:
	merge_list = ["None"]

if "corrected" in config["assemble"]["assemble_at_stage"]:
	correct_list = list(config["illumina_correction"])
	if not "None" in merge_list:
		merge_list.append("None")
else:
	if "corrected" in config["illumina_merge"]["merge_at_stage"] and "merged" in config["assemble"]["assemble_at_stage"]:
		correct_list = list(config["illumina_correction"])
	else:
		correct_list = ["None"]

#if trimmed is selected, None needs to be added to corrected and merged so we get a trimmed-None-None return
if "trimmed" in config["assemble"]["assemble_at_stage"]:
	if not "None" in merge_list:
		merge_list.append("None")
	if not "None" in correct_list:
		correct_list.append("None")

#print("FINAL - Illumina_trim for assembly: "+str(trim_list))
#print("FINAL - Illumina_correct for assembly: "+str(correct_list))
#print("FINAL - Illumina_merging for assembly: "+str(merge_list))

rule ail_kmergenie:
	input:
		reads = get_illumina_assembly_input,
	output:
		report = directory("results/{sample}/assembly/kmergenie/{trimmer}-{corrector}-{merger}/report"),
		bestk = "results/{sample}/assembly/kmergenie/{trimmer}-{corrector}-{merger}/{sample}.bestk",
		bestkcutoff = "results/{sample}/assembly/kmergenie/{trimmer}-{corrector}-{merger}/{sample}.bestk-cutoff",
	log:
		stdout = "results/{sample}/logs/kmergenie.{trimmer}-{corrector}-{merger}.stdout.txt",
		stderr = "results/{sample}/logs/kmergenie.{trimmer}-{corrector}-{merger}.stderr.txt"
	params:
		sample = "{sample}",
		mink = 21,
		stepk = 10,
		maxk = 121
	singularity: "docker://reslp/kmergenie:1.7051"
	threads: config["threads"]["kmergenie"]
	resources:
		mem_gb=config["max_mem_in_GB"]["kmergenie"]
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
				echo -e "\\n[ $(date) ]\\tBest k could not be found - moving on" 1>> {log.stdout}
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

rule ail_abyss:
	input:
		reads = get_illumina_assembly_input,
		bestk = rules.ail_kmergenie.output.bestk,
		bestkcutoff = rules.ail_kmergenie.output.bestkcutoff
	output:
		ok = "results/{sample}/assembly/abyss/{trimmer}-{corrector}-{merger}/bestk/abyss.ok",
		unitigs = "results/{sample}/assembly/abyss/{trimmer}-{corrector}-{merger}/bestk/{sample}-unitigs.fa",
		scaffolds = "results/{sample}/assembly/abyss/{trimmer}-{corrector}-{merger}/bestk/{sample}-scaffolds.fa",
	log:
		stdout = "results/{sample}/logs/abyss.{trimmer}-{corrector}-{merger}.bestk.stdout.txt",
		stderr = "results/{sample}/logs/abyss.{trimmer}-{corrector}-{merger}.bestk.stderr.txt"
	benchmark: "results/{sample}/benchmarks/abyss.{trimmer}-{corrector}-{merger}.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		dir = "results/{sample}/assembly/abyss/{trimmer}-{corrector}-{merger}/bestk",
		defaultk = k
	singularity: "docker://reslp/abyss:2.2.5"
#	shadow: "minimal"
	threads: config["threads"]["abyss"]
	resources:
		mem_gb=config["max_mem_in_GB"]["abyss"]
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
rule ail_minia:
	input:
		reads = get_illumina_assembly_input,
		bestk = rules.ail_kmergenie.output.bestk,
		bestkcutoff = rules.ail_kmergenie.output.bestkcutoff
	output:
		ok = "results/{sample}/assembly/minia/{trimmer}-{corrector}-{merger}/bestk/minia.ok",
		unitigs = "results/{sample}/assembly/minia/{trimmer}-{corrector}-{merger}/bestk/{sample}_bestk.unitigs.fa",
		contigs = "results/{sample}/assembly/minia/{trimmer}-{corrector}-{merger}/bestk/{sample}_bestk.contigs.fa",
	log:
		stdout = "results/{sample}/logs/minia.{trimmer}-{corrector}-{merger}.bestk.stdout.txt",
		stderr = "results/{sample}/logs/minia.{trimmer}-{corrector}-{merger}.bestk.stderr.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		dir = "results/{sample}/assembly/minia/{trimmer}-{corrector}-{merger}/bestk",
		defaultk = k
	singularity:
		"docker://chrishah/minia:3.2.4"
#	shadow: "minimal"
	threads: config["threads"]["minia"]
	resources:
		mem_gb=config["max_mem_in_GB"]["minia"]
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

rule ail_platanus_assemble:
	input:
		reads = get_illumina_assembly_input,
	output:
		ok = "results/{sample}/assembly/platanus/{trimmer}-{corrector}-{merger}/auto/platanus.assemble.ok"
	log:
		stdout = "results/{sample}/logs/platanus.assemble.{trimmer}-{corrector}-{merger}.stdout.txt",
		stderr = "results/{sample}/logs/platanus.assemble.{trimmer}-{corrector}-{merger}.stderr.txt"
	benchmark: "results/{sample}/benchmarks/platanus.assemble.{trimmer}-{corrector}-{merger}.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		dir = "results/{sample}/assembly/platanus/{trimmer}-{corrector}-{merger}/auto",
		min = 100,
		tmp = "."
	singularity:
		"docker://chrishah/platanus:v1.2.4"
#	shadow: "minimal"
	threads: config["threads"]["platanus_assemble"]
	resources:
		mem_gb=config["max_mem_in_GB"]["platanus_assemble"]
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
		-m {resources.mem_gb} -tmp {params.tmp} 1>> {params.wd}/{log.stdout} 2>> {params.wd}/{log.stderr}

		touch {params.wd}/{output}
		"""

rule ail_platanus_scaffold:
	input:
		reads = get_illumina_assembly_input,
		assemble = rules.ail_platanus_assemble.output
	output:
		ok = "results/{sample}/assembly/platanus/{trimmer}-{corrector}-{merger}/auto/platanus.scaffold.ok"
	log:
		stdout = "results/{sample}/logs/platanus.scaffold.{trimmer}-{corrector}-{merger}.stdout.txt",
		stderr = "results/{sample}/logs/platanus.scaffold.{trimmer}-{corrector}-{merger}.stderr.txt"
	benchmark: "results/{sample}/benchmarks/platanus.scaffold.{trimmer}-{corrector}-{merger}.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		dir = "results/{sample}/assembly/platanus/{trimmer}-{corrector}-{merger}/auto",
		min = 100,
		tmp = "."
	singularity:
		"docker://chrishah/platanus:v1.2.4"
#	shadow: "minimal"
	threads: config["threads"]["platanus_scaffold"]
	resources:
		mem_gb=config["max_mem_in_GB"]["platanus_scaffold"]
	shell:
		"""
		cd {params.dir}
		#filtering paired end reads by length (minimum 100)
		echo -e "\n$(date)\tFiltering paired end reads by length (mininum {params.min} bp) before scaffolding" 1>> {params.wd}/{log.stdout} 2>> {params.wd}/{log.stderr}
		paste <(zcat {params.wd}/{input.reads[0]}) <(zcat {params.wd}/{input.reads[1]}) | perl -ne 'chomp; $h=$_; $s=<>; chomp $s; $p=<>; $q=<>; chomp $p; chomp $q; @s=split("\\t",$s); if ((length($s[0]) >= {params.min}) && (length($s[1]) >= {params.min})){{@h=split("\\t",$h); $h[0] =~ s/^@/>/g; $h[1] =~ s/^@/>/g; @q=split("\\t",$q); print STDOUT "$h[0]\\n$s[0]\\n"; print STDERR "$h[1]\\n$s[1]\\n";}}' 1> {wildcards.sample}.min{params.min}.1.fasta 2> {wildcards.sample}.min{params.min}.2.fasta

		echo -e "\n$(date)\tRunning platanus scaffold" 1>> {params.wd}/{log.stdout} 2>> {params.wd}/{log.stderr}
		platanus scaffold \
		-o {params.sample} \
		-c {params.sample}_contig.fa \
		-b {params.sample}_contigBubble.fa \
		-IP1 {wildcards.sample}.min{params.min}.1.fasta {wildcards.sample}.min{params.min}.2.fasta \
		-t {threads} -tmp {params.tmp} 1>> {params.wd}/{log.stdout} 2>> {params.wd}/{log.stderr}
#		-IP1 temp.f1.fasta temp.f2.fasta \

		touch {params.wd}/{output}
		"""

rule ail_platanus_gapclose:
	input:
		reads = get_illumina_assembly_input,
		scaffold = rules.ail_platanus_scaffold.output
	output:
		ok = "results/{sample}/assembly/platanus/{trimmer}-{corrector}-{merger}/auto/platanus.gapclose.ok"
	log:
		stdout = "results/{sample}/logs/platanus.gapclose.{trimmer}-{corrector}-{merger}.stdout.txt",
		stderr = "results/{sample}/logs/platanus.gapclose.{trimmer}-{corrector}-{merger}.stderr.txt"
	benchmark: "results/{sample}/benchmarks/platanus.gapclose.{trimmer}-{corrector}-{merger}.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		dir = "results/{sample}/assembly/platanus/{trimmer}-{corrector}-{merger}/auto",
		min = 100,
		tmp = "."
	singularity:
		"docker://chrishah/platanus:v1.2.4"
#	shadow: "minimal"
	threads: config["threads"]["platanus_gapclose"]
	resources:
		mem_gb=config["max_mem_in_GB"]["platanus_gapclose"]
	shell:
		"""
		cd {params.dir}
		echo -e "\n$(date)\tRunning platanus gap close" 1>> {params.wd}/{log.stdout} 2>> {params.wd}/{log.stderr}
		platanus gap_close \
		-o {params.sample} \
		-c {params.sample}_scaffold.fa \
		-f <(zcat {params.wd}/{input.reads[2]}) \
		-IP1 {wildcards.sample}.min{params.min}.1.fasta {wildcards.sample}.min{params.min}.2.fasta \
		-t {threads} -tmp {params.tmp} 1>> {params.wd}/{log.stdout} 2>> {params.wd}/{log.stderr}

		rm {wildcards.sample}.min{params.min}.1.fasta {wildcards.sample}.min{params.min}.2.fasta
		touch {params.wd}/{output}
		"""
		#filtering paired end reads by length (minimum 100)
#		paste <(zcat {params.wd}/{input.f}) <(zcat {params.wd}/{input.r}) | perl -ne 'chomp; $h=$_; $s=<>; chomp $s; $p=<>; $q=<>; chomp $p; chomp $q; @s=split("\\t",$s); if ((length($s[0]) >= {params.min}) && (length($s[1]) >= {params.min})){{@h=split("\\t",$h); $h[0] =~ s/^@/>/g; $h[1] =~ s/^@/>/g; @q=split("\\t",$q); print STDOUT "$h[0]\\n$s[0]\\n"; print STDERR "$h[1]\\n$s[1]\\n";}}' 1> temp.f1.fasta 2> temp.f2.fasta


rule ail_spades:
	input:
		reads = get_illumina_assembly_input,
		bestk = rules.ail_kmergenie.output.bestk 
	output:
		ok = "results/{sample}/assembly/spades/{trimmer}-{corrector}-{merger}/{kmode}/spades.ok"
	log:
		stdout = "results/{sample}/logs/spades.{trimmer}-{corrector}-{merger}.{kmode}.stdout.txt",
		stderr = "results/{sample}/logs/spades.{trimmer}-{corrector}-{merger}.{kmode}.stderr.txt"
	benchmark: "results/{sample}/benchmarks/spades.{trimmer}-{corrector}-{merger}.{kmode}.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		dir = "results/{sample}/assembly/spades/{trimmer}-{corrector}-{merger}/{kmode}",
		defaultks = "21,33,55,77",
		longks = "21,33,55,77,99,127",
		kmode = "{kmode}",
		resume = config["assemble"]["spades_resume"],
		optional = config["assemble"]["spades_additional_options"],
		mode = "only-assembler", #could be careful, only-error-correction, only-assembler,
	singularity:
		"docker://reslp/spades:3.15.3"
#	shadow: "minimal"
	threads: config["threads"]["spades"]
	resources:
		mem_gb=config["max_mem_in_GB"]["spades"]
	shell:
		"""
		bestk=$(cat {input.bestk})
		if [ "{params.resume}" == "yes" ] && [ -d {params.dir}/{params.sample} ]
		then
			echo -e "\\n[$(date)]\\tContinuing previous run" 1>> {log.stdout} 2> {log.stderr}
			echo -e "[$(date)]\\tResuming SPAdes" 1>> {params.wd}/{log.stdout} 2>> {params.wd}/{log.stderr}
			cd {params.dir}
			spades.py -o {params.sample} --restart-from last -m $(( {resources.mem_gb} - 5 )) 1>> {params.wd}/{log.stdout} 2>> {params.wd}/{log.stderr}
			touch spades.ok
			exit 0
		elif [ $bestk -eq 0 ] && [ "{params.kmode}" == "bestk" ]
		then
			echo -e "[$(date)]\\tNo optimal k found" 1> {log.stdout} 2> {log.stderr}
			touch {output.ok}
			exit 0
		else
			if [[ "{params.kmode}" == "default" ]]
			then
				echo -e "[$(date)]\\tusing default kmer range" 1>> {log.stdout} 2>> {log.stderr}
				krange="-k auto"
			elif [[ "{params.kmode}" == "bestk" ]]
			then
				echo -e "[$(date)]\\tchecking read lengths" 1>> {log.stdout} 2>> {log.stderr}
				longest_read=$(for f in {input.reads}; do head -n 40000 <(zcat $f) | sed -n '2~4p'; done | perl -ne 'chomp; if (length($_) > $length){{$length=length($_)}}; if (eof()){{print "$length\\n"}}')
				echo -e "[$(date)]\\tlongest read was $longest_read" 1>> {log.stdout} 2>> {log.stderr}
				if [ "$longest_read" -ge 250 ]
				then
					ks=$(echo -e "{params.longks},$bestk" | tr ',' '\\n' | sort -n | uniq | tr '\\n' ',' | sed 's/,$//')
				else
					ks=$(echo -e "{params.defaultks},$bestk" | tr ',' '\\n' | sort -n | uniq | tr '\\n' ',' | sed 's/,$//')
				fi
				echo -e "[$(date)]\\tusing kmer range of: $ks" 1>> {log.stdout} 2>> {log.stderr}
				krange="-k $ks"
			fi
		fi

		if [ ! -d {params.dir} ]
		then
			mkdir -p {params.dir}
		else 
			rm -rf {params.dir}/* 
		fi

		cd {params.dir}

		echo -e "[$(date)]\\tStarting SPAdes" 1>> {params.wd}/{log.stdout} 2>> {params.wd}/{log.stderr}
		spades.py \
		-o ./{params.sample} \
		-1 {params.wd}/{input.reads[0]} -2 {params.wd}/{input.reads[1]} -s {params.wd}/{input.reads[2]} \
		$krange \
		--checkpoints last \
		--{params.mode} \
		-t $(( {threads} - 1 )) \
		-m $(( {resources.mem_gb} - 5 )) $(if [[ "{params.optional}" != "None" ]]; then echo "{params.optional}"; fi) 1>> {params.wd}/{log.stdout} 2>> {params.wd}/{log.stderr}

		touch spades.ok
		"""
rule ail_megahit:
	input:
		reads = get_illumina_assembly_input,
	output:
		ok = "results/{sample}/assembly/megahit/{trimmer}-{corrector}-{merger}/auto/megahit.ok",
	log:
		stdout = "results/{sample}/logs/megahit.{trimmer}-{corrector}-{merger}.bestk.stdout.txt",
		stderr = "results/{sample}/logs/megahit.{trimmer}-{corrector}-{merger}.bestk.stderr.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		dir = "results/{sample}/assembly/megahit/{trimmer}-{corrector}-{merger}/auto",
		optional = config["assemble"]["megahit_additional_options"],
	singularity:
		"docker://chrishah/megahit:v1.2.9-f8afe5d"
	threads: config["threads"]["megahit"]
	resources:
		mem_gb=config["max_mem_in_GB"]["megahit"]
	shell:
		"""
		cd {params.dir}

		megahit \
		-1 {params.wd}/{input.reads[0]} \
		-2 {params.wd}/{input.reads[1]} \
		-r {params.wd}/{input.reads[2]} \
		-m $(( ( {resources.mem_gb} - 5 ) * 1073741824 )) -t $(( {threads} - 1 )) $(if [[ "{params.optional}" != "None" ]]; then echo "{params.optional}"; fi) \
		--out-prefix {params.sample} 1>> {params.wd}/{log.stdout} 2> {params.wd}/{log.stderr}
		
		touch {params.wd}/{output.ok}
		"""
### this woudl be for restarting
#		if [ ! -d {params.dir} ]
#		then
#			mkdir -p {params.dir}
#		fi
#
#		cd {params.dir}
#
#		echo -e "[$(date)]\\tStarting SPAdes" 1>> {params.wd}/{log.stdout} 2>> {params.wd}/{log.stderr}
#		spades.py \
#		-o ./{params.sample} \
#		$krange \
#		--checkpoints last \
#		-t $(( {threads} - 1 )) \
#		--restart-from last \
#		-m $(( {resources.mem_gb} - 5 )) 1>> {params.wd}/{log.stdout} 2>> {params.wd}/{log.stderr}
###

rule ail_gather_illumina_assemblies:
	input:
		expand("results/{{sample}}/assembly/spades/{trimmer}-{corrector}-{merger}/{kmode}/spades.ok", sample=df["sample"],
                        trimmer=trim_list,
                        corrector=correct_list,
                        merger=merge_list,
			kmode=config["assemble"]["spades_kmode"])
	output:
		"results/{sample}/assembly/spades/spades.ok"
	shell:
		"""
		touch {output}
		"""
