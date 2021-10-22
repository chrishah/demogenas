rule flye_raw:
	input:
		long = get_long_assembly_input
	output:
		ok = "results/{sample}/assembly/flye/raw/{basecaller}/flye.ok",
		dir = directory("results/{sample}/assembly/flye/raw/{basecaller}/{sample}")
	log:
		stdout = "results/{sample}/logs/flye.raw.{basecaller}.stdout.txt",
		stderr = "results/{sample}/logs/flye.raw.{basecaller}.stderr.txt"
	benchmark:
		"results/{sample}/benchmarks/flye.raw.{basecaller}.benchmark.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		options = "--genome-size 350m --asm-coverage 20" #e.g.:--genome-size 200m --asm-coverage 25
	singularity: "docker://chrishah/flye:2.8.2-b1695"
	shadow: "minimal"
	threads: 90
	resources:
		mem_gb=750
	shell:
		"""
		#get all raw long reads (even if concurrency setting doesn't reflect the real number)
		echo -e "[$(date)]\\tConcatenating the following files into {wildcards.sample}.{wildcards.basecaller}.fastq.gz prior to assembly" > {log.stdout}
		echo -e "$(ls -1 $(echo -e "{input.long}" | tr ' ' '\\n' | cut -d "/" -f 1-6 | sort | uniq | sed 's?^?{params.wd}/?' | sed 's?$?/{wildcards.sample}.{wildcards.basecaller}.*.fastq.gz?g'))" >> {log.stdout}
		cat $(echo -e "{input.long}" | tr ' ' '\\n' | cut -d "/" -f 1-6 | sort | uniq | sed 's?^?{params.wd}/?' | sed 's?$?/{wildcards.sample}.{wildcards.basecaller}.*.fastq.gz?g') > {wildcards.sample}.{wildcards.basecaller}.fastq.gz

		flye \
		--out-dir {output.dir} \
		--threads {threads} {params.options} \
		--nano-raw {wildcards.sample}.{wildcards.basecaller}.fastq.gz 1>> {log.stdout} 2> {log.stderr}

		touch {output.ok}
		"""
		
rule flye_corrected:
	input:
#		long = "results/{sample}/errorcorrection/{longcorrection}/{trimmer}-{corrector}-{merger}-{basecaller}/{sample}.{basecaller}.{longcorrection}.fastq.gz"
		long = get_long_corrected_assembly_input
	output:
		ok = "results/{sample}/assembly/flye/{trimmer}-{corrector}-{merger}-{basecaller}-{longcorrection}/flye.ok",
		dir = directory("results/{sample}/assembly/flye/{trimmer}-{corrector}-{merger}-{basecaller}-{longcorrection}/{sample}")
	log:
		stdout = "results/{sample}/logs/flye.{trimmer}-{corrector}-{merger}-{basecaller}.{longcorrection}.stdout.txt",
		stderr = "results/{sample}/logs/flye.{trimmer}-{corrector}-{merger}-{basecaller}.{longcorrection}.stderr.txt"
	benchmark:
		"results/{sample}/benchmarks/flye.{trimmer}-{corrector}-{merger}-{basecaller}.{longcorrection}.benchmark.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		options = "" #e.g.:--genome-size 200m --asm-coverage 25
	singularity: "docker://chrishah/flye:2.8.2-b1695"
	shadow: "minimal"
	threads: 90
	resources:
		mem_gb=750
	shell:
		"""
		flye \
		--out-dir {output.dir} \
		--threads {threads} {params.options} \
		--nano-corr {input.long} 1> {log.stdout} 2> {log.stderr}

		touch {output.ok}
		"""

rule canu:
	input:
		long = get_long_assembly_input
	output:
		ok = "results/{sample}/assembly/canu/raw/{basecaller}/canu.ok",
		dir = directory("results/{sample}/assembly/canu/raw/{basecaller}/{sample}")
	log:
		stdout = "results/{sample}/logs/canu.raw.{basecaller}.stdout.txt",
		stderr = "results/{sample}/logs/canu.raw.{basecaller}.stderr.txt"
	benchmark:
		"results/{sample}/benchmarks/canu.raw.{basecaller}.benchmark.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		genome_size = "200m",
		dir = "results/{sample}/assembly/canu/raw/{basecaller}",
		options = "" #
	singularity: "docker://chrishah/canu:v2.1.1"
	shadow: "minimal"
	threads: 90
	resources:
		mem_gb=750
	shell:
		"""
		canu -assemble \
		-d {params.dir}/{wildcards.sample} \
		-p {sample} genomeSize={params.genome_size} \
		useGrid=false {params.options} \
		-nanopore {input.long} 1> {log.stdout} 2> {log.stderr}
		"""

rule canu_corrected:
	input:
#		long = "results/{sample}/errorcorrection/{longcorrection}/{trimmer}-{corrector}-{merger}-{basecaller}/{sample}.{basecaller}.{longcorrection}.fastq.gz"
		long = get_long_corrected_assembly_input
	output:
		ok = "results/{sample}/assembly/canu/{trimmer}-{corrector}-{merger}-{basecaller}-{longcorrection}/canu.ok",
		dir = directory("results/{sample}/assembly/canu/{trimmer}-{corrector}-{merger}-{basecaller}-{longcorrection}/{sample}")
	log:
		stdout = "results/{sample}/logs/canu.{trimmer}-{corrector}-{merger}-{basecaller}.{longcorrection}.stdout.txt",
		stderr = "results/{sample}/logs/canu.{trimmer}-{corrector}-{merger}-{basecaller}.{longcorrection}.stderr.txt"
	benchmark:
		"results/{sample}/benchmarks/canu.{trimmer}-{corrector}-{merger}-{basecaller}.{longcorrection}.benchmark.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		genome_size = "200m",
		dir = "results/{sample}/assembly/canu/{trimmer}-{corrector}-{merger}-{basecaller}-{longcorrection}",
		options = "" #
	singularity: "docker://chrishah/canu:v2.1.1"
	shadow: "minimal"
	threads: 90
	resources:
		mem_gb=750
	shell:
		"""
		canu -assemble \
		-d {params.dir}/{wildcards.sample} \
		-p {sample} genomeSize={params.genome_size} \
		useGrid=false {params.options} \
		-nanopore-corrected {input.long} 1> {log.stdout} 2> {log.stderr}
		"""


#	useGrid=

#rule marvel:
#	input:
#		do = "marvel_db/do.py",
#		long = get_long_assembly_input
#	output:
#		ok = "results/{sample}/assembly/marvel/raw/flappie/marvel.patching.ok",
##		dir = directory("results/{sample}/assembly/marvel/raw/flappie/{sample}")
#	log:
#		stdout = "results/{sample}/logs/marvel.raw.flappie.stdout.txt",
#		stderr = "results/{sample}/logs/marvel.raw.flappie.stderr.txt"
#	benchmark:
#		"results/{sample}/benchmarks/marvel.raw.flappie.benchmark.txt"
#	params:
#		wd = os.getcwd(),
#		sample = "{sample}",
#		genome_size = "200m",
#		minlen = "500",
#		dir = "results/{sample}/assembly/marvel/raw/flappie/{sample}",
#		options = "" #
#	singularity: "docker://chrishah/marvel:e3f3cae"
##	shadow: "minimal"
#	threads: 64
#	resources:
#		mem_gb=90
#	shell:
#		"""
#		home=$(pwd)
#		mkdir -p {params.dir}
#		zcat {input.long} | perl -ne 'chomp; $h=$_; $s=<>; $p=<>; $q=<>; $h =~ s/^@/>/; print "$h\\n$s"' > {params.dir}/long.fasta
##		cat {input.do} | \
##		sed 's/PARALLEL   =.*/PARALLEL   = {threads}/' | \
##		sed 's/DB         =.*/DB         = {wildcards.sample}/' > {params.dir}/do.py
#		cd {params.dir}
#		
#		DBprepare.py -j {threads} -x {params.minlen} {wildcards.sample} long.fasta 1> $home/{log.stdout} 2> $home/{log.stderr}
##		python3 {input.do} 1>> $home/{log.stdout} 2>> $home/{log.stderr}
#		#steps more or less from here: https://github.com/schloi/MARVEL/tree/master/examples/smed
#		
#		### Read patching phase
#		# create daligner and merge plans
#		HPCdaligner -v -t 100 -r1 -j{threads} --dal 32 --mrg 32 -o {wildcards.sample} {wildcards.sample}
#
#		# run the plans generated by HPCdaligner
#		echo -e "[$(date)]\tRunning and merging daligner jobs" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
#		bash ./{wildcards.sample}.dalign.plan 1>> $home/{log.stdout} 2>> $home/{log.stderr}
#		bash ./{wildcards.sample}.merge.plan 1>> $home/{log.stdout} 2>> $home/{log.stderr}
#
#		# compute quality and trim information for each block
#		echo -e "[$(date)]\tcompute quality and trim information for each block" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
#		for block in $(ls -1 {wildcards.sample}.*.las | sed 's/^{wildcards.sample}.//' | sed 's/.las$//')
#		do 
#			LAq -b $block {wildcards.sample} {wildcards.sample}.$block.las
#		done 1>> $home/{log.stdout} 2>> $home/{log.stderr}
#
#		# merge quality and trim tracks
#		echo -e "[$(date)]\tmerge quality and trim tracks" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
#		TKmerge -d {wildcards.sample} q 1>> $home/{log.stdout} 2>> $home/{log.stderr}
#		TKmerge -d {wildcards.sample} trim 1>> $home/{log.stdout} 2>> $home/{log.stderr}
#
#		echo -e "[$(date)]\tPatching sequence errors" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
#		for block in $(ls -1 {wildcards.sample}.*.las | sed 's/^{wildcards.sample}.//' | sed 's/.las$//')
#		do 
#			LAfix -x {params.minlen} {wildcards.sample} {wildcards.sample}.$block.las {wildcards.sample}.$block.fixed.fasta
#		done 1>> $home/{log.stdout} 2>> $home/{log.stderr}
#
#		touch ../marvel.patching.ok
#		"""

rule marvel_pat_init:
	input:
		long = get_long_assembly_input
	output:
		ok = "results/{sample}/assembly/marvel/raw/flappie/marvel.patching.initiate.ok",
	log:
		stdout = "results/{sample}/logs/marvel.raw.flappie.patching.initiate.stdout.txt",
		stderr = "results/{sample}/logs/marvel.raw.flappie.patching.initiate.stderr.txt"
	benchmark:
		"results/{sample}/benchmarks/marvel.raw.flappie.patching.initiate.benchmark.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		genome_size = "200m",
		minlen = "500",
		dir = "results/{sample}/assembly/marvel/raw/flappie/{sample}",
		threads = 32, #this is the number of threads set for subsequent daligner jobs, not the number of threads used by this rule
		blocksize = 200, #blocksize in megabases - with 200, daligner needs about 25G of RAM in the first testdataset
		options = "" #
	singularity: "docker://chrishah/marvel:e3f3cae"
#	shadow: "minimal"
	threads: 2
	resources:
		mem_gb=90
	shell:
		"""
		home=$(pwd)
		mkdir -p {params.dir}
		zcat {input.long} | perl -ne 'chomp; $h=$_; $s=<>; $p=<>; $q=<>; $h =~ s/^@/>/; print "$h\\n$s"' > {params.dir}/long.fasta
		cd {params.dir}
		
		DBprepare.py -j {threads} -x {params.minlen} -s {params.blocksize} {wildcards.sample} long.fasta 1> $home/{log.stdout} 2> $home/{log.stderr}
		#steps more or less from here: https://github.com/schloi/MARVEL/tree/master/examples/smed
		
		### Read patching phase
		# create daligner and merge plans
		HPCdaligner -v -t 100 -A -r1 -j{params.threads} --dal 32 --mrg 32 -o {wildcards.sample} {wildcards.sample}

		#extract base daligner command
		basecmd=$(cat {wildcards.sample}.dalign.plan | sed 's/{wildcards.sample}/\\t/' | cut -f 1 | sort -n | uniq)
		#find highest block
		last=$(cat {wildcards.sample}.dalign.plan | tail -n 1 | sed 's/{wildcards.sample}\./\\t/g' | cut -f 2)
		#create daligner cmds (considering the -A option, i.e. only do one way comparisons of blocks, appropriately)
		for i in $(seq 1 1 $last); do echo -en "$basecmd {wildcards.sample}.$i "; for j in $(seq 1 1 $last); do echo -en "{wildcards.sample}.$j "; done; echo; done > daligner.jobs
		#create quality trim cmds
		for i in $(seq 1 1 $last); do echo -e "LAq -b $i {wildcards.sample} {wildcards.sample}.$i.las"; done > LAq.jobs

		#create patching jobs
		for i in $(seq 1 1 $last); do echo -e "LAfix -x {params.minlen} {wildcards.sample} {wildcards.sample}.$i.las {wildcards.sample}.$i.fixed.fasta"; done > LAfix.jobs

		echo -e "[$(date)]\tDone!" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		touch $home/{output.ok}
		"""

rule marvel_pat_daligner_quality_trim:
	input:
		rules.marvel_pat_init.output.ok
	output:
		ok = "results/{sample}/assembly/marvel/raw/flappie/marvel.patching.daligner-quality-trim.{batch}.ok",
	log:
		stdout = "results/{sample}/logs/marvel.raw.flappie.patching.daligner-quality-trim.{batch}.stdout.txt",
		stderr = "results/{sample}/logs/marvel.raw.flappie.patching.daligner-quality-trim.{batch}.stderr.txt"
	benchmark:
		"results/{sample}/benchmarks/marvel.raw.flappie.patching.daligner-quality-trim.{batch}.benchmark.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		dir = "results/{sample}/assembly/marvel/raw/flappie/{sample}",
		n_batches = marvel_batches_n,
		batch = "{batch}",
		options = "" #
	singularity: "docker://chrishah/marvel:e3f3cae"
#	shadow: "minimal"
	threads: 32
	resources:
		mem_gb=30
	shell:
		"""
		home=$(pwd)
		cd {params.dir}
		# run the plans generated by HPCdaligner
		echo -e "[$(date)]\tRunning and merging daligner jobs" 1>> $home/{log.stdout} 2>> $home/{log.stderr}

		bash <(sed -n '{wildcards.batch}~{params.n_batches}p' daligner.jobs) 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		bash <(sed -n '{wildcards.batch}~{params.n_batches}p' {wildcards.sample}.merge.plan) 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		
		# compute quality and trim information for each block
		echo -e "[$(date)]\tcompute quality and trim information for each block" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		bash <(sed -n '{wildcards.batch}~{params.n_batches}p' LAq.jobs) 1>> $home/{log.stdout} 2>> $home/{log.stderr}

		echo -e "[$(date)]\tDone!" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		touch $home/{output.ok}	
		"""

rule marvel_pat_merge:
	input:
		lambda wildcards: expand(rules.marvel_pat_daligner_quality_trim.output.ok, sample="{sample}", batch=marvel_batch_list)
	output:
		ok = "results/{sample}/assembly/marvel/raw/flappie/marvel.patching.merge.ok",
	log:
		stdout = "results/{sample}/logs/marvel.raw.flappie.patching.merge.stdout.txt",
		stderr = "results/{sample}/logs/marvel.raw.flappie.patching.merge.stderr.txt"
	benchmark:
		"results/{sample}/benchmarks/marvel.raw.flappie.patching.merge.benchmark.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		dir = "results/{sample}/assembly/marvel/raw/flappie/{sample}",
		options = "" #
	singularity: "docker://chrishah/marvel:e3f3cae"
#	shadow: "minimal"
	threads: 2
	resources:
		mem_gb=30
	shell:
		"""
		home=$(pwd)
		cd {params.dir}
		# merge quality and trim tracks
		echo -e "[$(date)]\tmerge quality and trim tracks" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		TKmerge -d {wildcards.sample} q 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		TKmerge -d {wildcards.sample} trim 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		echo -e "[$(date)]\tDone!" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		touch $home/{output.ok}
		"""

rule marvel_pat_patch:
	input:
		rules.marvel_pat_merge.output
	output:
		ok = "results/{sample}/assembly/marvel/raw/flappie/marvel.patching.patching.{batch}.ok",
	log:
		stdout = "results/{sample}/logs/marvel.raw.flappie.patching.patching.{batch}.stdout.txt",
		stderr = "results/{sample}/logs/marvel.raw.flappie.patching.patching.{batch}.stderr.txt"
	benchmark:
		"results/{sample}/benchmarks/marvel.raw.flappie.patching.patching.{batch}.benchmark.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		genome_size = "200m",
		minlen = "500",
		dir = "results/{sample}/assembly/marvel/raw/flappie/{sample}",
		n_batches = marvel_batches_n,
		batch = "{batch}",
		options = "" #
	singularity: "docker://chrishah/marvel:e3f3cae"
#	shadow: "minimal"
	threads: 2
	resources:
		mem_gb=6
	shell:
		"""
		home=$(pwd)
		cd {params.dir}
		echo -e "[$(date)]\tPatching sequence errors" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		bash <(sed -n '{wildcards.batch}~{params.n_batches}p' LAfix.jobs) 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		echo -e "[$(date)]\tDone!" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		touch $home/{output.ok}
		"""

rule marvel_ass_init:
	input:
		lambda wildcards: expand(rules.marvel_pat_patch.output, sample="{sample}", batch=marvel_batch_list)
	output:
		ok = "results/{sample}/assembly/marvel/raw/flappie/marvel.assembly.initate.ok",
	log:
		stdout = "results/{sample}/logs/marvel.raw.flappie.assembly.initiate.stdout.txt",
		stderr = "results/{sample}/logs/marvel.raw.flappie.assembly.initiate.stderr.txt"
	benchmark:
		"results/{sample}/benchmarks/marvel.raw.flappie.assembly.initiate.benchmark.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		genome_size = "200m",
		minlen = "500",
		dir = "results/{sample}/assembly/marvel/raw/flappie/{sample}",
		threads = "32",
		blocksize = 200, #blocksize in megabases - with 200, daligner needs about 25G of RAM in the first testdataset
		track = "source",
		stitch_distance = 50, #LAstitch option -f (maximum stitch distance)
		expected_coverage = "-1", #-1 means 'auto-detect'
		above_multiple_repeat = "2.0", #LArepeat -h option
		below_multiple_repeat = "1.5", #LArepeat -l option
		exclude_stitch_dist = "300", #LAgap -s option (don't count alignments that would be stitchable with a maximum distance of n)
		options = "" #
	singularity: "docker://chrishah/marvel:e3f3cae"
#	shadow: "minimal"
	threads: 2
	resources:
		mem_gb=2
	shell:
		"""
		home=$(pwd)
		cd {params.dir}

		#cleanup
		rm -rf 001*
		rm {wildcards.sample}.*.las
		cat {wildcards.sample}.*.fixed.fasta > {wildcards.sample}.fixed.fasta
		rm {wildcards.sample}.*.fixed.fasta

		echo -e "[$(date)]\tcreate a new Database of fixed reads" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		# create a new Database of fixed reads
		FA2db -v -x {params.minlen} -c {params.track} {wildcards.sample}_FIX {wildcards.sample}.fixed.fasta 1>> $home/{log.stdout} 2>> $home/{log.stderr}

		# split database into blocks
		echo -e "[$(date)]\tsplit database into blocks" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		DBsplit -s {params.blocksize} {wildcards.sample}_FIX 1>> $home/{log.stdout} 2>> $home/{log.stderr}

		# create dust track
		echo -e "[$(date)]\tcreate dust track" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		DBdust {wildcards.sample}_FIX 1>> $home/{log.stdout} 2>> $home/{log.stderr}

		# create plans
		echo -e "[$(date)]\tcreate plans" 1>> $home/{log.stdout} 2>> $home/{log.stderr} 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		HPCdaligner -v -t 100 -A -r2 -j{params.threads} --dal 32 --mrg 32 -o {wildcards.sample}_FIX {wildcards.sample}_FIX 1>> $home/{log.stdout} 2>> $home/{log.stderr}

		#extract base daligner command
		basecmd=$(cat {wildcards.sample}_FIX.dalign.plan | sed 's/{wildcards.sample}_FIX/\\t/' | cut -f 1 | sort -n | uniq)
		#find highest block
		last=$(cat {wildcards.sample}_FIX.dalign.plan | tail -n 1 | sed 's/{wildcards.sample}_FIX\./\\t/g' | cut -f 2)
		#create daligner cmds (considering the -A option, i.e. only do one way comparisons of blocks, appropriately)
		for i in $(seq 1 1 $last); do echo -en "$basecmd {wildcards.sample}_FIX.$i "; for j in $(seq 1 1 $last); do echo -en "{wildcards.sample}_FIX.$j "; done; echo; done > daligner_FIX.jobs
		#create repair plan
		for i in $(seq 1 1 $last); do echo -e "LAstitch -f {params.stitch_distance} -L {wildcards.sample}_FIX {wildcards.sample}_FIX.$i.las {wildcards.sample}_FIX.$i.stitch.las"; done > LAstitch_FIX.jobs

		#create quality and trim plan
		for i in $(seq 1 1 $last); do echo -e "LAq -b $i -s 5 -T trim0 {wildcards.sample}_FIX {wildcards.sample}_FIX.$i.stitch.las"; done > LAq_FIX.jobs

		#LAfilter trim0
		for i in $(seq 1 1 $last); do echo -e "LAfilter -o 4000 -L -p -t trim0 {wildcards.sample}_FIX {wildcards.sample}_FIX.$i.stitch.las {wildcards.sample}_FIX.$i.ol.las"; done > LAfilter_FIX_trim0.jobs
	
		#create repeat annotation plan
		for i in $(seq 1 1 $last); do echo -e "LArepeat -c {params.expected_coverage} -l {params.below_multiple_repeat} -h {params.above_multiple_repeat} -t repeats -b $i {wildcards.sample}_FIX {wildcards.sample}_FIX.$i.ol.las"; done > LArepeat_FIX.jobs
		
		#create repeat homogenization plan
		for i in $(seq 1 1 $last); do echo -e "TKhomogenize -i repeats -I hrepeats -b $i {wildcards.sample}_FIX {wildcards.sample}_FIX.$i.ol.las"; done > TKhomogenize.jobs

		#create remove gaps plan
		for i in $(seq 1 1 $last); do echo -e "LAgap -s {params.exclude_stitch_dist} -L -t trim0 {wildcards.sample}_FIX {wildcards.sample}_FIX.$i.ol.las {wildcards.sample}_FIX.$i.gap.las"; done > LAgaps_FIX.jobs
		#create recalculate trim track based on cleaned up gaps plan
		for i in $(seq 1 1 $last); do echo -e "LAq -s 5 -u -t trim0 -T trim1 -b $i {wildcards.sample}_FIX {wildcards.sample}_FIX.$i.gap.las"; done > LAq_FIX_gap.jobs
		
		#create filter repeat induced alignments and resolve repeat modules plan
		for i in $(seq 1 1 $last); do echo -e "LAfilter -n 100 -u 0 -o 4200 -S 300 -L -p -T -t trim1 -r frepeats {wildcards.sample}_FIX {wildcards.sample}_FIX.$i.gap.las {wildcards.sample}_FIX.$i.filtered.las"; done > LAfilter_FIX.jobs
		for i in $(seq 1 1 $last); do echo -e "LAfilter -n 100 -u 0 -o 4200 -S 300 -M {params.expected_coverage} -L -p -T -t trim1 -r frepeats {wildcards.sample}_FIX {wildcards.sample}_FIX.$i.gap.las {wildcards.sample}_FIX.$i.filtered_M.las"; done > LAfilter_FIX_M.jobs
		for i in $(seq 1 1 $last); do echo -e "LAfilter -n 100 -u 0 -o 4200 -S 300 -MM {params.expected_coverage} -L -p -T -t trim1 -r frepeats {wildcards.sample}_FIX {wildcards.sample}_FIX.$i.gap.las {wildcards.sample}_FIX.$i.filtered_MM.las"; done > LAfilter_FIX_MM.jobs
		
		echo -e "[$(date)]\tDone!" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		touch $home/{output.ok}
		"""

rule marvel_ass_dali_rep_trim:
	input:
		rules.marvel_ass_init.output
	output:
		ok = "results/{sample}/assembly/marvel/raw/flappie/marvel.assembly.dalign-repair-trim.{batch}.ok",
	log:
		stdout = "results/{sample}/logs/marvel.raw.flappie.assembly.dalign-repair-trim.{batch}.stdout.txt",
		stderr = "results/{sample}/logs/marvel.raw.flappie.assembly.dalign-repair-trim.{batch}.stderr.txt"
	benchmark:
		"results/{sample}/benchmarks/marvel.raw.flappie.assembly.dalign-repair-trim.{batch}.benchmark.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		dir = "results/{sample}/assembly/marvel/raw/flappie/{sample}",
		n_batches = marvel_batches_n,
		batch = "{batch}",
		options = "" #
	singularity: "docker://chrishah/marvel:e3f3cae"
#	shadow: "minimal"
	threads: 32
	resources:
		mem_gb=30
	shell:
		"""
		home=$(pwd)
		cd {params.dir}

		echo -e "[$(date)]\tRunning daligner" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		bash <(sed -n '{wildcards.batch}~{params.n_batches}p' daligner_FIX.jobs) 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		
		echo -e "[$(date)]\tRunning merging" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		bash <(sed -n '{wildcards.batch}~{params.n_batches}p' {wildcards.sample}_FIX.merge.plan) 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		
		echo -e "[$(date)]\tRunning repair alignments jobs" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		bash <(sed -n '{wildcards.batch}~{params.n_batches}p' LAstitch_FIX.jobs) 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		
		echo -e "[$(date)]\tCompute quality and trim information for each block" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		bash <(sed -n '{wildcards.batch}~{params.n_batches}p' LAq_FIX.jobs) 1>> $home/{log.stdout} 2>> $home/{log.stderr}

		echo -e "[$(date)]\tDone!" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		touch $home/{output.ok}
		"""

rule marvel_ass_merge_qual_trim:
	input:
		lambda wildcards: expand(rules.marvel_ass_dali_rep_trim.output, sample="{sample}", batch=marvel_batch_list)
	output:
		ok = "results/{sample}/assembly/marvel/raw/flappie/marvel.assembly.merge-qual-trim.ok",
	log:
		stdout = "results/{sample}/logs/marvel.raw.flappie.assembly.merge-qual-trim.stdout.txt",
		stderr = "results/{sample}/logs/marvel.raw.flappie.assembly.merge-qual-trim.stderr.txt"
	benchmark:
		"results/{sample}/benchmarks/marvel.raw.flappie.assembly.merge-qual-trim.benchmark.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		dir = "results/{sample}/assembly/marvel/raw/flappie/{sample}",
		options = "" #
	singularity: "docker://chrishah/marvel:e3f3cae"
#	shadow: "minimal"
	threads: 2
	resources:
		mem_gb=4
	shell:
		"""
		home=$(pwd)
		cd {params.dir}

		echo -e "[$(date)]\tMerge quality and trim tracks" 1>> $home/{log.stdout} 2>> $home/{log.stderr}

		TKmerge -d {wildcards.sample}_FIX q
		TKmerge -d {wildcards.sample}_FIX trim0

		echo -e "[$(date)]\tDone!" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		touch $home/{output.ok}
		"""

rule marvel_ass_rep_annot:
	input:
		rules.marvel_ass_merge_qual_trim.output
	output:
		ok = "results/{sample}/assembly/marvel/raw/flappie/marvel.assembly.repeat-annotation.{batch}.ok",
	log:
		stdout = "results/{sample}/logs/marvel.raw.flappie.assembly.repeat-annotation.{batch}.stdout.txt",
		stderr = "results/{sample}/logs/marvel.raw.flappie.assembly.repeat-annotation.{batch}.stderr.txt"
	benchmark:
		"results/{sample}/benchmarks/marvel.raw.flappie.assembly.repeat-annotation.{batch}.benchmark.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		dir = "results/{sample}/assembly/marvel/raw/flappie/{sample}",
		n_batches = marvel_batches_n,
		batch = "{batch}",
		options = "" #
	singularity: "docker://chrishah/marvel:e3f3cae"
#	shadow: "minimal"
	threads: 2
	resources:
		mem_gb=8
	shell:
		"""
		home=$(pwd)
		cd {params.dir}

		echo -e "[$(date)]\tLAfilter trim0" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		bash <(sed -n '{wildcards.batch}~{params.n_batches}p' LAfilter_FIX_trim0.jobs) 1>> $home/{log.stdout} 2>> $home/{log.stderr}

		echo -e "[$(date)]\tCreate repeat annotatino for block {params.batch}" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		bash <(sed -n '{wildcards.batch}~{params.n_batches}p' LArepeat_FIX.jobs) 1>> $home/{log.stdout} 2>> $home/{log.stderr}

		echo -e "[$(date)]\tDone!" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		touch $home/{output.ok}
		"""
rule marvel_ass_merge_repeats:
	input:
		lambda wildcards: expand(rules.marvel_ass_rep_annot.output, sample="{sample}", batch=marvel_batch_list)
	output:
		ok = "results/{sample}/assembly/marvel/raw/flappie/marvel.assembly.merge-repeats.ok",
	log:
		stdout = "results/{sample}/logs/marvel.raw.flappie.assembly.merge-repeats.stdout.txt",
		stderr = "results/{sample}/logs/marvel.raw.flappie.assembly.merge-repeats.stderr.txt"
	benchmark:
		"results/{sample}/benchmarks/marvel.raw.flappie.assembly.merge-repeats.benchmark.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		dir = "results/{sample}/assembly/marvel/raw/flappie/{sample}",
		options = "" #
	singularity: "docker://chrishah/marvel:e3f3cae"
#	shadow: "minimal"
	threads: 2
	resources:
		mem_gb=6
	shell:
		"""
		home=$(pwd)
		cd {params.dir}

		echo -e "[$(date)]\tMerge quality and trim tracks" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		TKmerge -d {wildcards.sample}_FIX repeats
		
		echo -e "[$(date)]\tDone!" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		touch $home/{output.ok}
		"""

rule marvel_ass_hom_rep:
	input:
		rules.marvel_ass_merge_repeats.output
	output:
		ok = "results/{sample}/assembly/marvel/raw/flappie/marvel.assembly.homogenize-repeat-track.{batch}.ok",
	log:
		stdout = "results/{sample}/logs/marvel.raw.flappie.assembly.homogenize-repeat-track.{batch}.stdout.txt",
		stderr = "results/{sample}/logs/marvel.raw.flappie.assembly.homogenize-repeat-track.{batch}.stderr.txt"
	benchmark:
		"results/{sample}/benchmarks/marvel.raw.flappie.assembly.homogenize-repeat-track.{batch}.benchmark.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		dir = "results/{sample}/assembly/marvel/raw/flappie/{sample}",
		n_batches = marvel_batches_n,
		batch = "{batch}",
		options = "" #
	singularity: "docker://chrishah/marvel:e3f3cae"
#	shadow: "minimal"
	threads: 2
	resources:
		mem_gb=8
	shell:
		"""
		home=$(pwd)
		cd {params.dir}

		echo -e "[$(date)]\thomogenize repeat track for block {params.batch}" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		bash <(sed -n '{wildcards.batch}~{params.n_batches}p' TKhomogenize.jobs) 1>> $home/{log.stdout} 2>> $home/{log.stderr}

		echo -e "[$(date)]\tDone!" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		touch $home/{output.ok}
		"""

rule marvel_ass_combine_rem_dup_over:
	input:
		lambda wildcards: expand(rules.marvel_ass_hom_rep.output, sample="{sample}", batch=marvel_batch_list)
	output:
		ok = "results/{sample}/assembly/marvel/raw/flappie/marvel.assembly.combine-remove-duplicate-overlap-annotations.ok",
	log:
		stdout = "results/{sample}/logs/marvel.raw.flappie.assembly.combine-remove-duplicate-overlap-annotations.stdout.txt",
		stderr = "results/{sample}/logs/marvel.raw.flappie.assembly.combine-remove-duplicate-overlap-annotations.stderr.txt"
	benchmark:
		"results/{sample}/benchmarks/marvel.raw.flappie.assembly.combine-remove-duplicate-overlap-annotations.benchmark.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		dir = "results/{sample}/assembly/marvel/raw/flappie/{sample}",
		options = "" #
	singularity: "docker://chrishah/marvel:e3f3cae"
#	shadow: "minimal"
	threads: 2
	resources:
		mem_gb=6
	shell:
		"""
		home=$(pwd)
		cd {params.dir}

		echo -e "[$(date)]\tcombine (remove duplicate and overlapping annotations)" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		TKcombine -d {wildcards.sample}_FIX hrepeats \#.hrepeats 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		TKcombine {wildcards.sample}_FIX frepeats repeats hrepeats 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		
		echo -e "[$(date)]\tDone!" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		touch $home/{output.ok}
		"""

rule marvel_ass_gap_quality:
	input:
		rules.marvel_ass_combine_rem_dup_over.output
	output:
		ok = "results/{sample}/assembly/marvel/raw/flappie/marvel.assembly.remove-gaps-quality.{batch}.ok",
	log:
		stdout = "results/{sample}/logs/marvel.raw.flappie.assembly.remove-gaps-quality.{batch}.stdout.txt",
		stderr = "results/{sample}/logs/marvel.raw.flappie.assembly.remove-gaps-quality.{batch}.stderr.txt"
	benchmark:
		"results/{sample}/benchmarks/marvel.raw.flappie.assembly.remove-gaps-quality.{batch}.benchmark.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		dir = "results/{sample}/assembly/marvel/raw/flappie/{sample}",
		n_batches = marvel_batches_n,
		batch = "{batch}",
		options = "" #
	singularity: "docker://chrishah/marvel:e3f3cae"
#	shadow: "minimal"
	threads: 2
	resources:
		mem_gb=8
	shell:
		"""
		home=$(pwd)
		cd {params.dir}

		echo -e "[$(date)]\tremove gaps (ie. due to chimeric reads, ...) for block {params.batch}" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		bash <(sed -n '{wildcards.batch}~{params.n_batches}p' LAgaps_FIX.jobs) 1>> $home/{log.stdout} 2>> $home/{log.stderr}

		echo -e "[$(date)]\trecalculate the trim track based on cleanup gaps for block {params.batch}" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		bash <(sed -n '{wildcards.batch}~{params.n_batches}p' LAq_FIX_gap.jobs) 1>> $home/{log.stdout} 2>> $home/{log.stderr}

		echo -e "[$(date)]\tDone!" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		touch $home/{output.ok}
		"""

rule marvel_ass_merge_trim1:
	input:
		lambda wildcards: expand(rules.marvel_ass_gap_quality.output, sample="{sample}", batch=marvel_batch_list)
	output:
		ok = "results/{sample}/assembly/marvel/raw/flappie/marvel.assembly.merge-trim1.ok",
	log:
		stdout = "results/{sample}/logs/marvel.raw.flappie.assembly.merge-trim1.stdout.txt",
		stderr = "results/{sample}/logs/marvel.raw.flappie.assembly.merge-trim1.stderr.txt"
	benchmark:
		"results/{sample}/benchmarks/marvel.raw.flappie.assembly.merge-trim1.benchmark.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		dir = "results/{sample}/assembly/marvel/raw/flappie/{sample}",
		options = "" #
	singularity: "docker://chrishah/marvel:e3f3cae"
#	shadow: "minimal"
	threads: 2
	resources:
		mem_gb=6
	shell:
		"""
		home=$(pwd)
		cd {params.dir}

		echo -e "[$(date)]\tmerge trim1" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		TKmerge -d {wildcards.sample}_FIX trim1

		echo -e "[$(date)]\tDone!" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		touch $home/{output.ok}
		"""

rule marvel_ass_filter_rep_ali:
	input:
		rules.marvel_ass_merge_trim1.output
	output:
		ok = "results/{sample}/assembly/marvel/raw/flappie/marvel.assembly.filter-repeat-alignments.{batch}.ok",
	log:
		stdout = "results/{sample}/logs/marvel.raw.flappie.assembly.filter-repeat-alignments.{batch}.stdout.txt",
		stderr = "results/{sample}/logs/marvel.raw.flappie.assembly.filter-repeat-alignments.{batch}.stderr.txt"
	benchmark:
		"results/{sample}/benchmarks/marvel.raw.flappie.assembly.filter-repeat-alignments.{batch}.benchmark.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		dir = "results/{sample}/assembly/marvel/raw/flappie/{sample}",
		n_batches = marvel_batches_n,
		batch = "{batch}",
		options = "" #
	singularity: "docker://chrishah/marvel:e3f3cae"
#	shadow: "minimal"
	threads: 2
	resources:
		mem_gb=8
	shell:
		"""
		home=$(pwd)
		cd {params.dir}

		echo -e "[$(date)]\tfilter repeat induced alignments and try to resolve repeat modules for block {params.batch}" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		bash <(sed -n '{wildcards.batch}~{params.n_batches}p' LAfilter_FIX.jobs) 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		bash <(sed -n '{wildcards.batch}~{params.n_batches}p' LAfilter_FIX_M.jobs) 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		bash <(sed -n '{wildcards.batch}~{params.n_batches}p' LAfilter_FIX_MM.jobs) 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		echo -e "[$(date)]\tDone!" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		touch $home/{output.ok}
		"""

rule marvel_ass_merge_filtered:
	input:
		lambda wildcards: expand(rules.marvel_ass_filter_rep_ali.output, sample="{sample}", batch=marvel_batch_list)
	output:
		ok = "results/{sample}/assembly/marvel/raw/flappie/marvel.assembly.merge-filtered.ok",
	log:
		stdout = "results/{sample}/logs/marvel.raw.flappie.assembly.merge-filtered.stdout.txt",
		stderr = "results/{sample}/logs/marvel.raw.flappie.assembly.merge-filtered.stderr.txt"
	benchmark:
		"results/{sample}/benchmarks/marvel.raw.flappie.assembly.merge-filtered.benchmark.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		dir = "results/{sample}/assembly/marvel/raw/flappie/{sample}",
		options = "" #
	singularity: "docker://chrishah/marvel:e3f3cae"
#	shadow: "minimal"
	threads: 2
	resources:
		mem_gb=6
	shell:
		"""
		home=$(pwd)
		cd {params.dir}

		echo -e "[$(date)]\tmerge trim1" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		LAmerge -S filtered {wildcards.sample}_FIX {wildcards.sample}_FIX.filtered.las 1>> $home/{log.stdout} 2>> $home/{log.stderr}
#		LAmerge -S filtered_M {wildcards.sample}_FIX {wildcards.sample}_FIX.filtered_M.las 1>> $home/{log.stdout} 2>> $home/{log.stderr}
#		LAmerge -S filtered_MM {wildcards.sample}_FIX {wildcards.sample}_FIX.filtered_MM.las 1>> $home/{log.stdout} 2>> $home/{log.stderr}

		echo -e "[$(date)]\tCreate overlap graphs" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		if [ -d filtered ]; then rm -rf filtered; fi
		mkdir filtered
		OGbuild -t trim1 -s {wildcards.sample}_FIX {wildcards.sample}_FIX.filtered.las filtered/{wildcards.sample}_FIX 1>> $home/{log.stdout} 2>> $home/{log.stderr}
#		if [ -d filtered_M ]; then rm -rf filtered_M; fi
#		mkdir filtered_M
#		OGbuild -t trim1 -s {wildcards.sample}_FIX {wildcards.sample}_FIX.filtered_M.las filtered_M/{wildcards.sample}_FIX 1>> $home/{log.stdout} 2>> $home/{log.stderr}
#		if [ -d filtered_MM ]; then rm -rf filtered_MM; fi
#		mkdir filtered_MM
#		OGbuild -t trim1 -s {wildcards.sample}_FIX {wildcards.sample}_FIX.filtered_MM.las filtered_MM/{wildcards.sample}_FIX 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		
		echo -e "[$(date)]\tprepare graph touring jobs" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		for graphml in $(ls -1 filtered/*/*graphml)
		do
			echo -e "OGtour.py -c {wildcards.sample}_FIX $graphml"
		done > tour.jobs
#		for graphml in $(ls -1 filtered_M/*/*graphml)
#		do
#			echo -e "OGtour.py -c {wildcards.sample}_FIX $graphml"
#		done > tour_M.jobs
#		for graphml in $(ls -1 filtered_MM/*/*.graphml)
#		do
#			echo -e "OGtour.py -c {wildcards.sample}_FIX $graphml"
#		done > tour_MM.jobs

		echo -e "[$(date)]\tprepare convert the tours to sequence jobs" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
#		for graphml in $(ls -1 filtered/*/*.graphml)
		for graphml in $(ls -1 filtered/*/*.tour.graphml)
		do
#			path=$(echo -e "$graphml" | sed 's/graphml$/tour.paths/')
			path=$(echo -e "$graphml" | sed 's/graphml$/paths/')
			echo -e "if [ -f $path ]; then tour2fasta.py -t trim1 {wildcards.sample}_FIX $graphml $path; fi"
		done > convert.jobs


#		for graphml in $(ls -1 filtered/*/*.graphml)
#		do
#			path=$(echo -e "$graphml" | sed 's/graphml$/tour.paths/')
#			echo -e "tour2fasta.py -t trim1 {wildcards.sample}_FIX $graphml $path"
#		done > convert.jobs
			
#		for graphml in $(ls -1 filtered_M/*/*.graphml)
#		do
#			path=$(echo -e "$graphml" | sed 's/graphml$/tour.paths/')
#			echo -e "tour2fasta.py -t trim1 {wildcards.sample}_FIX $graphml $path"
#		done > convert_M.jobs
#
#		for graphml in $(ls -1 filtered_MM/*/*.graphml)
#		do
#			path=$(echo -e "$graphml" | sed 's/graphml$/tour.paths/')
#			echo -e "tour2fasta.py -t trim1 {wildcards.sample}_FIX $graphml $path"
#		done > convert_MM.jobs


		echo -e "[$(date)]\tDone!" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		touch $home/{output.ok}
		"""











rule marvel_ass_tour_convert:
	input:
		rules.marvel_ass_merge_filtered.output
	output:
		ok = "results/{sample}/assembly/marvel/raw/flappie/marvel.assembly.tour-convert.{batch}.ok",
	log:
		stdout = "results/{sample}/logs/marvel.raw.flappie.assembly.tour-convert.{batch}.stdout.txt",
		stderr = "results/{sample}/logs/marvel.raw.flappie.assembly.tour-convert.{batch}.stderr.txt"
	benchmark:
		"results/{sample}/benchmarks/marvel.raw.flappie.assembly.tour-convert.{batch}.benchmark.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		dir = "results/{sample}/assembly/marvel/raw/flappie/{sample}",
		n_batches = marvel_batches_n,
		batch = "{batch}",
		options = "" #
	singularity: "docker://chrishah/marvel:e3f3cae"
#	shadow: "minimal"
	threads: 2
	resources:
		mem_gb=8
	shell:
		"""
		home=$(pwd)
		cd {params.dir}

#		echo -e "[$(date)]\ttouring the graph for batch {params.batch}" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
#		bash <(sed -n '{wildcards.batch}~{params.n_batches}p' tour.jobs) 1>> $home/{log.stdout} 2>> $home/{log.stderr}

		echo -e "[$(date)]\tconverting tours to sequences for batch {params.batch}" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		bash <(sed -n '{wildcards.batch}~{params.n_batches}p' convert.jobs) 1>> $home/{log.stdout} 2>> $home/{log.stderr}

		echo -e "[$(date)]\tDone!" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		touch $home/{output.ok}
		"""

rule marvel_ass_finale:
	input:
		lambda wildcards: expand(rules.marvel_ass_tour_convert.output, sample="{sample}", batch=marvel_batch_list)
	output:
		ok = "results/{sample}/assembly/marvel/raw/flappie/marvel.assembly.finale.ok",
	log:
		stdout = "results/{sample}/logs/marvel.raw.flappie.assembly.finale.stdout.txt",
		stderr = "results/{sample}/logs/marvel.raw.flappie.assembly.finale.stderr.txt"
	benchmark:
		"results/{sample}/benchmarks/marvel.raw.flappie.assembly.finale.benchmark.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		dir = "results/{sample}/assembly/marvel/raw/flappie/{sample}",
		options = "" #
	singularity: "docker://chrishah/marvel:e3f3cae"
#	shadow: "minimal"
	threads: 2
	resources:
		mem_gb=6
	shell:
		"""
		home=$(pwd)
		cd {params.dir}

		echo -e "[$(date)]\tgather all sequences" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		cat filtered/*/*.fasta > filtered/{wildcards.sample}.marvel.fasta

		echo -e "[$(date)]\tDone!" 1>> $home/{log.stdout} 2>> $home/{log.stderr}
		touch $home/{output.ok}
		"""




#rule marvel:
#	input:
#		"marvel_db/E_coli.20x.fasta"
#	output:
#		ok = "results/{sample}/assembly/marvel/raw/flappie/marvel.ok",
##		dir = directory("results/{sample}/assembly/marvel/raw/marvel/{sample}")
##	log:
##		stdout = "results/{sample}/logs/marvel.raw.flappie.stdout.txt",
##		stderr = "results/{sample}/logs/marvel.raw.flappie.stderr.txt"
##	benchmark:
##		"results/{sample}/benchmarks/marvel.raw.flappie.benchmark.txt"
#	params:
#		wd = os.getcwd(),
#		sample = "{sample}",
#		genome_size = "200m",
#		dir = "results/{sample}/assembly/marvel/raw/marvel",
#		options = "" #
#	singularity: "docker://chrishah/marvel:e3f3cae"
##	shadow: "minimal"
#	threads: 4
#	resources:
#		mem_gb=5
#	run:
#		"""
#		
#		import multiprocessing
#		import marvel
#		import marvel.config
#		import marvel.queue
#		from pathlib import Path		
#		
#		### settings
#		
#		DB         = "marvel_db/ECOL"
#		COVERAGE   = 25
#		
#		DB_FIX     = DB + "_FIX"
#		PARALLEL   = {threads}
#		
#		q = marvel.queue.queue(DB, COVERAGE, PARALLEL)
#		
#		### patch raw reads
#		
#		### run daligner to create initial overlaps
#		q.plan("{db}.dalign.plan")
#		
#		### run LAmerge to merge overlap blocks
#		q.plan("{db}.merge.plan")
#		
#		# create quality and trim annotation (tracks) for each overlap block
#		q.block("{path}/LAq -d 35 -b {block} {db} {db}.{block}.las")
#		# merge quality and trim tracks
#		q.single("{path}/TKmerge -d {db} q")
#		q.single("{path}/TKmerge -d {db} trim")
#		
#		# run LAfix to patch reads based on overlaps
#		q.block("{path}/LAfix -Q 30 -g -1 {db} {db}.{block}.las {db}.{block}.fixed.fasta")
#		# join all fixed fasta files
#		q.single("!cat {db}.*.fixed.fasta > {db}.fixed.fasta")
#		
#		# create a new Database of fixed reads (-j numOfThreads, -g genome size)
#		q.single("{path_scripts}/DBprepare.py -c source -s 50 -r 2 -j 4 -g 4600000 {db_fixed} {db}.fixed.fasta", db_fixed = DB_FIX)
#		
#		q.process()
#		Path({output.ok}).touch()
#
#		"""
