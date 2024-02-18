#def def_or_data(string,sample):
#	print(string)
#	elements = string.split("|")
#	default = config[elements[0]]
#	print("the dict is: "+str(default))
#	for i in range(1,len(elements)):
#		print("current key is: "+elements[i])
#		current = elements[i]
#		print(current)
#		print(current,i,default[current])
#		default = default[current]
#	
##	print("after loop: "+default)
#	if "-".join(elements) in df:
#		print("-".join(elements)+" is there")
#		print(df.set_index(['sample']))

	
#def_or_data("ont_basecalling|guppy|flowcell","EcStj")
#./demogenas -m assemble -t slurm -s "-pr --until bca_guppy" --clusterconfig=data/testdata/test.cluster.SLURM --configfile=data/config.Eubo -i "-B $DATA -B $BINFL" --select="EcStj,EsKje_Sa_w" --dry

rule bca_guppy:
	input:
		dir = get_fast5_dir
	output:
		fastq = "results/{sample}/reads/ont/guppy/{lib}/{sample}.guppy.{unit}.fastq.gz"
	log:
		stdout = "results/{sample}/logs/guppy.{lib}.{unit}.stdout.txt",
		stderr = "results/{sample}/logs/guppy.{lib}.{unit}.stderr.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		unit = "{unit}",
		dir = "results/{sample}/reads/ont/guppy/{lib}",
		kit = config["ont_basecalling"]["guppy"]["kit"],
		flowcell = config["ont_basecalling"]["guppy"]["flowcell"],
		optional = config["ont_basecalling"]["guppy"]["optional_params"],
		nbatches = config["ont_basecalling"]["concurrency"],
		minlength = config["ont_basecalling"]["minlength"],
	singularity: "docker://chrishah/guppy:6.0.6-8a98bbcbd"
#	shadow: "minimal"
	threads: config["threads"]["guppy"]
	resources:
		mem_gb=10
	shell:
		"""
		cd {params.dir}
		unit={params.unit}
		step={params.nbatches}
		count=0
		echo "$(date)"
		for f in $(ls -1rt {params.wd}/{input.dir} | sed -n "$unit~$step p")
		do
			echo -e "\\n$(date) - Processing file $f (temp suffix: $unit-$count)"
			if [ ! -f $unit-$count.min-{params.minlength}.fastq.gz ]
			then
				if [ -d $unit-$count ]; then rm -rf $unit-$count/*; else mkdir $unit-$count; fi
				echo -e "$(date) - basecalling"
				guppy_basecaller --input_file_list <(echo {params.wd}/{input.dir}/$f) -s $unit-$count --cpu_threads_per_caller $(( {threads} - 1 )) --flowcell {params.flowcell} --kit {params.kit} {params.optional}
				cat $(find ./$unit-$count -name "*fastq" | grep -v "$unit-$count/fail/") | \
					perl -ne '$h=$_; $s=<>; $p=<>; $q=<>; if (length($s) > {params.minlength}){{print "$h$s$p$q"}}' | gzip > $unit-$count.min-{params.minlength}.fastq.gz.tmp
				mv $unit-$count.min-{params.minlength}.fastq.gz.tmp $unit-$count.min-{params.minlength}.fastq.gz
				cat $unit-$count/sequencing_summary.txt >> sequencing_summary-$unit.txt
#				cat $unit-$count/*.log >> {params.wd}/{log.stdout}
				rm -rf $unit-$count
			fi
			let "count+=1"
			echo -e "$(date) - Done"
		done 1> {params.wd}/{log.stdout} 2> {params.wd}/{log.stderr}
		cat $unit-*fastq.gz > $unit.min-{params.minlength}.fastq.gz
		rm $unit-*fastq.gz
		ln -s $unit.min-{params.minlength}.fastq.gz {params.wd}/{output.fastq}
		"""

rule bca_bonito_cpu:
	input:
		dir = get_fast5_dir
	output:
		fastq = "results/{sample}/reads/ont/bonito/{lib}/{sample}.bonito.{unit}.fastq.gz"
	log:
		stdout = "results/{sample}/logs/bonito.{lib}.{unit}.stdout.txt",
		stderr = "results/{sample}/logs/bonito.{lib}.{unit}.stderr.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		unit = "{unit}",
		dir = "results/{sample}/reads/ont/bonito",
		model = "dna_r10.3@v3", #"dna_r9.4.1@v2",
		nbatches = config["ont_basecalling"]["concurrency"],
		minlength = config["ont_basecalling"]["minlength"],
	singularity: "docker://reslp/bonito:0.3.2"
	shadow: "minimal"
	threads: 1
	resources:
		mem_gb=10
	shell:
		"""
		unit={params.unit}
		step={params.nbatches}

		export OPENBLAS_NUM_THREADS={threads}
		for f in $(ls -1rt {input.dir} | sed -n "$unit~$step p")
		do
			echo -e "\\n$(date) - Processing file $f (unit: {params.unit}) "
			mkdir tmp
			ln -s {input.dir}/$f tmp/
			bonito basecaller {params.model} --use_openvino --device=cpu tmp/ --fastq > temp.fastq
			cat temp.fastq | \
				perl -ne '$h=$_; $s=<>; $p=<>; $q=<>; if (length($s) > {params.minlength}){{print "$h$s$p$q"}}' | gzip > {output.fastq}
			echo -e "$(date) - Done"
		done 1> {log.stdout} 2> {log.stderr}

		"""
rule bca_flappie:
	input:
		dir = get_fast5_dir
	output:
		fastq = "results/{sample}/reads/ont/flappie/{lib}/{sample}.flappie.{unit}.fastq.gz"
	log:
		stdout = "results/{sample}/logs/flappie.{lib}.{unit}.stdout.txt",
		stderr = "results/{sample}/logs/flappie.{lib}.{unit}.stderr.txt"
	params:
		wd = os.getcwd(),
		sample = "{sample}",
		unit = "{unit}",
		dir = "results/{sample}/reads/ont/flappie/{lib}",
		model = "r941_native",
		nbatches = config["ont_basecalling"]["concurrency"],
		minlength = config["ont_basecalling"]["minlength"],
	singularity: "docker://reslp/flappie:4de542f"
#	shadow: "minimal"
	threads: 30
	resources:
		mem_gb=30
	shell:
		"""
		cd {params.dir}
		unit={params.unit}
		step={params.nbatches}
		model={params.model}
		export OPENBLAS_NUM_THREADS={threads}
		count=0

#		#this part only applies if the fast5 files already contains only one read per file (perhaps make an if later)
#		if [ -d $unit ]; then rm -rf $unit/*; else mkdir $unit; fi
#		for f in $(ls -1rt {input.dir} | sed -n "$unit~$step p")
#		do
#			ln -s {input.dir}/$f $unit/
#		done
#		flappie --model={params.model} $unit/* | \
#			perl -ne '$h=$_; $s=<>; $p=<>; $q=<>; if (length($s) > {params.minlength}){{print "$h$s$p$q"}}' | gzip > ../../../../../../{output.fastq} 1> ../../../../../../{log.stdout} 2> ../../../../../{log.stderr}
#		rm -rf $unit

		for f in $(ls -1rt {input.dir} | sed -n "$unit~$step p")
		do
			echo -e "\\n$(date) - Processing file $f (temp suffix: $unit-$count)"
			if [ ! -f $unit-$count.min-{params.minlength}.fastq.gz ]
			then
				if [ -d $unit-$count ]; then rm -rf $unit-$count/*; else mkdir $unit-$count; fi
				echo -e "$(date) - converting fast5"
				multi_to_single_fast5 -i {input.dir}/$f -t {threads} -s $unit-$count/
				echo -e "$(date) - basecalling"
				flappie --model={params.model} $unit-$count/*/ | \
					perl -ne '$h=$_; $s=<>; $p=<>; $q=<>; if (length($s) > {params.minlength}){{print "$h$s$p$q"}}' | gzip > $unit-$count.min-{params.minlength}.fastq.gz.tmp
				mv $unit-$count.min-{params.minlength}.fastq.gz.tmp $unit-$count.min-{params.minlength}.fastq.gz
				rm -rf $unit-$count
			fi
			let "count+=1"
			echo -e "$(date) - Done"
		done 1> {params.wd}/{log.stdout} 2> {params.wd}/{log.stderr}

		cat $unit-*fastq.gz > $unit.min-{params.minlength}.fastq.gz
		rm $unit-*fastq.gz
		ln -s $unit.min-{params.minlength}.fastq.gz {params.wd}/{output.fastq}
		"""

rule bca_x_gather_called_ont_reads:
	input: 
		expand("results/{{units.sample}}/reads/ont/{basecaller}/{units.lib}/{{units.sample}}.{basecaller}.{unit}.fastq.gz", basecaller=config["ont_basecalling"]["basecaller"], units=fast5_units.itertuples(), unit=flappie_unit_list)
	output: 
		"results/{sample}/reads/ont/ont-calling.done"
	shell:
		"""
		touch {output}
		"""
