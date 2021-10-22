wildcard_constraints:
	trimmer="trimgalore|None",
	corrector="bless|spades|None",
	merger="usearch|None",
	basecaller="guppy|flappie",
	longcorrection="canucorrect|consent|ratatosk"

#configfile from script

#get stuff from config file
k = config["assemble"]["default_k"]
#ec_concurrency = config["ectools"]["ec_concurrency"]
abyss_ass_to_scaffold = config["abyss"]["assembly_to_scaffold"]
longcor = config["long_correction"]

marvel_batches_n = 20
marvel_batch_list = range(1,marvel_batches_n+1)

#load snakefiles
include: "rules/functions.smk"
#reading in data for samples happens in 'rules/functions.smk'

include: "rules/illumina-process.smk"
include: "rules/illumina-correct.smk"
include: "rules/illumina-merge.smk"
include: "rules/illumina-assemble.smk"
include: "rules/ont-call.smk"
include: "rules/ont-correct.smk"
include: "rules/ont-assemble-hybrid.smk"
include: "rules/ont-assemble-denovo.smk"
include: "rules/eval_assembly.smk"






#different all rules that run the pipeline particular stages
rule eval_illumina:
	input:
		#fastqc
		lambda wildcards: expand("results/{test.sample}/read_qc/fastqc_raw/{test.lib}/{test.sample}.{test.lib}.status.ok", test=illumina_units.itertuples()),
		lambda wildcards: expand("results/{test.sample}/trimming/trimgalore/{test.lib}/{test.sample}.{test.lib}.fastqc.status.ok", test=illumina_units.itertuples()),
		#trim and gather
		lambda wildcards: expand("results/{test.sample}/trimming/trimgalore/{test.sample}-full/{test.sample}.cat.status.ok", test=illumina_units.itertuples()),
		lambda wildcards: expand("results/{test.sample}/plots/{test.sample}-k{k}-distribution-full.pdf", test=illumina_units.itertuples(), k=config["kmc"]["k"])
rule illumina_trim:
	input:
		#trim and gather
		lambda wildcards: expand("results/{test.sample}/trimming/trimgalore/{test.sample}-full/{test.sample}.cat.status.ok", test=illumina_units.itertuples()),

rule correct_illumina:
	input:
		expand(rules.cor_gather_illumina_corrected.output, sample=Illumina_process_df["sample"], trimmer=config["illumina_trimming"])

rule merge_illumina:
	input:
		expand(rules.mer_gather_illumina_merged.output, sample=Illumina_process_df["sample"], 
			trimmer=trimmed_list_for_merging,
			corrector=correct_list_for_merging,
			merger=merge_list_for_merging)

rule ont_call:
	input:
		expand("results/{units.sample}/reads/ont/{basecaller}/{units.lib}/{units.sample}.{basecaller}.{unit}.fastq.gz", basecaller=config["ont_basecalling"]["basecaller"], units=fast5_units.itertuples(), unit=flappie_unit_list),
		
rule long_correct:
	input:
		expand("results/{units.sample}/errorcorrection/{longcor}/{basecaller}/{units.sample}.{basecaller}.{longcor}.fastq.gz", units=fast5_units.itertuples(), longcor=config["long_correction"], basecaller=config["ont_basecalling"]["basecaller"])


#rule long_correct_dev:
#	input:
#		expand("results/{units.sample}/errorcorrection/dev/{longcor}/{basecaller}/{units.sample}.{basecaller}.{longcor}.fastq.gz", units=fast5_units.itertuples(), longcor=config["long_correction"], basecaller=config["basecaller"])

rule assemble:
	input:
		expand(rules.eva_quast.output, sample=df["sample"],
			trimmer=trim_list,
			corrector=correct_list,
			merger=merge_list)
		
#illumina_units = units[units["sample"].isin(df_fastq["sample"].tolist()+df_bam["sample"].tolist())]
#rule assemble_short:
#	input:
#		expand(rules.quast.output, sample=Illumina_process_df["sample"])

#rule assemble_long:
#	input:
#		expand("results/{test.sample}/assembly/quast/quast.ok", test=fast5_units.itertuples(), basecaller=config["basecaller"])


s_with_ass = find_samples_with_assemblies(all_samples)
#print(s_with_ass)
rule quastall:
	input:
		expand(rules.eva_just_quast.output, sample=s_with_ass)

rule test_marvel:
	input:
#		expand("results/{test.sample}/assembly/marvel/raw/flappie/marvel.patching.ok", test=fast5_units.itertuples())
		expand("results/{test.sample}/assembly/marvel/raw/flappie/marvel.assembly.finale.ok", test=fast5_units.itertuples(), batch=marvel_batch_list)

hybrid_units = units[units["sample"].isin(df_fast5["sample"].tolist()+df_fastq["sample"].tolist()+df_bam["sample"].tolist())]
rule assemble_hybrid:
	input:
		expand("results/{test.sample}/assembly/quast/quast.ok", test=hybrid_units.itertuples())
		

#rule illumina_correct:
#	input:
#		#fastqc
#		expand("results/{sample}/read_qc/fastqc_raw/{lib}/{sample}.{lib}.status.ok", zip, sample=Illumina_process_df["sample"].tolist(), lib=Illumina_process_df["lib"].tolist()),
#		expand("results/{sample}/trimming/trim_galore/{lib}/{sample}.{lib}.fastqc.status.ok", zip, sample=Illumina_process_df["sample"].tolist(), lib=Illumina_process_df["lib"].tolist()),
#		#correct and gather
#		expand("results/{sample}/errorcorrection/bless/{sample}-full/{sample}.cat.status.ok", zip, sample=Illumina_process_df["sample"].tolist(), lib=Illumina_process_df["lib"].tolist()),

#rule illumina_merge:
#	input:
#		#fastqc
#		expand("results/{sample}/read_qc/fastqc_raw/{lib}/{sample}.{lib}.status.ok", zip, sample=df["sample"].tolist(), lib=df["lib"].tolist()),
#		expand("results/{sample}/trimming/trim_galore/{lib}/{sample}.{lib}.fastqc.status.ok", zip, sample=df["sample"].tolist(), lib=df["lib"].tolist()),
#		#merge and gather
#		expand("results/{sample}/readmerging/usearch/{sample}-full/{sample}.cat.status.ok", zip, sample=df["sample"].tolist(), lib=df["lib"].tolist()),

#rule illumina_assemble_tricormer:
#	input:

		
#unit=df_fastq_unit.itertuples()), 
#		expand("results/{unit.sample}/Illumina/raw_reads/{unit.lib}/{unit.sample}.{unit.lib}.raw.1.fastq.gz", unit=df_fastq_unit.itertuples()) 


#rule prepro_trim:
#	input:
#		#fastqc
#		expand("results/{unit.sample}/read_qc/fastqc_raw/{unit.lib}/{unit.sample}.{unit.lib}.status.ok", unit=units.itertuples()),
#		expand("results/{unit.sample}/trimming/trim_galore/{unit.lib}/{unit.sample}.{unit.lib}.fastqc.status.ok", unit=units.itertuples()),
#		#trimgalore
#		expand("results/{unit.sample}/trimming/trim_galore/{unit.lib}/{unit.sample}.{unit.lib}.status.ok", unit=units.itertuples()),
#
#rule prepro:
#	input:
#		#fastqc
#		expand("results/{unit.sample}/read_qc/fastqc_raw/{unit.lib}/{unit.sample}.{unit.lib}.status.ok", unit=units.itertuples()),
#		expand("results/{unit.sample}/trimming/trim_galore/{unit.lib}/{unit.sample}.{unit.lib}.fastqc.status.ok", unit=units.itertuples()),
##		expand("results/{unit.sample}/trimming/trim_galore/{unit.lib}/{unit.sample}.{unit.lib}.status.ok", unit=units.itertuples()),
#		#kmc
#		expand("results/{unit.sample}/kmc/{unit.sample}.k{k}.histogram.txt", unit=units.itertuples(), k=config["kmc"]["k"]),
#		#plots
#		expand("results/{unit.sample}/plots/{unit.sample}-k{k}-distribution-full.pdf" , unit=units.itertuples(), k=config["kmc"]["k"]),
##		expand("results/{unit.sample}/kmc/{unit.sample}.k"+str(config["kmc"]["k"])+".histogram.txt", unit=units.itertuples()),
##		expand("results/{unit.sample}/plots/{unit.sample}-k"+str(config["kmc"]["k"])+"-distribution-full.pdf" , unit=units.itertuples()),
##		expand("results/{unit.sample}/errorcorrection/bless/{unit.lib}.{pe}.fastq.gz", unit=units.itertuples(), pe=["1","2"]),
##		expand("results/{sample}/errorcorrection/bless/{sample}.bestk", sample=unitdict.keys()),
##		expand("results/{unit.sample}/errorcorrection/bless/{unit.lib}.{pe}.corrected.fastq.gz", unit=units.itertuples(), pe=["1","2"]),
#		#ec se
#		expand("results/{unit.sample}/errorcorrection/bless/{unit.lib}/{unit.sample}.{unit.lib}.se.corrected.fastq.gz", unit=units.itertuples()),
#		#read merging
#		expand("results/{unit.sample}/readmerging/usearch/{unit.lib}/{unit.sample}.{unit.lib}.merged.fastq.gz", unit=units.itertuples()),
#		expand("results/{unit.sample}/readmerging/usearch/{unit.lib}/{unit.sample}.{unit.lib}_{pe}.nm.fastq.gz", unit=units.itertuples(), pe=["1","2"]),
#
#rule prepro_trim_clean:
#	input:
#		#trimming
#		expand("results/{unit.sample}/trimming/trim_galore/{unit.sample}-full/{unit.sample}.cat.status.ok", unit=units.itertuples())
#
#rule prepro_kmers:
#	input:
#		#trimgalore
#		expand("results/{unit.sample}/trimming/trim_galore/{unit.lib}/{unit.sample}.{unit.lib}.status.ok", unit=units.itertuples()),
#		#kmc
#		expand("results/{unit.sample}/kmc/{unit.sample}.k{k}.histogram.txt", unit=units.itertuples(), k=config["kmc"]["k"]),
#		#plots
#		expand("results/{unit.sample}/plots/{unit.sample}-k{k}-distribution-full.pdf" , unit=units.itertuples(), k=config["kmc"]["k"]),
#		
#rule prepro_correct:
#	input:
#		#trimgalore
#		expand("results/{unit.sample}/trimming/trim_galore/{unit.lib}/{unit.sample}.{unit.lib}.status.ok", unit=units.itertuples()),
#		#ec se
#		expand("results/{unit.sample}/errorcorrection/bless/{unit.lib}/{unit.sample}.{unit.lib}.se.corrected.fastq.gz", unit=units.itertuples()),
#		expand("results/{unit.sample}/errorcorrection/bless/{unit.lib}/{unit.sample}.{unit.lib}.{pe}.corrected.fastq.gz", unit=units.itertuples(), pe=["1","2"]),
#		
#rule prepro_correct_clean:
#	input:
#		expand("results/{unit.sample}/errorcorrection/bless/{unit.lib}/{unit.sample}.{unit.lib}.se.corrected.fastq.gz", unit=units.itertuples()),
#		expand("results/{unit.sample}/errorcorrection/bless/{unit.lib}/{unit.sample}.{unit.lib}.{pe}.corrected.fastq.gz", unit=units.itertuples(), pe=["1","2"]),
#		#corrected
#		expand("results/{unit.sample}/errorcorrection/blesss/{unit.sample}-full/{unit.sample}.cat.status.ok", unit=units.itertuples())
#
#rule prepro_merge:
#	input:
#		#read merging
#		expand("results/{unit.sample}/readmerging/usearch/{unit.lib}/{unit.sample}.{unit.lib}.merged.fastq.gz", unit=units.itertuples()),
#		expand("results/{unit.sample}/readmerging/usearch/{unit.lib}/{unit.sample}.{unit.lib}_{pe}.nm.fastq.gz", unit=units.itertuples(), pe=["1","2"]),
#
#rule prepro_merge_clean:
#	input:
#		#merging
#		expand("results/{unit.sample}/readmerging/usearch/{unit.sample}-full/{unit.sample}.cat.status.ok", unit=units.itertuples())
#
#rule assemble:
#	input:
##		expand("results/{sample}/assembly/kmergenie/{sample}.bestk", sample=df["sample"]),
#		expand("results/{sample}/assembly/abyss/tri-bestk/abyss.ok", sample=df["sample"]),
#		expand("results/{sample}/assembly/abyss/tricormer-bestk/abyss.ok", sample=df["sample"]),
##		expand("results/{sample}/assembly/abyss/bestk/abyss.ok", sample=df["sample"]),
#		expand("results/{sample}/assembly/spades/default-tricormer/spades.ok", sample=df["sample"]),
#		expand("results/{sample}/assembly/spades/bestk-tricormer/spades.ok", sample=df["sample"]),
#		expand("results/{sample}/assembly/platanus/auto/platanus.ok", sample=df["sample"]),
#		expand("results/{sample}/assembly/minia/bestk/minia.done", sample=df["sample"]),
##		expand("results/{sample}/assembly/abyss/{ass_to_scf}-bestk/{sample}-unitigs.nr-{similarity}.fa", sample=df_long["sample"], similarity=config["ectools"]["unitig_similarity"], ass_to_scf=abyss_ass_to_scaffold),
#		expand("results/{sample}/errorcorrection/ectools/pass1/partition-{similarity}.ok", sample=df_fast5["sample"], similarity=config["ectools"]["unitig_similarity"]),
#		expand("results/{sample}/errorcorrection/ectools/pass1/ectools.{similarity}.{unit}.ok", sample=df_fast5["sample"], unit=ec_unit_list, similarity=config["ectools"]["unitig_similarity"]),
#		expand("results/{sample}/errorcorrection/ectools/pass2/partition-{similarity}.ok", sample=df_fast5["sample"], similarity=config["ectools"]["unitig_similarity"]),
#		expand("results/{sample}/errorcorrection/ectools/pass2/ectools.{similarity}.{unit}.ok", sample=df_fast5["sample"], unit=ec_unit_list, similarity=config["ectools"]["unitig_similarity"]),
#		expand("results/{sample}/assembly/abyss/ectools-{similarity}-rescaffold-{ass_to_scf}/abyss.ok", sample=df["sample"], similarity=config["ectools"]["unitig_similarity"], ass_to_scf=abyss_ass_to_scaffold),
#		expand("results/{sample}/assembly/abyss/long-rescaffold/abyss.ok", sample=df_long["sample"]),
##		expand("results/{sample}/raw_reads/ont/flappie.{sample}.flappie.{unit}.fastq.gz", sample=df_fast5["sample"], unit=flappie_unit_list),
#		expand("results/{sample}/reads/ont/bonito/{sample}.flappie.{unit}.fastq.gz", sample=df_fast5["sample"], unit=flappie_unit_list),
#		expand("results/{sample}/assembly/abyss/flappie-rescaffold/abyss.ok", sample=df_fast5["sample"], unit=flappie_unit_list),
#		expand("results/{sample}/assembly/abyss/consent-rescaffold/abyss.ok", sample=df_fast5["sample"], unit=flappie_unit_list),
###		expand("results/{sample}/assembly/abyss/ectools-rescaffold/abyss.ok", sample=df_fast5["sample"], unit=flappie_unit_list),
#		expand("results/{sample}/assembly/spades/flappie-hybrid/spades.ok", sample=df_fast5["sample"], unit=flappie_unit_list),
#		expand("results/{sample}/assembly/spades/hybrid-ectools-{similarity}/spades.ok", sample=df_fast5["sample"], unit=ec_unit_list, similarity=config["ectools"]["unitig_similarity"]),
#		expand("results/{sample}/errorcorrection/consent/{sample}.consent.fastq.gz", sample=df_fast5["sample"], unit=flappie_unit_list),
#		expand("results/{sample}/errorcorrection/ratatosk/{sample}.ratatosk.fastq.gz", sample=df_fast5["sample"], unit=flappie_unit_list),
#		expand("results/{sample}/assembly/haslr/flappie/haslr.ok", sample=df_fast5["sample"], unit=flappie_unit_list),
#
