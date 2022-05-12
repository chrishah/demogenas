__default__:
   J: DEF
   ntasks: 1
   ntasks-per-core: 2
   hint: multithread
   qos: mem_0096
   partition: mem_0096
   mem: "4G"
   time: "00:10:00"
   output: $(pwd)/log/%x-%j.stdout.txt
   error: $(pwd)/log/%x-%j.stderr.txt
   mail-type: "BEGIN"
   mail-user: "dummy@gmail.at"
pre_prepare_fastq:
   J: PFQ
pre_sort_bam:
   J: SOB
   mem: "4G"
   time: "02:00:00"
pre_bam2fastq:
   J: B2F
   time: "02:00:00"
tri_trimgalore:
   J: TRI
   mem: "10G"
   time: "10:00:00"
tri_gather_short_trimmed_by_lib:
   J: GSL
   time: "05:30:00"
eva_fastqc_raw:
   J: FQR
   time: "05:00:00"
   mem: "6G"
eva_fastqc_trimmed:
   J: FQT
   time: "08:00:00"
   mem: "6G"
eva_stats:
   J: STA
   mem: "4G"
   time: "00:10:00"
eva_kmc:
   J: KMC
   time: "40:00:00"
eva_plot_k_hist:
   J: PKH
   mem: "2G"
   time: "00:01:00"

cor_bless_by_k:
   J: BBK
   time: "24:00:00"
   mem: 20G
cor_find_best_bless:
   J: FBK
   time: "00:00:10"
cor_bless_pe:
   J: BPE
   time: "24:00:00"
   mem: 20G
cor_bless_se:
   J: BSE
   time: "24:00:00"
   mem: 20G
cor_correct_spades:
   J: CSP
   N: 1
   qos: mem_0768
   partition: mem_0768
   time: "48:00:00"
cor_gather_illumina_corrected:
   J: GIC
   time: "00:00:10"

mer_setup_usearch:
   J: SUS
   mem: "2G"
   time: "00:10:00"
mer_usearch:
   J: MUS
   mem: "20G"
   time: "20:00:00"
mer_flash:
   J: FSH
   mem: "8G"
   time: "05:00:00"
mer_gather_illumina_merged:
   J: GIM
   mem: "4G"
   time: "00:00:10"

ail_kmergenie:
   J: KMG
   time: "30:00:00"
   mem: 20G
ail_abyss:
   J: ABY
   N: 1
   qos: mem_0768
   partition: mem_0768
   ntasks: "threads"
   time: "72:00:00"
   dependency: "Singleton"
ail_minia:
   J: MIN
   time: "30:00:00"
   mem: "20G"
ail_platanus:
   J: PLA
   N: 1
   qos: mem_0768
   partition: mem_0768
   time: "72:00:00"
   dependency: "Singleton"
ail_spades:
   J: SPA
   N: 1
   qos: mem_0768
   partition: mem_0768
   time: "72:00:00"
   dependency: "Singleton"
ail_gather_illumina_assemblies:
   J: GIA
   mem: "4G"
   time: "00:00:10"
ass_gather_assemblies:
   J: GAS
   mem: "4G"
   time: "00:00:10"
eva_prepare_assemblies:
   J: PRA
   mem: "4G"
   time: "00:01:00"
eva_quast:
   J: QUA
   time: "01:00:00"
eva_busco:
   J: BUS
   mem: "10G"
   time: "04:00:00"
eva_gather_busco:
   J: GBU
   time: "00:00:10"
eva_bowtie2_index:
   J: BOI
   time: "00:30:00"
   mem: "16G"
eva_bowtie2_mapping:
   J: BOM
   time: "03:00:00"
   mem: "8G"
eva_conversion_samtools:
   J: CST
   time: "02:00:00"
   mem: "8G"
eva_remove_duplicates:
   J: RMD
   time: "05:00:00"
   mem: "32G"
eva_index_duprmvd:
   J: SAI
   time: "03:00:00"
   mem: "8G"
eva_blobtools:
   J: BLO
   time: "03:00:00"
   mem: "8G"
eva_x_gather_evaluations:
   J: GEV
   time: "00:00:10"

alo_flye_raw:
   J: FLR
   N: 1
   qos: mem_0768
   partition: mem_0768
   time: "48:00:00"
alo_flye_corrected:
   J: FLC
   N: 1
   qos: mem_0768
   partition: mem_0768
   time: "48:00:00"
alo_canu:
   J: CAR
   N: 1
   qos: mem_0768
   partition: mem_0768
   time: "72:00:00"
alo_canu_corrected:
   J: CAC
   N: 1
   qos: mem_0768
   partition: mem_0768
   time: "72:00:00"
ahy_abyss_scaffold:
   J: ABS
   N: 1
   ntasks: "threads"
   qos: mem_0384
   partition: mem_0384
   time: "24:00:00"
ahy_spades_hybrid:
   J: SPH
   N: 1
   qos: mem_0768
   partition: mem_0768
   time: "48:00:00"
ahy_haslr:
   J: HSL
   N: 1
   qos: mem_0768
   partition: mem_0768
   time: "48:00:00"
bca_guppy:
   J: GUP
   time: "70:00:00"
bca_bonito_cpu:
   J: BOC
   time: "05:00:00"
bca_flappie:
   J: FLA
   time: "20:00:00"
bca_x_gather_called_ont_reads:
   J: GCO
   time: "00:00:05"
cor_ratatosk:
   J: RTT
   N: 1
   qos: mem_0768
   partition: mem_0768
   time: "20:00:00"
cor_gather_corrected_by_lib:
   J: GCL
   time: "00:00:10"
cor_consent:
   J: CON
   N: 1
   qos: mem_0768
   partition: mem_0768
   time: "70:00:00"
cor_correct_canu:
   N: CCA
   N: 1
   qos: mem_0768
   partition: mem_0768
   time: "70:00:00"

#old and future

eval_illumina:
   J: EVI
   time: "00:00:10"
assemble:
   J: assALL
assemble_short:
   J: assSHO
assemble_long:
   J: assLON
assemble_hybrid:
   J: assHYB
ont_call:
   J: F5C
abyss_tri:
   J: ABT
   N: 1
   qos: mem_0768
   partition: mem_0768
   time: "48:00:00"
abyss_rescaffold_long:
   J: ABL
   time: "24:00:00"
spades_hybrid_corrected:
   J: SPC
   N: 1
   qos: mem_0768
   partition: mem_0768
   time: "48:00:00"
cdhitest_unitigs:
   J: CDU
   time: "10:00:00"
   mem: "8G"
split:
   J: SPL
   time: "01:00:00"
ectools:
   J: EC1
   time: "48:00:00"
split_pass2:
   J: SP2
   time: "00:10:00"
ectools_pass2:
   J: EC2
   time: "10:00:00"
gather_ectools_corrected:
   J: ECG
   time: "05:00:00"
marvel:
   J: MAR
   time: "20:00:00"
marvel_pat_init:
   J: MPI
   time: "02:00:00"
marvel_pat_daligner_quality_trim:
   J: MPD
   time: "20:00:00"
marvel_pat_merge:
   J: MPM
   time: "02:00:00"
marvel_pat_patch:
   J: MPP
   time: "02:00:00"
marvel_ass_init:
   J: MAI
   time: "02:00:00"
marvel_ass_dali_rep_trim:
   J: MAD
   time: "10:00:00"
marvel_ass_merge_qual_trim:
   J: MAM
   time: "02:00:00"
marvel_ass_rep_annot:
   J: MRA
   time: "01:00:00"
marvel_ass_merge_repeats:
   J: MMR
   time: "01:00:00"
marvel_ass_hom_rep:
   J: MHR
   time: "01:00:00"
marvel_ass_combine_rem_dup_over:
   J: MCR
   time: "01:00:00"
marvel_ass_gap_quality:
   J: MGQ
   time: "01:00:00"
marvel_ass_merge_trim1:
   J: MMT
   time: "01:00:00"
marvel_ass_filter_rep_ali:
   J: MFR
   time: "01:00:00"
marvel_ass_merge_filtered:
   J: MMF
   time: "01:00:00"
marvel_ass_tour_convert:
   J: MTC
   time: "05:00:00"
marvel_ass_finale:
   J: MFI
   time: "00:10:00"
prepro:
   J: preALL
all_trim:
   J: preALT
all_trim_clean:
   J: preATC
all_kmers:
   J: preALK
all_correct:
   J: preALC
all_correct_clean:
   J: preACC
all_merge:
   J: preALM
all_merge_clean:
   J: preAMC
illumina_trim:
   J: ITR
   time: "00:01:00"
long_correct:
   J: LOC
reformat_read_headers:
   J: REF
   time: "24:00:00"
correct_illumina:
   J: CIL
   time: "00:00:10"
clean_merged_libs:
   J: CML
   time: "00:30:00"
