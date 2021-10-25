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
eval_illumina:
   J: EVI
   time: "00:00:10"
quast:
   J: QUA
   time: "01:00:00"
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
kmergenie:
   J: KMG
   time: "03:00:00"
   mem: 20G
abyss:
   J: ABY
   N: 1
   ntasks: "threads"
   qos: mem_0384
   partition: mem_0384
   time: "42:00:00"
abyss_tri:
   J: ABT
   N: 1
   qos: mem_0768
   partition: mem_0768
   time: "48:00:00"
abyss_rescaffold_long:
   J: ABL
   time: "24:00:00"
spades:
   J: SPA
   N: 1
   qos: mem_0768
   partition: mem_0768
   time: "48:00:00"
spades_hybrid:
   J: SPR
   N: 1
   qos: mem_0768
   partition: mem_0768
   time: "48:00:00"
spades_hybrid_corrected:
   J: SPC
   N: 1
   qos: mem_0768
   partition: mem_0768
   time: "48:00:00"
minia:
   J: MIN
   time: "05:00:00"
   mem: "40G"
platanus:
   J: PLA
   N: 1
   qos: mem_0384
   partition: mem_0384
   time: "24:00:00"
cdhitest_unitigs:
   J: CDU
   time: "10:00:00"
   mem: "8G"
flappie:
   J: FLA
   time: "20:00:00"
bonito_cpu:
   J: BOC
   time: "05:00:00"
guppy:
   J: GUP
   time: "15:00:00"
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
ratatosk:
   J: RTT
   N: 1
   qos: mem_0768
   partition: mem_0768
   time: "70:00:00"
gather_corrected_by_lib:
   J: GCL
   time: "01:00:00"
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
consent:
   J: CON
   time: "20:00:00"
abyss_rescaffold:
   J: ABR
   N: 1
   qos: mem_0768
   partition: mem_0768
   time: "10:00:00"
abyss_rescaffold_corrected:
   J: ABC
   N: 1
   qos: mem_0768
   partition: mem_0768
   time: "10:00:00"
spades_flappie:
   J: SPF
   N: 1
   qos: mem_0768
   partition: mem_0768
   time: "48:00:00"
spades_consent:
   J: SPC
   N: 1
   qos: mem_0768
   partition: mem_0768
   time: "48:00:00"
flye_raw:
   J: FLR
   N: 1
   qos: mem_0768
   partition: mem_0768
   time: "48:00:00"
flye_corrected:
   J: FLC
   N: 1
   qos: mem_0768
   partition: mem_0768
   time: "48:00:00"
canu:
   J: CAR
   N: 1
   qos: mem_0768
   partition: mem_0768
   time: "72:00:00"
canu_correct:
   J: CAC
   N: 1
   qos: mem_0768
   partition: mem_0768
   time: "72:00:00"
haslr:
   J: HSL
   N: 1
   qos: mem_0768
   partition: mem_0768
   time: "48:00:00"
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
prepare_fastq:
   J: PFQ
sort_bam:
   J: SOB
   mem: "4G"
   time: "02:00:00"
bam2fastq:
   J: B2F
   time: "02:00:00"
fastqc_raw:
   J: FQR
   time: "02:00:00"
   mem: "6G"
fastqc_trimmed:
   J: FQT
   time: "02:00:00"
   mem: "6G"
reformat_read_headers:
   J: REF
   time: "24:00:00"
trim_trimgalore:
   J: TRI
   mem: "10G"
   time: "10:00:00"
gather_short_trimmed_by_lib:
   J: GSL
   time: "00:00:10"
kmc:
   J: KMC
   time: "01:00:00"
bless_by_k:
   J: BBK
   time: "10:00:00"
   mem: 20G
bless_pe:
   J: BPE
   time: "24:00:00"
   mem: 20G
bless_se:
   J: BSE
   time: "24:00:00"
   mem: 20G
find_best_bless:
   J: FBK
   time: "00:00:10"
correct_spades:
   J: CSP
   N: 1
   qos: mem_0384
   partition: mem_0384
   time: "05:00:00"
gather_illumina_corrected:
   J: GIC
   time: "00:00:10"
correct_illumina:
   J: CIL
   time: "00:00:10"
setup_usearch:
   J: SUS
   mem: "2G"
   time: "00:10:00"
usearch:
   J: MUS
   mem: "20G"
   time: "05:00:00"
gather_illumina_merged:
   J: GIM
   mem: "4G"
   time: "00:00:10"
clean_merged_libs:
   J: CML
   time: "00:30:00"
stats:
   J: STA
   mem: "4G"
   time: "00:10:00"
plot_k_hist:
   J: PKH
   mem: "2G"
   time: "00:01:00"
