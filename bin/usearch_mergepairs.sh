#!/bin/bash

####
# This script merges paired end Illumina reads using the -fastq_mergepairs option from usearch (32-bit version). 
# Since the usearch 32-bit version has a memory limit the script breaks the fastq data into batches of 4,000,000 read pairs (worked for me) 
# and then merges the read pairs in each batch consecutively. Finally all files are combined to a final output gzipped fastq file.
# 
# Notes on requirements:
# The script requires:
# - a recent usearch (32-bit) version (tested with usearch 11 upwards)
# - GNU parallel
# 
# The script expects the fastq data to be gzipped (*.gz)
####

if [ $# -lt 3 ]
then
	echo -e "The script takes a minimum of three arguments (in this order): "
	echo -e "1) forward fastq file (gzipped)"
	echo -e "2) reverse fastq file (gzipped)"
	echo -e "3) prefix for output files"
	echo -e "optional:"
	echo -e "4) specify the number of processes to run in parallel (default=1)"
	echo -e "5) specify the number of read pairs per batch (default=4000000)"
	echo -e "Example: ./usearch_mergepairs.sh forward.fastq.gz reverse.fastq.gz prefix 5 4000000"
	exit 1
fi

if [ -z "$4" ]
then
	p=1
else
	p=$4
fi

if [ -z "$5" ]
then
	batchsize=4000000
else
	batchsize=$5
fi

usearch=usearch
prefix=$3

echo -e "\n[$(date)]\tbatchsize: $batchsize\tparallel: $p"

#count the number of reads in the forward read file
total=$(( $(zcat $1 | wc -l) / 4 ))

#calculate how many batches to do
batches=$(echo -e "$total\t$batchsize" | perl -ne 'chomp; @a=split("\t"); $ret=sprintf("%0.f",(($a[0]/$a[1])+0.5)); print "$ret\n"')

echo -e "\n[$(date)]\tDetected $total read pairs - Merging in $batches batches of $batchsize reads"

#function that splits up fastq files into batches
per_batch () {
        echo -en "[$(date)]\t\tPreparing batch\t$3 ...\t"

        paste <(zcat $1) <(zcat $2) <(echo "$3") <(echo "$4") | perl -ne 'chomp; @a=split("\t"); if ($. == 1){$batch=pop(@a); $factor=pop(@a);}; if ($. >= ($factor * ($batch*4) - ($batch*4) + 1)){ print "$a[0]\n"; print STDERR "$a[1]\n"}; if ($. >= ($factor * ($batch*4))){exit;}' 1> sub_$3\_1.fastq 2> sub_$3\_2.fastq

        echo -e "\n#####\nBatch $3\n#####\n" > merge_$3.log

        echo -e "[$(date)]\tDone!"
}
export -f per_batch

#using GNU parallel
echo -e "\n[$(date)]\tPreparing batches .."
parallel -j $p per_batch $1 $2 {1} $batchsize ::: $(echo $(seq $batches))

#merging reads in all batches
echo -e "\n[$(date)]\tMerging read pairs .."
for a in $(seq 1 $batches)
do
	echo -en "[$(date)]\t\tMerging reads from Batch\t$a ...\t"

	$usearch -threads $p -fastq_mergepairs sub_$a\_1.fastq -reverse sub_$a\_2.fastq -fastq_trimstagger -fastqout_notmerged_fwd sub_$a\_1.nm.fastq -fastqout_notmerged_rev sub_$a\_2.nm.fastq -fastqout sub_$a.merged.fastq &>> merge_$a.log

	echo -e "[$(date)]\tDone!"
done

echo -e "\n[$(date)]\tConcatenating and compressing .."
cat *merged.fastq | perl -ne 'if ((($.-1)%4) == 0){chomp; print substr($_,0,-2)."\n"}else{print $_}' | gzip > $prefix.merged.fastq.gz

for f in $(ls -1 *_1.nm.fastq | grep "^sub"); do for=$f; cat $for; done | gzip > $prefix\_1.nm.fastq.gz
#cat *1.nm.fastq | gzip > $prefix\_1.nm.fastq.gz

for f in $(ls -1 *_1.nm.fastq | grep "^sub"); do for=$f; rev=$(echo -e "$f" | sed 's/_1.nm.fastq/_2.nm.fastq/'); cat $rev; done | gzip > $prefix\_2.nm.fastq.gz
#cat *2.nm.fastq | gzip > $prefix\_2.nm.fastq.gz

cat merge_* > merging.log

#cleanup
echo -e "\n[$(date)]\tRemoving temporary files .."
rm *.merged.fastq
rm *1.nm.fastq
rm *2.nm.fastq
rm sub_*
rm merge_*

echo -e "\n[$(date)]\tAll Done!\n"
