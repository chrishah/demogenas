#!/bin/bash

sample=$1
shift

for run in "$@" 
do
	k=$(echo $run | awk '{n=split($0,A,"-"); print A[n]}' | sed 's/^k//' | sed 's/.done$//')
	kmers=$(cat results/$sample/logs/bless-k$k.log.txt | grep "Number of unique solid k-mers" | perl -ne 'chomp; @a=split(" "); print "$a[-1]\n"')
	echo "$k $kmers"
done | perl -ne 'chomp; @a=split(" "); $res=$a[1]/(4**$a[0]); if ($res < 0.0001){print "$a[0]\t$a[1]\t$res\n"}' | sort -n -k2
