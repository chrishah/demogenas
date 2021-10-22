#!/bin/bash

for log in "$@" 
do
	kmers=$(cat $log | grep "Number of unique solid k-mers" | perl -ne 'chomp; @a=split(" "); print "$a[-1]\n"')
	k=$(echo $log | awk '{n=split($0,A,"/"); print A[n]}' | sed 's/bless-k//' | sed 's/\..*//')
	echo "$k $kmers"
done | perl -ne 'chomp; @a=split(" "); $res=$a[1]/(4**$a[0]); if ($res < 0.0001){print "$a[0]\t$a[1]\t$res\n"}' | sort -n -k2
