#!/bin/bash

##define prefix
echo -ne "" > counts
for prefix in $(echo $1 | sed 's/,/\n/g')
do
	base=$(echo $prefix | awk -F "/" '{print $NF}')
	#forward
	startend=$(unzip -p $prefix.1_fastqc.zip $base.1_fastqc/fastqc_data.txt | grep -n -e "#" -e ">>" | grep -P "Length\t" -A 1 | sed 's/:.*//' | tr '\n' ',' | sed 's/,$/\n/')
	start=$(echo "$startend" | cut -d "," -f 1)
	start=$(( start + 1 ))
	end=$(echo "$startend" | cut -d "," -f 2)
	end=$(( end - 1 ))
	unzip -p $prefix.1_fastqc.zip $base.1_fastqc/fastqc_data.txt | head -n $end | tail -n +$start | perl -ne 'chomp; @a=split("\t"); $readlength=$a[0]; $readlength=~s/-.*//; $readcount=$readcount+$a[1]; $cum_length=$cum_length+($readlength*$a[1]); if (eof()){print "$cum_length,$readcount\n"}' >> counts
	#reverse
	startend=$(unzip -p $prefix.2_fastqc.zip $base.2_fastqc/fastqc_data.txt | grep -n -e "#" -e ">>" | grep -P "Length\t" -A 1 | sed 's/:.*//' | tr '\n' ',' | sed 's/,$/\n/')
	start=$(echo "$startend" | cut -d "," -f 1)
	start=$(( start + 1 ))
	end=$(echo "$startend" | cut -d "," -f 2)
	end=$(( end - 1 ))
	unzip -p $prefix.2_fastqc.zip $base.2_fastqc/fastqc_data.txt | head -n $end | tail -n +$start | perl -ne 'chomp; @a=split("\t"); $readlength=$a[0]; $readlength=~s/-.*//; $readcount=$readcount+$a[1]; $cum_length=$cum_length+($readlength*$a[1]); if (eof()){print "$cum_length,$readcount\n"}' >> counts
done

#sum and average
cat counts | perl -ne 'chomp; ($length,$readcount) = split(","); $slen = $slen + $length; $src = $src + $readcount; if (eof()){print "$slen ".($slen/$src)."\n"}'
