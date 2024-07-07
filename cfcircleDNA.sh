#!/bin/bash

input_dir=$1
output_dir=$2
hg_index=$3
files=$(ls "$input_dir"/*.paired.R1.fq.gz)


for filename in $files
do
name=`basename $filename .paired.R1.fq.gz`
read1=$input_dir/$name.paired.R1.fq.gz
read2=$input_dir/$name.paired.R2.fq.gz

mkdir temp_d
sam_file=temp_d/$name.sam

   bwa mem -M -t 16 $hg_index $read1 $read2 > $sam_file

    samtools view -@ 16 -bS "$sam_file" -o "${sam_file%.sam}.bam"

    samtools sort -@ 16 "${sam_file%.sam}.bam" -o "${sam_file%.sam}.sorted.bam"

    samtools index -@ 16 "${sam_file%.sam}.sorted.bam"

    samtools view -@ 16 -hf 0x100 "${sam_file%.sam}.sorted.bam" -bS > "${sam_file%.sam}_sec.bam"

    bedtools"${sam_file%.sam}_sec.bam" | sed -e s/\\//\ /g | awk '{printf ("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8)}' > "${sam_file%.sam}_sec.txt"

    awk '{print $4}' "${sam_file%.sam}_sec.txt" | sort | uniq > "${sam_file%.sam}_sec_uid.txt"

    samtools view -h  -F 0X100 -N "${sam_file%.sam}_sec_uid.txt" "${sam_file%.sam}.sorted.bam" > "${sam_file%.sam}_primary.bam"


    samtools view -h "${sam_file%.sam}_primary.bam" | awk -F'\t' 'BEGIN {OFS=FS} $0 ~ /^@/ {print} $0 !~ /^@/ {for (i=1; i<=NF; i++) {if ($i ~ /^SA:Z:/) {split($i, sa, ","); if (sa[5] >= 0) {print}; break}}}' | samtools view -bS - > "${sam_file%.sam}_pri.bam"
    samtools sort -@ 16 "${sam_file%.sam}_pri.bam" -o "${sam_file%.sam}_pri.sorted.bam"
    samtools index -@ 16 "${sam_file%.sam}_pri.sorted.bam"
    samtools sort -@ 16 "${sam_file%.sam}_sec.bam" -o "${sam_file%.sam}_sec.sorted.bam"
    samtools index -@ 16 "${sam_file%.sam}_sec.sorted.bam"
    samtools merge  "${sam_file%.sam}_split.bam" "${sam_file%.sam}_pri.sorted.bam" "${sam_file%.sam}_sec.sorted.bam" -f

    bedtools bamtobed -cigar -i "${sam_file%.sam}_split.bam" | sed -e s/\\//\ /g | awk '{printf ("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8)}' | awk '{print $4}' | sort | uniq -c | awk '$1=="2" {print $2}' > "${sam_file%.sam}_split.freq2_id.txt"

    samtools view -h -N "${sam_file%.sam}_split.freq2_id.txt" "${sam_file%.sam}_split.bam" > "${sam_file%.sam}_split.freq2.bam"

    bedtools bamtobed -cigar -i "${sam_file%.sam}_split.freq2.bam" | sed -e s/\\//\ /g | awk '{printf ("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8)}'| sort -k4 | sed 'N;s/\n/\t/' | awk '$6>19 || $14>19 {print $4}' > "${sam_file%.sam}_split.freq2.id.filter.txt"

    samtools view -h -N "${sam_file%.sam}_split.freq2.id.filter.txt" "${sam_file%.sam}_split.freq2.bam" > "${sam_file%.sam}_split.freq2.quality.bam"
    samtools sort -@ 16 "${sam_file%.sam}_split.freq2.quality.bam" -o "${sam_file%.sam}_split.freq2.quality.sorted.bam"  
    mv "${sam_file%.sam}_split.freq2.quality.sorted.bam" $output_dir
    mv "${sam_file%.sam}.sorted.bam" $output_dir
    mv "${sam_file%.sam}.sorted.bam.bai" $output_dir
    rm "$sam_file"
done

