## filter reads that are unique mapped and paired-end reads are mapped to different chromosomes
for i in `ls *bam`
do
    file=`basename $i .bam`
    samtools view -F 14 -q 5 -h  $i \
    | awk '$7 !~ /=/' \
    | samtools view -Sb \
    > filter/${file}_map_2_diff_chroms.bam
done


## make sliding windows for genome
ln -s ../../data/00ref/dmel-all-chromosome-r6.46.fasta ./
samtools faidx dmel-all-chromosome-r6.46.fasta
bedtools makewindows -g dmel-all-chromosome-r6.46.fasta.fai -w 100000 -s 10000 >w100000_s10000.bed


## count reads on sliding window and get paired-ends that are from 2L,2R,3L,3R,4,X,and Y chromosomes
for i in `ls *chroms.bam`
do
    file=`basename $i _map_2_diff_chroms.bam`
    bedtools bamtobed -i $i >$file.map_2_diff_chroms.bed
    awk '{if($1=="2L" || $1=="2R" || $1=="3L" || $1=="3R" || $1=="4" || $1=="X" || $1=="Y")print}' $file.map_2_diff_chroms.bed|cut -f4|sed 's#/.*$##g'|sort|uniq -c|awk '{if($1==2)print $2}' >$file.filter.reads.txt
    grep -f $file.filter.reads.txt $file.map_2_diff_chroms.bed > $file.filter.bed
    bedtools intersect -a ../w100000_s10000.bed -b $file.filter.bed -bed -c >$file.w100000_s10000.txt
    bedtools intersect -a ../w100000_s10000.bed -b $file.filter.bed -wo >$file.w100000_s10000.txt2
    awk 'BEGIN{OFS="\t"}{print $1":"$2"-"$3,$0}' $file.w100000_s10000.txt|sort -k1,1|join -t$'\t' -1 1 - -2 1 <(awk 'BEGIN{OFS="\t"}{print $1":"$2"-"$3,$0}' $file.w100000_s10000.txt2|sort -k1,1) >$file.w100000_s10000.join.txt
    awk 'BEGIN{OFS="\t"}{print $4,$1,$2,$3}' $file.filter.bed|sed 's#/1##g;s#/2##g' |sort -k1,1>$file.filter.txt
done

## running circos
for i in `ls data/*txt`
do
    file=`basename $i .txt`
    sed "s/test.txt/$file.txt/" test.conf > circos.conf
    circos -outputfile $file.png
done

awk '{if($2=="X" && $3<1.4e7 && $3>1.2e7 && $5>20)print}' BM-1.w100000_s10000.join.txt|sort -k5,5nr|cut -f12|sort -u|sed 's#/[1/2]##g'|grep -f - BM-1.filter.txt
awk '{if($2=="X" && $3<1.2e7 && $3>1.1e7 && $5>20)print}' BM-6.w100000_s10000.join.txt|sort -k5,5nr|cut -f12|sort -u|sed 's#/[1/2]##g'|grep -f - BM-6.filter.txt
awk '{if($2=="X" && $3<7e6 && $3>6.33e6 && $5>20)print}' BM-46.w100000_s10000.join.txt|sort -k5,5nr|cut -f12|sort -u|sed 's#/[1/2]##g'|grep -f - BM-46.filter.txt
