## swarm:
#spades assembly
mkdir BM-1;cd BM-1;spades.py -t $SLURM_CPUS_PER_TASK -o BM-1 --pe1-1 /data/wangy80/TK117/data/00reads/BM-1.R1.fastq.gz --pe1-2 /data/wangy80/TK117/data/00reads/BM-1.R2.fastq.gz
mkdir BM-21;cd BM-21;spades.py -t $SLURM_CPUS_PER_TASK -o BM-21 --pe1-1 /data/wangy80/TK117/data/00reads/BM-21.R1.fastq.gz --pe1-2 /data/wangy80/TK117/data/00reads/BM-21.R2.fastq.gz
mkdir BM-24;cd BM-24;spades.py -t $SLURM_CPUS_PER_TASK -o BM-24 --pe1-1 /data/wangy80/TK117/data/00reads/BM-24.R1.fastq.gz --pe1-2 /data/wangy80/TK117/data/00reads/BM-24.R2.fastq.gz
mkdir BM-3;cd BM-3;spades.py -t $SLURM_CPUS_PER_TASK -o BM-3 --pe1-1 /data/wangy80/TK117/data/00reads/BM-3.R1.fastq.gz --pe1-2 /data/wangy80/TK117/data/00reads/BM-3.R2.fastq.gz
mkdir BM-46;cd BM-46;spades.py -t $SLURM_CPUS_PER_TASK -o BM-46 --pe1-1 /data/wangy80/TK117/data/00reads/BM-46.R1.fastq.gz --pe1-2 /data/wangy80/TK117/data/00reads/BM-46.R2.fastq.gz
mkdir BM-4;cd BM-4;spades.py -t $SLURM_CPUS_PER_TASK -o BM-4 --pe1-1 /data/wangy80/TK117/data/00reads/BM-4.R1.fastq.gz --pe1-2 /data/wangy80/TK117/data/00reads/BM-4.R2.fastq.gz
mkdir BM-5;cd BM-5;spades.py -t $SLURM_CPUS_PER_TASK -o BM-5 --pe1-1 /data/wangy80/TK117/data/00reads/BM-5.R1.fastq.gz --pe1-2 /data/wangy80/TK117/data/00reads/BM-5.R2.fastq.gz
mkdir BM-6;cd BM-6;spades.py -t $SLURM_CPUS_PER_TASK -o BM-6 --pe1-1 /data/wangy80/TK117/data/00reads/BM-6.R1.fastq.gz --pe1-2 /data/wangy80/TK117/data/00reads/BM-6.R2.fastq.gz
mkdir BM-7;cd BM-7;spades.py -t $SLURM_CPUS_PER_TASK -o BM-7 --pe1-1 /data/wangy80/TK117/data/00reads/BM-7.R1.fastq.gz --pe1-2 /data/wangy80/TK117/data/00reads/BM-7.R2.fastq.gz
mkdir BM-S;cd BM-S;spades.py -t $SLURM_CPUS_PER_TASK -o BM-S --pe1-1 /data/wangy80/TK117/data/00reads/BM-S.R1.fastq.gz --pe1-2 /data/wangy80/TK117/data/00reads/BM-S.R2.fastq.gz


## get genome sequences from 2L,2R,3L,3R,4,X,and Y chromosomes
ln -s ../../data/00ref/Drosophila_melanogaster.BDGP6.46.dna_rm.toplevel.fa ./
sed 's/ .*$//g' Drosophila_melanogaster.BDGP6.46.dna_rm.toplevel.fa >dmel.rm.fa
samtools faidx dmel.rm.fa

cat chrs.txt|while read line
do 
	samtools faidx /gpfs/gsfs12/users/wangy80/TK117/spades/nucmer/dmel.rm.fa $line
done>dmel.7chrs.rm.fa

for i in 2L 2R 3L 3R 4 X Y
do
	samtools faidx /gpfs/gsfs12/users/wangy80/TK117/spades/nucmer/dmel.rm.fa $i >dmel.$i.rm.fa
done


## nucmer alignment
##swarm:
nucmer dmel.7chrs.rm.fa --mum -p BM-1 /gpfs/gsfs12/users/wangy80/TK117/spades/BM-1/BM-1/scaffolds.fasta
nucmer dmel.7chrs.rm.fa --mum -p BM-3 /gpfs/gsfs12/users/wangy80/TK117/spades/BM-3/BM-3/scaffolds.fasta
nucmer dmel.7chrs.rm.fa --mum -p BM-4 /gpfs/gsfs12/users/wangy80/TK117/spades/BM-4/BM-4/scaffolds.fasta
nucmer dmel.7chrs.rm.fa --mum -p BM-5 /gpfs/gsfs12/users/wangy80/TK117/spades/BM-5/BM-5/scaffolds.fasta
nucmer dmel.7chrs.rm.fa --mum -p BM-6 /gpfs/gsfs12/users/wangy80/TK117/spades/BM-6/BM-6/scaffolds.fasta
nucmer dmel.7chrs.rm.fa --mum -p BM-7 /gpfs/gsfs12/users/wangy80/TK117/spades/BM-7/BM-7/scaffolds.fasta
nucmer dmel.7chrs.rm.fa --mum -p BM-21 /gpfs/gsfs12/users/wangy80/TK117/spades/BM-21/BM-21/scaffolds.fasta
nucmer dmel.7chrs.rm.fa --mum -p BM-24 /gpfs/gsfs12/users/wangy80/TK117/spades/BM-24/BM-24/scaffolds.fasta
nucmer dmel.7chrs.rm.fa --mum -p BM-46 /gpfs/gsfs12/users/wangy80/TK117/spades/BM-46/BM-46/scaffolds.fasta
nucmer dmel.7chrs.rm.fa --mum -p BM-S /gpfs/gsfs12/users/wangy80/TK117/spades/BM-S/BM-S/scaffolds.fasta

## only keep 1-to-1 alignments
for i in `ls *delta`
do
	file=`basename $i .delta`
	delta-filter -1 $i >$file.filter.delta
done


## mummerplot
for i in `ls *delta`
do
	file=`basename $i .delta`
	mummerplot -l $i --large -prefix ${file} -t png --color
done

## identify translocation sites
for i in `ls *delta`
do
	file=`basename $i .delta`
	show-diff $i -q >$file.diff2
	show-diff $i -r >$file.diff
	show-coords $i -d -T > $file.coords
done

for i in `ls *coords`
do
	file=`basename $i .coords`
	awk 'BEGIN{OFS="\t"}{if($9=="-1"){$3=$3+1;$4=$4-1}else if($9=="1"){$3=$3-1;$4=$4+1}print}' $i|egrep -v "mapped|2110000|mitochondrion|rDNA"|sort -k11,11|join -t$'\t' -1 11 - -2 1 <(grep "SEQ" $file.diff2|sort -k1,1)|awk '{if(($11==$15 && $4==$13) || ($11==$15 && $4==$14) || ($11==$15 && $5==$14) || ($11==$15 && $5==$13) || ($11==$16 && $4==$13) || ($11==$16 && $4==$14) || ($11==$16 && $5==$13) || ($11==$16 && $5==$13))print $1":"$13"-"$14,$0}' >$file.translocation.txt
done


## circos ploting
for i in `ls *translocation.txt`
do
file=`basename $i .translocation.txt`
egrep -v "mapped|2110000|mitochondrion|rDNA" $i|awk 'BEGIN{OFS="\t"}{print $17,$3,$4,$18,$5,$6}'|grep "X" >circos/data/${file}.txt
done

cd circos
for i in `ls data/*txt`
do     
	file=`basename $i .txt`;     
	sed "s/test.txt/$file.txt/" test.conf > circos.conf;     
	circos -outputfile $file.png; 
done
