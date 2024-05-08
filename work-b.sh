#swarm 1:
mkdir BM-1; configManta.py --bam /data/wangy80/TK117/results/03map_reads2/BM-1.bam --referenceFasta /data/wangy80/TK117/data/00ref/dmel-all-chromosome-r6.46.fasta --runDir BM-1
mkdir BM-21; configManta.py --bam /data/wangy80/TK117/results/03map_reads2/BM-21.bam --referenceFasta /data/wangy80/TK117/data/00ref/dmel-all-chromosome-r6.46.fasta --runDir BM-21
mkdir BM-24; configManta.py --bam /data/wangy80/TK117/results/03map_reads2/BM-24.bam --referenceFasta /data/wangy80/TK117/data/00ref/dmel-all-chromosome-r6.46.fasta --runDir BM-24
mkdir BM-3; configManta.py --bam /data/wangy80/TK117/results/03map_reads2/BM-3.bam --referenceFasta /data/wangy80/TK117/data/00ref/dmel-all-chromosome-r6.46.fasta --runDir BM-3
mkdir BM-46; configManta.py --bam /data/wangy80/TK117/results/03map_reads2/BM-46.bam --referenceFasta /data/wangy80/TK117/data/00ref/dmel-all-chromosome-r6.46.fasta --runDir BM-46
mkdir BM-4; configManta.py --bam /data/wangy80/TK117/results/03map_reads2/BM-4.bam --referenceFasta /data/wangy80/TK117/data/00ref/dmel-all-chromosome-r6.46.fasta --runDir BM-4
mkdir BM-5; configManta.py --bam /data/wangy80/TK117/results/03map_reads2/BM-5.bam --referenceFasta /data/wangy80/TK117/data/00ref/dmel-all-chromosome-r6.46.fasta --runDir BM-5
mkdir BM-6; configManta.py --bam /data/wangy80/TK117/results/03map_reads2/BM-6.bam --referenceFasta /data/wangy80/TK117/data/00ref/dmel-all-chromosome-r6.46.fasta --runDir BM-6
mkdir BM-7; configManta.py --bam /data/wangy80/TK117/results/03map_reads2/BM-7.bam --referenceFasta /data/wangy80/TK117/data/00ref/dmel-all-chromosome-r6.46.fasta --runDir BM-7
mkdir BM-S; configManta.py --bam /data/wangy80/TK117/results/03map_reads2/BM-S.bam --referenceFasta /data/wangy80/TK117/data/00ref/dmel-all-chromosome-r6.46.fasta --runDir BM-S

#swarm 2:
BM-1/runWorkflow.py -m local -j $SLURM_CPUS_PER_TASK -g $((SLURM_MEM_PER_NODE / 1024))
BM-21/runWorkflow.py -m local -j $SLURM_CPUS_PER_TASK -g $((SLURM_MEM_PER_NODE / 1024))
BM-24/runWorkflow.py -m local -j $SLURM_CPUS_PER_TASK -g $((SLURM_MEM_PER_NODE / 1024))
BM-3/runWorkflow.py -m local -j $SLURM_CPUS_PER_TASK -g $((SLURM_MEM_PER_NODE / 1024))
BM-46/runWorkflow.py -m local -j $SLURM_CPUS_PER_TASK -g $((SLURM_MEM_PER_NODE / 1024))
BM-4/runWorkflow.py -m local -j $SLURM_CPUS_PER_TASK -g $((SLURM_MEM_PER_NODE / 1024))
BM-5/runWorkflow.py -m local -j $SLURM_CPUS_PER_TASK -g $((SLURM_MEM_PER_NODE / 1024))
BM-6/runWorkflow.py -m local -j $SLURM_CPUS_PER_TASK -g $((SLURM_MEM_PER_NODE / 1024))
BM-7/runWorkflow.py -m local -j $SLURM_CPUS_PER_TASK -g $((SLURM_MEM_PER_NODE / 1024))
BM-S/runWorkflow.py -m local -j $SLURM_CPUS_PER_TASK -g $((SLURM_MEM_PER_NODE / 1024))

## get translocation sites
for i in `ls /data/wangy80/TK117/manta/*/results/variants/candidateSV.vcf.gz`
do
	file=`echo $i|sed 's#^.*manta/##;s#/results/.*$##'`
	zcat $i|grep -v "#"|egrep -v "2110000|mapped"|grep "MantaBND"|cut -f1-3|sed 's/:0$//g;s/:1$//g'|awk 'BEGIN{OFS="\t"}{print $3,$1,$2,$2+1}'|sort -k1,1 > data/${file}.BND.txt
done

## circos ploting
for i in `ls data/*txt`
do     
	file=`basename $i .txt`;     
	sed "s/test.txt/$file.txt/" test.conf > circos.conf;     
	circos -outputfile $file.png; 
done