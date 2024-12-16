for i in `ls *.bam`
do
samtools index  ./${i} ./${i}.bai
done


for i in d5_1 d5_2 d5_3 d5_4 d5_5 d5_6 d5_7 d5_8
do
velocyto run -m /home/liyaru/public_Data/mm10_rmsk.gtf -@ 10 ${i}_possorted_genome_bam.bam /home/liyaru/public_Data/refdata-gex-mm10-2020-A/genes/genes.gtf
done


