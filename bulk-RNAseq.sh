/public-supool/home/zhangzh/soft_share/fastp -i $R1 -I $R2 -o ${name}_1P.fq.gz -O ${name}_2P.fq.gz -w 16 --json ${name}.json --html ${name}.html
/public-supool/home/zhangzh/soft_share/hisat2-2.2.1/hisat2 -x /public-supool/home/zhangzh/index/hisat.index/CS -p 4 -5 10 --min-intronlen 20 --max-intronlen 4000 -1 ${name}_1P.fq.gz -2 ${name}_2P.fq.gz 2> ${name}.hisat2.log | /public-supool/home/zhangzh/conda/bin/samtools view -bS -q 20 -@ 4 | /public-supool/home/zhangzh/conda/bin/samtools sort - -o ${name}.unique.bam
/public-supool/home/zhangzh/soft_share/subread-2.0.6-Linux-x86_64/bin/featureCounts -p -T 8 \
	-g gene_id \
	-t exon \
	-a ${ann-/public-supool/home/zhangzh/index/IWGSC_v1.1_HC_20170706.gtf} \
	-o ${name}.rawreadcount \
	${name}.unique.bam
cat ${name}.rawreadcount | sed '1,2d'  | gawk -v FS='\t' -v OFS='\t' '{print$1,$7}' > ${name}.counts
