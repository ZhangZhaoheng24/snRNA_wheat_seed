#BSUB -L /bin/bash 
#BSUB -J scRNAyinxing
#BSUB -n 1 
#BSUB -e %RNAseq_%J.err 
#BSUB -o %RNAseq_%J.out 
#BSUB -q standardB
#BSUB -R "rusage[mem=200]"
STARindex=wheat01V
dir=snRNAdata
dnbc4tools rna run \
    --name DAP4_1 \
    --cDNAfastq1 ${dir}/E250037243_L01_6_1.fq.gz \
    --cDNAfastq2 ${dir}/E250037243_L01_6_2.fq.gz \
    --oligofastq1 ${dir}/E250037243_L01_14_1.fq.gz \
    --oligofastq2 ${dir}/E250037243_L01_14_2.fq.gz \
    --genomeDir $STARindex --threads 16

dnbc4tools rna run \
    --name DAP4_2 \
    --cDNAfastq1 ${dir}/E250037243_L01_7_1.fq.gz \
    --cDNAfastq2 ${dir}/E250037243_L01_7_2.fq.gz \
    --oligofastq1 ${dir}/E250037243_L01_15_1.fq.gz \
    --oligofastq2 ${dir}/E250037243_L01_15_2.fq.gz \
    --genomeDir $STARindex --threads 16

dnbc4tools rna run \
    --name DAP8_1 \
    --cDNAfastq1 ${dir}/E250037243_L01_8_1.fq.gz \
    --cDNAfastq2 ${dir}/E250037243_L01_8_2.fq.gz \
    --oligofastq1 ${dir}/E250037243_L01_16_1.fq.gz \
    --oligofastq2 ${dir}/E250037243_L01_16_2.fq.gz \
    --genomeDir $STARindex --threads 16

dnbc4tools rna run \
    --name DAP8_2 \
    --cDNAfastq1 ${dir}/E250037243_L01_9_1.fq.gz \
    --cDNAfastq2 ${dir}/E250037243_L01_9_2.fq.gz \
    --oligofastq1 ${dir}/E250037243_L01_17_1.fq.gz \
    --oligofastq2 ${dir}/E250037243_L01_17_2.fq.gz \
    --genomeDir $STARindex --threads 16


