#############################################################################################
########             shell script on HPC:  SNP annotation using snpEff              #########
#############################################################################################

java -jar snpEff.jar

#Once SnpEff is installed, we will enter the following commands to download the pre-built human database (GRCh37.75) that will be used to annotate our data.
cd snpEff
ml java/20+36-2344
java -jar snpEff.jar download -v GRCh37.75
java -jar snpEff.jar download -v GRCh38.99

#!/bin/bash
#BSUB -J annotate_vcf_grch37 # Job name
#BSUB -P acc_psychgen # allocation account
#BSUB -q premium # queue
#BSUB -n 1 # number of compute cores
#BSUB -W 96:00 # walltime in HH:MM
#BSUB -R rusage[mem=120000] # 32GB
#BSUB -R himem
#BSUB -R span[hosts=1] # all cores from one node
#BSUB -o %J.stdout # output log (%J : JobID)
#BSUB -eo %J.stderr # error log
#BSUB -L /bin/bash # Initialize the execution environment

ml java/20+36-2344
java -Xmx8g -jar /hpc/users/xuj24/snpeff/snpEff/snpEff.jar -c /hpc/users/xuj24/snpeff/snpEff/snpEff.config -v GRCh37.75 all.chr.brief.final.vcf > all.chr.brief.ann.vcf

bsub < annotate_vcf_grch37.sh

#!/bin/bash
#BSUB -J annotate_vcf_grch38 # Job name
#BSUB -P acc_psychgen # allocation account
#BSUB -q premium # queue
#BSUB -n 1 # number of compute cores
#BSUB -W 96:00 # walltime in HH:MM
#BSUB -R rusage[mem=120000] # 32GB
#BSUB -R himem
#BSUB -R span[hosts=1] # all cores from one node
#BSUB -o %J.stdout # output log (%J : JobID)
#BSUB -eo %J.stderr # error log
#BSUB -L /bin/bash # Initialize the execution environment

ml java/20+36-2344
java -Xmx8g -jar /hpc/users/xuj24/snpeff/snpEff/snpEff.jar -c /hpc/users/xuj24/snpeff/snpEff/snpEff.config -v GRCh38.99 all.chr.brief.final.vcf > all.chr.brief.ann.vcf

bsub < annotate_vcf_grch38.sh
