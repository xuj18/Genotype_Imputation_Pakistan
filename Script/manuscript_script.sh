## Part I. Pakistani GSA QC ###

## 1. Convert .idat to .bcf using the MoChA pipeline on Google Cloud (script: google_cloud.txt)
#to start the cromwell server on google VM
(java -XX:MaxRAMPercentage=90 -Dconfig.file=cromwell.conf -jar cromwell-85.jar server &)
java -jar cromwell-85.jar submit mocha.wdl -i pakistan_gsa_2433_vcf.json -o options_2433_gsa_vcf.json

# create the pakistan_gsa_2433_vcf.json
{
  "mocha.sample_set_id": "pakistan_gsa",
  "mocha.mode": "idat",
  "mocha.target": "vcf",
  "mocha.realign": true,
  "mocha.max_win_size_cm": 300.0,
  "mocha.overlap_size_cm": 5.0,
  "mocha.ref_name": "GRCh38",
  "mocha.ref_path": "gs://bigdeli-working/pakistan/GRCh38",
  "mocha.manifest_path": "gs://bigdeli-working/pakistan/manifests",
  "mocha.batch_tsv_file": "gs://bigdeli-working/pakistan/samples/pakistan_gsa.batch.tsv",
  "mocha.sample_tsv_file": "gs://bigdeli-working/pakistan/samples/pakistan_gsa.sample.tsv",
  "mocha.data_path": "gs://bigdeli-working/pakistan/idats"
}


# options_2433_gsa_vcf.json
{
  "delete_intermediate_output_files": true,
  "final_workflow_outputs_dir": "gs://bigdeli-working/pakistan/cromwell/outputs/gsa_2433_vcf",
  "use_relative_output_paths": true,
  "final_workflow_log_dir": "gs://bigdeli-working/pakistan/cromwell/wf_logs",
  "final_call_logs_dir": "gs://bigdeli-working/pakistan/cromwell/call_logs"
}


## 2. Convert .bcf to .vcf at Minerva 
ml bcftools
bcftools view --no-version -Oz -o pakistan_gsa.vcf.gz pakistan_gsa.bcf
# view:  subset, filter and convert VCF and BCF files
# --no-version: Do not append version and command line information to the output VCF header.
## -Oz: Output compressed VCF (z)
# -o: output file name

zcat pakistan_gsa.vcf.gz|wc -l # 653948 lines, 237 header lines --> 653,711 SNPs

zcat pakistan_gsa.vcf.gz|awk '(NR==237){print NF}' # 2442 columns, 9 non-ID columns (#CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT) --> 2433 individuals

# 3 remove SNPs with cluster separation scores <0.3, Gentrain score <0.4

## 3.1 GSA
zcat pakistan_gsa.vcf.gz|awk '{print $3,$8}'|sed 's/\;/ /g'|awk '(NR>237){print $1,$12,$13,$14,$15,$16}' > pakistan_gsa_GenTrain_Cluster_Sep

# in R
gsa <- read.table("pakistan_gsa_GenTrain_Cluster_Sep",header=F)
# 653711 6

# how many rows in the 2nd column are for GenTrain
sum(substring(gsa[,2],0,8)=="GenTrain") # 650037 rows have 3rd column for GenTrain
sum(substring(gsa[,3],0,8)=="GenTrain") # 3674 rows have 3rd column for GenTrain

sum(substring(gsa[,4],0,7)=="Cluster") # 633748 rows have 4th column for GenTrain
sum(substring(gsa[,5],0,7)=="Cluster") # 16289 rows have 4th column for GenTrain
sum(substring(gsa[,6],0,7)=="Cluster") # 3674 rows have 4th column for GenTrain

gentrain_1 <- gsa[which(substring(gsa[,2],0,8)=="GenTrain"),c(1,2)]
gentrain_2 <- gsa[which(substring(gsa[,3],0,8)=="GenTrain"),c(1,3)]
colnames(gentrain_1) <- c("V1","V2")
colnames(gentrain_2) <- c("V1","V2")
gentrain_bind <- rbind(gentrain_1, gentrain_2)

colnames(gentrain_bind) <- c("SNP","GenTrain_Score")
gentrain_bind$GenTrain_Score <- gsub('GenTrain_Score=','',gentrain_bind$GenTrain_Score)

cluster_sep_1 <- gsa[which(substring(gsa[,4],0,7)=="Cluster"),c(1,4)]
cluster_sep_2 <- gsa[which(substring(gsa[,5],0,7)=="Cluster"),c(1,5)]
cluster_sep_3 <- gsa[which(substring(gsa[,6],0,7)=="Cluster"),c(1,6)]
colnames(cluster_sep_1) <- c("V1","V3")
colnames(cluster_sep_2) <- c("V1","V3")
colnames(cluster_sep_3) <- c("V1","V3")
cluster_sep_bind <- rbind(cluster_sep_1,cluster_sep_2,cluster_sep_3)
colnames(cluster_sep_bind) <- c("SNP","Cluster_Sep")
cluster_sep_bind$Cluster_Sep <- gsub('Cluster_Sep=','',cluster_sep_bind$Cluster_Sep)

gentrain_cluster <- merge(gentrain_bind,cluster_sep_bind,by="SNP")
#  653711      3
gentrain_cluster$GenTrain_Score <- as.numeric(gentrain_cluster$GenTrain_Score)
gentrain_cluster$Cluster_Sep <- as.numeric(gentrain_cluster$Cluster_Sep)

par(mfrow=c(1,1))
hist(gentrain_cluster$GenTrain_Score)
abline(v=0.4,col=2)
hist(gentrain_cluster$Cluster_Sep)
abline(v=0.3,col=2)

gentrain_cluster$exclude_snp <- ifelse(gentrain_cluster$GenTrain_Score >= 0.4 & gentrain_cluster$Cluster_Sep >= 0.3, "keep","exclude") # 653,711 SNPs
keep <- gentrain_cluster[which(gentrain_cluster$exclude_snp =="keep"),] # 646113 SNPs to keep
exclude <- gentrain_cluster[which(gentrain_cluster$exclude_snp =="exclude"),"SNP"] # 7598 (1.2% SNPs will be excluded)

# GenTrain score: The GenTrain score is computed from the GenTrain 2.0 clustering algorithm. It is a measurement of SNP calling quality, ranging from 0 to 1, with higher value meaning better quality. 
# Cluster Separation score: which measures how well the AA, AB and BB clusters are separated. The cluster separation score also ranges from 0 to 1, with higher meaning better (more separation). 

write.table(exclude,"exclude_snp_gentrain_cluster_sep",col=FALSE,row=FALSE,quote=FALSE)


# Output a new vcf.gz file from the input vcf.gz file that removes SNPs with low GenTrain and Cluster Sep scores
vcftools --gzvcf pakistan_gsa.vcf.gz --exclude exclude_snp_gentrain_cluster_sep --recode --recode-INFO-all  --stdout  | gzip -c > pakistan_gsa_filtered.vcf.gz 
# 646,350 rows (= 646,113 SNPs + 237 headers)
# --stdout: Direct the vcftools output to standard out so it can be piped into another program or written directly to a filename of choice. However, a select few output functions cannot be written to standard out.
# --recode: Generate a new file in either VCF or BCF from the input VCF or BCF file after applying the filtering options specified by the user. The output file has the suffix ".recode.vcf" or ".recode.bcf".
# --recode: By default, the INFO fields are removed from the output file, as the INFO values may be invalidated by the recoding (e.g. the total depth may need to be recalculated if individuals are removed).
# --recode-INFO-all: can be used with the above recode options to define an INFO key name to keep in the output file. This option is used to keep all INFO values in the original file.


# how many SNPs in non-autosomal chromosomes
zcat pakistan_gsa_filtered.vcf.gz |awk '(NR>237) {print $1}' |sort|uniq -c # 4178 SNPs on chrY + 1128 SNPs on chrM + 41 SNPs (total: 5347 SNPs) on chr11_KI270721v1_random, chr15_KI270727v1_random, chr16_KI270728v1_random, chr22_KI270734v1_random, chr4_GL000008v2_random


## 4 convert to PLINK and exclude chrY and chrM

# genotype files contain chrX, chrY and chrM (mitochondria) too

# some weird chromosome names (chr11_KI270721v1_random) 
# Add the "--allow-extra-chr" (or "--aec" for short) flag to permit these codes in the input. (--chr 1-22 X will then filter them out.)

plink --vcf pakistan_gsa_filtered.vcf.gz --allow-extra-chr  --chr 1-23 --make-bed --out pakistan_gsa  # 646350-237-5347 = 640,766 SNPs
# 526004 variants and 2433 people pass filters and QC.
# --chr: if there are n autosomes, n+1 is the X chromosome, n+2 is Y, n+3 is XY, and n+4 is MT.

## 5.1a GSA SNP QC

# 5.1a.1 SNP call rate and duplicated SNPs
# No duplicated SNP names
awk '{print $2}' pakistan_gsa.bim|sort|uniq -d

# but there are duplicated SNPs based on chr:pos:ref:alt
plink --bfile pakistan_gsa --missing --out stats

# in R 
miss <- read.table("stats.lmiss",head=T)
bim <- read.table("pakistan_gsa.bim",head=F,sep='\t')

colnames(bim)[colnames(bim)=="V2"] <- "SNP"
miss_merge <- merge(miss,bim,by="SNP") # 640766 SNPs

# remove white space: trimws

miss_merge$chr_position_alleles      <- paste0(miss_merge$CHR,"_",miss_merge$V4,"_",miss_merge$V5,"_",miss_merge$V6)
# paste0: will not add any padding space introduced in the concatenated field.
# paste will add the space
miss_merge$chr_position_alleles_flip <- paste0(miss_merge$CHR,"_",miss_merge$V4,"_",miss_merge$V6,"_",miss_merge$V5)


# some SNPs could be duplicated (chr-pos-ref-alt matches chr-pos-ref-alt)
# some SNPs could be duplicated (chr-pos-ref-alt matches chr-pos-alt-ref)

forward_match <- miss_merge$chr_position_alleles[duplicated(miss_merge$chr_position_alleles)]  # 1609 SNPs (there are also duplicated SNPs within these 1609 SNPs)

reverse_match <- miss_merge$chr_position_alleles[which(miss_merge$chr_position_alleles %in% miss_merge$chr_position_alleles_flip)] # 374 SNPs

dim(miss_merge[which(miss_merge$chr_position_alleles  %in% forward_match),]) # 3175 lines to check forward duplicated pairs
dim(miss_merge[which(miss_merge$chr_position_alleles  %in% reverse_match),]) # 374 lines to check reverse duplicated pairs
dim(miss_merge[which(miss_merge$chr_position_alleles_flip  %in% reverse_match),]) # 374 lines to check reverse duplicated pairs

fwd_rvs_match <- forward_match[which(forward_match %in% reverse_match)] # 8 duplicated SNPs in both forward and reverse lines
# [1] "11_2570664_G_GTGGTCCGCCTC" "3_36996627_CA_C"          
# [3] "16_2088278_T_TC"           "5_112835090_C_CT"         
# [5] "3_37028890_GT_G"           "17_2680205_TG_T"          
# [7] "19_11105459_A_AG"          "13_32337943_G_GA"    

forward_match_new <- forward_match[!forward_match %in% fwd_rvs_match] # 1601 SNPs
reverse_match_new <- reverse_match[!reverse_match %in% fwd_rvs_match] # 358 SNPs

reverse_match[reverse_match %in% fwd_rvs_match]

dat <- miss_merge[miss_merge$chr_position_alleles %in% fwd_rvs_match |miss_merge$chr_position_alleles_flip %in% fwd_rvs_match ,  ]
dat[order(dat$V4),] # 24 rows (8+8+8)
# to keep: rs397515139, 11:2570674, rs113994201, 19:11216136, rs397507677, rs63750867, rs63749916, 5:112170790

forward_match_keep <- array(NA, dim=c(1601,1))
for (i in 1:length(forward_match_new)) {
   id <- forward_match_new[i]
   dat <- miss_merge[which(miss_merge$chr_position_alleles==id),]
   forward_match_keep[i,1] <- dat$SNP[which.min(dat$N_MISS)]
}
# if there are 3 lines, then the same SNP will be picked 3 times (as there will be 3 dats in the loop)

# need to remove duplicate
forward_match_keep[!duplicated(forward_match_keep)] # 43 duplication, 1558 non-duplicated SNPs


reverse_match_keep <- array(NA, dim=c(358,1))
for (i in 1:length(reverse_match_new)) {
   id <- reverse_match_new[i]
   dat <- miss_merge[which(miss_merge$chr_position_alleles==id| miss_merge$chr_position_alleles_flip==id),]
   reverse_match_keep[i,1] <- dat$SNP[which.min(dat$N_MISS)] # if the missing rate is the same, the first one is picked
}

reverse_match_keep[!duplicated(reverse_match_keep)] # 175 duplication, 183 non-duplicated SNPs (175+183 = 358)

# duplicated SNPs to keep: 1558 + 183 = 1741
forward_match_keep_final <- as.vector(forward_match_keep[!duplicated(forward_match_keep)])
reverse_match_keep_final <- as.vector(reverse_match_keep[!duplicated(reverse_match_keep)])
forward_reverse_match_keep <- c("rs397515139", "11:2570674", "rs113994201", "19:11216136", "rs397507677", "rs63750867", "rs63749916", "5:112170790") 
# these 8 SNPs are duplicated in reverse_match_keep_final

other_non_dup <- miss_merge[which(! (miss_merge$chr_position_alleles  %in% forward_match | miss_merge$chr_position_alleles  %in% reverse_match) ), "SNP"] # 637,233 SNPs 
length(miss_merge[which(miss_merge$chr_position_alleles  %in% forward_match | miss_merge$chr_position_alleles  %in% reverse_match ), "SNP"]) # 3533 SNPs (+637,233 = 640,766 SNPs)

snp_to_keep  <- c(forward_match_keep_final, reverse_match_keep_final, other_non_dup) # 638,974 SNPs

write.table(snp_to_keep,"snp_to_keep_no_dup.snplist",col=FALSE,row=FALSE,quote=FALSE) #638,974 SNPs

plink --bfile pakistan_gsa --extract snp_to_keep_no_dup.snplist --make-bed --out pakistan_gsa_QC1
# 638974 variants and 2433 people pass filters and QC.

# 5.1a.1.1 load case/control status (SCZ) (skip for now, have to figure out the sample swap issues)
# there are 2 different versions of keys, but neither of these seem to be correct, supposed duplicated samples are unrelated)

scz_status <- ayub[,c("Sema4.MatrixTubeBarcode","Sema4.MatrixTubeBarcode","Origin.CaseControlStatus")] # 2636 individuals

library(dplyr)
scz_status$pheno <- case_when(scz_status$Origin.CaseControlStatus =="case" ~ 2, scz_status$Origin.CaseControlStatus =="Case" ~ 2,scz_status$Origin.CaseControlStatus =="ctrl" ~ 1)
# 1583 cases (2), 1053 controls (1)

scz_status <- scz_status[,c("Sema4.MatrixTubeBarcode","Sema4.MatrixTubeBarcode","pheno")]

write.table(scz_status,"scz_status",col=FALSE,row=FALSE,quote=FALSE)

# update the phenotype 
plink --bfile pakistan_gsa_QC1 --make-bed --pheno scz_status --out pakistan_gsa_QC1_pheno
# 638974 variants and 2433 people pass filters and QC.

### 5.1a.1.2 check duplicate samples and remove (check for their call rate and phenotype) (skip for now)
# first prune common SNPs on autosomes 
plink --bfile pakistan_gsa_QC1_pheno --indep-pairwise 50 10 0.2 --out paut --maf 0.01  --autosome
# 183,902 SNPs in paut.prune.in

# calculate identify by descent IBD on these pruned SNPs to check for duplicated samples
plink --bfile pakistan_gsa_QC1_pheno --extract paut.prune.in --genome --out ibd
# duplicated samples: z0=0, z1=0, z2=1
# the probability of pairs sharing zero, one, and two IBD alleles
# 2,958,528 paris (number of possible pairs = 2433*2432/2 = 2,958,528)

genome <- read.table("ibd.genome",head=T)

dup_pairs <- genome[which(genome$Z0==0 & genome$Z1==0 & genome$Z2==1),]
plot(genome$Z1, genome$Z2)
plot(genome$Z0, genome$Z2)
dup_pairs_test <- genome[which(genome$Z0==0 & genome$Z1==0 ),]

# check the 39 pairs with same content after sample swap
duplicate_sample_pairs

for (i in 1:88) {
    id_1 <- duplicate_sample_pairs[i,1]
    id_2 <- duplicate_sample_pairs[i,2]
    dat <- genome[which(genome$FID1 == id_1 & genome$FID2 == id_2),]
    if (dim(dat)[1] ==0) next
    else {
    print(paste0(dat$FID1," ", dat$FID2, " ", dat$Z0," ", dat$Z1," ", dat$Z2))
    }
}

genome[which(genome$FID1 %in% dup_samples_old_id$old_sample_id & genome$FID2 %in% dup_samples_old_id$old_sample_id),]

# 5.1b.1.3 exclude SNPs with different missingness across cases and controls (SNP call rate by case control status)


# 5.1a.2 SNP call rate and MAF 
plink --bfile pakistan_gsa_QC1 --missing --out stats
plink --bfile pakistan_gsa_QC1 --freq --out stats #  638,974 SNPs (+1 header)



# in R
miss <- as.matrix(read.table("stats.lmiss",head=T))
freq <- as.matrix(read.table("stats.frq",head=T)) # variant-based missing data

missf <- as.numeric(miss[,5])
freqf <- as.numeric(freq[,5])

hist(1-missf,n=100,xlab="CallRate",main="") # n=100, more bars where each bar represents a smaller range
abline(v=.95,col=2,lwd=3) # lwd: line width
hist(freqf,n=100,xlab="MAF",main="")
abline(v=0.001,col=2,lwd=2)


# 5.1a.3 SNP HWE
# identify individuals with low autozygosity

# first prune common SNPs on autosomes (r2>0.8, same as Huang 2022, pick 0.8 over 0.2 because for psych chip, if we pick 0.2, no ROH will be found)
plink --bfile pakistan_gsa_QC1 --indep-pairwise 50 10 0.8 --out paut --maf 0.01  --autosome
# 376820 lines in paut.prune.in 
# --indep-pairwise: <window size>['kb'] <step size (variant ct)> <r^2 threshold>

# next, called runs of homozygosity (plink with default parameters) only in autosomes
plink --bfile pakistan_gsa_QC1 --extract paut.prune.in --chr 1-22 --homozyg

# calculate RoHs%, FROH: ROH-based inbreeding
# From runs of homozygosity (ROH), individual inbreeding/consanguinity coefficients can be calculated as: FROH=∑LROH/Lgenome
#where ∑LROHis the sum of the length of all ROH detected in an individual, and Lgenome is the total length of the genome that was used.

# To calculate FROH, a percentage of homozygosity was calculated by summing ROHs >1 Mb (1mb = 1000kb) across the covered autosomal genome 
# and dividing by the total autosomal base pairs represented in the SNP data. Specifically, the summed length of identified ROHs were divided by a factor of 2,772.7 
# and subsequently converted to a percent by multiplying the dividend by 100. A factor of 2,772.7 is the number of megabases covered by SNPs after imputation
# which was calculated by summing the distance between the first and the last available consecutive SNP on each chromosomal arm for each of the 22 autosomes. 
# Christofidou 2015 https://pubmed.ncbi.nlm.nih.gov/26166477/

# .hom.ind
# NSEG	Number of runs of homozygosity
# KB	Total length of runs (kb)
# KBAVG	Average length of runs (kb)

long_roh <- read.table("plink.hom",head=T)
ind_roh <- read.table("plink.hom.indiv",head=T)

sum(test_1_ind$KB) #  210425.9, consistent with the first line in plink.hom.indiv: NSEG=28, KB= 210426
head(ind_roh)

ind_roh$MB <- ind_roh$KB/1000
ind_roh$FROH <- ind_roh$MB/2772.7  *100 # 1-100%

# low autozygosity: individuals who had a fraction of the genome in RoHs <0.5%
hist(ind_roh$FROH,n=100) # FROH ranges from 0 to 43
abline(v=0.5,col=2,lwd=3)

dim(ind_roh[which(ind_roh$FROH<0.5),]) # 459 out of 2433 individuals (19%) have RoHs < 0.5% # 826 individuals if the r2 cutoff changes to 0.2

keep <- which(ind_roh$FROH<0.5)

write.table(ind_roh[keep,1:2],"keepind_low_ROH.fam",row=F,col=F,quote=F)

# Run the HWE test among those with low autozygosity

wc -l keepind_low_ROH.fam
plink --bfile pakistan_gsa_QC1 --hardy --keep keepind_low_ROH.fam  --out stats

hardy <- as.matrix(read.table("stats.hwe",head=T))
hardyp <- as.numeric(hardy[,9])

hist(-log10(hardyp),n=100,xlab="-log10(P) from HWE Test") # ranges from 0 to 69
abline(v=6,col=2,lwd=3)
summary(-log10(hardyp))


## Extract SNPs that have low call rate (<0.01), common (MAF>0.01), and in HWE (P>1e-6)
length(which(missf>=0.05)) # 1,593 SNPs
length(which(freqf <= 0.001)) # 46,288 SNPs
length(which(-log10(hardyp) >= 6)) # 68,735 SNPs  # if the r2 cutoff changes to 0.8, the number of SNPs excluded is 67,527 SNPs

length(which(missf>=0.05 |freqf <= 0.001 | -log10(hardyp) >= 6) ) # 116,268 SNPs
extract <- which(missf<0.05 & freqf > 0.001 & -log10(hardyp)< 6) #  522,706 SNPs  (+ 116,268 = 638,974)
summary(missf[extract])
summary(freqf[extract])
summary(-log10(hardyp[extract]))
write.table(hardy[extract,2],"extract.snplist",row=F,col=F,quote=F) # extract 522,706 SNPs

## Exclude SNPS that fail call rate (missingness), MAF and HWE
plink --bfile pakistan_gsa_QC1 --make-bed --extract extract.snplist --out pakistan_gsa_QC2
# 522,706 variants and 2433 people pass filters and QC.

########## 5.2a GSA individual QC
# 5.2a.1 update the phenotype 
plink --bfile pakistan_gsa_QC2 --make-bed --pheno scz_status --out pakistan_gsa_QC2_pheno
# 522706 variants and 2433 people pass filters and QC.

# 5.2a.2 individual heterozygosity and missingness
# first prune common SNPs on autosomes
plink --bfile pakistan_gsa_QC2_pheno --indep-pairwise 50 10 .2 --out prune_aut --maf 0.01  --autosome
# 181,820 lines in prune_aut.prune.in
plink --bfile pakistan_gsa_QC2_pheno --missing --out stats --extract prune_aut.prune.in
plink --bfile pakistan_gsa_QC2_pheno --het --out stats --extract prune_aut.prune.in

# in R
miss <- as.matrix(read.table("stats.imiss",header=T))
het <- as.matrix(read.table("stats.het",header=T))
missf <- as.numeric(miss[,6])
hetf <- 1-as.numeric(het[,3])/as.numeric(het[,5])

hist(missf,xlab="Missingness")
hist(hetf,xlab="Heterozygosity")

plot(missf,hetf,xlab="Missingness",ylab="Heterozygosity")
abline(v=0.05,col=2);abline(h=c(0.15,0.23),col=2) # decide by eye
keep <- which(missf<0.05&hetf<0.23&hetf>0.15) # 72 excluded, 2361 kept
lines(missf[keep],hetf[keep],col=2,pch=19,type="p") 
# type="p": point
# col: color
# pch: modify the symbol of the points
# https://r-coder.com/plot-r/
write.table(miss[keep,1:2],"keepind.fam",row=F,col=F,quote=F)

# remove those with high missing rate (5%) or are heterozygosity outliers, + only keep autosomes 
plink --bfile pakistan_gsa_QC2_pheno --make-bed --keep keepind.fam --chr 1-22 --out pakistan_gsa_QC3
# 520234 variants and 2361 people pass filters and QC.

# 5.2a.3 check duplicates and relatedness (select the ones with highest call rate)
# first prune common SNPs on autosomes 
plink --bfile pakistan_gsa_QC3 --indep-pairwise 50 10 0.2 --out paut --maf 0.01  --autosome
# 174,557 SNPs in paut.prune.in

## Round 1 
## use the kinship coef
# use plink2 to get KING coefficient
/hpc/packages/minerva-common/plink2/2.0a/plink2 --bfile pakistan_gsa_QC3 --extract paut.prune.in  --make-king-table --king-table-filter 0.088  --out related
# 198 pairs
# which plink (to check for path for plink2)
# 2,785,980 paris (number of possible pairs = 2361*2360/2 = 2,785,980)
# KINSHIP column has the kinship coefficient 
# mono twins (kinship coef > 0.707), 1st degree (0.354), 2nd degree (0.177), 3rd degree (0.088), unrelated (<0.088)
# get the same results using pruned set of SNPs or all SNPs
# get the same results by manually restricting to 0.088 or above
# /hpc/packages/minerva-common/plink2/2.0a/plink2 --bfile pakistan_psych_chip_QC3 --extract paut.prune.in  --make-king-table  --out related_test 

# remove the # at the front of the column headers in related.kin0
king <- read.table("related.kin0",head=T) # [1] 197   8
head(king)

miss <- read.table("stats.imiss",header=T)

miss_fid1 <- miss[,c("FID","F_MISS")]
miss_fid2 <- miss[,c("FID","F_MISS")]
colnames(miss_fid1) <- c("FID1","F_MISS1")
colnames(miss_fid2) <- c("FID2","F_MISS2")

king <- merge(king,miss_fid1,by="FID1")
king <- merge(king,miss_fid2,by="FID2") # 197 10

length(unique(king$FID1)) # 140 unique IDs
length(unique(king$FID2)) # 137 unique IDs

related_id <- unique(c(king$FID1,king$FID2)) # 246 related individuals

all_related_id <- c(king$FID1,king$FID2) # 394 related individuals

all_related_id[duplicated(all_related_id)] # 148 individuals occur more than 1 time

# first round (for each pair, remove the one with lower call rate)

ID_remove <- array(NA,dim=c(197,1))
for (i in 1: 197) {
    if (king$F_MISS1[i] <= king$F_MISS2 [i]) {
      ID_remove[i,1] <- king$FID2[i]
    }
    else if (king$F_MISS1[i] > king$F_MISS2 [i]) {
      ID_remove[i,1] <- king$FID1[i]
    }
}

ID_remove_unique <- unique(ID_remove) # 137 subjects
ID_remove_unique <- ID_remove_unique[,c(1,1)]

write.table(ID_remove_unique,"related_ID_remove_round1.fam",col=F,row=F,quote=F)


plink --bfile pakistan_gsa_QC3 --make-bed --remove related_ID_remove_round1.fam --out pakistan_gsa_QC4
# 2361 - 137 = 2224 individuals,  520234  variants

## Round 2
## use the kinship coef
# use plink2 to get KING coefficient
/hpc/packages/minerva-common/plink2/2.0a/plink2 --bfile pakistan_gsa_QC4 --extract paut.prune.in  --make-king-table --king-table-filter 0.088  --out related_check
# all 3rd degree or closer relatives are excluded!! 

# calculate identify by descent IBD on these pruned SNPs to check for duplicated samples
plink --bfile pakistan_gsa_QC4 --extract paut.prune.in --genome --out ibd
# duplicated samples: z0=0, z1=0, z2=1
# the probability of pairs sharing zero, one, and two IBD alleles
# 2,471,976 paris (number of possible pairs = 2224*2223/2 = 2,471,976)

# in R
ibd <- read.table("ibd.genome",head=T)

plot(ibd$Z0,ibd$Z1,xlim=c(0,1),ylim=c(0,1))
plot(ibd$Z0,ibd$Z2,xlim=c(0,1),ylim=c(0,1))
plot(ibd$Z1,ibd$Z2,xlim=c(0,1),ylim=c(0,1))

# 5.2a.3 infer sex
# currently in the .fam file, the column for sex has all 0 values
# prune common SNPs on chr X 

# create the reported sex file with the file
ayub <- read.table("../Ayub_sample",head=T) # 2636 individuals

summary(as.factor(ayub$Gender))
#     F Female      M   Male   NA's 
#   352     84    979    168   1053 
library(dplyr)
ayub$reported_gender = case_when(ayub$Gender =="M" ~ 1, ayub$Gender =="Male" ~ 1,ayub$Gender =="F" ~ 2, ayub$Gender =="Female" ~ 2,is.na(ayub$Gender) ~ 0)
table(ayub$reported_gender)
#   0    1    2 
# 1053 1147  436 

rep_sex <- ayub[,c("Sema4.MatrixTubeBarcode","Sema4.MatrixTubeBarcode","reported_gender")] # 2636 IDs

# remove 2 duplicated IDs
rep_sex[duplicated(rep_sex[,1]),]

dedup_rep_sex <- rep_sex[!duplicated(rep_sex[,1]),] # 2634 IDs

write.table(dedup_rep_sex,"gsa_reported_gender",row=F,col=F,quote=F)
# 2634 lines

# update sex column with reported gender (need to use the original plink data, coz the current one doesn't have X chr included, also after QC, very few SNPs in chrX left)
plink --bfile pakistan_gsa --make-bed --update-sex gsa_reported_gender --out pakistan_gsa_updatesex
# 640766 variants and 2433 people pass filters and QC.
# include the individuals in pakistan_gsa_QC4.fam
awk '{print $1,$2}' pakistan_gsa_QC4.fam > keepind.fam 
plink --bfile pakistan_gsa_updatesex --make-bed --keep keepind.fam --split-x hg38 --out pakistan_gsa_sexcheck
# --split-x: 727 chromosome codes changed.
# 2224 individuals
plink --bfile pakistan_gsa_sexcheck --indep-pairwise 50 10 0.2 --out pX --maf 0.01 --chr 23 
# 4,849 SNPs on chrX are kept in pX.prune.in
plink --bfile pakistan_gsa_sexcheck --extract pX.prune.in --check-sex --out sex
# 2224 individuals

# in R
sex <- read.table("sex.sexcheck",header=T) # 2224 individuals
# check how many mismatch sex among those with reported sex
rep_snp_sex <- sex[which(sex$PEDSEX !=0),] # 1332 individuals
# 1015 reported males, 385 reported females
table(rep_snp_sex$PEDSEX, rep_snp_sex$SNPSEX) 

#       0   1   2 (SNPSEX)
#  1    5 924  37
#  2   29  41 296

# among individuals with reported gender, 34 have ambiguous sex 
# among individuals with reported gender, 37 reported males are biologically females, 41 reported females are biologically males (balanced)
# sex mismatches did not favor a particular direction, indicating that sample swaps, not contamination, was the likely problem.
 

plot(sex$F,col=as.numeric(sex$SNPSEX)+1, pch=19,ylab="Estimated Heterozygosity",xlab="Individual") 
# Females (XX) will have (expected) heterozygosity 0, males (XY) will have expected heterozygosity of 1
table(sex$SNPSEX)
# sex[,4]= SNPSEX (values: 0, 1, 2)
# To perform table(sex$column_name), it has to be read as a dataframe, not a matrix
# values:    0   1    2 
# counts:   36 1728  460 
male_assign <- sex[which(sex$SNPSEX==1),] # 1728 male sex
female_assign <- sex[which(sex$SNPSEX==2),] # 460 female sex

summary(as.numeric(sex$F)) 
summary(as.numeric(male_assign$F)) # ranges from 0.9742 to  0.9992
summary(as.numeric(female_assign$F)) # ranges from -0.070890 to 0.1975


male=which(as.numeric(sex[,6])>0.8) # 1728 males
female=which(as.numeric(sex[,6])<0.2) # 460 females 
ambiguous =which(as.numeric(sex[,6])>=0.2 & as.numeric(sex[,6])<=0.8 ) # 36 ambiguous

newsex=array(0,nrow(sex))
newsex[male]=1
newsex[female]=2
table(sex$SNPSEX,newsex) # 0 - ambiguous, 1 - males, 2 - females
#    newsex
#        0   1   2
#  0   36    0    0
#  1    0 1728    0
#  2    0    0  460
write.table(cbind(sex[,1:2],newsex),"all.sex",row=F,col=F,quote=F) # 2224 individuals
write.table(sex[c(male,female),1:2],"keepind2.fam",row=F,col=F,quote=F) # 2188 individuals to keep (36 ambiguous sex to be excluded)


het <- read.table("stats.het",header=T)

sex_het <- merge(sex,het,by=c("FID","IID")) # 2224

sex_het$new_hetf <- 1-as.numeric(sex_het$"O.HOM.")/as.numeric(sex_het$"N.NM.")
plot(sex_het$F.x,sex_het$new_hetf,xlab="ChrX Heterozygosity",ylab="Autosome Heterozygosity", col=as.numeric(sex_het$SNPSEX)+1, pch=19)

# exclude ambiguous X heterozygosity
plink --bfile pakistan_gsa_QC4 -make-bed --keep keepind2.fam --update-sex all.sex --out pakistan_gsa_final
# 520234 variants and 2188 people pass filters and QC.
# 1728 males, 460 females, 1298 cases, 890 controls

# exclude all potential mismatched IDs (n=374)
plink --bfile pakistan_gsa_final --make-bed --remove /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/incorrect_gsa_ID_union --out pakistan_gsa_final_match_only
# 520234 variants and 1814 people pass filters and QC.
# 1456 males, 358 females, 1079 cases, 735 controls.

# check if the compute_sex file from mocha matches with my plink-assigned sex: YES, they match! 
compute_sex <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_2433_unphased/pakistan_gsa.sample.tsv",header=T)
plink_sex <- read.table("pakistan_gsa_final_match_only.fam",header=F) # 1814 row
colnames(compute_sex)[colnames(compute_sex)=="sample_id"] <- "V1"

dat <- merge(plink_sex,compute_sex,by="V1") #1814 rows
table(dat$V5,dat$computed_gender) # all

# 5.2a.4 check for Pakistani ancestry
#Always use both pruning and removing of long-range LD regions when computing the PCs

#remove high LD regions
# high_LD_region_grch38 has the starting and ending positions of SNPS in high LD regions

plink --bfile pakistan_gsa_final --make-set /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/high_LD_region_grch38 --write-set --out high_ld
#With --make-set, you can define sets from a list of named bp ranges instead. 
# Each line of the --make-set input file is expected to have the following 4-5 fields in front:
# 1. Chromosome code, 2. Start of range (base-pair units), 3. End of range (this position is included in the interval), 4. Set ID
# 15807 lines in high_ld.set (there are lines for 1, 2, 3, ... 43 sets, and lines for END, and empty lines)
# example
# rs4494819
# GSA-rs75758696
# END

# 42
# END

# 43
# END
plink --bfile pakistan_gsa_final --exclude high_ld.set  --write-snplist  --make-just-fam --out pakistan_gsa_final_lowLD
# 504570 variants and 2188 people pass filters and QC.

# prune SNPs before PCA 
plink --bfile pakistan_gsa_final --keep pakistan_gsa_final_lowLD.fam --extract pakistan_gsa_final_lowLD.snplist --indep-pairwise 50 10 0.2 --out paut --maf 0.01  --autosome
# 172,305 SNPs in paut.prune.in

#run the PCA analysis
plink --bfile pakistan_gsa_final --extract  paut.prune.in --pca 20 --out pca
# 172305 variants and 2188 people pass filters and QC.
# either --pca 20 or --pca 80 have the same results for PC1 ~ PC 20

# in R 
eigenval <- read.table("pca.eigenval", header=FALSE)
eigenvec <- read.table("pca.eigenvec", header=FALSE) # 2188 individuals

# variance per PC
plot(eigenval$V1,type="p")
plot(eigenval$V1[c(1:20)],type="p") # seems to plateau after first 5 PCs


plot(eigenvec$V3,eigenvec$V4,xlim=c(-0.15,0.15),ylim=c(-0.15,0.15))
# V3 is PC1, V4 = PC2
plot(eigenvec$V5,eigenvec$V6)
plot(eigenvec$V7,eigenvec$V8)
plot(eigenvec$V9,eigenvec$V10)
plot(eigenvec$V11,eigenvec$V12)

ub_pc1 <- mean(eigenvec$V3)+3*sd(eigenvec$V3) # 0.06414992
lb_pc1 <- mean(eigenvec$V3)-3*sd(eigenvec$V3) # -0.06415014
ub_pc2 <- mean(eigenvec$V4)+3*sd(eigenvec$V4) # 0.06414958
lb_pc2 <- mean(eigenvec$V4)-3*sd(eigenvec$V4) # -0.06415049

dim(eigenvec) # 2188 individuals

eigenvec_new <- eigenvec[which(eigenvec$V3 >= lb_pc1 & eigenvec$V3 <= ub_pc1),] # 2187 individuals
keep <- which(eigenvec$V3 >= lb_pc1 & eigenvec$V3 <= ub_pc1)
lines(eigenvec$V3[keep],eigenvec$V4[keep],col=2,pch=19,type="p") 

plot(eigenvec$V3,eigenvec$V4,xlim=c(-0.15,0.15),ylim=c(-0.15,0.15))
test <- which(eigenvec$V3 >= lb_pc1 & eigenvec$V3 <= ub_pc1 & eigenvec$V4 >= lb_pc2 & eigenvec$V4 <= ub_pc2) # 15 excluded, 1566 kept
lines(eigenvec$V3[test],eigenvec$V4[test],col=2,pch=19,type="p") 

# decided not to exclude any PC outliers

### 5.2a.5 merge with 1000 Genome data and run PCA plot

## use the hg38 1000G that Emily merged
# path: /sc/arion/projects/psychgen/resources/genotype_ref_panel/1000g/1000g_hg38_merged_plink

## Update SNP names to chr:position in 1000G

awk '{print $2"\t"$1":"$4}' /sc/arion/projects/psychgen/resources/genotype_ref_panel/1000g/1000g_hg38_merged_plink/full.autosomes.rsid.snpdedup.hg38.dedup.bim | sort -u > 1000genomes_updateSNPs 
# 39192947 SNPs 
# To sort and remove duplicates pass the -u option to sort. This will write a sorted list to standard output and remove duplicates. 

# in the folder of  genotype_call_gsa_2219
plink --bfile /sc/arion/projects/psychgen/resources/genotype_ref_panel/1000g/1000g_hg38_merged_plink/full.autosomes.rsid.snpdedup.hg38.dedup --update-name 1000genomes_updateSNPs --make-bed --out 1000genomesUpdateNames
# 39192947 variants and 2504 people pass filters and QC.

## Update SNP names to chr:position in GSA
awk '{print $2"\t"$1":"$4}' pakistan_gsa_final_match_only.bim  | sort -u > gsa_updateSNPs # 520234 SNPs 
plink --bfile pakistan_gsa_final_match_only --update-name gsa_updateSNPs --make-bed --out pakistan_gsa_final_match_only_chrbp
# 520234 variants and 1814 people pass filters and QC.

awk '{print $2}'  pakistan_gsa_final_match_only_chrbp.bim |sort > pakistan_gsa_final_match_only_chrbp.snp # 520234 SNPs
awk '{print $2}'  1000genomesUpdateNames.bim  |sort > 1000genomesUpdateNames.snp # 39192947 SNPs

comm -12 pakistan_gsa_final_match_only_chrbp.snp 1000genomesUpdateNames.snp > overlapsnps # 508,812 out of 520,234 SNPs (~98% are overlapped)


### Merge 1000genomes and Pakistani plink files using only overlapping SNPs
#newmergelist includes Pakistani bed,bim and fam
plink --bfile 1000genomesUpdateNames --merge-list newmergelist --extract overlapsnps --make-bed --out pakistan_1000g_merged
# Error: 12007 variants with 3+ alleles present.

#multiallelic SNPs
grep -wFf pakistan_1000g_merged-merge.missnp 1000genomesUpdateNames.bim|wc -l 
# 23,457 SNPs
grep -wFf pakistan_1000g_merged-merge.missnp  pakistan_gsa_final_match_only_chrbp.bim|wc -l
# 871 SNPs

plink --bfile 1000genomesUpdateNames --exclude pakistan_1000g_merged-merge.missnp    --make-bed --out 1000genomesUpdateNames.wo.missnp # 39192947 SNPs --> 39169472 SNPs
plink --bfile pakistan_gsa_final_match_only_chrbp --exclude pakistan_1000g_merged-merge.missnp    --make-bed --out pakistan_gsa_final_match_only_chrbp.wo.missnp # 520234 SNPs --> 519363 SNPs

#excluded multiallelic SNPs, try merging again
awk '{print $2}' 1000genomesUpdateNames.wo.missnp.bim |sort >  1000genomesUpdateNames.wo.missnp.snplist # 39,169,472 SNPs
awk '{print $2}' pakistan_gsa_final_match_only_chrbp.wo.missnp.bim |sort >  pakistan_gsa_final_match_only_chrbp.wo.missnp.snplist # 519,363 SNPs

comm -12 pakistan_gsa_final_match_only_chrbp.wo.missnp.snplist 1000genomesUpdateNames.wo.missnp.snplist >  overlapsnps # 508419 SNPs out of 520,234 SNPs (~98% are overlapped)

#newmergelist2 includes pakistan_gsa_final_match_only_chrbp.wo.missnp bed,bim and fam
plink --bfile 1000genomesUpdateNames.wo.missnp --merge-list newmergelist2 --extract overlapsnps --make-bed --out pakistan_1000g_merged2
# 508419 variants and 4318 people pass filters and QC.

### prepare for PCA

#remove high LD regions
plink --bfile pakistan_1000g_merged2 --make-set /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/high_LD_region_grch38 --write-set --out high_ld

plink --bfile pakistan_1000g_merged2 --exclude high_ld.set  --write-snplist  --make-just-fam --out pakistan_1000g_merged2_lowLD
# 493256 variants and 4318 people pass filters and QC.

# prune SNPs before PCA 
plink --bfile pakistan_1000g_merged2 --keep pakistan_1000g_merged2_lowLD.fam --extract pakistan_1000g_merged2_lowLD.snplist --indep-pairwise 50 10 0.2 --out paut --maf 0.01  --autosome
# 193,330 SNPs in paut.prune.in
# another paper used 146,671 SNPs for PCA (We applied ProPCA to compute the top five PCs on genotype data from the UK Biobank,consisting of 488,363 individuals and 146,671 SNPs, https://doi.org/10.1371/journal. pgen.1008773)

#run the PCA analysis
plink --bfile pakistan_1000g_merged2 --extract  paut.prune.in --pca 20 --out gsa_1000g_pca
# 193330 variants and 4318 people pass filters and QC.

# sample ancestry information for 1000G ppts
# path /sc/arion/projects/psychgen/resources/1000Genome_phase3_ref_panel
# ancestry info all together for 1000G: /sc/arion/projects/psychgen/resources/genotype_ref_panel/1000g/integrated_call_samples_v3.20130502.ALL.panel

## in R
ancestry_1000G <- read.table("/sc/arion/projects/psychgen/resources/genotype_ref_panel/1000g/integrated_call_samples_v3.20130502.ALL.panel",header=TRUE) # super-population category for 2504 individuals in 1000G
eigenval <- read.table("gsa_1000g_pca.eigenval", header=FALSE)
eigenvec <- read.table("gsa_1000g_pca.eigenvec", header=FALSE) 

colnames(eigenvec)[colnames(eigenvec)=="V2"] <- "sample"
eigenvec_new <- merge(eigenvec,ancestry_1000G,by="sample",all.x=TRUE) # 4318 individuals

table(eigenvec_new$super_pop) # 2504 individuals in 1000G
# AFR AMR EAS EUR SAS 
# 661 347 504 503 489 
tmp <- as.data.frame(eigenval[c(1:10),])
plot(tmp[,1],type="l")
#keep 6 PCs

# put Pakistani samples as PAK in super_pop 
eigenvec_new$super_pop[is.na(eigenvec_new$super_pop)] <- "PAK"

table(eigenvec_new$super_pop) # 1814 PAK
# AFR  AMR  EAS  EUR  PAK  SAS 
# 661  347  504  503 2219  489 

# save the plot 

pdf("PCA_plot_Pakistani_samples_with_1000G.pdf", width=10, height=10)
library(ggplot2)
p <- ggplot(eigenvec_new, aes(x = V3, y = V4)) + geom_point(aes(color = factor(super_pop))) + theme_minimal()
p + scale_color_manual(values=c("#FF595E","#1982C4","#0C7C59","#8AC926","#FF924C","#6A4C93")) + labs( x="PC1", y = "PC2")
dev.off()

# only to plot SAS vs PAK
eigenvec_sas <- eigenvec_new[which(eigenvec_new$super_pop %in% c("PAK","SAS")),]
ggplot(eigenvec_sas, aes(x = V3, y = V4)) + geom_point(aes(color = factor(super_pop))) + theme_minimal()


## 6. Create a separate vcf.gz file for each chromosome, and variations must be sorted by genomic position.

# 6.1 Recode plink to vcf by each chromosome
/hpc/packages/minerva-common/plink2/2.0a/plink2 --bfile pakistan_gsa_final_match_only --recode vcf   --out pakistan_gsa_final_match_only # only autosomes
# 28 lines of headers for pakistan_gsa_final_match_only.vcf
# should be 9+ 1814 = 1823 columns
# 28+520,234 = 520,262 lines
# 1823 columns, 520,262 rows (yes!) 
# all hard call, no dosage as this is genotyped data, no imputed data yet with probabilities!

# for the imputation website, if it is grch38 SNPs, the chromosome format should be chr1, instead of 1
/hpc/packages/minerva-common/plink2/2.0a/plink2 --bfile pakistan_gsa_final_match_only --recode vcf  --chr 1 --out pakistan_gsa_final_match_only_chr1 --output-chr chrM
/hpc/packages/minerva-common/plink2/2.0a/plink2 --bfile pakistan_gsa_final_match_only --recode vcf  --chr 2 --out pakistan_gsa_final_match_only_chr2 --output-chr chrM
/hpc/packages/minerva-common/plink2/2.0a/plink2 --bfile pakistan_gsa_final_match_only --recode vcf  --chr 3 --out pakistan_gsa_final_match_only_chr3 --output-chr chrM
/hpc/packages/minerva-common/plink2/2.0a/plink2 --bfile pakistan_gsa_final_match_only --recode vcf  --chr 4 --out pakistan_gsa_final_match_only_chr4 --output-chr chrM
/hpc/packages/minerva-common/plink2/2.0a/plink2 --bfile pakistan_gsa_final_match_only --recode vcf  --chr 5 --out pakistan_gsa_final_match_only_chr5 --output-chr chrM
/hpc/packages/minerva-common/plink2/2.0a/plink2 --bfile pakistan_gsa_final_match_only --recode vcf  --chr 6 --out pakistan_gsa_final_match_only_chr6 --output-chr chrM
/hpc/packages/minerva-common/plink2/2.0a/plink2 --bfile pakistan_gsa_final_match_only --recode vcf  --chr 7 --out pakistan_gsa_final_match_only_chr7 --output-chr chrM
/hpc/packages/minerva-common/plink2/2.0a/plink2 --bfile pakistan_gsa_final_match_only --recode vcf  --chr 8 --out pakistan_gsa_final_match_only_chr8 --output-chr chrM
/hpc/packages/minerva-common/plink2/2.0a/plink2 --bfile pakistan_gsa_final_match_only --recode vcf  --chr 9 --out pakistan_gsa_final_match_only_chr9 --output-chr chrM
/hpc/packages/minerva-common/plink2/2.0a/plink2 --bfile pakistan_gsa_final_match_only --recode vcf  --chr 10 --out pakistan_gsa_final_match_only_chr10 --output-chr chrM
/hpc/packages/minerva-common/plink2/2.0a/plink2 --bfile pakistan_gsa_final_match_only --recode vcf  --chr 11 --out pakistan_gsa_final_match_only_chr11 --output-chr chrM
/hpc/packages/minerva-common/plink2/2.0a/plink2 --bfile pakistan_gsa_final_match_only --recode vcf  --chr 12 --out pakistan_gsa_final_match_only_chr12 --output-chr chrM
/hpc/packages/minerva-common/plink2/2.0a/plink2 --bfile pakistan_gsa_final_match_only --recode vcf  --chr 13 --out pakistan_gsa_final_match_only_chr13 --output-chr chrM
/hpc/packages/minerva-common/plink2/2.0a/plink2 --bfile pakistan_gsa_final_match_only --recode vcf  --chr 14 --out pakistan_gsa_final_match_only_chr14 --output-chr chrM
/hpc/packages/minerva-common/plink2/2.0a/plink2 --bfile pakistan_gsa_final_match_only --recode vcf  --chr 15 --out pakistan_gsa_final_match_only_chr15 --output-chr chrM
/hpc/packages/minerva-common/plink2/2.0a/plink2 --bfile pakistan_gsa_final_match_only --recode vcf  --chr 16 --out pakistan_gsa_final_match_only_chr16 --output-chr chrM
/hpc/packages/minerva-common/plink2/2.0a/plink2 --bfile pakistan_gsa_final_match_only --recode vcf  --chr 17 --out pakistan_gsa_final_match_only_chr17 --output-chr chrM
/hpc/packages/minerva-common/plink2/2.0a/plink2 --bfile pakistan_gsa_final_match_only --recode vcf  --chr 18 --out pakistan_gsa_final_match_only_chr18 --output-chr chrM
/hpc/packages/minerva-common/plink2/2.0a/plink2 --bfile pakistan_gsa_final_match_only --recode vcf  --chr 19 --out pakistan_gsa_final_match_only_chr19 --output-chr chrM
/hpc/packages/minerva-common/plink2/2.0a/plink2 --bfile pakistan_gsa_final_match_only --recode vcf  --chr 20 --out pakistan_gsa_final_match_only_chr20 --output-chr chrM
/hpc/packages/minerva-common/plink2/2.0a/plink2 --bfile pakistan_gsa_final_match_only --recode vcf  --chr 21 --out pakistan_gsa_final_match_only_chr21 --output-chr chrM
/hpc/packages/minerva-common/plink2/2.0a/plink2 --bfile pakistan_gsa_final_match_only --recode vcf  --chr 22 --out pakistan_gsa_final_match_only_chr22 --output-chr chrM # 7833 lines (7 headers + 7826 SNPs), 1823 columns (9 header columns + 1814 participants)

# 6.2 Create a sorted vcf.gz file using BCFtools:

ml bcftools
bcftools sort pakistan_gsa_final_match_only_chr1.vcf -Oz -o pakistan_gsa_final_match_only_chr1.vcf.gz
bcftools sort pakistan_gsa_final_match_only_chr2.vcf -Oz -o pakistan_gsa_final_match_only_chr2.vcf.gz
bcftools sort pakistan_gsa_final_match_only_chr3.vcf -Oz -o pakistan_gsa_final_match_only_chr3.vcf.gz
bcftools sort pakistan_gsa_final_match_only_chr4.vcf -Oz -o pakistan_gsa_final_match_only_chr4.vcf.gz
bcftools sort pakistan_gsa_final_match_only_chr5.vcf -Oz -o pakistan_gsa_final_match_only_chr5.vcf.gz
bcftools sort pakistan_gsa_final_match_only_chr6.vcf -Oz -o pakistan_gsa_final_match_only_chr6.vcf.gz
bcftools sort pakistan_gsa_final_match_only_chr7.vcf -Oz -o pakistan_gsa_final_match_only_chr7.vcf.gz
bcftools sort pakistan_gsa_final_match_only_chr8.vcf -Oz -o pakistan_gsa_final_match_only_chr8.vcf.gz
bcftools sort pakistan_gsa_final_match_only_chr9.vcf -Oz -o pakistan_gsa_final_match_only_chr9.vcf.gz
bcftools sort pakistan_gsa_final_match_only_chr10.vcf -Oz -o pakistan_gsa_final_match_only_chr10.vcf.gz
bcftools sort pakistan_gsa_final_match_only_chr11.vcf -Oz -o pakistan_gsa_final_match_only_chr11.vcf.gz
bcftools sort pakistan_gsa_final_match_only_chr12.vcf -Oz -o pakistan_gsa_final_match_only_chr12.vcf.gz
bcftools sort pakistan_gsa_final_match_only_chr13.vcf -Oz -o pakistan_gsa_final_match_only_chr13.vcf.gz
bcftools sort pakistan_gsa_final_match_only_chr14.vcf -Oz -o pakistan_gsa_final_match_only_chr14.vcf.gz
bcftools sort pakistan_gsa_final_match_only_chr15.vcf -Oz -o pakistan_gsa_final_match_only_chr15.vcf.gz
bcftools sort pakistan_gsa_final_match_only_chr16.vcf -Oz -o pakistan_gsa_final_match_only_chr16.vcf.gz
bcftools sort pakistan_gsa_final_match_only_chr17.vcf -Oz -o pakistan_gsa_final_match_only_chr17.vcf.gz
bcftools sort pakistan_gsa_final_match_only_chr18.vcf -Oz -o pakistan_gsa_final_match_only_chr18.vcf.gz
bcftools sort pakistan_gsa_final_match_only_chr19.vcf -Oz -o pakistan_gsa_final_match_only_chr19.vcf.gz
bcftools sort pakistan_gsa_final_match_only_chr20.vcf -Oz -o pakistan_gsa_final_match_only_chr20.vcf.gz
bcftools sort pakistan_gsa_final_match_only_chr21.vcf -Oz -o pakistan_gsa_final_match_only_chr21.vcf.gz
bcftools sort pakistan_gsa_final_match_only_chr22.vcf -Oz -o pakistan_gsa_final_match_only_chr22.vcf.gz


# 7. transfer the files to Google Cloud Storage 
# go to transfer node
# ssh li03c01
# this has to loaded using home node, not computing node such as sklar1
ml google_cloud_sdk
ml firefox proxies
gcloud auth login 
#https://cloud.google.com/sdk/gcloud/reference/auth/login

#gsutil to copy
# /hpc/users/xuj24/google_cloud_sdk_345.0.0/google-cloud-sdk/bin/
gsutil cp pakistan_gsa_final_match_only.vcf gs://bigdeli-working/pakistan/cromwell/outputs/vcf_for_topmed/gsa_1814/
gsutil cp pakistan_gsa.batch.tsv gs://bigdeli-working/pakistan/samples/gsa_1814/
gsutil cp gs://bigdeli-working/pakistan/cromwell/outputs/gsa_2433_vcf/pakistan_gsa.sample.tsv /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_2433_unphased/
gsutil cp pakistan_gsa_removed_619_samples.tsv gs://bigdeli-working/pakistan/samples/gsa_1814/
gsutil cp pakistan_gsa_removed_SNP.recode.vcf  gs://bigdeli-working/pakistan/samples/gsa_1814/
gsutil cp gs://bigdeli-working/pakistan/cromwell/outputs/gsa_1814_phasing/pakistan_gsa_match.pgt.bcf .
gsutil cp gs://bigdeli-working/pakistan/cromwell/outputs/gsa_1814_phasing/pakistan_gsa_match.pgt.bcf.csi .
###

# get the list of samples removed (2433-1814 = 619 individuals removed)
og_fam <- read.table("pakistan_gsa.fam",header=F)
final_fam <- read.table("pakistan_gsa_final_match_only.fam",header=F)

removed_sample <- og_fam[which(!og_fam$V1 %in% final_fam$V1),"V1"] # 619 individuals
include_sample <- og_fam[which(og_fam$V1 %in% final_fam$V1),"V1"] # 1814 individuals

write.table(removed_sample,"pakistan_gsa_removed_619_samples.tsv",row=F,col=F,quote=F)
write.table(include_sample,"pakistan_gsa_included_1814_samples.tsv",row=F,col=F,quote=F)

# get the list of variants removed (653,711 - 520,234 = 133,477 SNPs)
zcat pakistan_gsa.vcf.gz|awk '(NR>237){print $1,$3}' > pakistan_gsa.bcf.snplist # 653,711 SNPs
awk '{print $2}' pakistan_gsa_final_match_only.bim > pakistan_gsa.match.snplist # 520,234 SNPs

og_snp <- read.table("pakistan_gsa.bcf.snplist",header=F)
final_snp <- read.table("pakistan_gsa.match.snplist",header=F)

removed_snp <- og_snp[which(!og_snp$V2 %in% final_snp$V1),] # 133,477 SNPs
removed_snp_wo_chrX <- removed_snp[which(removed_snp$V1 !="chrX"),"V2"] # 105,522 SNPs

write.table(removed_snp_wo_chrX,"pakistan_gsa_removed_SNPs_wo_chrX.tsv",row=F,col=F,quote=F,sep="\t")

# create a vcf file with list of additional variants to exclude from analysis
ml vcftools
vcftools --gzvcf pakistan_gsa.vcf.gz --snps pakistan_gsa_removed_SNPs_wo_chrX.tsv --recode --recode-INFO-all --out pakistan_gsa_removed_SNP
# pakistan_gsa_removed_SNP.recode.vcf  kept 105,522 out of a possible 653,711 Sites
# 105,759 lines (237 headers + 105,522 SNPs )
### phase using MoChA


# create the pakistan_gsa_1814_phasing.json
{
  "mocha.sample_set_id": "pakistan_gsa_match",
  "mocha.mode": "vcf",
  "mocha.max_win_size_cm": 50.0,
  "mocha.overlap_size_cm": 5.0,
  "mocha.ref_name": "GRCh38",
  "mocha.ref_path": "gs://bigdeli-working/pakistan/GRCh38",
  "mocha.batch_tsv_file": "gs://bigdeli-working/pakistan/samples/gsa_1814/pakistan_gsa.batch.tsv",
  "mocha.sample_tsv_file": "gs://bigdeli-working/pakistan/cromwell/outputs/gsa_2433_vcf/pakistan_gsa.sample.tsv",
  "mocha.data_path": "gs://bigdeli-working/pakistan/cromwell/outputs/gsa_2433_vcf/",
  "mocha.extra_xcl_vcf_file": "gs://bigdeli-working/pakistan/samples/gsa_1814/pakistan_gsa_removed_SNP.recode.vcf",
  "mocha.duplicate_samples_file": "gs://bigdeli-working/pakistan/samples/gsa_1814/pakistan_gsa_removed_619_samples.tsv"
}

# batch.tsv 
# batch_id	n_smpls	vcf	vcf_index
# 1	1814	pakistan_gsa.bcf	pakistan_gsa.bcf.csi

# sample.tsv 
# sample_id	computed_gender	call_rate

# options_1814_gsa_phasing.json
{
  "delete_intermediate_output_files": true,
  "final_workflow_outputs_dir": "gs://bigdeli-working/pakistan/cromwell/outputs/gsa_1814_phasing",
  "use_relative_output_paths": true,
  "final_workflow_log_dir": "gs://bigdeli-working/pakistan/cromwell/wf_logs",
  "final_call_logs_dir": "gs://bigdeli-working/pakistan/cromwell/call_logs"
}

(java -XX:MaxRAMPercentage=90 -Dconfig.file=cromwell.conf -jar cromwell-85.jar server &)
java -jar cromwell-85.jar submit mocha.wdl -i pakistan_gsa_1814_phasing.json -o options_1814_gsa_phasing.json
# Workflow 0f61d643-f420-4fb0-9a7f-097ad2c0ef95 submitted
curl -X GET http://localhost:8000/api/workflows/v1/0f61d643-f420-4fb0-9a7f-097ad2c0ef95/metadata | jq


## 8. exclude the participants in bcf before uploading to imputation server
# path: /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_phased_mocha
ml bcftools
bcftools view --no-version pakistan_gsa_match.pgt.bcf|head -n206 # 206 lines of header
bcftools view --no-version pakistan_gsa_match.pgt.bcf|wc -l
# 527,142 lines (526,936 SNPs)

bcftools view --no-version pakistan_gsa_match.pgt.bcf|awk '(NR==206){print NF}'  
# 2442 columns (9 columns of header + 2433 participants)
# 2433 participants with 526,936 SNPs

# convert bcf to vcf.gz
bcftools view --no-version -Oz -o pakistan_gsa_match.pgt.vcf.gz pakistan_gsa_match.pgt.bcf

zcat pakistan_gsa_match.pgt.vcf.gz|awk 'NR==206{print NF}' # 2442 columns
zcat pakistan_gsa_match.pgt.vcf.gz|wc -l # 527,142 lines

# only include 1814 participants
bcftools view -S /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_2433_unphased/pakistan_gsa_included_1814_samples.tsv  -Oz -o pakistan_gsa_match_1814.pgt.vcf.gz  pakistan_gsa_match.pgt.vcf.gz

zcat pakistan_gsa_match_1814.pgt.vcf.gz|head -n209 # 209 lines of header
zcat pakistan_gsa_match_1814.pgt.vcf.gz|awk 'NR==209{print NF}' # 1823 columns (=1814+9)
zcat pakistan_gsa_match_1814.pgt.vcf.gz|wc -l # 527,145 lines (209+ 526,936 SNPs)
# 1814 participants with 526,936 SNPs

# create index for vcf.gz
bcftools index -f -t pakistan_gsa_match_1814.pgt.vcf.gz
# -f: overwrite index if it already exists
# -t: generate TBI-format index for VCF files

# break down into different chromosomes
bcftools view --no-version -Oz -o pakistan_gsa_match_1814.chr1.pgt.vcf.gz -r chr1 pakistan_gsa_match_1814.pgt.vcf.gz
bcftools view --no-version -Oz -o pakistan_gsa_match_1814.chr2.pgt.vcf.gz -r chr2 pakistan_gsa_match_1814.pgt.vcf.gz
bcftools view --no-version -Oz -o pakistan_gsa_match_1814.chr3.pgt.vcf.gz -r chr3 pakistan_gsa_match_1814.pgt.vcf.gz
bcftools view --no-version -Oz -o pakistan_gsa_match_1814.chr4.pgt.vcf.gz -r chr4 pakistan_gsa_match_1814.pgt.vcf.gz
bcftools view --no-version -Oz -o pakistan_gsa_match_1814.chr5.pgt.vcf.gz -r chr5 pakistan_gsa_match_1814.pgt.vcf.gz
bcftools view --no-version -Oz -o pakistan_gsa_match_1814.chr6.pgt.vcf.gz -r chr6 pakistan_gsa_match_1814.pgt.vcf.gz
bcftools view --no-version -Oz -o pakistan_gsa_match_1814.chr7.pgt.vcf.gz -r chr7 pakistan_gsa_match_1814.pgt.vcf.gz
bcftools view --no-version -Oz -o pakistan_gsa_match_1814.chr8.pgt.vcf.gz -r chr8 pakistan_gsa_match_1814.pgt.vcf.gz
bcftools view --no-version -Oz -o pakistan_gsa_match_1814.chr9.pgt.vcf.gz -r chr9 pakistan_gsa_match_1814.pgt.vcf.gz
bcftools view --no-version -Oz -o pakistan_gsa_match_1814.chr10.pgt.vcf.gz -r chr10 pakistan_gsa_match_1814.pgt.vcf.gz
bcftools view --no-version -Oz -o pakistan_gsa_match_1814.chr11.pgt.vcf.gz -r chr11 pakistan_gsa_match_1814.pgt.vcf.gz
bcftools view --no-version -Oz -o pakistan_gsa_match_1814.chr12.pgt.vcf.gz -r chr12 pakistan_gsa_match_1814.pgt.vcf.gz
bcftools view --no-version -Oz -o pakistan_gsa_match_1814.chr13.pgt.vcf.gz -r chr13 pakistan_gsa_match_1814.pgt.vcf.gz
bcftools view --no-version -Oz -o pakistan_gsa_match_1814.chr14.pgt.vcf.gz -r chr14 pakistan_gsa_match_1814.pgt.vcf.gz
bcftools view --no-version -Oz -o pakistan_gsa_match_1814.chr15.pgt.vcf.gz -r chr15 pakistan_gsa_match_1814.pgt.vcf.gz
bcftools view --no-version -Oz -o pakistan_gsa_match_1814.chr16.pgt.vcf.gz -r chr16 pakistan_gsa_match_1814.pgt.vcf.gz
bcftools view --no-version -Oz -o pakistan_gsa_match_1814.chr17.pgt.vcf.gz -r chr17 pakistan_gsa_match_1814.pgt.vcf.gz
bcftools view --no-version -Oz -o pakistan_gsa_match_1814.chr18.pgt.vcf.gz -r chr18 pakistan_gsa_match_1814.pgt.vcf.gz
bcftools view --no-version -Oz -o pakistan_gsa_match_1814.chr19.pgt.vcf.gz -r chr19 pakistan_gsa_match_1814.pgt.vcf.gz
bcftools view --no-version -Oz -o pakistan_gsa_match_1814.chr20.pgt.vcf.gz -r chr20 pakistan_gsa_match_1814.pgt.vcf.gz
bcftools view --no-version -Oz -o pakistan_gsa_match_1814.chr21.pgt.vcf.gz -r chr21 pakistan_gsa_match_1814.pgt.vcf.gz
bcftools view --no-version -Oz -o pakistan_gsa_match_1814.chr22.pgt.vcf.gz -r chr22 pakistan_gsa_match_1814.pgt.vcf.gz
bcftools view --no-version -Oz -o pakistan_gsa_match_1814.chr22.pgt.vcf.gz -r chrX pakistan_gsa_match_1814.pgt.vcf.gz # ChrX: 24116-209 = 23,907 SNPs
# all SNPs without ChrX =  526,936 - 23,907 =503,029 SNPs

# sort the vcf.gz files
bcftools sort  pakistan_gsa_match_1814.chr1.pgt.vcf.gz -Oz -o pakistan_gsa_match_1814.sort.chr1.pgt.vcf.gz
bcftools sort  pakistan_gsa_match_1814.chr2.pgt.vcf.gz -Oz -o pakistan_gsa_match_1814.sort.chr2.pgt.vcf.gz
bcftools sort  pakistan_gsa_match_1814.chr3.pgt.vcf.gz -Oz -o pakistan_gsa_match_1814.sort.chr3.pgt.vcf.gz
bcftools sort  pakistan_gsa_match_1814.chr4.pgt.vcf.gz -Oz -o pakistan_gsa_match_1814.sort.chr4.pgt.vcf.gz
bcftools sort  pakistan_gsa_match_1814.chr5.pgt.vcf.gz -Oz -o pakistan_gsa_match_1814.sort.chr5.pgt.vcf.gz
bcftools sort  pakistan_gsa_match_1814.chr6.pgt.vcf.gz -Oz -o pakistan_gsa_match_1814.sort.chr6.pgt.vcf.gz
bcftools sort  pakistan_gsa_match_1814.chr7.pgt.vcf.gz -Oz -o pakistan_gsa_match_1814.sort.chr7.pgt.vcf.gz
bcftools sort  pakistan_gsa_match_1814.chr8.pgt.vcf.gz -Oz -o pakistan_gsa_match_1814.sort.chr8.pgt.vcf.gz
bcftools sort  pakistan_gsa_match_1814.chr9.pgt.vcf.gz -Oz -o pakistan_gsa_match_1814.sort.chr9.pgt.vcf.gz
bcftools sort  pakistan_gsa_match_1814.chr10.pgt.vcf.gz -Oz -o pakistan_gsa_match_1814.sort.chr10.pgt.vcf.gz
bcftools sort  pakistan_gsa_match_1814.chr11.pgt.vcf.gz -Oz -o pakistan_gsa_match_1814.sort.chr11.pgt.vcf.gz
bcftools sort  pakistan_gsa_match_1814.chr12.pgt.vcf.gz -Oz -o pakistan_gsa_match_1814.sort.chr12.pgt.vcf.gz
bcftools sort  pakistan_gsa_match_1814.chr13.pgt.vcf.gz -Oz -o pakistan_gsa_match_1814.sort.chr13.pgt.vcf.gz
bcftools sort  pakistan_gsa_match_1814.chr14.pgt.vcf.gz -Oz -o pakistan_gsa_match_1814.sort.chr14.pgt.vcf.gz
bcftools sort  pakistan_gsa_match_1814.chr15.pgt.vcf.gz -Oz -o pakistan_gsa_match_1814.sort.chr15.pgt.vcf.gz
bcftools sort  pakistan_gsa_match_1814.chr16.pgt.vcf.gz -Oz -o pakistan_gsa_match_1814.sort.chr16.pgt.vcf.gz
bcftools sort  pakistan_gsa_match_1814.chr17.pgt.vcf.gz -Oz -o pakistan_gsa_match_1814.sort.chr17.pgt.vcf.gz
bcftools sort  pakistan_gsa_match_1814.chr18.pgt.vcf.gz -Oz -o pakistan_gsa_match_1814.sort.chr18.pgt.vcf.gz
bcftools sort  pakistan_gsa_match_1814.chr19.pgt.vcf.gz -Oz -o pakistan_gsa_match_1814.sort.chr19.pgt.vcf.gz
bcftools sort  pakistan_gsa_match_1814.chr20.pgt.vcf.gz -Oz -o pakistan_gsa_match_1814.sort.chr20.pgt.vcf.gz
bcftools sort  pakistan_gsa_match_1814.chr21.pgt.vcf.gz -Oz -o pakistan_gsa_match_1814.sort.chr21.pgt.vcf.gz
bcftools sort  pakistan_gsa_match_1814.chr22.pgt.vcf.gz -Oz -o pakistan_gsa_match_1814.sort.chr22.pgt.vcf.gz


gsutil cp /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_phased_mocha/pakistan_gsa_match_1814.sort.*.pgt.vcf.gz gs://bigdeli-working/pakistan/cromwell/outputs/vcf_for_topmed/gsa_1814/


# 9. imputation on TOPMed and Michigan website

# Expand 1000G options (minimac4 1.7.3 imputation): 1000G phase 3 30X (grch38/hg38, BETA), GRCh38, rsq filter off, no phasing, population: ALL (to get the QC report), quality control and imputation, generate meta-imputation file (job name: gsa_1814_phased_1000G)
# Genome Asia options (minimac4 1.7.3 imputation): Genome Asia Pilot - GAsP (GRCH37/hg19), GRCh38, rsq filter off, no phasing, population: ASN (to get the QC report, there is no SAS option), quality control and imputation, generate meta-imputation file (job name: gsa_1814_phased_GA)
# TOPMed options (minimac4 1.7.3 imputation): TOPMed r2, GRCh38, rsq filter off, no phasing, population: versus TOPMed panel (to get the QC report), quality control and imputation, generate meta-imputation file (job name: gsa_1814_phased_topmed)
# 1000G (grch38) options (minimac4 1.7.3 imputation): 1000G phase 3 (grch38/hg38), GRCh38, rsq filter off, no phasing, population: ALL (to get the QC report), quality control and imputation, generate meta-imputation file (job name: gsa_1814_phased_1000G_OG)
# 1000G (grch37) options (minimac4 1.7.3 imputation): 1000G phase 3 (grch37/hg19), GRCh38, rsq filter off, no phasing, population: SAS (to get the QC report), quality control and imputation, generate meta-imputation file (job name: gsa_1814_phased_1000G_OG_grch37)


# TOPMed: gsa_1814_phased_topmed (password: Dkb-xKMujHx6G5 ) (excluded SNPs: 11609 , typed-only SNPs: 4165) (some SNPs in typed-only.txt (16,398 lines) are shown as Genotyped in the INFO files (e.g., chr9:138121868:G:A))
# Expanded 1000G: gsa_1814_phased_1000G (password: 8nZzKiLbh9erXF ) (excluded SNPs: 14041, typed-only SNPs: 4594) (typed-only.txt has 4489 lines)
# GenomeAsia: gsa_1814_phased_GA (password: mQRcb0KkvS8$M4 )  (excluded SNPs: 11522, typed-only SNPs: 30072)
# 1000G grch38: gsa_1814_phased_1000G_OG (password: TazSkFN2Vhpui7 )(excluded SNPs: 12245, typed-only SNPs: 2949)
# 1000G grch37: gsa_1814_phased_1000G_OG_grch37 (password: D59UmeXC73Pgtk ) (excluded SNPs: 11577, typed-only SNPs: 200)

#### Imputation DONE! #####

## inflate all the zip from Michigan/TopMed server on minerva
for z in *.zip; do unzip $z; done
 
# 53 participants with genotyping data do not have sequencing data


#### 9. meta-imputation #####

#Please do load python before cmake.
ml git python
ml cmake
git clone https://github.com/yukt/MetaMinimac2.git
cd MetaMinimac2
bash install.sh


#!/bin/bash
#BSUB -J meta_imputation_chr1_test # Job name
#BSUB -P acc_psychgen # allocation account
#BSUB -q premium # queue
#BSUB -n 2 # number of compute cores
#BSUB -W 48:00 # walltime in HH:MM
#BSUB -R rusage[mem=64000] # 80 GB of memory (8GB per core)
#BSUB -R span[hosts=1] # all cores from one node
#BSUB -o %J.stdout # output log (%J : JobID)
#BSUB -eo %J.stderr # error log
#BSUB -L /bin/bash # Initialize the execution environment

ml git python
ml cmake

/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/MetaMinimac2/release-build/MetaMinimac2 -i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/chr1:/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/chr1 -o chr1_topmed_1000g.meta.testrun

bsub < meta_imputation_chr1_test.sh

## for the remaining 22 chromosomes

#!/bin/bash
#BSUB -J meta_imputation_topmed_1000g # Job name
#BSUB -P acc_psychgen # allocation account
#BSUB -q premium # queue
#BSUB -n 1 # number of compute cores
#BSUB -W 120:00 # walltime in HH:MM
#BSUB -R himem
#BSUB -R rusage[mem=120000] # 80 GB of memory (8GB per core)
#BSUB -R span[hosts=1] # all cores from one node
#BSUB -o %J.stdout # output log (%J : JobID)
#BSUB -eo %J.stderr # error log
#BSUB -L /bin/bash # Initialize the execution environment

ml git python
ml cmake

/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/MetaMinimac2/release-build/MetaMinimac2 -i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/chr2:/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/chr2 -o chr2_topmed_1000g.meta
/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/MetaMinimac2/release-build/MetaMinimac2 -i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/chr3:/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/chr3 -o chr3_topmed_1000g.meta
/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/MetaMinimac2/release-build/MetaMinimac2 -i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/chr4:/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/chr4 -o chr4_topmed_1000g.meta
/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/MetaMinimac2/release-build/MetaMinimac2 -i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/chr5:/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/chr5 -o chr5_topmed_1000g.meta
/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/MetaMinimac2/release-build/MetaMinimac2 -i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/chr6:/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/chr6 -o chr6_topmed_1000g.meta
/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/MetaMinimac2/release-build/MetaMinimac2 -i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/chr7:/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/chr7 -o chr7_topmed_1000g.meta
/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/MetaMinimac2/release-build/MetaMinimac2 -i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/chr8:/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/chr8 -o chr8_topmed_1000g.meta
/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/MetaMinimac2/release-build/MetaMinimac2 -i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/chr9:/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/chr9 -o chr9_topmed_1000g.meta
/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/MetaMinimac2/release-build/MetaMinimac2 -i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/chr10:/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/chr10 -o chr10_topmed_1000g.meta
/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/MetaMinimac2/release-build/MetaMinimac2 -i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/chr11:/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/chr11 -o chr11_topmed_1000g.meta
/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/MetaMinimac2/release-build/MetaMinimac2 -i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/chr12:/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/chr12 -o chr12_topmed_1000g.meta
/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/MetaMinimac2/release-build/MetaMinimac2 -i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/chr13:/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/chr13 -o chr13_topmed_1000g.meta
/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/MetaMinimac2/release-build/MetaMinimac2 -i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/chr14:/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/chr14 -o chr14_topmed_1000g.meta
/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/MetaMinimac2/release-build/MetaMinimac2 -i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/chr15:/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/chr15 -o chr15_topmed_1000g.meta
/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/MetaMinimac2/release-build/MetaMinimac2 -i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/chr16:/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/chr16 -o chr16_topmed_1000g.meta
/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/MetaMinimac2/release-build/MetaMinimac2 -i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/chr17:/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/chr17 -o chr17_topmed_1000g.meta
/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/MetaMinimac2/release-build/MetaMinimac2 -i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/chr18:/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/chr18 -o chr18_topmed_1000g.meta
/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/MetaMinimac2/release-build/MetaMinimac2 -i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/chr19:/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/chr19 -o chr19_topmed_1000g.meta
/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/MetaMinimac2/release-build/MetaMinimac2 -i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/chr20:/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/chr20 -o chr20_topmed_1000g.meta
/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/MetaMinimac2/release-build/MetaMinimac2 -i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/chr21:/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/chr21 -o chr21_topmed_1000g.meta
/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/MetaMinimac2/release-build/MetaMinimac2 -i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/chr22:/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/chr22 -o chr22_topmed_1000g.meta

bsub < meta_imputation_topmed_1000g.sh


##### 10. checking the imputation INFO files

for z in *.info.gz; do gunzip $z; done

grep '^SNP' chr1.info > all.chr.info
grep -v '^SNP' chr1.info  chr2.info  chr3.info chr4.info chr5.info chr6.info chr7.info chr8.info chr9.info chr10.info chr11.info chr12.info chr13.info chr14.info chr15.info chr16.info chr17.info chr18.info chr19.info chr20.info chr21.info chr22.info >> all.chr.info

# check Rsq and MAF columns
library(data.table)
info <- fread("all.chr.info",header=T)
table(info$Genotyped)
# GA:
# Genotyped    Imputed Typed_Only 
#    461435   21033213      30072 
# 1000G:
# Genotyped    Imputed Typed_Only 
#    484357   48406882       4594 
# TOPMed:
# Genotyped    Imputed Typed_Only 
#    487217  291649245       4165 
# 1000G GRCh38 original
# Genotyped    Imputed Typed_Only 
#   487835    44719577       2949 
# 1000G GRCh37 original
#  Genotyped    Imputed Typed_Only 
#    491089   46608462        201 

summary(as.numeric(info$Rsq))
# GA: ranges from 0 to 1, 30072 NAs
# 1000G: ranges from 0 to 1, 4594 NAs
# TOPMed: ranges from 0 to 1, 4165 NAs
# 1000G GRCh38 original: ranges from 0 to 1,  2949 NAs
# 1000G GRCh37 original: ranges from 0 to 1,  201 NAs

summary(as.numeric(info$AvgCall))
# GA: ranges from 0.506 to 1, 30072 NAs
# 1000G: ranges from 0.572 to 1, 4594 NAs
# TOPMed: ranges from 0.530 to 1, 4165 NAs
# 1000G GRCh38 original: ranges from 0.5172 to 1, 2949 NAs
# 1000G GRCh37 original: ranges from 0.5006 to 1, 201 NAs

summary(as.numeric(info$MAF))
# GA: ranged from 0 to 0.5
# 1000G: ranged from 0 to 0.5
# TOPMed: ranged from 0 to 0.5
# 1000G GRCh38 original: ranges from 0 to 0.5
# 1000G GRCh37 original: ranges from 0 to 0.5

summary(as.numeric(info$ALT_Frq))
# GA: ranged from 0 to 1
# 1000G: ranged from 0 to 1
# TOPMed: ranged from 0 to 1
# 1000G GRCh38 original: ranges from 0 to 1
# 1000G GRCh37 original: ranges from 0 to 1

## number of imputed SNPs and well-imputed SNPs (Rsq >= 0.8)
dim(info[which(info$Genotyped=="Imputed"),])
dim(info[which(info$Genotyped=="Imputed" & info$Rsq >= 0.8),])

# 1000G GRCh38 original: 44719577 imputed SNPs, 9880149 well-imputed SNPs
# 1000G GRCh37 original: 46608462 imputed SNPs, 8013431 well-imputed SNPs

dim(info[which(info$Genotyped=="Genotyped"),])
dim(info[which(info$Genotyped=="Genotyped" & info$Rsq >= 0.8),])
dim(info[which(info$Genotyped=="Genotyped" & info$EmpRsq >= 0.8),])



###########################################################
####                 R script: R2 by MAF               ####
###########################################################

#!/bin/bash
#BSUB -J r2_by_maf_GA # Job name
#BSUB -n 1
#BSUB -P acc_psychgen # allocation account
#BSUB -q premium # queue
#BSUB -W 48:00 # walltime in HH:MM
#BSUB -R himem
#BSUB -R rusage[mem=96G] # 96 GB of memory 
#BSUB -R span[hosts=1] # all cores from one node
#BSUB -o %J.stdout # output log (%J : JobID)
#BSUB -eo %J.stderr # error log
#BSUB -L /bin/bash # Initialize the execution environment


module load R

Rscript r2_by_maf.r

#!/usr/bin/env Rscript

options(echo=TRUE)
print("In R now!")

## SCRIPT STARTS BELOW

library(data.table)
info <- fread("all.chr.info",header=T)
info$Rsq <- as.character(info$Rsq)
info$Rsq <- as.numeric(info$Rsq)

#explain how Rsq is calculated: https://genome.sph.umich.edu/wiki/Minimac3_Info_File

imputed <- info[which(info$Genotyped=="Imputed"),]

dim(info)
dim(imputed)

# calculate R2 by MAF
mean_r2 <- array(NA, dim=c(5000,2))
colnames(mean_r2) <- c("MAF","Rsq")

for (i in 1:5000) {
	if (i==1 ) { 
   dat <- subset(imputed, MAF >= (i/10000 - 0.0001) & MAF <= i/10000)
   mean_r2[i,1] <- i/10000
   mean_r2[i,2] <- mean(dat$Rsq)
	 }
	else  {
   dat <- subset(imputed, MAF > (i/10000 - 0.0001) & MAF <= i/10000)
   mean_r2[i,1] <- i/10000
   mean_r2[i,2] <- mean(dat$Rsq)
   }
}

write.table(mean_r2,'mean_r2',sep='\t', row.names=FALSE,quote=FALSE)

bsub < r2_by_maf.sh

####################################################################
####         R script: R2 by MAF (using bigger MAF bins)        ####
####################################################################

# what about calculating MAF per MAF bins (0, 0.0005, 0.001, 0.002, 0.005, 0.010, 0.015, 0.020, 0.035, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)

#!/bin/bash
#BSUB -J r2_by_maf_new # Job name
#BSUB -n 1
#BSUB -P acc_psychgen # allocation account
#BSUB -q premium # queue
#BSUB -W 48:00 # walltime in HH:MM
#BSUB -R himem
#BSUB -R rusage[mem=96G] # 120 GB of memory 
#BSUB -R span[hosts=1] # all cores from one node
#BSUB -o %J.stdout # output log (%J : JobID)
#BSUB -eo %J.stderr # error log
#BSUB -L /bin/bash # Initialize the execution environment


module load R

Rscript r2_by_maf_new.r

#!/usr/bin/env Rscript

options(echo=TRUE)
print("In R now!")

## SCRIPT STARTS BELOW

library(data.table)
info <- fread("all.chr.info",header=T)
info$Rsq <- as.character(info$Rsq)
info$Rsq <- as.numeric(info$Rsq)

#explain how Rsq is calculated: https://genome.sph.umich.edu/wiki/Minimac3_Info_File

imputed <- info[which(info$Genotyped!="Typed_Only"),]

dim(info)
dim(imputed)

# calculate R2 by MAF
mean_r2 <- array(NA, dim=c(14,2))
colnames(mean_r2) <- c("MAF_bin","Rsq")

for (i in 1:14) {
	if (i==1 ) { 
   dat <- subset(imputed, MAF >= 0 & MAF <= 0.0005)
   mean_r2[i,1] <- "[0,0.0005]"
   mean_r2[i,2] <- mean(dat$Rsq)
	 }
	else if (i==2 ) { 
   dat <- subset(imputed, MAF > 0.0005 & MAF <= 0.001)
   mean_r2[i,1] <- "(0.0005,0.001]"
   mean_r2[i,2] <- mean(dat$Rsq)
	 }
	else if (i==3 ) { 
   dat <- subset(imputed, MAF > 0.001 & MAF <= 0.002)
   mean_r2[i,1] <- "(0.001,0.002]"
   mean_r2[i,2] <- mean(dat$Rsq)
	 }
	else if (i==4 ) { 
   dat <- subset(imputed, MAF > 0.002 & MAF <= 0.005)
   mean_r2[i,1] <- "(0.002,0.005]"
   mean_r2[i,2] <- mean(dat$Rsq)
	 }
	else if (i==5 ) { 
   dat <- subset(imputed, MAF > 0.005 & MAF <= 0.01)
   mean_r2[i,1] <- "(0.005,0.01]"
   mean_r2[i,2] <- mean(dat$Rsq)
	 }
	else if (i==6 ) { 
   dat <- subset(imputed, MAF > 0.01 & MAF <= 0.015)
   mean_r2[i,1] <- "(0.01,0.015]"
   mean_r2[i,2] <- mean(dat$Rsq)
	 }
	else if (i==7 ) { 
   dat <- subset(imputed, MAF > 0.015 & MAF <= 0.02)
   mean_r2[i,1] <- "(0.015,0.02]"
   mean_r2[i,2] <- mean(dat$Rsq)
	 }
	else if (i==8 ) { 
   dat <- subset(imputed, MAF > 0.02 & MAF <= 0.035)
   mean_r2[i,1] <- "(0.02,0.035]"
   mean_r2[i,2] <- mean(dat$Rsq)
	 }
	else if (i==9 ) { 
   dat <- subset(imputed, MAF > 0.035 & MAF <= 0.05)
   mean_r2[i,1] <- "(0.035,0.05]"
   mean_r2[i,2] <- mean(dat$Rsq)
	 }
	else if (i==10 ) { 
   dat <- subset(imputed, MAF > 0.05 & MAF <= 0.1)
   mean_r2[i,1] <- "(0.05,0.1]"
   mean_r2[i,2] <- mean(dat$Rsq)
	 }
	else if (i==11 ) { 
   dat <- subset(imputed, MAF > 0.1 & MAF <= 0.2)
   mean_r2[i,1] <- "(0.1,0.2]"
   mean_r2[i,2] <- mean(dat$Rsq)
	 }
	else if (i==12 ) { 
   dat <- subset(imputed, MAF > 0.2 & MAF <= 0.3)
   mean_r2[i,1] <- "(0.2,0.3]"
   mean_r2[i,2] <- mean(dat$Rsq)
	 }
	else if (i==13 ) { 
   dat <- subset(imputed, MAF > 0.3 & MAF <= 0.4)
   mean_r2[i,1] <- "(0.3,0.4]"
   mean_r2[i,2] <- mean(dat$Rsq)
	 }
	else if (i==14 ) { 
   dat <- subset(imputed, MAF > 0.4 & MAF <= 0.5)
   mean_r2[i,1] <- "(0.4,0.5]"
   mean_r2[i,2] <- mean(dat$Rsq)
	 }
}

mean_r2 <- as.data.frame(mean_r2)

write.table(mean_r2,'mean_r2_new_MAF_bin',sep='\t', row.names=FALSE,quote=FALSE)

bsub < r2_by_maf_new.sh

######## for meta-imputation #############
#!/usr/bin/env Rscript

options(echo=TRUE)
print("In R now!")

## SCRIPT STARTS BELOW

library(data.table)
imputed <- fread("all.chr.info",header=T)
imputed$Rsq <- as.character(imputed$Rsq)
imputed$Rsq <- as.numeric(imputed$Rsq)

# calculate R2 by MAF
mean_r2 <- array(NA, dim=c(14,2))
colnames(mean_r2) <- c("MAF_bin","Rsq")

for (i in 1:14) {
	if (i==1 ) { 
   dat <- subset(imputed, MAF >= 0 & MAF <= 0.0005)
   mean_r2[i,1] <- "[0,0.0005]"
   mean_r2[i,2] <- mean(dat$Rsq)
	 }
	else if (i==2 ) { 
   dat <- subset(imputed, MAF > 0.0005 & MAF <= 0.001)
   mean_r2[i,1] <- "(0.0005,0.001]"
   mean_r2[i,2] <- mean(dat$Rsq)
	 }
	else if (i==3 ) { 
   dat <- subset(imputed, MAF > 0.001 & MAF <= 0.002)
   mean_r2[i,1] <- "(0.001,0.002]"
   mean_r2[i,2] <- mean(dat$Rsq)
	 }
	else if (i==4 ) { 
   dat <- subset(imputed, MAF > 0.002 & MAF <= 0.005)
   mean_r2[i,1] <- "(0.002,0.005]"
   mean_r2[i,2] <- mean(dat$Rsq)
	 }
	else if (i==5 ) { 
   dat <- subset(imputed, MAF > 0.005 & MAF <= 0.01)
   mean_r2[i,1] <- "(0.005,0.01]"
   mean_r2[i,2] <- mean(dat$Rsq)
	 }
	else if (i==6 ) { 
   dat <- subset(imputed, MAF > 0.01 & MAF <= 0.015)
   mean_r2[i,1] <- "(0.01,0.015]"
   mean_r2[i,2] <- mean(dat$Rsq)
	 }
	else if (i==7 ) { 
   dat <- subset(imputed, MAF > 0.015 & MAF <= 0.02)
   mean_r2[i,1] <- "(0.015,0.02]"
   mean_r2[i,2] <- mean(dat$Rsq)
	 }
	else if (i==8 ) { 
   dat <- subset(imputed, MAF > 0.02 & MAF <= 0.035)
   mean_r2[i,1] <- "(0.02,0.035]"
   mean_r2[i,2] <- mean(dat$Rsq)
	 }
	else if (i==9 ) { 
   dat <- subset(imputed, MAF > 0.035 & MAF <= 0.05)
   mean_r2[i,1] <- "(0.035,0.05]"
   mean_r2[i,2] <- mean(dat$Rsq)
	 }
	else if (i==10 ) { 
   dat <- subset(imputed, MAF > 0.05 & MAF <= 0.1)
   mean_r2[i,1] <- "(0.05,0.1]"
   mean_r2[i,2] <- mean(dat$Rsq)
	 }
	else if (i==11 ) { 
   dat <- subset(imputed, MAF > 0.1 & MAF <= 0.2)
   mean_r2[i,1] <- "(0.1,0.2]"
   mean_r2[i,2] <- mean(dat$Rsq)
	 }
	else if (i==12 ) { 
   dat <- subset(imputed, MAF > 0.2 & MAF <= 0.3)
   mean_r2[i,1] <- "(0.2,0.3]"
   mean_r2[i,2] <- mean(dat$Rsq)
	 }
	else if (i==13 ) { 
   dat <- subset(imputed, MAF > 0.3 & MAF <= 0.4)
   mean_r2[i,1] <- "(0.3,0.4]"
   mean_r2[i,2] <- mean(dat$Rsq)
	 }
	else if (i==14 ) { 
   dat <- subset(imputed, MAF > 0.4 & MAF <= 0.5)
   mean_r2[i,1] <- "(0.4,0.5]"
   mean_r2[i,2] <- mean(dat$Rsq)
	 }
}

mean_r2 <- as.data.frame(mean_r2)

write.table(mean_r2,'mean_r2_new_MAF_bin',sep='\t', row.names=FALSE,quote=FALSE)


#####################################################################################
####                 R script: R2 by MAF  (for meta-imputation only)             ####
#####################################################################################

#!/bin/bash
#BSUB -J info_meta_imputation # Job name
#BSUB -n 1
#BSUB -P acc_psychgen # allocation account
#BSUB -q premium # queue
#BSUB -W 48:00 # walltime in HH:MM
#BSUB -R himem
#BSUB -R rusage[mem=150G] # 96 GB of memory 
#BSUB -R span[hosts=1] # all cores from one node
#BSUB -o %J.stdout # output log (%J : JobID)
#BSUB -eo %J.stderr # error log
#BSUB -L /bin/bash # Initialize the execution environment

ml bcftools
bcftools annotate -x ^INFO/MAF,INFO/R2 chr1_topmed_1000g.meta.metaDose.vcf.gz |grep -v '^#' |awk '{print $3"\t"$8}' |sed 's/\;/\t/' |sed 's/MAF=//' |sed 's/R2=//' |sed '1i SNP\tMAF\tR2'  > chr1_topmed_1000g.info 
bcftools annotate -x ^INFO/MAF,INFO/R2 chr2_topmed_1000g.meta.metaDose.vcf.gz |grep -v '^#' |awk '{print $3"\t"$8}' |sed 's/\;/\t/' |sed 's/MAF=//' |sed 's/R2=//' |sed '1i SNP\tMAF\tR2'  > chr2_topmed_1000g.info 
bcftools annotate -x ^INFO/MAF,INFO/R2 chr3_topmed_1000g.meta.metaDose.vcf.gz |grep -v '^#' |awk '{print $3"\t"$8}' |sed 's/\;/\t/' |sed 's/MAF=//' |sed 's/R2=//' |sed '1i SNP\tMAF\tR2'  > chr3_topmed_1000g.info 
bcftools annotate -x ^INFO/MAF,INFO/R2 chr4_topmed_1000g.meta.metaDose.vcf.gz |grep -v '^#' |awk '{print $3"\t"$8}' |sed 's/\;/\t/' |sed 's/MAF=//' |sed 's/R2=//' |sed '1i SNP\tMAF\tR2'  > chr4_topmed_1000g.info 
bcftools annotate -x ^INFO/MAF,INFO/R2 chr5_topmed_1000g.meta.metaDose.vcf.gz |grep -v '^#' |awk '{print $3"\t"$8}' |sed 's/\;/\t/' |sed 's/MAF=//' |sed 's/R2=//' |sed '1i SNP\tMAF\tR2'  > chr5_topmed_1000g.info 
bcftools annotate -x ^INFO/MAF,INFO/R2 chr6_topmed_1000g.meta.metaDose.vcf.gz |grep -v '^#' |awk '{print $3"\t"$8}' |sed 's/\;/\t/' |sed 's/MAF=//' |sed 's/R2=//' |sed '1i SNP\tMAF\tR2'  > chr6_topmed_1000g.info 
bcftools annotate -x ^INFO/MAF,INFO/R2 chr7_topmed_1000g.meta.metaDose.vcf.gz |grep -v '^#' |awk '{print $3"\t"$8}' |sed 's/\;/\t/' |sed 's/MAF=//' |sed 's/R2=//' |sed '1i SNP\tMAF\tR2'  > chr7_topmed_1000g.info 
bcftools annotate -x ^INFO/MAF,INFO/R2 chr8_topmed_1000g.meta.metaDose.vcf.gz |grep -v '^#' |awk '{print $3"\t"$8}' |sed 's/\;/\t/' |sed 's/MAF=//' |sed 's/R2=//' |sed '1i SNP\tMAF\tR2'  > chr8_topmed_1000g.info 
bcftools annotate -x ^INFO/MAF,INFO/R2 chr9_topmed_1000g.meta.metaDose.vcf.gz |grep -v '^#' |awk '{print $3"\t"$8}' |sed 's/\;/\t/' |sed 's/MAF=//' |sed 's/R2=//' |sed '1i SNP\tMAF\tR2'  > chr9_topmed_1000g.info 
bcftools annotate -x ^INFO/MAF,INFO/R2 chr10_topmed_1000g.meta.metaDose.vcf.gz |grep -v '^#' |awk '{print $3"\t"$8}' |sed 's/\;/\t/' |sed 's/MAF=//' |sed 's/R2=//' |sed '1i SNP\tMAF\tR2'  > chr10_topmed_1000g.info 
bcftools annotate -x ^INFO/MAF,INFO/R2 chr11_topmed_1000g.meta.metaDose.vcf.gz |grep -v '^#' |awk '{print $3"\t"$8}' |sed 's/\;/\t/' |sed 's/MAF=//' |sed 's/R2=//' |sed '1i SNP\tMAF\tR2'  > chr11_topmed_1000g.info 
bcftools annotate -x ^INFO/MAF,INFO/R2 chr12_topmed_1000g.meta.metaDose.vcf.gz |grep -v '^#' |awk '{print $3"\t"$8}' |sed 's/\;/\t/' |sed 's/MAF=//' |sed 's/R2=//' |sed '1i SNP\tMAF\tR2'  > chr12_topmed_1000g.info 
bcftools annotate -x ^INFO/MAF,INFO/R2 chr13_topmed_1000g.meta.metaDose.vcf.gz |grep -v '^#' |awk '{print $3"\t"$8}' |sed 's/\;/\t/' |sed 's/MAF=//' |sed 's/R2=//' |sed '1i SNP\tMAF\tR2'  > chr13_topmed_1000g.info 
bcftools annotate -x ^INFO/MAF,INFO/R2 chr14_topmed_1000g.meta.metaDose.vcf.gz |grep -v '^#' |awk '{print $3"\t"$8}' |sed 's/\;/\t/' |sed 's/MAF=//' |sed 's/R2=//' |sed '1i SNP\tMAF\tR2'  > chr14_topmed_1000g.info 
bcftools annotate -x ^INFO/MAF,INFO/R2 chr15_topmed_1000g.meta.metaDose.vcf.gz |grep -v '^#' |awk '{print $3"\t"$8}' |sed 's/\;/\t/' |sed 's/MAF=//' |sed 's/R2=//' |sed '1i SNP\tMAF\tR2'  > chr15_topmed_1000g.info 
bcftools annotate -x ^INFO/MAF,INFO/R2 chr16_topmed_1000g.meta.metaDose.vcf.gz |grep -v '^#' |awk '{print $3"\t"$8}' |sed 's/\;/\t/' |sed 's/MAF=//' |sed 's/R2=//' |sed '1i SNP\tMAF\tR2'  > chr16_topmed_1000g.info 
bcftools annotate -x ^INFO/MAF,INFO/R2 chr17_topmed_1000g.meta.metaDose.vcf.gz |grep -v '^#' |awk '{print $3"\t"$8}' |sed 's/\;/\t/' |sed 's/MAF=//' |sed 's/R2=//' |sed '1i SNP\tMAF\tR2'  > chr17_topmed_1000g.info 
bcftools annotate -x ^INFO/MAF,INFO/R2 chr18_topmed_1000g.meta.metaDose.vcf.gz |grep -v '^#' |awk '{print $3"\t"$8}' |sed 's/\;/\t/' |sed 's/MAF=//' |sed 's/R2=//' |sed '1i SNP\tMAF\tR2'  > chr18_topmed_1000g.info 
bcftools annotate -x ^INFO/MAF,INFO/R2 chr19_topmed_1000g.meta.metaDose.vcf.gz |grep -v '^#' |awk '{print $3"\t"$8}' |sed 's/\;/\t/' |sed 's/MAF=//' |sed 's/R2=//' |sed '1i SNP\tMAF\tR2'  > chr19_topmed_1000g.info 
bcftools annotate -x ^INFO/MAF,INFO/R2 chr20_topmed_1000g.meta.metaDose.vcf.gz |grep -v '^#' |awk '{print $3"\t"$8}' |sed 's/\;/\t/' |sed 's/MAF=//' |sed 's/R2=//' |sed '1i SNP\tMAF\tR2'  > chr20_topmed_1000g.info 
bcftools annotate -x ^INFO/MAF,INFO/R2 chr21_topmed_1000g.meta.metaDose.vcf.gz |grep -v '^#' |awk '{print $3"\t"$8}' |sed 's/\;/\t/' |sed 's/MAF=//' |sed 's/R2=//' |sed '1i SNP\tMAF\tR2'  > chr21_topmed_1000g.info 
bcftools annotate -x ^INFO/MAF,INFO/R2 chr22_topmed_1000g.meta.metaDose.vcf.gz |grep -v '^#' |awk '{print $3"\t"$8}' |sed 's/\;/\t/' |sed 's/MAF=//' |sed 's/R2=//' |sed '1i SNP\tMAF\tR2'  > chr22_topmed_1000g.info 

cat chr*_topmed_1000g.info  > all.topmed_1000g.info 

bsub < info_meta_imputation.sh

grep -v SNP all.topmed_1000g.info |sed '1i SNP\tMAF\tRsq' > all.chr.info # 296353413
# there are 22 header lines of "SNP	MAF	R2", need to remove 21 lines

#!/bin/bash
#BSUB -J r2_by_maf_meta_imputation # Job name
#BSUB -n 1
#BSUB -P acc_psychgen # allocation account
#BSUB -q premium # queue
#BSUB -W 48:00 # walltime in HH:MM
#BSUB -R himem
#BSUB -R rusage[mem=96G] # 96 GB of memory 
#BSUB -R span[hosts=1] # all cores from one node
#BSUB -o %J.stdout # output log (%J : JobID)
#BSUB -eo %J.stderr # error log
#BSUB -L /bin/bash # Initialize the execution environment


module load R

Rscript r2_by_maf.r

#!/usr/bin/env Rscript

options(echo=TRUE)
print("In R now!")

## SCRIPT STARTS BELOW

library(data.table)
info <- fread("all.chr.info",header=T)

#explain how Rsq is calculated: https://genome.sph.umich.edu/wiki/Minimac3_Info_File

# calculate R2 by MAF
mean_r2 <- array(NA, dim=c(5000,2))
colnames(mean_r2) <- c("MAF","Rsq")

for (i in 1:5000) {
	if (i==1 ) { 
   dat <- subset(info, MAF >= (i/10000 - 0.0001) & MAF <= i/10000)
   mean_r2[i,1] <- i/10000
   mean_r2[i,2] <- mean(info$Rsq)
	 }
	else  {
   dat <- subset(info, MAF > (i/10000 - 0.0001) & MAF <= i/10000)
   mean_r2[i,1] <- i/10000
   mean_r2[i,2] <- mean(info$Rsq)
   }
}

write.table(mean_r2,'mean_r2',sep='\t', row.names=FALSE,quote=FALSE)

bsub < r2_by_maf.sh

#################################################################################################################
####               Supplemental Table 2: count numbers of imputed SNPs and number of well-imputed SNPs       ####
#################################################################################################################


#############################################################################################################
####               testing: R2 by MAF in chr1 for meta-imputed results (TOPMed and expanded 1000G)       ####
#############################################################################################################

zcat ../chr1_topmed_1000g.meta.testrun.metaDose.vcf.gz|awk '(NR>16){print $3,$8}' |sed 's/\;/ /g'|sed 's/R2=//1'|sed 's/NST=//'|sed 's/AF=//'| sed 's/MAF=//'|sed '1i ID NST S2 AF MAF R2' > chr1_topmed_1000g.meta.info

# SNPs in only 1 reference panel
# chr1:10481:A:G 1 S2 0.00001 0.00001 0.00161

# SNPs in 2 reference panels
# chr1:13917:T:TC 2 S1 S2 0.00008 0.00008 0.00842

# SNPs used as training (all of them are in both reference panels)
# chr1:858952:G:A 2 S1 S2 TRAINING 0.09773 0.09773 0.97241

awk '(NF==6){print $1,$4,$5,$6}'  chr1_topmed_1000g.meta.info | awk '{print $0,"Imputed"}' > chr1_topmed_1000g.meta.1panel.info # 20,508,689 lines
awk '(NF==7){print $1,$5,$6,$7}'  chr1_topmed_1000g.meta.info | awk '{print $0,"Imputed"}' >  chr1_topmed_1000g.meta.2panel.info # 3,523,710 lines 
awk '(NF==8){print $1,$6,$7,$8}'  chr1_topmed_1000g.meta.info | awk '{print $0,"Genotyped"}' >  chr1_topmed_1000g.meta.2panel.training.info # 38,483 lines (genotyped SNPs)
# 38712 genotyped SNPs on chr1 in TOPMed, 38515 genotyped SNPs on chr1 in expanded 1000G

cat chr1_topmed_1000g.meta.1panel.info chr1_topmed_1000g.meta.2panel.info  chr1_topmed_1000g.meta.2panel.training.info  > chr1_topmed_1000g.meta.final.info 
# 24,070,882 lines 
awk '{print NF}' chr1_topmed_1000g.meta.final.info |sort |uniq -c
# each line has 5 columns

### meta-imputed chr1

library(data.table)
info <- fread("chr1_topmed_1000g.meta.final.info",header=T)
info$R2 <- as.numeric(info$R2)
colnames(info)[colnames(info)=="R2"] <- "Rsq"
colnames(info)[colnames(info)=="Imputed"] <- "Genotyped"

#explain how Rsq is calculated: https://genome.sph.umich.edu/wiki/Minimac3_Info_File

imputed <- info[which(info$Genotyped=="Imputed"),]

dim(info)
dim(imputed)

# calculate R2 by MAF
mean_r2 <- array(NA, dim=c(5000,2))
colnames(mean_r2) <- c("MAF","Rsq")

for (i in 1:5000) {
	if (i==1 ) { 
   dat <- subset(imputed, MAF >= (i/10000 - 0.0001) & MAF <= i/10000)
   mean_r2[i,1] <- i/10000
   mean_r2[i,2] <- mean(dat$Rsq)
	 }
	else  {
   dat <- subset(imputed, MAF > (i/10000 - 0.0001) & MAF <= i/10000)
   mean_r2[i,1] <- i/10000
   mean_r2[i,2] <- mean(dat$Rsq)
   }
}

write.table(mean_r2,'mean_r2_chr1',sep='\t', row.names=FALSE,quote=FALSE)

### TOPMed and expanded 1000G chr1

library(data.table)
info <- fread("chr1.info",header=T)
info$Rsq <- as.numeric(info$Rsq)

#explain how Rsq is calculated: https://genome.sph.umich.edu/wiki/Minimac3_Info_File

imputed <- info[which(info$Genotyped=="Imputed"),]

dim(info)
dim(imputed)

# calculate R2 by MAF
mean_r2 <- array(NA, dim=c(5000,2))
colnames(mean_r2) <- c("MAF","Rsq")

for (i in 1:5000) {
	if (i==1 ) { 
   dat <- subset(imputed, MAF >= (i/10000 - 0.0001) & MAF <= i/10000)
   mean_r2[i,1] <- i/10000
   mean_r2[i,2] <- mean(dat$Rsq)
	 }
	else  {
   dat <- subset(imputed, MAF > (i/10000 - 0.0001) & MAF <= i/10000)
   mean_r2[i,1] <- i/10000
   mean_r2[i,2] <- mean(dat$Rsq)
   }
}

write.table(mean_r2,'mean_r2_chr1',sep='\t', row.names=FALSE,quote=FALSE)


################################################################
####         R script: empirical R2 by MAF (masking)        ####
################################################################

#!/bin/bash
#BSUB -J emp_r2_by_maf_GA # Job name
#BSUB -n 1
#BSUB -P acc_psychgen # allocation account
#BSUB -q premium # queue
#BSUB -W 48:00 # walltime in HH:MM
#BSUB -R himem
#BSUB -R rusage[mem=96G] # 120 GB of memory 
#BSUB -R span[hosts=1] # all cores from one node
#BSUB -o %J.stdout # output log (%J : JobID)
#BSUB -eo %J.stderr # error log
#BSUB -L /bin/bash # Initialize the execution environment


module load R

Rscript emp_r2_by_maf.r

#!/usr/bin/env Rscript

options(echo=TRUE)
print("In R now!")

## SCRIPT STARTS BELOW

library(data.table)
info <- fread("all.chr.info",header=T)
info$EmpRsq <- as.character(info$EmpRsq)
info$EmpRsq <- as.numeric(info$EmpRsq)

# Typed_Only SNPs don't have EmpRsq

#explain how Rsq is calculated: https://genome.sph.umich.edu/wiki/Minimac3_Info_File

# dat <- info[which(!is.na(info$EmpRsq)),]
# table(dat$Genotyped) # All SNPs with EmpRsq are the genotyped SNPs

genotyped <- info[which(info$Genotyped=="Genotyped"),]

dim(info)
dim(genotyped)

# calculate R2 by MAF
mean_r2 <- array(NA, dim=c(5000,2))
colnames(mean_r2) <- c("MAF","EmpRsq")

for (i in 1:5000) {
	if (i==1 ) { 
   dat <- subset(genotyped, MAF >= (i/10000 - 0.0001) & MAF <= i/10000)
   mean_r2[i,1] <- i/10000
   mean_r2[i,2] <- mean(dat$EmpRsq)
	 }
	else  {
   dat <- subset(genotyped, MAF > (i/10000 - 0.0001) & MAF <= i/10000)
   mean_r2[i,1] <- i/10000
   mean_r2[i,2] <- mean(dat$EmpRsq)
   }
}

write.table(mean_r2,'mean_emp_r2',sep='\t', row.names=FALSE,quote=FALSE)

bsub < emp_r2_by_maf.sh

## the genotyped SNP MAF are concentrated at the rare end than the common end
# among those genotyped SNPs with MAF < 1%, the number of SNPs per 0.0005 bin are sort of evenly distributed (more towards the lower end, 0 - 0.002)

## examine why the empRsq zigzag

dat1 <- genotyped[which(genotyped$MAF >= 0 & genotyped$MAF <=0.0001 ),]
dat2 <- genotyped[which(genotyped$MAF > 0.0001 & genotyped$MAF <=0.0002),]
dat3 <- genotyped[which(genotyped$MAF > 0.0002 & genotyped$MAF <=0.0003),]
dat4 <- genotyped[which(genotyped$MAF > 0.0003 & genotyped$MAF <=0.0004),]
dat5 <- genotyped[which(genotyped$MAF > 0.0004 & genotyped$MAF <=0.0005),]
dat6 <- genotyped[which(genotyped$MAF > 0.0005 & genotyped$MAF <=0.0006),]
dat7 <- genotyped[which(genotyped$MAF > 0.0006 & genotyped$MAF <=0.0007),]
dat8 <- genotyped[which(genotyped$MAF > 0.0007 & genotyped$MAF <=0.0008),]
dat9 <- genotyped[which(genotyped$MAF > 0.0008 & genotyped$MAF <=0.0009),]
dat10 <- genotyped[which(genotyped$MAF > 0.0009 & genotyped$MAF <=0.0010),]
dat11 <- genotyped[which(genotyped$MAF > 0.0010 & genotyped$MAF <=0.0011),]
dat12 <- genotyped[which(genotyped$MAF > 0.0011 & genotyped$MAF <=0.0012),]
dat13 <- genotyped[which(genotyped$MAF > 0.0012 & genotyped$MAF <=0.0013),]

dim(dat1) # 430 (0 to 0.0001) 
dim(dat2) # 96 (0.0001 to 0.0002)
dim(dat3) # 2203 (0.0002 to 0.0003) (Minor allele count =1 out of 3628 alleles at a locus (1814*2))
dim(dat4) # 985 (0.0003 to 0.0004)
dim(dat5) # 241 (0.0004 to 0.0005)
dim(dat6) # 1946 (0.0005 to 0.0006) (Minor allele count =2 out of 3628 alleles at a locus (1814*2))
dim(dat7) # 560 (0.0006 to 0.0007)
dim(dat8) # 283 (0.0007 to 0.0008)
dim(dat9) # 1963 (0.0008 to 0.0009) (Minor allele count =3 out of 3628 alleles at a locus (1814*2))
dim(dat10) # 474 (0.0009 to 0.0010)
dim(dat11) # 608 (0.0010 to 0.0011) 
dim(dat12) # 1829 (0.0011 to 0.0012)  (Minor allele count =4 out of 3628 alleles at a locus (1814*2))
dim(dat13) # 424 (0.0012 to 0.0013)


# the more SNPs in each bin, the higher the average empR2
par(mfrow = c(4, 3))
hist(dat1$EmpRsq)
hist(dat2$EmpRsq)
hist(dat3$EmpRsq)
hist(dat4$EmpRsq)
hist(dat5$EmpRsq)
hist(dat6$EmpRsq)
hist(dat7$EmpRsq)
hist(dat8$EmpRsq)
hist(dat9$EmpRsq)
hist(dat10$EmpRsq)

##############################################################################
####         R script: empirical R2 by MAF (using bigger MAF bins)        ####
##############################################################################

# what about calculating MAF per MAF bins (0, 0.0005, 0.001, 0.002, 0.005, 0.010, 0.015, 0.020, 0.035, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)

#!/bin/bash
#BSUB -J emp_r2_by_maf_new # Job name
#BSUB -n 1
#BSUB -P acc_psychgen # allocation account
#BSUB -q premium # queue
#BSUB -W 24:00 # walltime in HH:MM
#BSUB -R rusage[mem=32G] # 120 GB of memory 
#BSUB -R span[hosts=1] # all cores from one node
#BSUB -o %J.stdout # output log (%J : JobID)
#BSUB -eo %J.stderr # error log
#BSUB -L /bin/bash # Initialize the execution environment


module load R

Rscript emp_r2_by_maf_new.r

#!/usr/bin/env Rscript

options(echo=TRUE)
print("In R now!")

## SCRIPT STARTS BELOW

library(data.table)
info <- fread("all.chr.info",header=T)
info$EmpRsq <- as.character(info$EmpRsq)
info$EmpRsq <- as.numeric(info$EmpRsq)

# Typed_Only SNPs don't have EmpRsq

#explain how Rsq is calculated: https://genome.sph.umich.edu/wiki/Minimac3_Info_File

# dat <- info[which(!is.na(info$EmpRsq)),]
# table(dat$Genotyped) # All SNPs with EmpRsq are the genotyped SNPs

genotyped <- info[which(info$Genotyped=="Genotyped"),]

dim(info)
dim(genotyped)

# calculate R2 by MAF
mean_r2 <- array(NA, dim=c(14,2))
colnames(mean_r2) <- c("MAF_bin","EmpRsq")


for (i in 1:14) {
	if (i==1 ) { 
   dat <- subset(genotyped, MAF >= 0 & MAF <= 0.0005)
   mean_r2[i,1] <- "[0,0.0005]"
   mean_r2[i,2] <- mean(dat$EmpRsq)
	 }
	else if (i==2 ) { 
   dat <- subset(genotyped, MAF > 0.0005 & MAF <= 0.001)
   mean_r2[i,1] <- "(0.0005,0.001]"
   mean_r2[i,2] <- mean(dat$EmpRsq)
	 }
	else if (i==3 ) { 
   dat <- subset(genotyped, MAF > 0.001 & MAF <= 0.002)
   mean_r2[i,1] <- "(0.001,0.002]"
   mean_r2[i,2] <- mean(dat$EmpRsq)
	 }
	else if (i==4 ) { 
   dat <- subset(genotyped, MAF > 0.002 & MAF <= 0.005)
   mean_r2[i,1] <- "(0.002,0.005]"
   mean_r2[i,2] <- mean(dat$EmpRsq)
	 }
	else if (i==5 ) { 
   dat <- subset(genotyped, MAF > 0.005 & MAF <= 0.01)
   mean_r2[i,1] <- "(0.005,0.01]"
   mean_r2[i,2] <- mean(dat$EmpRsq)
	 }
	else if (i==6 ) { 
   dat <- subset(genotyped, MAF > 0.01 & MAF <= 0.015)
   mean_r2[i,1] <- "(0.01,0.015]"
   mean_r2[i,2] <- mean(dat$EmpRsq)
	 }
	else if (i==7 ) { 
   dat <- subset(genotyped, MAF > 0.015 & MAF <= 0.02)
   mean_r2[i,1] <- "(0.015,0.02]"
   mean_r2[i,2] <- mean(dat$EmpRsq)
	 }
	else if (i==8 ) { 
   dat <- subset(genotyped, MAF > 0.02 & MAF <= 0.035)
   mean_r2[i,1] <- "(0.02,0.035]"
   mean_r2[i,2] <- mean(dat$EmpRsq)
	 }
	else if (i==9 ) { 
   dat <- subset(genotyped, MAF > 0.035 & MAF <= 0.05)
   mean_r2[i,1] <- "(0.035,0.05]"
   mean_r2[i,2] <- mean(dat$EmpRsq)
	 }
	else if (i==10 ) { 
   dat <- subset(genotyped, MAF > 0.05 & MAF <= 0.1)
   mean_r2[i,1] <- "(0.05,0.1]"
   mean_r2[i,2] <- mean(dat$EmpRsq)
	 }
	else if (i==11 ) { 
   dat <- subset(genotyped, MAF > 0.1 & MAF <= 0.2)
   mean_r2[i,1] <- "(0.1,0.2]"
   mean_r2[i,2] <- mean(dat$EmpRsq)
	 }
	else if (i==12 ) { 
   dat <- subset(genotyped, MAF > 0.2 & MAF <= 0.3)
   mean_r2[i,1] <- "(0.2,0.3]"
   mean_r2[i,2] <- mean(dat$EmpRsq)
	 }
	else if (i==13 ) { 
   dat <- subset(genotyped, MAF > 0.3 & MAF <= 0.4)
   mean_r2[i,1] <- "(0.3,0.4]"
   mean_r2[i,2] <- mean(dat$EmpRsq)
	 }
	else if (i==14 ) { 
   dat <- subset(genotyped, MAF > 0.4 & MAF <= 0.5)
   mean_r2[i,1] <- "(0.4,0.5]"
   mean_r2[i,2] <- mean(dat$EmpRsq)
	 }
}

mean_r2 <- as.data.frame(mean_r2)

write.table(mean_r2,'mean_emp_r2_new_MAF_bin',sep='\t', row.names=FALSE,quote=FALSE)

bsub < emp_r2_by_maf_new.sh
##################################################
####             plot the Rsq by MAF          ####
##################################################

# in the folder: /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_compare 

# plot the Rsq by MAF
all_topmed <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/mean_r2", header=T) # 5000 rows
all_1000g <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/mean_r2", header=T) # 5000 rows
all_ga <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_GA/mean_r2", header=T) # 5000 rows

all_topmed$ref_panel <- "TOPMed"
all_1000g$ref_panel <- "1000G"
all_ga$ref_panel <- "GenomeAsia"

merged <- rbind(all_topmed,all_1000g,all_ga)
# 15,000 rows 
merged$ref_panel <- factor(merged$ref_panel,  levels=c("TOPMed","1000G","GenomeAsia"))

# color blind friendly
colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
scales::show_col(colorBlindGrey8)

library(ggplot2)
library(dplyr)
p <- ggplot(merged, aes(x = MAF, y = Rsq,group=ref_panel, color=ref_panel))  + geom_line() + geom_point(shape=1)+ theme_classic(base_size=20) + coord_cartesian(ylim = c(0, 1)) + theme(text = element_text(size = 15,colour="black"))   #all font size
p + scale_color_manual(values=c( "#E68F00", "#56B4E9","#CC79A7"),name="Reference Panel") + labs(x="Minor Allele Frequency", y="Estimated R2")

# ulta rare
ultra_rare <- merged[which(merged$MAF <= 0.001),]
p <- ggplot(ultra_rare, aes(x = MAF, y = Rsq,group=ref_panel, color=ref_panel))  + geom_line() + geom_point(shape=1)+ theme_classic(base_size=20) + coord_cartesian(ylim = c(0, 1)) + theme(text = element_text(size = 15,colour="black"))   #all font size
p + scale_color_manual(values=c( "#E68F00", "#56B4E9","#CC79A7"),name="Reference Panel") + labs(x="Minor Allele Frequency", y="Estimated R2")

# rare 
rare <- merged[which(merged$MAF > 0.001 & merged$MAF <= 0.01),]
p <- ggplot(rare, aes(x = MAF, y = Rsq,group=ref_panel, color=ref_panel))  + geom_line() + geom_point(shape=1)+theme_classic(base_size=20) + coord_cartesian(ylim = c(0, 1)) + theme(text = element_text(size = 15 ,colour="black"))   #all font size
p + scale_color_manual(values=c( "#E68F00", "#56B4E9","#CC79A7"),name="Reference Panel") + labs(x="Minor Allele Frequency", y="Estimated R2")

# common
common <- merged[which(merged$MAF > 0.01 ),]
p <- ggplot(common, aes(x = MAF, y = Rsq,group=ref_panel, color=ref_panel))  + geom_line() + geom_point(shape=1)+ theme_classic(base_size=20) + coord_cartesian(ylim = c(0, 1)) + theme(text = element_text(size = 15 ,colour="black"))   #all font size
p + scale_color_manual(values=c( "#E68F00", "#56B4E9","#CC79A7"),name="Reference Panel") + labs(x="Minor Allele Frequency", y="Estimated R2")

## all 5 panels

# in the folder: /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_compare 


# plot the Rsq by MAF
all_topmed <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/mean_r2", header=T) # 5000 rows
all_1000g <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/mean_r2", header=T) # 5000 rows
all_ga <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_GA/mean_r2", header=T) # 5000 rows
all_1000g_37 <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_37/mean_r2", header=T) # 5000 rows
all_1000g_38 <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_38/mean_r2", header=T) # 5000 rows

all_topmed$ref_panel <- "TOPMed"
all_1000g$ref_panel <- "Expanded_1000G"
all_ga$ref_panel <- "GenomeAsia_Pilot"
all_1000g_37$ref_panel <- "1000G_GRCh37_SAS"
all_1000g_38$ref_panel <- "1000G_GRCh38"

merged <- rbind(all_topmed,all_ga,all_1000g,all_1000g_38,all_1000g_37)
merged$ref_panel <- factor(merged$ref_panel,  levels=c("TOPMed","GenomeAsia_Pilot","Expanded_1000G","1000G_GRCh38","1000G_GRCh37_SAS"))
# 25,000

# color blind friendly
colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
scales::show_col(colorBlindGrey8)

library(ggplot2)
library(dplyr)
p <- ggplot(merged, aes(x = MAF, y = Rsq,group=ref_panel, color=ref_panel))  + geom_line() + geom_point(shape=1)+ theme_classic(base_size=20) + coord_cartesian(ylim = c(0, 1)) + theme(text = element_text(size = 15,colour="black"))   #all font size
p + scale_color_manual(values=c( "#E68F00", "#56B4E9","#CC79A7","#009E73","#F0E442"),name="Reference Panel") + labs(x="Minor Allele Frequency", y="Estimated R2")

# ulta rare
ultra_rare <- merged[which(merged$MAF <= 0.001),]
p <- ggplot(ultra_rare, aes(x = MAF, y = Rsq,group=ref_panel, color=ref_panel))  + geom_line() + geom_point(shape=1)+ theme_classic(base_size=20) + coord_cartesian(ylim = c(0, 1)) + theme(text = element_text(size = 15,colour="black"))   #all font size
p + scale_color_manual(values=c( "#E68F00", "#56B4E9","#CC79A7","#009E73","#F0E442"),name="Reference Panel") + labs(x="Minor Allele Frequency", y="Estimated R2")

# rare 
rare <- merged[which(merged$MAF > 0.001 & merged$MAF <= 0.01),]
p <- ggplot(rare, aes(x = MAF, y = Rsq,group=ref_panel, color=ref_panel))  + geom_line() + geom_point(shape=1)+theme_classic(base_size=20) + coord_cartesian(ylim = c(0, 1)) + theme(text = element_text(size = 15 ,colour="black"))   #all font size
p + scale_color_manual(values=c( "#E68F00", "#56B4E9","#CC79A7","#009E73","#F0E442"),name="Reference Panel") + labs(x="Minor Allele Frequency", y="Estimated R2")

# common
common <- merged[which(merged$MAF > 0.01 ),]
p <- ggplot(common, aes(x = MAF, y = Rsq,group=ref_panel, color=ref_panel))  + geom_line() + geom_point(shape=1)+ theme_classic(base_size=20) + coord_cartesian(ylim = c(0, 1)) + theme(text = element_text(size = 15 ,colour="black"))   #all font size
p + scale_color_manual(values=c( "#E68F00", "#56B4E9","#CC79A7","#009E73","#F0E442"),name="Reference Panel") + labs(x="Minor Allele Frequency", y="Estimated R2")



########## plot the Rsq by MAF using the AggRSquare MAF bins
# calculate the mean of R2 for meta-imputation in between MAF 0.0005 and 0.001
awk '($2>0.0005 && $2 <= 0.001){print $0}' all.chr.info| awk '{ total += $3 } END { print total/NR } '
awk '($2>0.1&& $2 <= 0.2){print $0}' all.chr.info| awk '{ total += $3 } END { print total/NR } '

# in R
metaimp <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_topmed_1000G/mean_r2_new_MAF_bin",header=T)
topmed <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/mean_r2_new_MAF_bin",header=T)
ex1000G <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/mean_r2_new_MAF_bin",header=T)
all_1000G_38 <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_38/mean_r2_new_MAF_bin",header=T)
sas_1000G_37 <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_37/mean_r2_new_MAF_bin",header=T)
GA <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_GA/mean_r2_new_MAF_bin",header=T)

# update the SNP ID to the same as tarseq ID (check forward and reverse)

metaimp$ref_panel <- "TOPMed_1000G"
topmed$ref_panel <- "TOPMed"
ex1000G$ref_panel <- "Expanded_1000G"
all_1000G_38$ref_panel <- "1000G_GRCh38"
sas_1000G_37$ref_panel <- "1000G_GRCh37_SAS"
GA$ref_panel <- "GenomeAsia_Pilot"

merged <- rbind(metaimp,topmed,ex1000G,all_1000G_38,sas_1000G_37,GA) #  84 rows, 7 columns
merged$ref_panel <- factor(merged$ref_panel,  levels=c("TOPMed_1000G","TOPMed","Expanded_1000G","1000G_GRCh38", "1000G_GRCh37_SAS", "GenomeAsia_Pilot"))
merged$MAF_bin <- factor(merged$MAF_bin,  levels=c("[0,0.0005]","(0.0005,0.001]","(0.001,0.002]","(0.002,0.005]","(0.005,0.01]","(0.01,0.015]","(0.015,0.02]","(0.02,0.035]","(0.035,0.05]","(0.05,0.1]","(0.1,0.2]","(0.2,0.3]","(0.3,0.4]","(0.4,0.5]"))

# color blind friendly
colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
scales::show_col(colorBlindGrey8)

library(ggplot2)
library(dplyr)

pdf('r2_maf_all_snp.pdf',height=10,width=15)
p <- ggplot(merged, aes(x = MAF_bin, y = Rsq,group=ref_panel, color=ref_panel))  + geom_line(size=1) + geom_point(shape=1)+ theme_classic(base_size=30) + coord_cartesian(ylim = c(0, 1)) + scale_x_discrete(guide = guide_axis(angle = 45))  + theme(text = element_text(size = 30,colour="black"))  + geom_hline(yintercept=0.8,  color = "dark grey") + geom_hline(yintercept=0.6, linetype="dashed", color = "dark grey")  + geom_vline(xintercept="(0.01,0.015]",  color = "dark grey",linetype="dashed") + geom_vline(xintercept="(0.05,0.1]",  color = "dark grey") #all font size
p + scale_color_manual(values=c("#661100", "#D55E00", "#E69F00","#0072B2","#332288","#999999"),name="Reference Panel") + labs(x="Minor Allele Frequency", y="Average Estimated R2") +
  annotate('text', x = "[0,0.0005]", y = 0.2, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.0005,0.001]", y = 0.5, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.001,0.002]", y = 0.6, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.002,0.005]", y = 0.7, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.005,0.01]", y = 0.75, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.01,0.015]", y = 0.8, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.015,0.02]", y = 0.85, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.02,0.035]", y = 0.9, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.035,0.05]", y = 0.94, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.05,0.1]", y = 0.95, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.1,0.2]", y = 0.96, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.2,0.3]", y = 0.97, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.3,0.4]", y = 0.98, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.4,0.5]", y = 0.99, label='"**"', parse=TRUE, color = "red",size=10)  + theme(plot.title = element_text(hjust = 0.5)) + labs(title="All SNPs")
dev.off()

merged[which(merged$MAF_bin=="[0,0.0005]"),]
## what about significant difference?

# 5 imputation panels
awk '($8!="Typed_Only" && $5 >= 0 && $5 <= 0.0005) {print $0}' all.chr.info > all.maf.0.0005.info
awk '($8!="Typed_Only" && $5 > 0.0005 && $5 <= 0.001) {print $0}' all.chr.info > all.maf.0.001.info
awk '($8!="Typed_Only" && $5 > 0.001 && $5 <= 0.002) {print $0}' all.chr.info > all.maf.0.002.info
awk '($8!="Typed_Only" && $5 > 0.002 && $5 <= 0.005) {print $0}' all.chr.info > all.maf.0.005.info
awk '($8!="Typed_Only" && $5 > 0.005 && $5 <= 0.01) {print $0}' all.chr.info > all.maf.0.01.info
awk '($8!="Typed_Only" && $5 > 0.01 && $5 <= 0.015) {print $0}' all.chr.info > all.maf.0.015.info
awk '($8!="Typed_Only" && $5 > 0.015 && $5 <= 0.02) {print $0}' all.chr.info > all.maf.0.02.info
awk '($8!="Typed_Only" && $5 > 0.02 && $5 <= 0.035) {print $0}' all.chr.info > all.maf.0.035.info
awk '($8!="Typed_Only" && $5 > 0.035 && $5 <= 0.05) {print $0}' all.chr.info > all.maf.0.05.info
awk '($8!="Typed_Only" && $5 > 0.05 && $5 <= 0.1) {print $0}' all.chr.info > all.maf.0.1.info
awk '($8!="Typed_Only" && $5 > 0.1 && $5 <= 0.2) {print $0}' all.chr.info > all.maf.0.2.info
awk '($8!="Typed_Only" && $5 > 0.2 && $5 <= 0.3) {print $0}' all.chr.info > all.maf.0.3.info
awk '($8!="Typed_Only" && $5 > 0.3 && $5 <= 0.4) {print $0}' all.chr.info > all.maf.0.4.info
awk '($8!="Typed_Only" && $5 > 0.4 && $5 <= 0.5) {print $0}' all.chr.info > all.maf.0.5.info

# meta-imputation 
awk '( $2 >= 0 && $2 <= 0.0005) {print $0}' all.topmed_1000g.info > all.maf.0.0005.info
awk '( $2 > 0.0005 && $2 <= 0.001) {print $0}' all.topmed_1000g.info > all.maf.0.001.info
awk '( $2 > 0.001 && $2 <= 0.002) {print $0}' all.topmed_1000g.info > all.maf.0.002.info
awk '( $2 > 0.002 && $2 <= 0.005) {print $0}' all.topmed_1000g.info > all.maf.0.005.info
awk '( $2 > 0.005 && $2 <= 0.01) {print $0}' all.topmed_1000g.info > all.maf.0.01.info
awk '( $2 > 0.01 && $2 <= 0.015) {print $0}' all.topmed_1000g.info > all.maf.0.015.info
awk '( $2 > 0.015 && $2 <= 0.02) {print $0}' all.topmed_1000g.info > all.maf.0.02.info
awk '( $2 > 0.02 && $2 <= 0.035) {print $0}' all.topmed_1000g.info > all.maf.0.035.info
awk '( $2 > 0.035 && $2 <= 0.05) {print $0}' all.topmed_1000g.info > all.maf.0.05.info
awk '( $2 > 0.05 && $2 <= 0.1) {print $0}' all.topmed_1000g.info > all.maf.0.1.info
awk '( $2 > 0.1 && $2 <= 0.2) {print $0}' all.topmed_1000g.info > all.maf.0.2.info
awk '( $2 > 0.2 && $2 <= 0.3) {print $0}' all.topmed_1000g.info > all.maf.0.3.info
awk '( $2 > 0.3 && $2 <= 0.4) {print $0}' all.topmed_1000g.info > all.maf.0.4.info
awk '( $2 > 0.4 && $2 <= 0.5) {print $0}' all.topmed_1000g.info > all.maf.0.5.info


# in R

#!/bin/bash
#BSUB -J get_p_value_for_rsq_by_panel # Job name
#BSUB -n 1
#BSUB -P acc_psychgen # allocation account
#BSUB -q premium # queue
#BSUB -W 48:00 # walltime in HH:MM
#BSUB -R himem
#BSUB -R rusage[mem=200G] # 120 GB of memory 
#BSUB -R span[hosts=1] # all cores from one node
#BSUB -o %J.stdout # output log (%J : JobID)
#BSUB -eo %J.stderr # error log
#BSUB -L /bin/bash # Initialize the execution environment


module load R

Rscript get_p_value_for_rsq_by_panel.r

bsub < get_p_value_for_rsq_by_panel.sh


#!/usr/bin/env Rscript

options(echo=TRUE)
print("In R now!")

## SCRIPT STARTS BELOW

library(data.table)
metaimp <- fread("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_topmed_1000G/all.maf.0.5.info",header=F)
topmed <- fread("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/all.maf.0.5.info",header=F)
ex1000G <- fread("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/all.maf.0.5.info",header=F)
all_1000G_38 <- fread("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_38/all.maf.0.5.info",header=F)
sas_1000G_37 <- fread("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_37/all.maf.0.5.info",header=F)
GA <- fread("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_GA/all.maf.0.5.info",header=F)

metaimp <- as.data.frame(metaimp)
topmed <- as.data.frame(topmed)
ex1000G <- as.data.frame(ex1000G)
all_1000G_38 <- as.data.frame(all_1000G_38)
sas_1000G_37 <- as.data.frame(sas_1000G_37)
GA <- as.data.frame(GA)

topmed <- topmed[,c(1,5,7)]
ex1000G <- ex1000G[,c(1,5,7)]
all_1000G_38 <- all_1000G_38[,c(1,5,7)]
sas_1000G_37 <- sas_1000G_37[,c(1,5,7)]
GA <- GA[,c(1,5,7)]

colnames(metaimp) <- c("SNP","MAF","Rsq")
colnames(topmed) <- c("SNP","MAF","Rsq")
colnames(ex1000G) <- c("SNP","MAF","Rsq")
colnames(all_1000G_38) <- c("SNP","MAF","Rsq")
colnames(sas_1000G_37) <- c("SNP","MAF","Rsq")
colnames(GA) <- c("SNP","MAF","Rsq")

metaimp$ref_panel <- "TOPMed_1000G"
topmed$ref_panel <- "TOPMed"
ex1000G$ref_panel <- "Expanded_1000G"
all_1000G_38$ref_panel <- "1000G_GRCh38"
sas_1000G_37$ref_panel <- "1000G_GRCh37_SAS"
GA$ref_panel <- "GenomeAsia_Pilot"

dim(metaimp) # [1] 1034699       4
dim(topmed)  # [1] 1027957       4
dim(ex1000G)  # [1] 1020480       4
dim(all_1000G_38)  # [1] 1107046       4
dim(sas_1000G_37)  # [1] 1195460       4
dim(GA)   # [1] 859417      4

merged <- rbind(metaimp,topmed,ex1000G,all_1000G_38,sas_1000G_37,GA)  # [1] 6245059       4

kruskal.test(Rsq ~ ref_panel, data = merged)
pairwise.t.test(merged$Rsq, merged$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) 

# Non-parametric alternative to one-way ANOVA test

# "[0,0.0005]"
kruskal.test(Rsq ~ ref_panel, data = merged) # p-value < 2.2e-16
# "(0.0005,0.001]"
kruskal.test(Rsq ~ ref_panel, data = merged) # p-value < 2.2e-16
# "(0.001,0.002]"
kruskal.test(Rsq ~ ref_panel, data = merged) # p-value < 2.2e-16
# "(0.002,0.005]"
kruskal.test(Rsq ~ ref_panel, data = merged) # p-value < 2.2e-16
# "(0.005,0.01]"
kruskal.test(Rsq ~ ref_panel, data = merged) # p-value < 2.2e-16
# "(0.01,0.015]"
kruskal.test(Rsq ~ ref_panel, data = merged) # p-value < 2.2e-16
# "(0.015,0.02]"
kruskal.test(Rsq ~ ref_panel, data = merged) # p-value < 2.2e-16
# "(0.02,0.035]"
kruskal.test(Rsq ~ ref_panel, data = merged) # p-value < 2.2e-16
# "(0.035,0.05]"
kruskal.test(Rsq ~ ref_panel, data = merged) # p-value < 2.2e-16
# "(0.05,0.1]"
kruskal.test(Rsq ~ ref_panel, data = merged) # p-value < 2.2e-16
# "(0.1,0.2]"
kruskal.test(Rsq ~ ref_panel, data = merged) # p-value < 2.2e-16
# "(0.2,0.3]"
kruskal.test(Rsq ~ ref_panel, data = merged) # p-value < 2.2e-16
# "(0.3,0.4]"
kruskal.test(Rsq ~ ref_panel, data = merged) # p-value < 2.2e-16
# "(0.4,0.5]"
kruskal.test(Rsq ~ ref_panel, data = merged) # p-value < 2.2e-16
# 0.05/14 = 0.003571 (** Bonferroni significant )
# 0.05 (* nominal significant)

# check for pairwise difference
# Pairwise t-tests with no assumption of equal variances

# "[0,0.0005]"
pairwise.t.test(merged$Rsq, merged$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) # 1000G GRCh38 (highest) better than all other panels 
# if we only take Bonferroni significant results (0.05/15 = 0.0033), then expanded 1000G is as good as 1000G GRCh38

#                 1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     <2e-16           -            -              -               
# Expanded_1000G   <2e-16           0.0068       -              -               
# GenomeAsia_Pilot <2e-16           <2e-16       <2e-16         -               
# TOPMed           <2e-16           <2e-16       <2e-16         <2e-16          
# TOPMed_1000G     <2e-16           <2e-16       <2e-16         <2e-16          
#                  TOPMed
# 1000G_GRCh38     -     
# Expanded_1000G   -     
# GenomeAsia_Pilot -     
# TOPMed           -     
# TOPMed_1000G     <2e-16

# "(0.0005,0.001]"
pairwise.t.test(merged$Rsq, merged$ref_panel,  p.adjust.method = "none", pool.sd = FALSE)  # meta-imputation better than all other panels

#                  1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     <2e-16           -            -              -               
# Expanded_1000G   <2e-16           <2e-16       -              -               
# GenomeAsia_Pilot <2e-16           <2e-16       <2e-16         -               
# TOPMed           <2e-16           <2e-16       <2e-16         <2e-16          
# TOPMed_1000G     <2e-16           <2e-16       <2e-16         <2e-16          
#                  TOPMed
# 1000G_GRCh38     -     
# Expanded_1000G   -     
# GenomeAsia_Pilot -     
# TOPMed           -     
# TOPMed_1000G     <2e-16

# "(0.001,0.002]"
pairwise.t.test(merged$Rsq, merged$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) # meta-imputation better than all other panels

#                  1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     <2e-16           -            -              -               
# Expanded_1000G   <2e-16           <2e-16       -              -               
# GenomeAsia_Pilot <2e-16           <2e-16       <2e-16         -               
# TOPMed           <2e-16           <2e-16       <2e-16         <2e-16          
# TOPMed_1000G     <2e-16           <2e-16       <2e-16         <2e-16          
#                  TOPMed
# 1000G_GRCh38     -     
# Expanded_1000G   -     
# GenomeAsia_Pilot -     
# TOPMed           -     
# TOPMed_1000G     <2e-16

# "(0.002,0.005]"
pairwise.t.test(merged$Rsq, merged$ref_panel,  p.adjust.method = "none", pool.sd = FALSE)  # meta-imputation better than all other panels

#                  1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     < 2e-16          -            -              -               
# Expanded_1000G   < 2e-16          < 2e-16      -              -               
# GenomeAsia_Pilot < 2e-16          < 2e-16      < 2e-16        -               
# TOPMed           < 2e-16          < 2e-16      2.7e-14        < 2e-16         
# TOPMed_1000G     < 2e-16          < 2e-16      < 2e-16        < 2e-16         
#                  TOPMed 
# 1000G_GRCh38     -      
# Expanded_1000G   -      
# GenomeAsia_Pilot -      
# TOPMed           -      
# TOPMed_1000G     < 2e-16

# "(0.005,0.01]"
pairwise.t.test(merged$Rsq, merged$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) # meta-imputation better than all other panels

#                  1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     <2e-16           -            -              -               
# Expanded_1000G   <2e-16           <2e-16       -              -               
# GenomeAsia_Pilot <2e-16           <2e-16       <2e-16         -               
# TOPMed           <2e-16           <2e-16       <2e-16         <2e-16          
# TOPMed_1000G     <2e-16           <2e-16       <2e-16         <2e-16          
#                  TOPMed
# 1000G_GRCh38     -     
# Expanded_1000G   -     
# GenomeAsia_Pilot -     
# TOPMed           -     
# TOPMed_1000G     <2e-16

# "(0.01,0.015]"
pairwise.t.test(merged$Rsq, merged$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) # expanded 1000G better than all other panels, second place: meta-imputation better than all remaining panels except TOPMed

#                  1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     <2e-16           -            -              -               
# Expanded_1000G   <2e-16           <2e-16       -              -               
# GenomeAsia_Pilot <2e-16           <2e-16       <2e-16         -               
# TOPMed           <2e-16           <2e-16       <2e-16         <2e-16          
# TOPMed_1000G     <2e-16           <2e-16       <2e-16         <2e-16          
#                  TOPMed
# 1000G_GRCh38     -     
# Expanded_1000G   -     
# GenomeAsia_Pilot -     
# TOPMed           -     
# TOPMed_1000G     0.17  

# "(0.015,0.02]"
pairwise.t.test(merged$Rsq, merged$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) # expanded 1000G better than all other panels, second place: meta-imputation better than all remaining panels

#                 1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     < 2e-16          -            -              -               
# Expanded_1000G   < 2e-16          < 2e-16      -              -               
# GenomeAsia_Pilot 1.8e-12          < 2e-16      < 2e-16        -               
# TOPMed           < 2e-16          1.0e-06      < 2e-16        < 2e-16         
# TOPMed_1000G     < 2e-16          < 2e-16      < 2e-16        < 2e-16         
#                  TOPMed 
# 1000G_GRCh38     -      
# Expanded_1000G   -      
# GenomeAsia_Pilot -      
# TOPMed           -      
# TOPMed_1000G     < 2e-16

# "(0.02,0.035]"
pairwise.t.test(merged$Rsq, merged$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) # expanded 1000G better than all other panels, second place: meta-imputation better than all remaining panels

#                  1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
#  1000G_GRCh38     <2e-16           -            -              -               
#  Expanded_1000G   <2e-16           <2e-16       -              -               
#  GenomeAsia_Pilot <2e-16           <2e-16       <2e-16         -               
#  TOPMed           <2e-16           <2e-16       <2e-16         <2e-16          
#  TOPMed_1000G     <2e-16           <2e-16       <2e-16         <2e-16          
#                   TOPMed
#  1000G_GRCh38     -     
#  Expanded_1000G   -     
#  GenomeAsia_Pilot -     
#  TOPMed           -     
#  TOPMed_1000G     <2e-16

# "(0.035,0.05]"
pairwise.t.test(merged$Rsq, merged$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) # expanded 1000G better than all other panels, second place: meta-imputation better than all remaining panels

#                  1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
#  1000G_GRCh38     <2e-16           -            -              -               
#  Expanded_1000G   <2e-16           <2e-16       -              -               
#  GenomeAsia_Pilot <2e-16           <2e-16       <2e-16         -               
#  TOPMed           <2e-16           <2e-16       <2e-16         <2e-16          
#  TOPMed_1000G     <2e-16           <2e-16       <2e-16         <2e-16          
#                   TOPMed
#  1000G_GRCh38     -     
#  Expanded_1000G   -     
#  GenomeAsia_Pilot -     
#  TOPMed           -     
#  TOPMed_1000G     <2e-16

# "(0.05,0.1]"
pairwise.t.test(merged$Rsq, merged$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) # expanded 1000G better than all other panels, second place: meta-imputation better than all remaining panels

#                  1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     <2e-16           -            -              -               
# Expanded_1000G   <2e-16           <2e-16       -              -               
# GenomeAsia_Pilot 0.49             <2e-16       <2e-16         -               
# TOPMed           <2e-16           <2e-16       <2e-16         <2e-16          
# TOPMed_1000G     <2e-16           <2e-16       <2e-16         <2e-16          
#                  TOPMed
# 1000G_GRCh38     -     
# Expanded_1000G   -     
# GenomeAsia_Pilot -     
# TOPMed           -     
# TOPMed_1000G     <2e-16

# "(0.1,0.2]"
pairwise.t.test(merged$Rsq, merged$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) # expanded 1000G better than all other panels, second place: meta-imputation better than all remaining panels

#                  1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     <2e-16           -            -              -               
# Expanded_1000G   <2e-16           <2e-16       -              -               
# GenomeAsia_Pilot <2e-16           <2e-16       <2e-16         -               
# TOPMed           <2e-16           <2e-16       <2e-16         <2e-16          
# TOPMed_1000G     <2e-16           <2e-16       <2e-16         <2e-16          
#                  TOPMed
# 1000G_GRCh38     -     
# Expanded_1000G   -     
# GenomeAsia_Pilot -     
# TOPMed           -     
# TOPMed_1000G     <2e-16

# "(0.2,0.3]"
pairwise.t.test(merged$Rsq, merged$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) # expanded 1000G better than all other panels, second place: meta-imputation better than all remaining panels

#                  1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     < 2e-16          -            -              -               
# Expanded_1000G   < 2e-16          < 2e-16      -              -               
# GenomeAsia_Pilot 1.1e-08          < 2e-16      < 2e-16        -               
# TOPMed           < 2e-16          < 2e-16      < 2e-16        < 2e-16         
# TOPMed_1000G     < 2e-16          < 2e-16      < 2e-16        < 2e-16         
#                  TOPMed 
# 1000G_GRCh38     -      
# Expanded_1000G   -      
# GenomeAsia_Pilot -      
# TOPMed           -      
# TOPMed_1000G     < 2e-16

# "(0.3,0.4]"
pairwise.t.test(merged$Rsq, merged$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) # expanded 1000G better than all other panels, second place: meta-imputation better than all remaining panels

#                  1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     <2e-16           -            -              -               
# Expanded_1000G   <2e-16           <2e-16       -              -               
# GenomeAsia_Pilot <2e-16           <2e-16       <2e-16         -               
# TOPMed           <2e-16           <2e-16       <2e-16         <2e-16          
# TOPMed_1000G     <2e-16           <2e-16       <2e-16         <2e-16          
#                  TOPMed
# 1000G_GRCh38     -     
# Expanded_1000G   -     
# GenomeAsia_Pilot -     
# TOPMed           -     
# TOPMed_1000G     <2e-16

# "(0.4,0.5]
pairwise.t.test(merged$Rsq, merged$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) # expanded 1000G better than all other panels, second place: meta-imputation better than all remaining panels

#                 1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     <2e-16           -            -              -               
# Expanded_1000G   <2e-16           <2e-16       -              -               
# GenomeAsia_Pilot <2e-16           <2e-16       <2e-16         -               
# TOPMed           <2e-16           <2e-16       <2e-16         <2e-16          
# TOPMed_1000G     <2e-16           <2e-16       <2e-16         <2e-16          
#                  TOPMed
# 1000G_GRCh38     -     
# Expanded_1000G   -     
# GenomeAsia_Pilot -     
# TOPMed           -     
# TOPMed_1000G     <2e-16

##############################################################################################################
####             plot the Rsq by MAF   (chr1 testing in TOPMed, expanded 1000G and meta-imputation)       ####
##############################################################################################################

# in the folder: /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_compare 

# plot the Rsq by MAF
all_topmed <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/mean_r2_chr1", header=T) # 5000 rows
all_1000g <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/mean_r2_chr1", header=T) # 5000 rows
all_meta <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_topmed_1000G/scratch/mean_r2_chr1", header=T) # 5000 rows

all_topmed$ref_panel <- "TOPMed"
all_1000g$ref_panel <- "Expanded_1000G"
all_meta$ref_panel <- "TOPMed_1000G"

merged <- rbind(all_topmed,all_1000g,all_meta)
# 15,000 rows 
merged$ref_panel <- factor(merged$ref_panel,  levels=c("TOPMed","Expanded_1000G","TOPMed_1000G"))

# color blind friendly
colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
scales::show_col(colorBlindGrey8)

library(ggplot2)
library(dplyr)
p <- ggplot(merged, aes(x = MAF, y = Rsq,group=ref_panel, color=ref_panel))  + geom_line() + geom_point(shape=1)+ theme_classic(base_size=20) + coord_cartesian(ylim = c(0, 1)) + theme(text = element_text(size = 15,colour="black"))   #all font size
p + scale_color_manual(values=c( "#E68F00", "#56B4E9","#CC79A7"),name="Reference Panel") + labs(x="Minor Allele Frequency", y="Estimated R2")

# ulta rare
ultra_rare <- merged[which(merged$MAF <= 0.001),]
p <- ggplot(ultra_rare, aes(x = MAF, y = Rsq,group=ref_panel, color=ref_panel))  + geom_line() + geom_point(shape=1)+ theme_classic(base_size=20) + coord_cartesian(ylim = c(0, 1)) + theme(text = element_text(size = 15,colour="black"))   #all font size
p + scale_color_manual(values=c( "#E68F00", "#56B4E9","#CC79A7"),name="Reference Panel") + labs(x="Minor Allele Frequency", y="Estimated R2")

# rare 
rare <- merged[which(merged$MAF > 0.001 & merged$MAF <= 0.01),]
p <- ggplot(rare, aes(x = MAF, y = Rsq,group=ref_panel, color=ref_panel))  + geom_line() + geom_point(shape=1)+theme_classic(base_size=20) + coord_cartesian(ylim = c(0, 1)) + theme(text = element_text(size = 15 ,colour="black"))   #all font size
p + scale_color_manual(values=c( "#E68F00", "#56B4E9","#CC79A7"),name="Reference Panel") + labs(x="Minor Allele Frequency", y="Estimated R2")

# common
common <- merged[which(merged$MAF > 0.01 ),]
p <- ggplot(common, aes(x = MAF, y = Rsq,group=ref_panel, color=ref_panel))  + geom_line() + geom_point(shape=1)+ theme_classic(base_size=20) + coord_cartesian(ylim = c(0, 1)) + theme(text = element_text(size = 15 ,colour="black"))   #all font size
p + scale_color_manual(values=c( "#E68F00", "#56B4E9","#CC79A7"),name="Reference Panel") + labs(x="Minor Allele Frequency", y="Estimated R2")


################################################################
####             plot the emp Rsq by MAF (masking)          ####
################################################################

########  MAF bin of 0.0001
# in the folder: /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_compare 

all_topmed <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/mean_emp_r2", header=T) # 5000 rows
all_1000g <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/mean_emp_r2", header=T) # 5000 rows
all_ga <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_GA/mean_emp_r2", header=T) # 5000 rows
all_1000g_37 <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_37/mean_emp_r2", header=T) # 5000 rows
all_1000g_38 <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_38/mean_emp_r2", header=T) # 5000 rows

all_topmed$ref_panel <- "TOPMed"
all_1000g$ref_panel <- "Expanded_1000G"
all_ga$ref_panel <- "GenomeAsia_Pilot"
all_1000g_37$ref_panel <- "1000G_GRCh37_SAS"
all_1000g_38$ref_panel <- "1000G_GRCh38"

merged <- rbind(all_topmed,all_ga,all_1000g,all_1000g_38,all_1000g_37)
merged$ref_panel <- factor(merged$ref_panel,  levels=c("TOPMed","GenomeAsia_Pilot","Expanded_1000G","1000G_GRCh38","1000G_GRCh37_SAS"))
# [1] 25000     3

# color blind friendly
library(ggplot2)
library(dplyr)
p <- ggplot(merged, aes(x = MAF, y = EmpRsq,group=ref_panel, color=ref_panel))  + geom_line() + geom_point(shape=1)+ theme_classic(base_size=20)  + coord_cartesian(ylim = c(0, 1)) + theme(text = element_text(size = 15))   #all font size
p + scale_color_manual(values=c( "#E68F00", "#56B4E9","#CC79A7","#009E73","#F0E442"),name="Reference Panel") + labs(x="Minor Allele Frequency", y="Empirical R2")
# original 3 colors: "#E68F00", "#56B4E9","#CC79A7"

# ulta rare
ultra_rare <- merged[which(merged$MAF <= 0.001),]
p <- ggplot(ultra_rare, aes(x = MAF, y = EmpRsq,group=ref_panel, color=ref_panel))  + geom_line() + geom_point(shape=1)+  theme_classic(base_size=20) + coord_cartesian(ylim = c(0, 1)) + theme(text = element_text(size = 20))   #all font size
p + scale_color_manual(values=c( "#E68F00", "#56B4E9","#CC79A7","#009E73","#F0E442"),name="Reference Panel") + labs(x="Minor Allele Frequency", y="Empirical R2")

# rare 
rare <- merged[which(merged$MAF > 0.001 & merged$MAF <= 0.01),]
p <- ggplot(rare, aes(x = MAF, y = EmpRsq,group=ref_panel, color=ref_panel))  + geom_line() + geom_point(shape=1)+ theme_classic(base_size=20) + coord_cartesian(ylim = c(0, 1)) + theme(text = element_text(size = 20))   #all font size
p + scale_color_manual(values=c( "#E68F00", "#56B4E9","#CC79A7","#009E73","#F0E442"),name="Reference Panel") + labs(x="Minor Allele Frequency", y="Empirical R2")

# common
common <- merged[which(merged$MAF > 0.01 ),]
p <- ggplot(common, aes(x = MAF, y = EmpRsq,group=ref_panel, color=ref_panel))  + geom_line() + geom_point(shape=1)+ theme_classic(base_size=20) + coord_cartesian(ylim = c(0, 1)) + theme(text = element_text(size = 20))   #all font size
p + scale_color_manual(values=c( "#E68F00", "#56B4E9","#CC79A7","#009E73","#F0E442"),name="Reference Panel") + labs(x="Minor Allele Frequency", y="Empirical R2")

# for GWAS, keep GSA and psych chip separate and meta-analyze 
# Given different genotyping arrays, genotyping locations, there are different factors contributing to heterogeneity, and have to be controlled for if merging these 2 datasets together)


########  MAF bin same as aggRSquare default

all_topmed <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/mean_emp_r2_new_MAF_bin", header=T) # 5000 rows
all_1000g <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/mean_emp_r2_new_MAF_bin", header=T) # 5000 rows
all_ga <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_GA/mean_emp_r2_new_MAF_bin", header=T) # 5000 rows
all_1000g_37 <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_37/mean_emp_r2_new_MAF_bin", header=T) # 5000 rows
all_1000g_38 <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_38/mean_emp_r2_new_MAF_bin", header=T) # 5000 rows

all_topmed$ref_panel <- "TOPMed"
all_1000g$ref_panel <- "Expanded_1000G"
all_ga$ref_panel <- "GenomeAsia_Pilot"
all_1000g_37$ref_panel <- "1000G_GRCh37_SAS"
all_1000g_38$ref_panel <- "1000G_GRCh38"

merged <- rbind(all_topmed,all_ga,all_1000g,all_1000g_38,all_1000g_37)
merged$ref_panel <- factor(merged$ref_panel,  levels=c("TOPMed","Expanded_1000G","1000G_GRCh38","1000G_GRCh37_SAS","GenomeAsia_Pilot"))
merged$MAF_bin <- factor(merged$MAF_bin,  levels=c("[0,0.0005]","(0.0005,0.001]","(0.001,0.002]","(0.002,0.005]","(0.005,0.01]","(0.01,0.015]","(0.015,0.02]","(0.02,0.035]","(0.035,0.05]","(0.05,0.1]","(0.1,0.2]","(0.2,0.3]","(0.3,0.4]","(0.4,0.5]"))
# [1] 70     3

# color blind friendly
library(ggplot2)
library(dplyr)
p <- ggplot(merged, aes(x = MAF_bin, y = EmpRsq,group=ref_panel, color=ref_panel))  + geom_line() + geom_point(shape=1)+ theme_classic(base_size=20)  + coord_cartesian(ylim = c(0, 1)) + scale_x_discrete(guide = guide_axis(angle = 45)) + theme(text = element_text(size = 15,,colour="black"))  + geom_hline(yintercept=0.8,  color = "red") + geom_hline(yintercept=0.6, linetype="dashed", color = "red")  + geom_vline(xintercept="(0.01,0.015]",  color = "black") + geom_vline(xintercept="(0.05,0.1]",  color = "black", linetype="dashed")  #all font size
p + scale_color_manual(values=c("#E68F00", "#56B4E9","#CC79A7","#999999","#F0E442"),name="Reference Panel") + labs(x="Minor Allele Frequency", y="Average Empirical R2")
# original 3 colors: "#E68F00", "#56B4E9","#CC79A7"

#### test for significance  ####
awk '($8=="Genotyped") {print $0}' all.chr.info > all.genotyped.info

library(data.table)
all_topmed <- fread("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/all.genotyped.info", header=T) # 5000 rows
all_1000g <- fread("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/all.genotyped.info", header=T) # 5000 rows
all_1000g_38 <- fread("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_38/all.genotyped.info", header=T) # 5000 rows
all_1000g_37 <- fread("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_37/all.genotyped.info", header=T) # 5000 rows
all_ga <- fread("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_GA/all.genotyped.info", header=T) # 5000 rows

all_topmed <- as.data.frame(all_topmed)
all_1000g <- as.data.frame(all_1000g)
all_1000g_38 <- as.data.frame(all_1000g_38)
all_1000g_37 <- as.data.frame(all_1000g_37)
all_ga <- as.data.frame(all_ga)

all_topmed$ref_panel <- "TOPMed"
all_1000g$ref_panel <- "Expanded_1000G"
all_1000g_38$ref_panel <- "1000G_GRCh38"
all_1000g_37$ref_panel <- "1000G_GRCh37_SAS"
all_ga$ref_panel <- "GenomeAsia_Pilot"

dim(all_topmed) # [1] 487217     14
dim(all_1000g) # [1] 484357     14
dim(all_1000g_38) # [1] 487835     14
dim(all_1000g_37) # [1] 491089     14
dim(all_ga) # [1] 461435     14

merged <- rbind(all_topmed,all_1000g,all_1000g_38,all_1000g_37,all_ga) # [1] 2411933      14

table(merged$ref_panel)
summary(merged$MAF)

library(dplyr)
merged$MAF_bin <- case_when(merged$MAF >= 0 & merged$MAF <= 0.0005 ~ "[0,0.0005]", merged$MAF > 0.0005 & merged$MAF <= 0.001 ~ "(0.0005,0.001]", merged$MAF > 0.001 & merged$MAF <= 0.002 ~ "(0.001,0.002]", merged$MAF > 0.002 & merged$MAF <= 0.005 ~ "(0.002,0.005]",merged$MAF > 0.005 & merged$MAF <= 0.01 ~ "(0.005,0.01]",merged$MAF > 0.01 & merged$MAF <= 0.015 ~ "(0.01,0.015]",merged$MAF > 0.015 & merged$MAF <= 0.02 ~ "(0.015,0.02]", merged$MAF > 0.02 & merged$MAF <= 0.035 ~ "(0.02,0.035]", merged$MAF > 0.035 & merged$MAF <= 0.05 ~ "(0.035,0.05]", merged$MAF > 0.05 & merged$MAF <= 0.1 ~ "(0.05,0.1]", merged$MAF > 0.1 & merged$MAF <= 0.2  ~ "(0.1,0.2]", merged$MAF > 0.2 & merged$MAF <= 0.3  ~ "(0.2,0.3]", merged$MAF > 0.3 & merged$MAF <= 0.4  ~ "(0.3,0.4]", merged$MAF > 0.4 & merged$MAF <= 0.5  ~ "(0.4,0.5]")
table(merged$MAF_bin) # add up to 2411933

# install.packages("plotrix")
# library(plotrix)
emp_r2_maf <- aggregate(EmpRsq~ref_panel+MAF_bin, data=merged, mean)
# FUN = function(x) c(mean = mean(x), se = std.error(x))
emp_r2_maf$ref_panel <- factor(emp_r2_maf$ref_panel,  levels=c("TOPMed","Expanded_1000G","1000G_GRCh38", "1000G_GRCh37_SAS", "GenomeAsia_Pilot"))
emp_r2_maf$MAF_bin <- factor(emp_r2_maf$MAF_bin,  levels=c("[0,0.0005]","(0.0005,0.001]","(0.001,0.002]","(0.002,0.005]","(0.005,0.01]","(0.01,0.015]","(0.015,0.02]","(0.02,0.035]","(0.035,0.05]","(0.05,0.1]","(0.1,0.2]","(0.2,0.3]","(0.3,0.4]","(0.4,0.5]"))

library(ggplot2)
library(dplyr)
pdf('emp_r2_maf_all_snp.pdf',height=10,width=15)
p <- ggplot(emp_r2_maf, aes(x = MAF_bin, y = EmpRsq,group=ref_panel, color=ref_panel))  + geom_line(size=1) + geom_point(shape=1)+ theme_classic(base_size=30) + coord_cartesian(ylim = c(0, 1)) + scale_x_discrete(guide = guide_axis(angle = 45))  + theme(text = element_text(size = 30,colour="black"))  + geom_hline(yintercept=0.8,  color = "dark grey") + geom_hline(yintercept=0.6, linetype="dashed", color = "dark grey")  + geom_vline(xintercept="(0.01,0.015]",  color = "dark grey",linetype="dashed") + geom_vline(xintercept="(0.05,0.1]",  color = "dark grey") #all font size
p + scale_color_manual(values=c("#D55E00", "#E69F00","#0072B2","#332288","#999999"),name="Reference Panel") + labs(x="Minor Allele Frequency", y="Average Empirical R2") +
  annotate('text', x = "[0,0.0005]", y = 0.45, label='"*"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.0005,0.001]", y = 0.5, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.001,0.002]", y = 0.6, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.002,0.005]", y = 0.7, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.005,0.01]", y = 0.75, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.01,0.015]", y = 0.8, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.015,0.02]", y = 0.85, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.02,0.035]", y = 0.85, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.035,0.05]", y = 0.9, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.05,0.1]", y = 0.9, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.1,0.2]", y = 0.95, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.2,0.3]", y = 0.95, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.3,0.4]", y = 0.95, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.4,0.5]", y = 0.95, label='"**"', parse=TRUE, color = "red",size=10) 
dev.off()
# expanded 1000G is the best or equivalently the best at all MAF bins

# one way ANOVA
merged_1 <- merged[which(merged$MAF_bin == "[0,0.0005]"),] # 18950 rows
merged_2 <- merged[which(merged$MAF_bin == "(0.0005,0.001]"),] # 25292 rows
merged_3 <- merged[which(merged$MAF_bin == "(0.001,0.002]"),] # 51051 rows
merged_4 <- merged[which(merged$MAF_bin == "(0.002,0.005]"),] # 123591 rows
merged_5 <- merged[which(merged$MAF_bin == "(0.005,0.01]"),] # 176304 rows
merged_6 <- merged[which(merged$MAF_bin == "(0.01,0.015]"),] # 130665 rows
merged_7 <- merged[which(merged$MAF_bin == "(0.015,0.02]"),] # 94074 rows
merged_8 <- merged[which(merged$MAF_bin == "(0.02,0.035]"),] # 195160 rows
merged_9 <- merged[which(merged$MAF_bin == "(0.035,0.05]"),] # 141582 rows
merged_10 <- merged[which(merged$MAF_bin == "(0.05,0.1]"),] # 306597 rows
merged_11 <- merged[which(merged$MAF_bin == "(0.1,0.2]"),] # 379952 rows
merged_12 <- merged[which(merged$MAF_bin == "(0.2,0.3]"),] # 291348 rows
merged_13 <- merged[which(merged$MAF_bin == "(0.3,0.4]"),] # 249210 rows
merged_14 <- merged[which(merged$MAF_bin == "(0.4,0.5]"),] # 228157 rows

dim(merged_1)
dim(merged_2)
dim(merged_3)
dim(merged_4)
dim(merged_5)
dim(merged_6)
dim(merged_7)
dim(merged_8)
dim(merged_9)
dim(merged_10)
dim(merged_11)
dim(merged_12)
dim(merged_13)
dim(merged_14)

table(merged$MAF_bin)

# Non-parametric alternative to one-way ANOVA test
kruskal.test(EmpRsq ~ ref_panel, data = merged_1) # P 0.07233 (*)
kruskal.test(EmpRsq ~ ref_panel, data = merged_2) # P < 2.2e-16 (**) 
kruskal.test(EmpRsq ~ ref_panel, data = merged_3) # P < 2.2e-16 (**)
kruskal.test(EmpRsq ~ ref_panel, data = merged_4) # P < 2.2e-16 (**)
kruskal.test(EmpRsq ~ ref_panel, data = merged_5) # P < 2.2e-16 (**)
kruskal.test(EmpRsq ~ ref_panel, data = merged_6) # P < 2.2e-16 (**)
kruskal.test(EmpRsq ~ ref_panel, data = merged_7) # P < 2.2e-16 (**)
kruskal.test(EmpRsq ~ ref_panel, data = merged_8) # P < 2.2e-16 (**)
kruskal.test(EmpRsq ~ ref_panel, data = merged_9) # P < 2.2e-16 (**)
kruskal.test(EmpRsq ~ ref_panel, data = merged_10) # P < 2.2e-16 (**)
kruskal.test(EmpRsq ~ ref_panel, data = merged_11) # P < 2.2e-16 (**)
kruskal.test(EmpRsq ~ ref_panel, data = merged_12) # P < 2.2e-16 (**)
kruskal.test(EmpRsq ~ ref_panel, data = merged_13) # P < 2.2e-16 (**)
kruskal.test(EmpRsq ~ ref_panel, data = merged_14) # P < 2.2e-16 (**)
# 0.05/14 = 0.003571 (** Bonferroni significant )
# 0.05 (* nominal significant)

# check for pairwise difference
# Pairwise t-tests with no assumption of equal variances
pairwise.t.test(merged_1$EmpRsq, merged_1$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) # TOPMed better than 1000G GRCh38, 1000G SAS GRCh37 and GenomeAsia

#                 1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     0.55313          -            -              -               
# Expanded_1000G   0.05669          0.18106      -              -               
# GenomeAsia_Pilot 0.07885          0.02009      0.00043        -               
# TOPMed           0.00064          0.00386      0.11184        1.1e-06  

pairwise.t.test(merged_2$EmpRsq, merged_2$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) # Expanded 1000G better than 1000G GRCh38, 1000G SAS GRCh37 and GenomeAsia

#                  1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     3.1e-06          -            -              -               
# Expanded_1000G   2.3e-12          0.018        -              -               
# GenomeAsia_Pilot 1.6e-06          < 2e-16      < 2e-16        -               
# TOPMed           4.8e-09          0.203        0.288          < 2e-16   

pairwise.t.test(merged_3$EmpRsq, merged_3$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) # Expanded 1000G better than all

#                  1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     1.1e-06          -            -              -               
# Expanded_1000G   < 2e-16          1.2e-10      -              -               
# GenomeAsia_Pilot < 2e-16          < 2e-16      < 2e-16        -               
# TOPMed           0.00062          0.14224      2.1e-15        < 2e-16  

pairwise.t.test(merged_4$EmpRsq, merged_4$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) # Expanded 1000G better than all

#                  1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     9.5e-16          -            -              -               
# Expanded_1000G   < 2e-16          4.6e-16      -              -               
# GenomeAsia_Pilot < 2e-16          < 2e-16      < 2e-16        -               
# TOPMed           0.94             < 2e-16      < 2e-16        < 2e-16      

pairwise.t.test(merged_5$EmpRsq, merged_5$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) # Expanded 1000G better than all

#                  1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     <2e-16           -            -              -               
# Expanded_1000G   <2e-16           <2e-16       -              -               
# GenomeAsia_Pilot <2e-16           <2e-16       <2e-16         -               
# TOPMed           0.15             <2e-16       <2e-16         <2e-16   

pairwise.t.test(merged_6$EmpRsq, merged_6$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) # Expanded 1000G better than all

#                  1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     < 2e-16          -            -              -               
# Expanded_1000G   < 2e-16          < 2e-16      -              -               
# GenomeAsia_Pilot < 2e-16          < 2e-16      < 2e-16        -               
# TOPMed           0.00039          5.9e-16      < 2e-16        < 2e-16     

pairwise.t.test(merged_7$EmpRsq, merged_7$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) # Expanded 1000G better than all


#                  1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     < 2e-16          -            -              -               
# Expanded_1000G   < 2e-16          4.2e-15      -              -               
# GenomeAsia_Pilot < 2e-16          < 2e-16      < 2e-16        -               
# TOPMed           0.001            1.5e-13      < 2e-16        < 2e-16    

pairwise.t.test(merged_8$EmpRsq, merged_8$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) # Expanded 1000G better than all

#                  1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     < 2e-16          -            -              -               
# Expanded_1000G   < 2e-16          < 2e-16      -              -               
# GenomeAsia_Pilot < 2e-16          < 2e-16      < 2e-16        -               
# TOPMed           2.6e-05          < 2e-16      < 2e-16        < 2e-16      

pairwise.t.test(merged_9$EmpRsq, merged_9$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) # Expanded 1000G better than all

#                  1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     < 2e-16          -            -              -               
# Expanded_1000G   < 2e-16          4.2e-14      -              -               
# GenomeAsia_Pilot 7.2e-12          < 2e-16      < 2e-16        -               
# TOPMed           0.0013           8.6e-09      < 2e-16        < 2e-16   

pairwise.t.test(merged_10$EmpRsq, merged_10$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) # Expanded 1000G better than all

#                  1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     < 2e-16          -            -              -               
# Expanded_1000G   < 2e-16          < 2e-16      -              -               
# GenomeAsia_Pilot 2.4e-14          < 2e-16      < 2e-16        -               
# TOPMed           0.012            6.8e-10      < 2e-16        < 2e-16         

pairwise.t.test(merged_11$EmpRsq, merged_11$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) # Expanded 1000G better than all

#                  1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     2.1e-12          -            -              -               
# Expanded_1000G   < 2e-16          < 2e-16      -              -               
# GenomeAsia_Pilot 1.4e-15          < 2e-16      < 2e-16        -               
# TOPMed           3.3e-05          0.0043       < 2e-16        < 2e-16      

pairwise.t.test(merged_12$EmpRsq, merged_12$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) # Expanded 1000G better than all

#                  1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     2.1e-06          -            -              -               
# Expanded_1000G   < 2e-16          9.1e-16      -              -               
# GenomeAsia_Pilot 3.4e-15          < 2e-16      < 2e-16        -               
# TOPMed           3.5e-06          0.9          < 2e-16        < 2e-16         

pairwise.t.test(merged_13$EmpRsq, merged_13$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) # Expanded 1000G better than all


#                  1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     4.3e-05          -            -              -               
# Expanded_1000G   < 2e-16          4.4e-14      -              -               
# GenomeAsia_Pilot 1.1e-10          < 2e-16      < 2e-16        -               
# TOPMed           7.6e-07          0.41         1.3e-11        < 2e-16     

pairwise.t.test(merged_14$EmpRsq, merged_14$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) # Expanded 1000G better than all

#                  1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     0.00014          -            -              -               
# Expanded_1000G   < 2e-16          4.9e-13      -              -               
# GenomeAsia_Pilot 1.6e-10          < 2e-16      < 2e-16        -               
# TOPMed           2.4e-06          0.38442      1.6e-10        < 2e-16       


#######   Empirical Rsq by MAF for well-imputed SNPs ###########

good_topmed <- all_topmed[which(all_topmed$Rsq >= 0.8),]
good_1000g <- all_1000g[which(all_1000g$Rsq >= 0.8),]
good_1000g_38 <- all_1000g_38[which(all_1000g_38$Rsq >= 0.8),]
good_1000g_37 <- all_1000g_37[which(all_1000g_37$Rsq >= 0.8),]
good_ga <- all_ga[which(all_ga$Rsq >= 0.8),]

dim(good_topmed) # [1] 485045     14
dim(good_1000g) # [1] 482451     14
dim(good_1000g_38) # [1] 485734     14
dim(good_1000g_37) # [1] 489612     14
dim(good_ga) # [1] 460411     14

good_merged <- rbind(good_topmed,good_1000g,good_1000g_38,good_1000g_37,good_ga) # [1] 2403253      14

table(good_merged$ref_panel)
summary(good_merged$MAF)

library(dplyr)
good_merged$MAF_bin <- case_when(good_merged$MAF >= 0 & good_merged$MAF <= 0.0005 ~ "[0,0.0005]", good_merged$MAF > 0.0005 & good_merged$MAF <= 0.001 ~ "(0.0005,0.001]", good_merged$MAF > 0.001 & good_merged$MAF <= 0.002 ~ "(0.001,0.002]", good_merged$MAF > 0.002 & good_merged$MAF <= 0.005 ~ "(0.002,0.005]",good_merged$MAF > 0.005 & good_merged$MAF <= 0.01 ~ "(0.005,0.01]",good_merged$MAF > 0.01 & good_merged$MAF <= 0.015 ~ "(0.01,0.015]",good_merged$MAF > 0.015 & good_merged$MAF <= 0.02 ~ "(0.015,0.02]", good_merged$MAF > 0.02 & good_merged$MAF <= 0.035 ~ "(0.02,0.035]", good_merged$MAF > 0.035 & good_merged$MAF <= 0.05 ~ "(0.035,0.05]", good_merged$MAF > 0.05 & good_merged$MAF <= 0.1 ~ "(0.05,0.1]", good_merged$MAF > 0.1 & good_merged$MAF <= 0.2  ~ "(0.1,0.2]", good_merged$MAF > 0.2 & good_merged$MAF <= 0.3  ~ "(0.2,0.3]", good_merged$MAF > 0.3 & good_merged$MAF <= 0.4  ~ "(0.3,0.4]", good_merged$MAF > 0.4 & good_merged$MAF <= 0.5  ~ "(0.4,0.5]")
table(good_merged$MAF_bin) # add up to 2403253

# install.packages("plotrix")
# library(plotrix)
good_emp_r2_maf <- aggregate(EmpRsq~ref_panel+MAF_bin, data=good_merged, mean)
# FUN = function(x) c(mean = mean(x), se = std.error(x))
good_emp_r2_maf$ref_panel <- factor(good_emp_r2_maf$ref_panel,  levels=c("TOPMed","Expanded_1000G","1000G_GRCh38", "1000G_GRCh37_SAS", "GenomeAsia_Pilot"))
good_emp_r2_maf$MAF_bin <- factor(good_emp_r2_maf$MAF_bin,  levels=c("[0,0.0005]","(0.0005,0.001]","(0.001,0.002]","(0.002,0.005]","(0.005,0.01]","(0.01,0.015]","(0.015,0.02]","(0.02,0.035]","(0.035,0.05]","(0.05,0.1]","(0.1,0.2]","(0.2,0.3]","(0.3,0.4]","(0.4,0.5]"))

library(ggplot2)
library(dplyr)
pdf('emp_r2_maf_well_imputed_snp.pdf',height=10,width=15)
p <- ggplot(good_emp_r2_maf, aes(x = MAF_bin, y = EmpRsq,group=ref_panel, color=ref_panel))  + geom_line(size=1) + geom_point(shape=1)+ theme_classic(base_size=30) + coord_cartesian(ylim = c(0, 1)) + scale_x_discrete(guide = guide_axis(angle = 45))  + theme(text = element_text(size = 30,colour="black"))  + geom_hline(yintercept=0.8,  color = "dark grey") + geom_hline(yintercept=0.6, linetype="dashed", color = "dark grey")  + geom_vline(xintercept="(0.01,0.015]",  color = "dark grey",linetype="dashed") + geom_vline(xintercept="(0.05,0.1]",  color = "dark grey") #all font size
p + scale_color_manual(values=c("#D55E00", "#E69F00","#0072B2","#332288","#999999"),name="Reference Panel") + labs(x="Minor Allele Frequency", y="Average Empirical R2") +
  annotate('text', x = "[0,0.0005]", y = 0.53, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.0005,0.001]", y = 0.55, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.001,0.002]", y = 0.65, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.002,0.005]", y = 0.7, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.005,0.01]", y = 0.75, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.01,0.015]", y = 0.8, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.015,0.02]", y = 0.85, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.02,0.035]", y = 0.87, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.035,0.05]", y = 0.9, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.05,0.1]", y = 0.92, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.1,0.2]", y = 0.95, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.2,0.3]", y = 0.95, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.3,0.4]", y = 0.95, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.4,0.5]", y = 0.95, label='"**"', parse=TRUE, color = "red",size=10) 
dev.off()

# one way ANOVA
merged_1 <- good_merged[which(good_merged$MAF_bin == "[0,0.0005]"),] # 14010 rows
merged_2 <- good_merged[which(good_merged$MAF_bin == "(0.0005,0.001]"),] # 23749 rows
merged_3 <- good_merged[which(good_merged$MAF_bin == "(0.001,0.002]"),] # 49851 rows
merged_4 <- good_merged[which(good_merged$MAF_bin == "(0.002,0.005]"),] # 122978 rows
merged_5 <- good_merged[which(good_merged$MAF_bin == "(0.005,0.01]"),] # 176192 rows
merged_6 <- good_merged[which(good_merged$MAF_bin == "(0.01,0.015]"),] # 130640 rows
merged_7 <- good_merged[which(good_merged$MAF_bin == "(0.015,0.02]"),] # 94057 rows
merged_8 <- good_merged[which(good_merged$MAF_bin == "(0.02,0.035]"),] # 195125 rows
merged_9 <- good_merged[which(good_merged$MAF_bin == "(0.035,0.05]"),] # 141565 rows
merged_10 <- good_merged[which(good_merged$MAF_bin == "(0.05,0.1]"),] # 306562 rows
merged_11 <- good_merged[which(good_merged$MAF_bin == "(0.1,0.2]"),] # 379898 rows
merged_12 <- good_merged[which(good_merged$MAF_bin == "(0.2,0.3]"),] # 291318 rows
merged_13 <- good_merged[which(good_merged$MAF_bin == "(0.3,0.4]"),] # 249186 rows
merged_14 <- good_merged[which(good_merged$MAF_bin == "(0.4,0.5]"),] # 228122 rows

dim(merged_1)
dim(merged_2)
dim(merged_3)
dim(merged_4)
dim(merged_5)
dim(merged_6)
dim(merged_7)
dim(merged_8)
dim(merged_9)
dim(merged_10)
dim(merged_11)
dim(merged_12)
dim(merged_13)
dim(merged_14)

table(good_merged$MAF_bin)

# Non-parametric alternative to one-way ANOVA test
kruskal.test(EmpRsq ~ ref_panel, data = merged_1) # P 1.169e-15 (**)
kruskal.test(EmpRsq ~ ref_panel, data = merged_2) # P < 2.2e-16 (**) 
kruskal.test(EmpRsq ~ ref_panel, data = merged_3) # P < 2.2e-16 (**)
kruskal.test(EmpRsq ~ ref_panel, data = merged_4) # P < 2.2e-16 (**)
kruskal.test(EmpRsq ~ ref_panel, data = merged_5) # P < 2.2e-16 (**)
kruskal.test(EmpRsq ~ ref_panel, data = merged_6) # P < 2.2e-16 (**)
kruskal.test(EmpRsq ~ ref_panel, data = merged_7) # P < 2.2e-16 (**)
kruskal.test(EmpRsq ~ ref_panel, data = merged_8) # P < 2.2e-16 (**)
kruskal.test(EmpRsq ~ ref_panel, data = merged_9) # P < 2.2e-16 (**)
kruskal.test(EmpRsq ~ ref_panel, data = merged_10) # P < 2.2e-16 (**)
kruskal.test(EmpRsq ~ ref_panel, data = merged_11) # P < 2.2e-16 (**)
kruskal.test(EmpRsq ~ ref_panel, data = merged_12) # P < 2.2e-16 (**)
kruskal.test(EmpRsq ~ ref_panel, data = merged_13) # P < 2.2e-16 (**)
kruskal.test(EmpRsq ~ ref_panel, data = merged_14) # P < 2.2e-16 (**)
# 0.05/14 = 0.003571 (** Bonferroni significant )
# 0.05 (* nominal significant)


#################################################################
######          Empirical R2 vs. Estimated R2         ###########
#################################################################


# in R
cor(all_topmed$Rsq, all_topmed$EmpRsq) # 0.5240815
cor(all_1000g$Rsq, all_1000g$EmpRsq) # 0.5387407
cor(all_1000g_38$Rsq, all_1000g_38$EmpRsq) # 0.5370734
cor(all_1000g_37$Rsq, all_1000g_37$EmpRsq) # 0.4908013
cor(all_ga$Rsq, all_ga$EmpRsq) # 0.4337825


plot(all_topmed$Rsq, all_topmed$EmpRsq , xlab = "Estimated R Square", ylab = "Empirical R Square", pch=20, frame = FALSE)
# Add regression line
abline(lm(EmpRsq ~ Rsq , data = all_topmed), col = "red", lwd=2)

plot(all_1000g$Rsq, all_1000g$EmpRsq , xlab = "Estimated R Square", ylab = "Empirical R Square", pch=20, frame = FALSE)
# Add regression line
abline(lm(EmpRsq ~ Rsq , data = all_1000g), col = "red", lwd=2)

plot(all_1000g_38$Rsq, all_1000g_38$EmpRsq , xlab = "Estimated R Square", ylab = "Empirical R Square", pch=20, frame = FALSE)
# Add regression line
abline(lm(EmpRsq ~ Rsq , data = all_1000g_38), col = "red", lwd=2)

plot(all_1000g_37$Rsq, all_1000g_37$EmpRsq , xlab = "Estimated R Square", ylab = "Empirical R Square", pch=20, frame = FALSE)
# Add regression line
abline(lm(EmpRsq ~ Rsq , data = all_1000g_37), col = "red", lwd=2)

plot(all_ga$Rsq, all_ga$EmpRsq , xlab = "Estimated R Square", ylab = "Empirical R Square", pch=20, frame = FALSE)
# Add regression line
abline(lm(EmpRsq ~ Rsq , data = all_ga), col = "red", lwd=2)



####################################################################################
####             plot number of SNPs and well-imputed SNPs (rare, common)       ####
####################################################################################

# parallel running in different folders to save time 


#!/bin/bash
#BSUB -J snp_count_rare_common # Job name
#BSUB -n 1
#BSUB -P acc_psychgen # allocation account
#BSUB -q premium # queue
#BSUB -W 24:00 # walltime in HH:MM
#BSUB -R himem
#BSUB -R rusage[mem=120G] # 120 GB of memory 
#BSUB -R span[hosts=1] # all cores from one node
#BSUB -o snp_count_rare_common.stdout # output log (%J : JobID)
#BSUB -eo snp_count_rare_common.stderr # error log
#BSUB -L /bin/bash # Initialize the execution environment


module load R

Rscript snp_count_rare_common.r

#!/usr/bin/env Rscript

options(echo=TRUE)
print("In R now!")

## SCRIPT STARTS BELOW

library(data.table)
dat <- as.data.frame(fread("all.chr.info", select = c("SNP","MAF","Rsq","Genotyped")))
dat <- dat[which(dat$Genotyped!="Typed_Only"),]

rare <- dat[which(dat$MAF < 0.01),]
common <- dat[which(dat$MAF >= 0.01 ),]

dim(dat)
dim(rare)
dim(common)
# rare + common add up to total

dim(rare[which(rare$Rsq >=0.8),])
dim(common[which(common$Rsq >=0.8),])

bsub < snp_count_rare_common.sh

## for meta-imputation

# All SNPs
awk '($2>=0.01){print $0}' all.chr.info|wc -l
awk '($2<0.01){print $0}' all.chr.info|wc -l
wc -l all.chr.info

# meta imputation: rare, common: 286871988 9481424
# TOPMed  rare, common:  282753704 9382758
# 1000G  rare, common: 39617700 9273539
# 1000G GRCh38: rare, common: 35213874 9993538
# 1000G GRCh37 SAS:  rare, common: 36691953 10407598
# GA  rare, common: 13803012 7691636

# high quality SNPs

awk '($3>=0.8 && $2>=0.01){print $0}' all.chr.info|wc -l
awk '($3>=0.8 && $2<0.01){print $0}' all.chr.info|wc -l
awk '($3>=0.8){print $0}' all.chr.info|wc -l

# meta imputation well-imputed SNPs: rare, common well-imputed SNPs: 4747618 7932142
# TOPMed rare, common well-imputed SNPs: 4202423 7370972
# 1000G rare, common well-imputed SNPs: 2854147 7851403
# 1000G GRCh38: rare, common well-imputed SNPs: 2346512 8019371
# 1000G GRCh37 SAS: rare, common well-imputed SNPs: 1426862 7076181
# GA rare, common well-imputed SNPs: 668162 5227987

#### merge results from different panels together
# for example, the # of common variants in GA is obtained from: awk '($8!="Typed_Only" && $5 >=0.01 && NR>1){print $1}' all.chr.info|wc -l
snp_count <- read.table("snp_count_new",header=T)
table(snp_count$ref_panel)
table(snp_count$maf_cat)
table(snp_count$snp_cat)

snp_count$ref_panel <- factor(snp_count$ref_panel,  levels=c("TOPMed_1000G","TOPMed","Expanded_1000G","1000G_GRCh38", "1000G_GRCh37_SAS", "GenomeAsia_Pilot"))
snp_count$maf_cat <- factor(snp_count$maf_cat,  levels=c("common","rare"))

snp_count_all <- snp_count[which(snp_count$snp_cat=="all_snps"),]
snp_count_well <- snp_count[which(snp_count$snp_cat=="well_imputed_snps"),]

# options(scipen = 999)

# pdf('~/www/test.pdf',width=10,height=5)

library(ggplot2)
pdf('all_snp_count.pdf',height=10,width=15)
p <- ggplot(data=snp_count_all, aes(x=maf_cat, y=count, fill=ref_panel)) + xlab("Variant Type") + ylab("Count") + geom_bar(stat="identity", position=position_dodge()) +  geom_text(aes(label=scales::label_number_si(accuracy=0.01)(count)),size=4, position=position_dodge(width=0.9), vjust=-0.25) + theme_classic(base_size=30)  + scale_x_discrete(labels=c("Common","Rare")) + scale_y_continuous(labels=scales::label_number_si())+theme(axis.text = element_text(size = 30,colour="black")) 
p + scale_fill_manual(values=c("#661100", "#D55E00", "#E69F00","#0072B2","#332288","#999999"),name="Reference Panel",labels=c("TOPMed_1000G","TOPMed","Expanded_1000G","1000G_GRCh38", "1000G_GRCh37_SAS", "GenomeAsia_Pilot")) + theme(legend.position = "none" ,plot.title = element_text(hjust = 0.5)) + labs(title="All SNPs")
dev.off()

pdf('well_imputed_snp_count.pdf',height=10,width=15)
p <- ggplot(data=snp_count_well, aes(x=maf_cat, y=count, fill=ref_panel)) + xlab("Variant Type") + ylab("Count") + geom_bar(stat="identity", position=position_dodge()) +  geom_text(aes(label=scales::label_number_si(accuracy=0.01)(count)),size=4, position=position_dodge(width=0.9), vjust=-0.25) + theme_classic(base_size=30)  + scale_x_discrete(labels=c("Common","Rare")) + scale_y_continuous(labels=scales::label_number_si())+theme(axis.text = element_text(size = 30,colour="black")) 
p + scale_fill_manual(values=c("#661100", "#D55E00", "#E69F00","#0072B2","#332288","#999999"),name="Reference Panel",labels=c("TOPMed_1000G","TOPMed","Expanded_1000G","1000G_GRCh38", "1000G_GRCh37_SAS", "GenomeAsia_Pilot")) + theme(legend.position = "none",plot.title = element_text(hjust = 0.5)) + labs(title="Well-imputed SNPs")
dev.off()


########################################################################################
####             plot number of imputed SNPs with MAF <= 0.0005 across panels       ####
########################################################################################

# meta-imputation
awk '( $2 <= 0.0005 && NR>1 && $3 < 0.8){print $1}' all.chr.info|wc -l # 269,143,511

# TOPMed
awk '($8!="Typed_Only" && $5 <= 0.0005 && NR>1 && $7 < 0.8){print $1}' all.chr.info|wc -l # 267,025,690

# ex1KG
awk '($8!="Typed_Only" && $5 <= 0.0005 && NR>1 && $7 < 0.8){print $1}' all.chr.info|wc -l # 25,584,975

# 1KG38
awk '($8!="Typed_Only" && $5 <= 0.0005 && NR>1 && $7 < 0.8){print $1}' all.chr.info|wc -l # 22,211,916

# 1KG37
awk '($8!="Typed_Only" && $5 <= 0.0005 && NR>1 && $7 < 0.8){print $1}' all.chr.info|wc -l # 21,866,257

# GA
awk '($8!="Typed_Only" && $5 <= 0.0005 && NR>1 && $7 < 0.8){print $1}' all.chr.info|wc -l # 5,370,126

snp_count <- read.table("snp_count_maf_0.0005",header=T)

snp_count$ref_panel <- factor(snp_count$ref_panel,  levels=c("TOPMed_1000G","TOPMed","Expanded_1000G","1000G_GRCh38", "1000G_GRCh37_SAS", "GenomeAsia_Pilot"))

library(ggplot2)
pdf('poorly_imputed_snp_count_MAF_le_0.0005.pdf',height=10,width=20)
p <- ggplot(data=snp_count, aes(x=ref_panel, y=count, fill=ref_panel)) + xlab("Reference Panel") + ylab("Count")  + geom_bar(stat="identity", position=position_dodge()) +  geom_text(aes(label=scales::label_number_si(accuracy=0.01)(count)),size=7, position=position_dodge(width=0.9), vjust=-0.25) + theme_classic(base_size=20)  + scale_x_discrete() + scale_y_continuous(labels=scales::label_number_si())+theme(axis.text = element_text(size = 20,colour="black")) 
p + scale_fill_manual(values=c("#661100", "#D55E00", "#E69F00","#0072B2","#332288","#999999"))  + theme(legend.position = "none",plot.title = element_text(hjust = 0.5)) + labs(title="Poorly imputed SNPs")
dev.off()

###################################################################################################################
####             plot number of SNPs and well-imputed SNPs by SNVs and indels (ultra rare, rare, common)       ####
###################################################################################################################

#!/bin/bash
#BSUB -J snv_indel_count # Job name
#BSUB -n 1
#BSUB -P acc_psychgen # allocation account
#BSUB -q premium # queue
#BSUB -W 48:00 # walltime in HH:MM
#BSUB -R himem
#BSUB -R rusage[mem=96G] # 120 GB of memory 
#BSUB -R span[hosts=1] # all cores from one node
#BSUB -o %J.stdout  # output log (%J : JobID)
#BSUB -eo %J.stderr # error log
#BSUB -L /bin/bash # Initialize the execution environment


module load R

Rscript snv_indel_count.r

bsub < snv_indel_count.sh

#!/usr/bin/env Rscript

options(echo=TRUE)
print("In R now!")

## SCRIPT STARTS BELOW

library(data.table)
info <- fread("all.chr.info", header = T)

require(stringr)
info$ref_allele_number <- str_count(info$"REF(0)")
info$alt_allele_number <- str_count(info$"ALT(1)")

summary(info$ref_allele_number)
summary(info$alt_allele_number)

info$indel <- ifelse(info$ref_allele_number==1 & info$alt_allele_number==1, "SNV","indel")

table(info$indel)

info$Rsq <- as.numeric(info$Rsq)

library(dplyr)
info$MAF_cat <- case_when(info$MAF <= 0.001 ~ "ultra_rare",info$MAF > 0.001 & info$MAF <= 0.01 ~ "rare",  info$MAF >0.01 ~ "common")

imputed <- info[which(info$Genotyped=="Imputed"),]
summary(imputed$Rsq)

# imputed$high_quality_SNP <- ifelse(imputed$Rsq >= 0.8,1,0)

snp_type <- data.frame(matrix(ncol = 4, nrow = 8))
colnames(snp_type) <- c('indel', 'MAF_cat','counts', 'counts_high_r2')

snp_type$indel[c(1,3,5,7)] <- "SNV"
snp_type$indel[c(2,4,6,8)] <- "indel"

snp_type$MAF_cat[c(1,2)] <- "all"
snp_type$MAF_cat[c(3,4)] <- "ultra_rare"
snp_type$MAF_cat[c(5,6)] <- "rare"
snp_type$MAF_cat[c(7,8)] <- "common"

snp_type$counts[1] <- dim(imputed[which(imputed$indel=="SNV"),])[1]
snp_type$counts[2] <- dim(imputed[which(imputed$indel=="indel"),])[1]
snp_type$counts[3] <- dim(imputed[which(imputed$indel=="SNV" & imputed$MAF_cat =="ultra_rare"),])[1]
snp_type$counts[4] <- dim(imputed[which(imputed$indel=="indel" & imputed$MAF_cat =="ultra_rare"),])[1]
snp_type$counts[5] <- dim(imputed[which(imputed$indel=="SNV" & imputed$MAF_cat =="rare"),])[1]
snp_type$counts[6] <- dim(imputed[which(imputed$indel=="indel" & imputed$MAF_cat =="rare"),])[1]
snp_type$counts[7] <- dim(imputed[which(imputed$indel=="SNV" & imputed$MAF_cat =="common"),])[1]
snp_type$counts[8] <- dim(imputed[which(imputed$indel=="indel" & imputed$MAF_cat =="common"),])[1]


snp_type$counts_high_r2[1] <- dim(imputed[which(imputed$indel=="SNV" & imputed$Rsq >= 0.8),])[1]
snp_type$counts_high_r2[2] <- dim(imputed[which(imputed$indel=="indel" &  imputed$Rsq >= 0.8),])[1]
snp_type$counts_high_r2[3] <- dim(imputed[which(imputed$indel=="SNV" & imputed$MAF_cat =="ultra_rare" & imputed$Rsq >= 0.8),])[1]
snp_type$counts_high_r2[4] <- dim(imputed[which(imputed$indel=="indel" & imputed$MAF_cat =="ultra_rare" & imputed$Rsq >= 0.8),])[1]
snp_type$counts_high_r2[5] <- dim(imputed[which(imputed$indel=="SNV" & imputed$MAF_cat =="rare" &  imputed$Rsq >= 0.8),])[1]
snp_type$counts_high_r2[6] <- dim(imputed[which(imputed$indel=="indel" & imputed$MAF_cat =="rare" &  imputed$Rsq >= 0.8),])[1]
snp_type$counts_high_r2[7] <- dim(imputed[which(imputed$indel=="SNV" & imputed$MAF_cat =="common" &  imputed$Rsq >= 0.8),])[1]
snp_type$counts_high_r2[8] <- dim(imputed[which(imputed$indel=="indel" & imputed$MAF_cat =="common" &  imputed$Rsq >= 0.8),])[1]

dim(imputed[which( imputed$high_quality_SNP ==1),])

snp_type$prop_high_r2 <- snp_type$counts_high_r2 / snp_type$counts * 100

write.table(snp_type,"snp_type_indel.txt", quote=F, row.names=F)


#######################################################################################################
####             plot number of SNPs and well-imputed SNPs by SNVs and indels (rare, common)       ####
#######################################################################################################

#!/bin/bash
#BSUB -J snv_indel_count_new # Job name
#BSUB -n 1
#BSUB -P acc_psychgen # allocation account
#BSUB -q premium # queue
#BSUB -W 48:00 # walltime in HH:MM
#BSUB -R himem
#BSUB -R rusage[mem=96G] # 120 GB of memory 
#BSUB -R span[hosts=1] # all cores from one node
#BSUB -o %J.stdout  # output log (%J : JobID)
#BSUB -eo %J.stderr # error log
#BSUB -L /bin/bash # Initialize the execution environment


module load R

Rscript snv_indel_count_new.r

bsub < snv_indel_count_new.sh

#!/usr/bin/env Rscript

options(echo=TRUE)
print("In R now!")

## SCRIPT STARTS BELOW

library(data.table)
info <- fread("all.chr.info", header = T)

require(stringr)
info$ref_allele_number <- str_count(info$"REF(0)")
info$alt_allele_number <- str_count(info$"ALT(1)")

summary(info$ref_allele_number)
summary(info$alt_allele_number)

info$indel <- ifelse(info$ref_allele_number==1 & info$alt_allele_number==1, "SNV","indel")

table(info$indel)

info$Rsq <- as.numeric(info$Rsq)

library(dplyr)
info$MAF_cat <- case_when( info$MAF < 0.01 ~ "rare",  info$MAF >= 0.01 ~ "common")

imputed <- info[which(info$Genotyped!="Typed_Only"),]
summary(imputed$Rsq)

# imputed$high_quality_SNP <- ifelse(imputed$Rsq >= 0.8,1,0)

snp_type <- data.frame(matrix(ncol = 4, nrow = 6))
colnames(snp_type) <- c('indel', 'MAF_cat','counts', 'counts_high_r2')

snp_type$indel[c(1,3,5)] <- "SNV"
snp_type$indel[c(2,4,6)] <- "indel"

snp_type$MAF_cat[c(1,2)] <- "all"
snp_type$MAF_cat[c(3,4)] <- "rare"
snp_type$MAF_cat[c(5,6)] <- "common"

snp_type$counts[1] <- dim(imputed[which(imputed$indel=="SNV"),])[1]
snp_type$counts[2] <- dim(imputed[which(imputed$indel=="indel"),])[1]
snp_type$counts[3] <- dim(imputed[which(imputed$indel=="SNV" & imputed$MAF_cat =="rare"),])[1]
snp_type$counts[4] <- dim(imputed[which(imputed$indel=="indel" & imputed$MAF_cat =="rare"),])[1]
snp_type$counts[5] <- dim(imputed[which(imputed$indel=="SNV" & imputed$MAF_cat =="common"),])[1]
snp_type$counts[6] <- dim(imputed[which(imputed$indel=="indel" & imputed$MAF_cat =="common"),])[1]


snp_type$counts_high_r2[1] <- dim(imputed[which(imputed$indel=="SNV" & imputed$Rsq >= 0.8),])[1]
snp_type$counts_high_r2[2] <- dim(imputed[which(imputed$indel=="indel" &  imputed$Rsq >= 0.8),])[1]
snp_type$counts_high_r2[3] <- dim(imputed[which(imputed$indel=="SNV" & imputed$MAF_cat =="rare" &  imputed$Rsq >= 0.8),])[1]
snp_type$counts_high_r2[4] <- dim(imputed[which(imputed$indel=="indel" & imputed$MAF_cat =="rare" &  imputed$Rsq >= 0.8),])[1]
snp_type$counts_high_r2[5] <- dim(imputed[which(imputed$indel=="SNV" & imputed$MAF_cat =="common" &  imputed$Rsq >= 0.8),])[1]
snp_type$counts_high_r2[6] <- dim(imputed[which(imputed$indel=="indel" & imputed$MAF_cat =="common" &  imputed$Rsq >= 0.8),])[1]


snp_type$prop_high_r2 <- snp_type$counts_high_r2 / snp_type$counts * 100

write.table(snp_type,"snp_type_indel_new.txt", quote=F, row.names=F)

########       for meta-imputation       ############

#!/bin/bash
#BSUB -J snv_indel_count_new # Job name
#BSUB -n 1
#BSUB -P acc_psychgen # allocation account
#BSUB -q premium # queue
#BSUB -W 48:00 # walltime in HH:MM
#BSUB -R himem
#BSUB -R rusage[mem=150G] # 120 GB of memory 
#BSUB -R span[hosts=1] # all cores from one node
#BSUB -o %J.stdout  # output log (%J : JobID)
#BSUB -eo %J.stderr # error log
#BSUB -L /bin/bash # Initialize the execution environment


module load R

Rscript snv_indel_count_new.r

bsub < snv_indel_count_new.sh

#!/usr/bin/env Rscript

options(echo=TRUE)
print("In R now!")

## SCRIPT STARTS BELOW

library(data.table)
info <- fread("all.chr.info", header = T)

info <- as.data.frame(info)
library(stringr)
#split 'player' column using '_' as the separator
info[c('chr', 'pos','a1','a2')] <- str_split_fixed(info$SNP, ':', 4)

info$ref_allele_number <- str_count(info$"a1")
info$alt_allele_number <- str_count(info$"a2")

summary(info$ref_allele_number)
summary(info$alt_allele_number)

info$indel <- ifelse(info$ref_allele_number==1 & info$alt_allele_number==1, "SNV","indel")

table(info$indel)

info$Rsq <- as.numeric(info$Rsq)

library(dplyr)
info$MAF_cat <- case_when( info$MAF < 0.01 ~ "rare",  info$MAF >= 0.01 ~ "common")

snp_type <- data.frame(matrix(ncol = 4, nrow = 6))
colnames(snp_type) <- c('indel', 'MAF_cat','counts', 'counts_high_r2')

snp_type$indel[c(1,3,5)] <- "SNV"
snp_type$indel[c(2,4,6)] <- "indel"

snp_type$MAF_cat[c(1,2)] <- "all"
snp_type$MAF_cat[c(3,4)] <- "rare"
snp_type$MAF_cat[c(5,6)] <- "common"

snp_type$counts[1] <- dim(info[which(info$indel=="SNV"),])[1]
snp_type$counts[2] <- dim(info[which(info$indel=="indel"),])[1]
snp_type$counts[3] <- dim(info[which(info$indel=="SNV" & info$MAF_cat =="rare"),])[1]
snp_type$counts[4] <- dim(info[which(info$indel=="indel" & info$MAF_cat =="rare"),])[1]
snp_type$counts[5] <- dim(info[which(info$indel=="SNV" & info$MAF_cat =="common"),])[1]
snp_type$counts[6] <- dim(info[which(info$indel=="indel" & info$MAF_cat =="common"),])[1]


snp_type$counts_high_r2[1] <- dim(info[which(info$indel=="SNV" & info$Rsq >= 0.8),])[1]
snp_type$counts_high_r2[2] <- dim(info[which(info$indel=="indel" &  info$Rsq >= 0.8),])[1]
snp_type$counts_high_r2[3] <- dim(info[which(info$indel=="SNV" & info$MAF_cat =="rare" &  info$Rsq >= 0.8),])[1]
snp_type$counts_high_r2[4] <- dim(info[which(info$indel=="indel" & info$MAF_cat =="rare" &  info$Rsq >= 0.8),])[1]
snp_type$counts_high_r2[5] <- dim(info[which(info$indel=="SNV" & info$MAF_cat =="common" &  info$Rsq >= 0.8),])[1]
snp_type$counts_high_r2[6] <- dim(info[which(info$indel=="indel" & info$MAF_cat =="common" &  info$Rsq >= 0.8),])[1]

snp_type$prop_high_r2 <- snp_type$counts_high_r2 / snp_type$counts * 100

write.table(snp_type,"snp_type_indel_new.txt", quote=F, row.names=F)

# change rare to Rare, common to Common, all to All

### now it's time to plot!

all_topmed <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/snp_type_indel_new.txt", header=T)
all_1000g <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/snp_type_indel_new.txt", header=T) 
all_ga <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_GA/snp_type_indel_new.txt", header=T) 
all_1000g_37 <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_37/snp_type_indel_new.txt", header=T) 
all_1000g_38 <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_38/snp_type_indel_new.txt", header=T) 
all_meta_imp <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_topmed_1000G/snp_type_indel_new.txt", header=T) 

dim(all_topmed) # [1] 6 5
dim(all_1000g) # [1] 6 5
dim(all_ga) # [1] 6 5
dim(all_1000g_37) # [1] 6 5
dim(all_1000g_38) # [1] 6 5
dim(all_meta_imp) # [1] 6 5

all_meta_imp$ref_panel <- "TOPMed_1000G"
all_topmed$ref_panel <- "TOPMed"
all_1000g$ref_panel <- "Expanded_1000G"
all_1000g_38$ref_panel <- "1000G_GRCh38"
all_1000g_37$ref_panel <- "1000G_GRCh37_SAS"
all_ga$ref_panel <- "GenomeAsia_Pilot"

snp_count <- rbind(all_meta_imp, all_topmed, all_1000g, all_1000g_38, all_1000g_37, all_ga) # [1] 36  6

snp_count_all <- snp_count[,c(1,2,3,6)]
snp_count_well <- snp_count[,c(1,2,4,6)]

colnames(snp_count_well)[colnames(snp_count_well)=="counts_high_r2"] <- "counts"

snp_count_all$snp_cat <- "all"
snp_count_well$snp_cat <- "well_imp"

snp_count_merge <- rbind(snp_count_all,snp_count_well)
snp_count_merge$ref_panel <- factor(snp_count_merge$ref_panel,  levels=c("TOPMed_1000G","TOPMed","Expanded_1000G","1000G_GRCh38", "1000G_GRCh37_SAS", "GenomeAsia_Pilot"))
snp_count_merge$MAF_cat <- factor(snp_count_merge$MAF_cat,  levels=c("All","Common","Rare"))
snp_count_merge$snp_cat <- factor(snp_count_merge$snp_cat,  levels=c("all","well_imp"))
snp_count_merge$indel <-  factor(snp_count_merge$indel,  levels=c("SNV","indel"))

snp_count_all <- snp_count_merge[which(snp_count_merge$MAF_cat=="All"),]
snp_count_maf <- snp_count_merge[which(snp_count_merge$MAF_cat !="All"),]

# Separated by all SNPs and well-imputed SNPs for rare and common variants
snp_count_all_maf <- snp_count_maf[which(snp_count_maf$snp_cat =="all"),]
snp_count_well_maf <- snp_count_maf[which(snp_count_maf$snp_cat =="well_imp"),]

library(ggplot2)
pdf('snv_indel_all_snp_count.pdf',height=10,width=20)
p <- ggplot(data=snp_count_all_maf, aes(x=interaction(MAF_cat,indel), y=counts, fill=ref_panel)) + xlab("Variant Type") + ylab("Count") + geom_bar(stat="identity", position=position_dodge()) +  geom_text(aes(label=scales::label_number_si(accuracy=0.1)(counts)),size=3, position=position_dodge(width=0.9), vjust=-0.25) + theme_classic(base_size=25)  + scale_x_discrete() + scale_y_continuous(labels=scales::label_number_si())+theme(axis.text = element_text(size = 25,colour="black"))
p + scale_fill_manual(values=c("#661100", "#D55E00", "#E69F00","#0072B2","#332288","#999999"),name="Reference Panel",labels=c("TOPMed_1000G","TOPMed","Expanded_1000G","1000G_GRCh38", "1000G_GRCh37_SAS", "GenomeAsia_Pilot")) + theme(plot.title = element_text(hjust = 0.5)) + labs(title="All SNPs")
dev.off()

pdf('snv_indel_well_imputed_snp_count.pdf',height=10,width=20)
p <- ggplot(data=snp_count_well_maf, aes(x=interaction(MAF_cat,indel), y=counts, fill=ref_panel)) + xlab("Variant Type") + ylab("Count") + geom_bar(stat="identity", position=position_dodge()) +  geom_text(aes(label=scales::label_number_si(accuracy=0.1)(counts)),size=4, position=position_dodge(width=0.9), vjust=-0.25) + theme_classic(base_size=25)  + scale_x_discrete() + scale_y_continuous(labels=scales::label_number_si())+theme(axis.text = element_text(size = 25,colour="black"))
p + scale_fill_manual(values=c("#661100", "#D55E00", "#E69F00","#0072B2","#332288","#999999"),name="Reference Panel",labels=c("TOPMed_1000G","TOPMed","Expanded_1000G","1000G_GRCh38", "1000G_GRCh37_SAS", "GenomeAsia_Pilot")) + theme(legend.position="none",plot.title = element_text(hjust = 0.5)) + labs(title="Well-imputed SNPs")
dev.off()



####################################################################################################################
####             plot number of SNPs and well-imputed SNPs by multi-allelic vs. bi-allelic (rare, common)       ####
####################################################################################################################

# create an info file with chr:pos
sed 's/\:/\t/g' all.chr.info |awk '(NR>1){print $2":"$3"\t"$1":"$2":"$3":"$4":"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17}' |sed '1i pos\tSNP\tREF(0)\tALT(1)\tALT_Frq\tMAF\tAvgCall\tRsq\tGenotyped\tLooRsq\tEmpR\tEmpRsq\tDose0\tDose1' > all.chr.pos.info 


#!/bin/bash
#BSUB -J multi_allele_bi_allele_count # Job name
#BSUB -n 1
#BSUB -P acc_psychgen # allocation account
#BSUB -q premium # queue
#bsub -n 1
#BSUB -W 48:00 # walltime in HH:MM
#BSUB -R himem
#BSUB -R rusage[mem=120G] # 120 GB of memory 
#BSUB -R span[hosts=1] # all cores from one node
#BSUB -o %J.stdout # output log (%J : JobID)
#BSUB -eo %J.stderr # error log
#BSUB -L /bin/bash # Initialize the execution environment


module load R

Rscript multi_allele_bi_allele_count.r

#!/usr/bin/env Rscript

options(echo=TRUE)
print("In R now!")

## SCRIPT STARTS BELOW

library(data.table)
info <- fread("all.chr.pos.info",header=T)
info <- as.data.frame(info)
dup_pos <- info[duplicated(info$pos),"pos"]
# if the SNP ID appeared 3 times, then the 2nd and 3rd will be both included in the dup_pos file

info$dup <- ifelse(info$pos %in% dup_pos, 1,0)

info$Rsq <- as.numeric(info$Rsq)

imputed <- info[which(info$Genotyped!="Typed_Only"),]


snp_type <- data.frame(matrix(ncol = 4, nrow = 6))
colnames(snp_type) <- c('multi_allele', 'MAF_cat','counts', 'counts_high_r2')

snp_type$multi_allele[c(1,3,5)] <- "multi"
snp_type$multi_allele[c(2,4,6)] <- "bi"

snp_type$MAF_cat[c(1,2)] <- "all"
snp_type$MAF_cat[c(3,4)] <- "rare"
snp_type$MAF_cat[c(5,6)] <- "common"

snp_type$counts[1] <- dim(imputed[which(imputed$dup==1),])[1]
snp_type$counts[2] <- dim(imputed[which(imputed$dup==0),])[1]
snp_type$counts[3] <- dim(imputed[which(imputed$dup==1 & imputed$MAF < 0.01),])[1]
snp_type$counts[4] <- dim(imputed[which(imputed$dup==0 & imputed$MAF < 0.01),])[1]
snp_type$counts[5] <- dim(imputed[which(imputed$dup==1 & imputed$MAF >= 0.01),])[1]
snp_type$counts[6] <- dim(imputed[which(imputed$dup==0 & imputed$MAF >= 0.01),])[1]


snp_type$counts_high_r2[1] <- dim(imputed[which(imputed$dup==1 & imputed$Rsq >= 0.8),])[1]
snp_type$counts_high_r2[2] <- dim(imputed[which(imputed$dup==0 & imputed$Rsq >= 0.8),])[1]
snp_type$counts_high_r2[3] <- dim(imputed[which(imputed$dup==1 & imputed$MAF < 0.01 & imputed$Rsq >= 0.8),])[1]
snp_type$counts_high_r2[4] <- dim(imputed[which(imputed$dup==0 & imputed$MAF < 0.01 & imputed$Rsq >= 0.8),])[1]
snp_type$counts_high_r2[5] <- dim(imputed[which(imputed$dup==1 & imputed$MAF >= 0.01 & imputed$Rsq >= 0.8),])[1]
snp_type$counts_high_r2[6] <- dim(imputed[which(imputed$dup==0 & imputed$MAF >= 0.01 & imputed$Rsq >= 0.8),])[1]

snp_type$prop_high_r2 <- snp_type$counts_high_r2 / snp_type$counts * 100

write.table(snp_type,"snp_type_multi_allele_bi_allele.txt", quote=F, row.names=F)

bsub < multi_allele_bi_allele_count.sh


## for meta-imputation

# create an info file with chr:pos
sed 's/\:/\t/g' all.chr.info |awk '(NR>1){print $1":"$2"\t"$1":"$2":"$3":"$4"\t"$5"\t"$6}' | sed '1i pos\tSNP\tMAF\tRsq' > all.chr.pos.info 


#!/usr/bin/env Rscript

options(echo=TRUE)
print("In R now!")

## SCRIPT STARTS BELOW

library(data.table)
info <- fread("all.chr.pos.info",header=T)
info <- as.data.frame(info)
dup_pos <- info[duplicated(info$pos),"pos"]
# if the SNP ID appeared 3 times, then the 2nd and 3rd will be both included in the dup_pos file

info$dup <- ifelse(info$pos %in% dup_pos, 1,0)

info$Rsq <- as.numeric(info$Rsq)


snp_type <- data.frame(matrix(ncol = 4, nrow = 6))
colnames(snp_type) <- c('multi_allele', 'MAF_cat','counts', 'counts_high_r2')

snp_type$multi_allele[c(1,3,5)] <- "multi"
snp_type$multi_allele[c(2,4,6)] <- "bi"

snp_type$MAF_cat[c(1,2)] <- "all"
snp_type$MAF_cat[c(3,4)] <- "rare"
snp_type$MAF_cat[c(5,6)] <- "common"

imputed <- info 

snp_type$counts[1] <- dim(imputed[which(imputed$dup==1),])[1]
snp_type$counts[2] <- dim(imputed[which(imputed$dup==0),])[1]
snp_type$counts[3] <- dim(imputed[which(imputed$dup==1 & imputed$MAF < 0.01),])[1]
snp_type$counts[4] <- dim(imputed[which(imputed$dup==0 & imputed$MAF < 0.01),])[1]
snp_type$counts[5] <- dim(imputed[which(imputed$dup==1 & imputed$MAF >= 0.01),])[1]
snp_type$counts[6] <- dim(imputed[which(imputed$dup==0 & imputed$MAF >= 0.01),])[1]


snp_type$counts_high_r2[1] <- dim(imputed[which(imputed$dup==1 & imputed$Rsq >= 0.8),])[1]
snp_type$counts_high_r2[2] <- dim(imputed[which(imputed$dup==0 & imputed$Rsq >= 0.8),])[1]
snp_type$counts_high_r2[3] <- dim(imputed[which(imputed$dup==1 & imputed$MAF < 0.01 & imputed$Rsq >= 0.8),])[1]
snp_type$counts_high_r2[4] <- dim(imputed[which(imputed$dup==0 & imputed$MAF < 0.01 & imputed$Rsq >= 0.8),])[1]
snp_type$counts_high_r2[5] <- dim(imputed[which(imputed$dup==1 & imputed$MAF >= 0.01 & imputed$Rsq >= 0.8),])[1]
snp_type$counts_high_r2[6] <- dim(imputed[which(imputed$dup==0 & imputed$MAF >= 0.01 & imputed$Rsq >= 0.8),])[1]

snp_type$prop_high_r2 <- snp_type$counts_high_r2 / snp_type$counts * 100

write.table(snp_type,"snp_type_multi_allele_bi_allele.txt", quote=F, row.names=F)

# change rare to Rare, common to Common, all to All
# change multi to multi_allelic, bi to bi_allelic
#########################################################################################################
####             plot number of SNPs and well-imputed SNPs by multi-allelic and bi-allelic SNPs       ###
#########################################################################################################

all_meta_imp <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_topmed_1000G/snp_type_multi_allele_bi_allele.txt", header=T) 
all_topmed <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/snp_type_multi_allele_bi_allele.txt", header=T)
all_1000g <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/snp_type_multi_allele_bi_allele.txt", header=T) 
all_1000g_38 <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_38/snp_type_multi_allele_bi_allele.txt", header=T) 
all_1000g_37 <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_37/snp_type_multi_allele_bi_allele.txt", header=T) 
all_ga <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_GA/snp_type_multi_allele_bi_allele.txt", header=T) 

dim(all_topmed) # [1] 6 5
dim(all_1000g) # [1] 6 5
dim(all_ga) # [1] 6 5
dim(all_1000g_37) # [1] 6 5
dim(all_1000g_38) # [1] 6 5
dim(all_meta_imp) # [1] 6 5

all_meta_imp$ref_panel <- "TOPMed_1000G"
all_topmed$ref_panel <- "TOPMed"
all_1000g$ref_panel <- "Expanded_1000G"
all_1000g_38$ref_panel <- "1000G_GRCh38"
all_1000g_37$ref_panel <- "1000G_GRCh37_SAS"
all_ga$ref_panel <- "GenomeAsia_Pilot"

snp_count <- rbind(all_meta_imp, all_topmed, all_1000g, all_1000g_38, all_1000g_37, all_ga) # [1] 36  6

snp_count_all <- snp_count[,c(1,2,3,6)]
snp_count_well <- snp_count[,c(1,2,4,6)]

colnames(snp_count_well)[colnames(snp_count_well)=="counts_high_r2"] <- "counts"

snp_count_all$snp_cat <- "all"
snp_count_well$snp_cat <- "well_imp"

snp_count_merge <- rbind(snp_count_all,snp_count_well)
snp_count_merge$ref_panel <- factor(snp_count_merge$ref_panel,  levels=c("TOPMed_1000G","TOPMed","Expanded_1000G","1000G_GRCh38", "1000G_GRCh37_SAS", "GenomeAsia_Pilot"))
snp_count_merge$MAF_cat <- factor(snp_count_merge$MAF_cat,  levels=c("All","Common","Rare"))
snp_count_merge$snp_cat <- factor(snp_count_merge$snp_cat,  levels=c("all","well_imp"))

snp_count_all <- snp_count_merge[which(snp_count_merge$MAF_cat=="All"),]
snp_count_maf <- snp_count_merge[which(snp_count_merge$MAF_cat !="All"),]

# Separated by all SNPs and well-imputed SNPs for rare and common variants
snp_count_all_maf <- snp_count_maf[which(snp_count_maf$snp_cat =="all"),]
snp_count_well_maf <- snp_count_maf[which(snp_count_maf$snp_cat =="well_imp"),]

library(ggplot2)
pdf('multi_bi_allele_all_snp_count.pdf',height=10,width=20)
p <- ggplot(data=snp_count_all_maf, aes(x=interaction(MAF_cat,multi_allele), y=counts, fill=ref_panel)) + xlab("Variant Type") + ylab("Count") + geom_bar(stat="identity", position=position_dodge()) +  geom_text(aes(label=scales::label_number_si(accuracy=0.1)(counts)),size=3, position=position_dodge(width=0.9), vjust=-0.25) + theme_classic(base_size=20)  + scale_x_discrete() + scale_y_continuous(labels=scales::label_number_si())+theme(axis.text = element_text(size = 20,colour="black"))
p + scale_fill_manual(values=c("#661100", "#D55E00", "#E69F00","#0072B2","#332288","#999999"),name="Reference Panel",labels=c("TOPMed_1000G","TOPMed","Expanded_1000G","1000G_GRCh38", "1000G_GRCh37_SAS", "GenomeAsia_Pilot")) + theme(plot.title = element_text(hjust = 0.5)) + labs(title="All SNPs")
dev.off()

pdf('multi_bi_allele_well_imputed_snp_count.pdf',height=10,width=20)
p <- ggplot(data=snp_count_well_maf, aes(x=interaction(MAF_cat,multi_allele), y=counts, fill=ref_panel)) + xlab("Variant Type") + ylab("Count") + geom_bar(stat="identity", position=position_dodge()) +  geom_text(aes(label=scales::label_number_si(accuracy=0.1)(counts)),size=4, position=position_dodge(width=0.9), vjust=-0.25) + theme_classic(base_size=20)  + scale_x_discrete() + scale_y_continuous(labels=scales::label_number_si())+theme(axis.text = element_text(size = 20,colour="black"))
p + scale_fill_manual(values=c("#661100", "#D55E00", "#E69F00","#0072B2","#332288","#999999"),name="Reference Panel",labels=c("TOPMed_1000G","TOPMed","Expanded_1000G","1000G_GRCh38", "1000G_GRCh37_SAS", "GenomeAsia_Pilot"))+ theme(legend.position="none",plot.title = element_text(hjust = 0.5)) + labs(title="Well-imputed SNPs")
dev.off()



##################################################################################
####             plot number of SNPs and well-imputed SNPs by SNP type         ###
##################################################################################

### create the short version of vcf for snp annotation


#!/bin/bash
#BSUB -J brief_vcf # Job name
#BSUB -P acc_psychgen # allocation account
#BSUB -q premium # queue
#BSUB -n 1 # number of compute cores
#BSUB -W 96:00 # walltime in HH:MM
#BSUB -R himem
#BSUB -R rusage[mem=150G] # 32GB
#BSUB -R span[hosts=1] # all cores from one node
#BSUB -o %J.stdout # output log (%J : JobID)
#BSUB -eo %J.stderr # error log
#BSUB -L /bin/bash # Initialize the execution environment

for file in chr*.dose.vcf.gz
do

output_file=$(echo ${file%%.*}".brief.vcf")
zcat $file | awk '(NR>19){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' > $output_file

done

bsub < brief_vcf.sh

# meta-imputation only 

#!/bin/bash
#BSUB -J brief_vcf # Job name
#BSUB -P acc_psychgen # allocation account
#BSUB -q premium # queue
#BSUB -n 1 # number of compute cores
#BSUB -W 96:00 # walltime in HH:MM
#BSUB -R himem
#BSUB -R rusage[mem=150G] # 32GB
#BSUB -R span[hosts=1] # all cores from one node
#BSUB -o %J.stdout # output log (%J : JobID)
#BSUB -eo %J.stderr # error log
#BSUB -L /bin/bash # Initialize the execution environment

for file in chr*vcf.gz
do

output_file=$(echo ${file%%.*}".brief.vcf")
zcat $file | awk '(NR>15){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' > $output_file

done

bsub < brief_vcf.sh

# meta-imputation
grep '^#CHROM' chr1_topmed_1000g.brief.vcf > all.chr.brief.vcf
grep -v '^#CHROM' chr1_topmed_1000g.brief.vcf  chr2_topmed_1000g.brief.vcf  chr3_topmed_1000g.brief.vcf chr4_topmed_1000g.brief.vcf chr5_topmed_1000g.brief.vcf chr6_topmed_1000g.brief.vcf chr7_topmed_1000g.brief.vcf chr8_topmed_1000g.brief.vcf chr9_topmed_1000g.brief.vcf chr10_topmed_1000g.brief.vcf chr11_topmed_1000g.brief.vcf chr12_topmed_1000g.brief.vcf chr13_topmed_1000g.brief.vcf chr14_topmed_1000g.brief.vcf chr15_topmed_1000g.brief.vcf chr16_topmed_1000g.brief.vcf chr17_topmed_1000g.brief.vcf chr18_topmed_1000g.brief.vcf chr19_topmed_1000g.brief.vcf chr20_topmed_1000g.brief.vcf chr21_topmed_1000g.brief.vcf chr22_topmed_1000g.brief.vcf >> all.chr.brief.vcf

# meta-imputation

#!/bin/bash
#BSUB -J brief_vcf_merge # Job name
#BSUB -P acc_psychgen # allocation account
#BSUB -q premium # queue
#BSUB -n 1 # number of compute cores
#BSUB -W 24:00 # walltime in HH:MM
#BSUB -R himem
#BSUB -R rusage[mem=96G] # 32GB
#BSUB -R span[hosts=1] # all cores from one node
#BSUB -o %J.stdout # output log (%J : JobID)
#BSUB -eo %J.stderr # error log
#BSUB -L /bin/bash # Initialize the execution environment

grep '^#CHROM' chr1.brief.vcf > all.chr.brief.vcf
grep -v '^#CHROM' chr1.brief.vcf  chr2.brief.vcf  chr3.brief.vcf chr4.brief.vcf chr5.brief.vcf chr6.brief.vcf chr7.brief.vcf chr8.brief.vcf chr9.brief.vcf chr10.brief.vcf chr11.brief.vcf chr12.brief.vcf chr13.brief.vcf chr14.brief.vcf chr15.brief.vcf chr16.brief.vcf chr17.brief.vcf chr18.brief.vcf chr19.brief.vcf chr20.brief.vcf chr21.brief.vcf chr22.brief.vcf >> all.chr.brief.vcf

awk '(NR>1){print $0}' all.chr.brief.vcf  | sed 's/\:/\t/1' | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' |sed '1i CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' > all.chr.brief.final.vcf 

bsub < brief_vcf_merge.sh

wc -l  all.chr.brief.vcf
# 296353412 SNPs (meta-imputation)
# 292140627 SNPs (TopMed)
# 48895833 SNPs (expanded 1000G)
# 45210361 SNPs (1000G GRCh38)
# 47099752 SNPs (1000G GRCh37 SAS)
# 21524720 SNPs (GenomeAsia)

awk '(NR>6) {print $0}' all.chr.brief.ann.vcf | wc -l
# 21524720 SNPs (GenomeAsia)

#snpeff located at /hpc/users/xuj24/snpeff/snpEff/snpEff.jar

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

#for file in chr*.brief.vcf
#do
#output_file=$(echo ${file%%.*}".brief.ann.vcf")
#java -Xmx8g -jar /hpc/users/xuj18/snpeff/snpEff/snpEff.jar -c /hpc/users/xuj18/snpeff/snpEff/snpEff.config -v GRCh37.75 $file > $output_file
#done

# -Xmx8g: set the amount of memory in your java virtual machine to 8g
# when you run SnpEff from a different directory than your install directory, you have to specify where the config file is located using the '-c' command line option.
# -v: verbose, this makes SnpEff to show a lot of information which can be useful for debugging.

# grep different SNPs by their effect type

grep 3_prime_UTR_variant  all.chr.brief.ann.vcf|wc -l # GenomeAsia: 261929; expanded 1000G: 732889
grep 5_prime_UTR_premature_start_codon_gain_variant  all.chr.brief.ann.vcf|wc -l # GenomeAsia: 14085; expanded 1000G: 30789
grep 5_prime_UTR_variant  all.chr.brief.ann.vcf|wc -l # GenomeAsia: 78213; expanded 1000G: 182085
grep downstream_gene_variant  all.chr.brief.ann.vcf|wc -l # GenomeAsia: 3141614; expanded 1000G: 7934803
grep initiator_codon_variant  all.chr.brief.ann.vcf|wc -l # GenomeAsia: 32; expanded 1000G: 82
grep intergenic_region  all.chr.brief.ann.vcf|wc -l # GenomeAsia: 9867172; expanded 1000G: 19180753
grep intragenic_variant  all.chr.brief.ann.vcf|wc -l # GenomeAsia: 4226; expanded 1000G: 9781
grep intron_variant  all.chr.brief.ann.vcf|wc -l # GenomeAsia: 11228836; expanded 1000G: 28620565
grep missense_variant  all.chr.brief.ann.vcf|wc -l # GenomeAsia: 103406; expanded 1000G: 248158
grep non_coding_transcript_exon_variant  all.chr.brief.ann.vcf|wc -l # GenomeAsia: 514453; expanded 1000G: 1540544
grep splice_acceptor_variant  all.chr.brief.ann.vcf|wc -l # GenomeAsia: 1189; expanded 1000G: 3719
grep splice_donor_variant  all.chr.brief.ann.vcf|wc -l # GenomeAsia: 1790; expanded 1000G: 5257
grep splice_region_variant  all.chr.brief.ann.vcf|wc -l # GenomeAsia: 31028; expanded 1000G: 82073
grep start_lost  all.chr.brief.ann.vcf|wc -l # GenomeAsia: 326; expanded 1000G: 749
grep start_retained_variant  all.chr.brief.ann.vcf|wc -l # GenomeAsia: 2; expanded 1000G: 22
grep stop_gained  all.chr.brief.ann.vcf|wc -l # GenomeAsia: 1418; expanded 1000G: 4131
grep stop_lost  all.chr.brief.ann.vcf|wc -l # GenomeAsia: 181; expanded 1000G: 473
grep stop_retained_variant  all.chr.brief.ann.vcf|wc -l # GenomeAsia: 108; expanded 1000G: 247
grep synonymous_variant all.chr.brief.ann.vcf|wc -l # GenomeAsia: 88749; expanded 1000G: 187416
grep upstream_gene_variant  all.chr.brief.ann.vcf|wc -l # GenomeAsia: 3043457; expanded 1000G: 7564128


## merge the R2 info with annotation info

####   for individual imputation panels

#!/bin/bash
#BSUB -J merge_info_annotation # Job name
#BSUB -n 1
#BSUB -P acc_psychgen # allocation account
#BSUB -q premium # queue
#BSUB -W 4:00 # walltime in HH:MM
#BSUB -R himem
#BSUB -R rusage[mem=200G] # 120 GB of memory 
#BSUB -R span[hosts=1] # all cores from one node
#BSUB -o merge_info_annotation.stdout # output log (%J : JobID)
#BSUB -eo merge_info_annotation.stderr # error log
#BSUB -L /bin/bash # Initialize the execution environment

# get the brief annotation file
awk '(NR>6) {print $0}'  all.chr.brief.ann.vcf  | awk '{print $3"\t"$8}'  | sort -k1 >  all.chr.brief.ann.sort.vcf

# get the INFO file
sed 's/^[^:]*://g' all.chr.info  | awk '(NR>1) {print $0} ' |sort -k1  >  all.chr.sort.info
# In other words, select the beginning of the line, any number of things that aren't a colon, and the first colon.
# The first ^ means beginning of the line. The [^:] is just the only way I know how to write not a colon.
# The * after the colon means any number of the things right before me (in this case the not-colon). Finally, the : selects the colon.

# check to see if they can be directly pasted together
paste  all.chr.sort.info  all.chr.brief.ann.sort.vcf > all.chr.ann.sort.info

bsub < merge_info_annotation.sh

####   for meta-imputation

#!/bin/bash
#BSUB -J merge_info_annotation # Job name
#BSUB -n 1
#BSUB -P acc_psychgen # allocation account
#BSUB -q premium # queue
#BSUB -W 4:00 # walltime in HH:MM
#BSUB -R himem
#BSUB -R rusage[mem=200G] # 120 GB of memory 
#BSUB -R span[hosts=1] # all cores from one node
#BSUB -o merge_info_annotation.stdout # output log (%J : JobID)
#BSUB -eo merge_info_annotation.stderr # error log
#BSUB -L /bin/bash # Initialize the execution environment

# get the brief annotation file
awk '(NR>6) {print $0}'  all.chr.brief.ann.vcf  | awk '{print $3"\t"$8}'  | sort -k1 >  all.chr.brief.ann.sort.vcf

# get the INFO file
awk '(NR>1) {print $0} '  all.chr.info  |sort -k1  >  all.chr.sort.info

# check to see if they can be directly pasted together
paste  all.chr.sort.info  all.chr.brief.ann.sort.vcf > all.chr.ann.sort.info

bsub < merge_info_annotation.sh


#!/bin/bash
#BSUB -J quality_by_snp_type # Job name
#BSUB -n 1
#BSUB -P acc_psychgen # allocation account
#BSUB -q premium # queue
#BSUB -W 96:00 # walltime in HH:MM
#BSUB -R himem
#BSUB -R rusage[mem=150G] # 120 GB of memory 
#BSUB -R span[hosts=1] # all cores from one node
#BSUB -o quality_by_snp_type.stdout # output log (%J : JobID)
#BSUB -eo quality_by_snp_type.stderr # error log
#BSUB -L /bin/bash # Initialize the execution environment


module load R

Rscript quality_by_snp_type.r

bsub < quality_by_snp_type.sh

#!/usr/bin/env Rscript

options(echo=TRUE)
print("In R now!")

## SCRIPT STARTS BELOW

 # in R
library(data.table)
all <- as.data.frame(fread("all.chr.ann.sort.info", header=F))

colnames(all) <- c("SNP","REF.0.","ALT.1.","ALT_Frq","MAF","AvgCall","Rsq","Genotyped","LooRsq","EmpR","EmpRsq","Dose0","Dose1","dup_col_SNP","info")
all$Rsq <- as.numeric(all$Rsq)


imputed <- all[which(all$Genotyped!="Typed_Only"),]

variant_3_prime_UTR <- imputed[which(grepl("3_prime_UTR_variant", imputed$info)),]
variant_5_prime_UTR_premature_start_codon_gain <- imputed[which(grepl("5_prime_UTR_premature_start_codon_gain_variant", imputed$info)),]
variant_5_prime_UTR <- imputed[which(grepl("5_prime_UTR_variant", imputed$info)),]
variant_downstream <- imputed[which(grepl("downstream_gene_variant", imputed$info)),]
variant_initiator_codon <- imputed[which(grepl("initiator_codon_variant", imputed$info)),]
variant_intergenic_region <- imputed[which(grepl("intergenic_region",imputed$info)),]
variant_intragenic <- imputed[which(grepl("intragenic_variant", imputed$info)),]
variant_intron <- imputed[which(grepl("intron_variant", imputed$info)),]
variant_missense <- imputed[which(grepl("missense_variant", imputed$info)),]
variant_non_coding_transcript_exon <- imputed[which(grepl("non_coding_transcript_exon_variant", imputed$info)),]
variant_splice_acceptor <- imputed[which(grepl("splice_acceptor_variant", imputed$info)),]
variant_splice_donor <- imputed[which(grepl("splice_donor_variant", imputed$info)),]
variant_splice_region <- imputed[which(grepl("splice_region_variant", imputed$info)),]
variant_start_lost <- imputed[which(grepl("start_lost", imputed$info)),]
variant_start_retained <- imputed[which(grepl("start_retained_variant", imputed$info)),]
variant_stop_gained <- imputed[which(grepl("stop_gained", imputed$info)),]
variant_stop_lost <- imputed[which(grepl("stop_lost",imputed$info)),]
variant_stop_retained <- imputed[which(grepl("stop_retained_variant", imputed$info)),]
variant_synonymous <- imputed[which(grepl("synonymous_variant", imputed$info)),]
variant_upstream <- imputed[which(grepl("upstream_gene_variant", imputed$info)),]

good_variant_3_prime_UTR <- variant_3_prime_UTR[which(variant_3_prime_UTR$Rsq >= 0.8),]
good_variant_5_prime_UTR_premature_start_codon_gain  <- variant_5_prime_UTR_premature_start_codon_gain[which(variant_5_prime_UTR_premature_start_codon_gain$Rsq >= 0.8),]
good_variant_5_prime_UTR  <- variant_5_prime_UTR[which(variant_5_prime_UTR$Rsq >= 0.8),]
good_variant_downstream  <- variant_downstream[which(variant_downstream$Rsq >= 0.8),]
good_variant_initiator_codon <- variant_initiator_codon[which(variant_initiator_codon$Rsq >= 0.8),]
good_variant_intergenic_region  <- variant_intergenic_region[which(variant_intergenic_region$Rsq >= 0.8),]
good_variant_intragenic <- variant_intragenic[which(variant_intragenic$Rsq >= 0.8),]
good_variant_intron <- variant_intron[which(variant_intron$Rsq >= 0.8),]
good_variant_missense  <- variant_missense[which(variant_missense$Rsq >= 0.8),]
good_variant_non_coding_transcript_exon  <- variant_non_coding_transcript_exon[which(variant_non_coding_transcript_exon$Rsq >= 0.8),]
good_variant_splice_acceptor <- variant_splice_acceptor[which(variant_splice_acceptor$Rsq >= 0.8),]
good_variant_splice_donor  <- variant_splice_donor[which(variant_splice_donor$Rsq >= 0.8),]
good_variant_splice_region  <- variant_splice_region[which(variant_splice_region$Rsq >= 0.8),]
good_variant_start_lost  <- variant_start_lost[which(variant_start_lost$Rsq >= 0.8),]
good_variant_start_retained  <- variant_start_retained[which(variant_start_retained$Rsq >= 0.8),]
good_variant_stop_gained <- variant_stop_gained[which(variant_stop_gained$Rsq >= 0.8),]
good_variant_stop_lost  <- variant_stop_lost[which(variant_stop_lost$Rsq >= 0.8),]
good_variant_stop_retained  <- variant_stop_retained[which(variant_stop_retained$Rsq >= 0.8),]
good_variant_synonymous  <- variant_synonymous[which(variant_synonymous$Rsq >= 0.8),]
good_variant_upstream  <- variant_upstream[which(variant_upstream$Rsq >= 0.8),]

snp_count <- data.frame(matrix(ncol = 5, nrow = 40))
colnames(snp_count) <- c("ref_panel","maf_cat", 'snp_type', 'counts', 'counts_high_r2')

snp_count$ref_panel <- "GenomeAsia_Pilot"
snp_count$snp_type <-  rep(c("3_prime_UTR_variant","5_prime_UTR_premature_start_codon_gain_variant","5_prime_UTR_variant","downstream_gene_variant","initiator_codon_variant","intergenic_region",
"intragenic_variant","intron_variant","missense_variant","non_coding_transcript_exon_variant","splice_acceptor_variant","splice_donor_variant","splice_region_variant","start_lost",
"start_retained_variant","stop_gained","stop_lost","stop_retained_variant","synonymous_variant","upstream_gene_variant"),2)
snp_count$maf_cat[1:20] <- "Common"
snp_count$maf_cat[21:40] <- "Rare"

dim(variant_3_prime_UTR) 
dim(variant_3_prime_UTR[which(variant_3_prime_UTR$MAF >=0.01),])
dim(variant_3_prime_UTR[which(variant_3_prime_UTR$MAF <0.01),])


snp_count$counts[1]  <- dim(variant_3_prime_UTR[which(variant_3_prime_UTR$MAF >=0.01),])[1]
snp_count$counts[2]  <- dim(variant_5_prime_UTR_premature_start_codon_gain[which(variant_5_prime_UTR_premature_start_codon_gain$MAF >=0.01),])[1]
snp_count$counts[3]  <- dim(variant_5_prime_UTR[which(variant_5_prime_UTR$MAF >=0.01),])[1]
snp_count$counts[4]  <- dim(variant_downstream[which(variant_downstream$MAF >=0.01),])[1]
snp_count$counts[5]  <- dim(variant_initiator_codon[which(variant_initiator_codon$MAF >=0.01),])[1]
snp_count$counts[6]  <- dim(variant_intergenic_region[which(variant_intergenic_region$MAF >=0.01),])[1]
snp_count$counts[7]  <- dim(variant_intragenic[which(variant_intragenic$MAF >=0.01),])[1]
snp_count$counts[8]  <- dim(variant_intron[which(variant_intron$MAF >=0.01),])[1]
snp_count$counts[9]  <- dim(variant_missense[which(variant_missense$MAF >=0.01),])[1]
snp_count$counts[10] <- dim(variant_non_coding_transcript_exon[which(variant_non_coding_transcript_exon$MAF >=0.01),])[1]
snp_count$counts[11] <- dim(variant_splice_acceptor[which(variant_splice_acceptor$MAF >=0.01),])[1]
snp_count$counts[12] <- dim(variant_splice_donor[which(variant_splice_donor$MAF >=0.01),])[1]
snp_count$counts[13] <- dim(variant_splice_region[which(variant_splice_region$MAF >=0.01),])[1]
snp_count$counts[14] <- dim(variant_start_lost[which(variant_start_lost$MAF >=0.01),])[1]
snp_count$counts[15] <- dim(variant_start_retained[which(variant_start_retained$MAF >=0.01),])[1]
snp_count$counts[16] <- dim(variant_stop_gained[which(variant_stop_gained$MAF >=0.01),])[1]
snp_count$counts[17] <- dim(variant_stop_lost[which(variant_stop_lost$MAF >=0.01),])[1]
snp_count$counts[18] <- dim(variant_stop_retained[which(variant_stop_retained$MAF >=0.01),])[1]
snp_count$counts[19] <- dim(variant_synonymous[which(variant_synonymous$MAF >=0.01),])[1]
snp_count$counts[20] <- dim(variant_upstream[which(variant_upstream$MAF >=0.01),])[1]

snp_count$counts[21]  <- dim(variant_3_prime_UTR[which(variant_3_prime_UTR$MAF <0.01),])[1]
snp_count$counts[22]  <- dim(variant_5_prime_UTR_premature_start_codon_gain[which(variant_5_prime_UTR_premature_start_codon_gain$MAF <0.01),])[1]
snp_count$counts[23]  <- dim(variant_5_prime_UTR[which(variant_5_prime_UTR$MAF <0.01),])[1]
snp_count$counts[24]  <- dim(variant_downstream[which(variant_downstream$MAF <0.01),])[1]
snp_count$counts[25]  <- dim(variant_initiator_codon[which(variant_initiator_codon$MAF <0.01),])[1]
snp_count$counts[26]  <- dim(variant_intergenic_region[which(variant_intergenic_region$MAF <0.01),])[1]
snp_count$counts[27]  <- dim(variant_intragenic[which(variant_intragenic$MAF <0.01),])[1]
snp_count$counts[28]  <- dim(variant_intron[which(variant_intron$MAF <0.01),])[1]
snp_count$counts[29]  <- dim(variant_missense[which(variant_missense$MAF <0.01),])[1]
snp_count$counts[30] <- dim(variant_non_coding_transcript_exon[which(variant_non_coding_transcript_exon$MAF <0.01),])[1]
snp_count$counts[31] <- dim(variant_splice_acceptor[which(variant_splice_acceptor$MAF <0.01),])[1]
snp_count$counts[32] <- dim(variant_splice_donor[which(variant_splice_donor$MAF <0.01),])[1]
snp_count$counts[33] <- dim(variant_splice_region[which(variant_splice_region$MAF <0.01),])[1]
snp_count$counts[34] <- dim(variant_start_lost[which(variant_start_lost$MAF <0.01),])[1]
snp_count$counts[35] <- dim(variant_start_retained[which(variant_start_retained$MAF <0.01),])[1]
snp_count$counts[36] <- dim(variant_stop_gained[which(variant_stop_gained$MAF <0.01),])[1]
snp_count$counts[37] <- dim(variant_stop_lost[which(variant_stop_lost$MAF <0.01),])[1]
snp_count$counts[38] <- dim(variant_stop_retained[which(variant_stop_retained$MAF <0.01),])[1]
snp_count$counts[39] <- dim(variant_synonymous[which(variant_synonymous$MAF <0.01),])[1]
snp_count$counts[40] <- dim(variant_upstream[which(variant_upstream$MAF <0.01),])[1]



snp_count$counts_high_r2[1]  <- dim(good_variant_3_prime_UTR[which(good_variant_3_prime_UTR$MAF >=0.01),])[1]
snp_count$counts_high_r2[2]  <- dim(good_variant_5_prime_UTR_premature_start_codon_gain[which(good_variant_5_prime_UTR_premature_start_codon_gain$MAF >=0.01),])[1]
snp_count$counts_high_r2[3]  <- dim(good_variant_5_prime_UTR[which(good_variant_5_prime_UTR$MAF >=0.01),])[1]
snp_count$counts_high_r2[4]  <- dim(good_variant_downstream[which(good_variant_downstream$MAF >=0.01),])[1]
snp_count$counts_high_r2[5]  <- dim(good_variant_initiator_codon[which(good_variant_initiator_codon$MAF >=0.01),])[1]
snp_count$counts_high_r2[6]  <- dim(good_variant_intergenic_region[which(good_variant_intergenic_region$MAF >=0.01),])[1]
snp_count$counts_high_r2[7]  <- dim(good_variant_intragenic[which(good_variant_intragenic$MAF >=0.01),])[1]
snp_count$counts_high_r2[8]  <- dim(good_variant_intron[which(good_variant_intron$MAF >=0.01),])[1]
snp_count$counts_high_r2[9]  <- dim(good_variant_missense[which(good_variant_missense$MAF >=0.01),])[1]
snp_count$counts_high_r2[10] <- dim(good_variant_non_coding_transcript_exon[which(good_variant_non_coding_transcript_exon$MAF >=0.01),])[1]
snp_count$counts_high_r2[11] <- dim(good_variant_splice_acceptor[which(good_variant_splice_acceptor$MAF >=0.01),])[1]
snp_count$counts_high_r2[12] <- dim(good_variant_splice_donor[which(good_variant_splice_donor$MAF >=0.01),])[1]
snp_count$counts_high_r2[13] <- dim(good_variant_splice_region[which(good_variant_splice_region$MAF >=0.01),])[1]
snp_count$counts_high_r2[14] <- dim(good_variant_start_lost[which(good_variant_start_lost$MAF >=0.01),])[1]
snp_count$counts_high_r2[15] <- dim(good_variant_start_retained[which(good_variant_start_retained$MAF >=0.01),])[1]
snp_count$counts_high_r2[16] <- dim(good_variant_stop_gained[which(good_variant_stop_gained$MAF >=0.01),])[1]
snp_count$counts_high_r2[17] <- dim(good_variant_stop_lost[which(good_variant_stop_lost$MAF >=0.01),])[1]
snp_count$counts_high_r2[18] <- dim(good_variant_stop_retained[which(good_variant_stop_retained$MAF >=0.01),])[1]
snp_count$counts_high_r2[19] <- dim(good_variant_synonymous[which(good_variant_synonymous$MAF >=0.01),])[1]
snp_count$counts_high_r2[20] <- dim(good_variant_upstream[which(good_variant_upstream$MAF >=0.01),])[1]

snp_count$counts_high_r2[21]  <- dim(good_variant_3_prime_UTR[which(good_variant_3_prime_UTR$MAF <0.01),])[1]
snp_count$counts_high_r2[22]  <- dim(good_variant_5_prime_UTR_premature_start_codon_gain[which(good_variant_5_prime_UTR_premature_start_codon_gain$MAF <0.01),])[1]
snp_count$counts_high_r2[23]  <- dim(good_variant_5_prime_UTR[which(good_variant_5_prime_UTR$MAF <0.01),])[1]
snp_count$counts_high_r2[24]  <- dim(good_variant_downstream[which(good_variant_downstream$MAF <0.01),])[1]
snp_count$counts_high_r2[25]  <- dim(good_variant_initiator_codon[which(good_variant_initiator_codon$MAF <0.01),])[1]
snp_count$counts_high_r2[26]  <- dim(good_variant_intergenic_region[which(good_variant_intergenic_region$MAF <0.01),])[1]
snp_count$counts_high_r2[27]  <- dim(good_variant_intragenic[which(good_variant_intragenic$MAF <0.01),])[1]
snp_count$counts_high_r2[28]  <- dim(good_variant_intron[which(good_variant_intron$MAF <0.01),])[1]
snp_count$counts_high_r2[29]  <- dim(good_variant_missense[which(good_variant_missense$MAF <0.01),])[1]
snp_count$counts_high_r2[30] <- dim(good_variant_non_coding_transcript_exon[which(good_variant_non_coding_transcript_exon$MAF <0.01),])[1]
snp_count$counts_high_r2[31] <- dim(good_variant_splice_acceptor[which(good_variant_splice_acceptor$MAF <0.01),])[1]
snp_count$counts_high_r2[32] <- dim(good_variant_splice_donor[which(good_variant_splice_donor$MAF <0.01),])[1]
snp_count$counts_high_r2[33] <- dim(good_variant_splice_region[which(good_variant_splice_region$MAF <0.01),])[1]
snp_count$counts_high_r2[34] <- dim(good_variant_start_lost[which(good_variant_start_lost$MAF <0.01),])[1]
snp_count$counts_high_r2[35] <- dim(good_variant_start_retained[which(good_variant_start_retained$MAF <0.01),])[1]
snp_count$counts_high_r2[36] <- dim(good_variant_stop_gained[which(good_variant_stop_gained$MAF <0.01),])[1]
snp_count$counts_high_r2[37] <- dim(good_variant_stop_lost[which(good_variant_stop_lost$MAF <0.01),])[1]
snp_count$counts_high_r2[38] <- dim(good_variant_stop_retained[which(good_variant_stop_retained$MAF <0.01),])[1]
snp_count$counts_high_r2[39] <- dim(good_variant_synonymous[which(good_variant_synonymous$MAF <0.01),])[1]
snp_count$counts_high_r2[40] <- dim(good_variant_upstream[which(good_variant_upstream$MAF <0.01),])[1]

dim(variant_upstream)
dim(variant_upstream[which(variant_upstream$MAF >=0.01),]) 
dim(variant_upstream[which(variant_upstream$MAF <0.01),]) 


dim(good_variant_upstream)
dim(good_variant_upstream[which(good_variant_upstream$MAF >=0.01),]) 
dim(good_variant_upstream[which(good_variant_upstream$MAF <0.01),]) 


snp_count$prop_high_r2 <- snp_count$counts_high_r2 / snp_count$counts * 100

write.table(snp_count,"snp_type_counts.txt", quote=F, row.names=F)



###############         For meta-imputation        ###############

#!/usr/bin/env Rscript

options(echo=TRUE)
print("In R now!")

## SCRIPT STARTS BELOW

 # in R
library(data.table)
all <- as.data.frame(fread("all.chr.ann.sort.info", header=F))

colnames(all) <- c("SNP","MAF","Rsq","dup_col_SNP","info")
all$Rsq <- as.numeric(all$Rsq)

imputed <- all 

variant_3_prime_UTR <- imputed[which(grepl("3_prime_UTR_variant", imputed$info)),]
variant_5_prime_UTR_premature_start_codon_gain <- imputed[which(grepl("5_prime_UTR_premature_start_codon_gain_variant", imputed$info)),]
variant_5_prime_UTR <- imputed[which(grepl("5_prime_UTR_variant", imputed$info)),]
variant_downstream <- imputed[which(grepl("downstream_gene_variant", imputed$info)),]
variant_initiator_codon <- imputed[which(grepl("initiator_codon_variant", imputed$info)),]
variant_intergenic_region <- imputed[which(grepl("intergenic_region",imputed$info)),]
variant_intragenic <- imputed[which(grepl("intragenic_variant", imputed$info)),]
variant_intron <- imputed[which(grepl("intron_variant", imputed$info)),]
variant_missense <- imputed[which(grepl("missense_variant", imputed$info)),]
variant_non_coding_transcript_exon <- imputed[which(grepl("non_coding_transcript_exon_variant", imputed$info)),]
variant_splice_acceptor <- imputed[which(grepl("splice_acceptor_variant", imputed$info)),]
variant_splice_donor <- imputed[which(grepl("splice_donor_variant", imputed$info)),]
variant_splice_region <- imputed[which(grepl("splice_region_variant", imputed$info)),]
variant_start_lost <- imputed[which(grepl("start_lost", imputed$info)),]
variant_start_retained <- imputed[which(grepl("start_retained_variant", imputed$info)),]
variant_stop_gained <- imputed[which(grepl("stop_gained", imputed$info)),]
variant_stop_lost <- imputed[which(grepl("stop_lost",imputed$info)),]
variant_stop_retained <- imputed[which(grepl("stop_retained_variant", imputed$info)),]
variant_synonymous <- imputed[which(grepl("synonymous_variant", imputed$info)),]
variant_upstream <- imputed[which(grepl("upstream_gene_variant", imputed$info)),]

good_variant_3_prime_UTR <- variant_3_prime_UTR[which(variant_3_prime_UTR$Rsq >= 0.8),]
good_variant_5_prime_UTR_premature_start_codon_gain  <- variant_5_prime_UTR_premature_start_codon_gain[which(variant_5_prime_UTR_premature_start_codon_gain$Rsq >= 0.8),]
good_variant_5_prime_UTR  <- variant_5_prime_UTR[which(variant_5_prime_UTR$Rsq >= 0.8),]
good_variant_downstream  <- variant_downstream[which(variant_downstream$Rsq >= 0.8),]
good_variant_initiator_codon <- variant_initiator_codon[which(variant_initiator_codon$Rsq >= 0.8),]
good_variant_intergenic_region  <- variant_intergenic_region[which(variant_intergenic_region$Rsq >= 0.8),]
good_variant_intragenic <- variant_intragenic[which(variant_intragenic$Rsq >= 0.8),]
good_variant_intron <- variant_intron[which(variant_intron$Rsq >= 0.8),]
good_variant_missense  <- variant_missense[which(variant_missense$Rsq >= 0.8),]
good_variant_non_coding_transcript_exon  <- variant_non_coding_transcript_exon[which(variant_non_coding_transcript_exon$Rsq >= 0.8),]
good_variant_splice_acceptor <- variant_splice_acceptor[which(variant_splice_acceptor$Rsq >= 0.8),]
good_variant_splice_donor  <- variant_splice_donor[which(variant_splice_donor$Rsq >= 0.8),]
good_variant_splice_region  <- variant_splice_region[which(variant_splice_region$Rsq >= 0.8),]
good_variant_start_lost  <- variant_start_lost[which(variant_start_lost$Rsq >= 0.8),]
good_variant_start_retained  <- variant_start_retained[which(variant_start_retained$Rsq >= 0.8),]
good_variant_stop_gained <- variant_stop_gained[which(variant_stop_gained$Rsq >= 0.8),]
good_variant_stop_lost  <- variant_stop_lost[which(variant_stop_lost$Rsq >= 0.8),]
good_variant_stop_retained  <- variant_stop_retained[which(variant_stop_retained$Rsq >= 0.8),]
good_variant_synonymous  <- variant_synonymous[which(variant_synonymous$Rsq >= 0.8),]
good_variant_upstream  <- variant_upstream[which(variant_upstream$Rsq >= 0.8),]

snp_count <- data.frame(matrix(ncol = 5, nrow = 40))
colnames(snp_count) <- c("ref_panel","maf_cat", 'snp_type', 'counts', 'counts_high_r2')

snp_count$ref_panel <- "TOPMed_1000G"
snp_count$snp_type <-  rep(c("3_prime_UTR_variant","5_prime_UTR_premature_start_codon_gain_variant","5_prime_UTR_variant","downstream_gene_variant","initiator_codon_variant","intergenic_region",
"intragenic_variant","intron_variant","missense_variant","non_coding_transcript_exon_variant","splice_acceptor_variant","splice_donor_variant","splice_region_variant","start_lost",
"start_retained_variant","stop_gained","stop_lost","stop_retained_variant","synonymous_variant","upstream_gene_variant"),2)
snp_count$maf_cat[1:20] <- "Common"
snp_count$maf_cat[21:40] <- "Rare"

dim(variant_3_prime_UTR) 
dim(variant_3_prime_UTR[which(variant_3_prime_UTR$MAF >=0.01),])
dim(variant_3_prime_UTR[which(variant_3_prime_UTR$MAF <0.01),])


snp_count$counts[1]  <- dim(variant_3_prime_UTR[which(variant_3_prime_UTR$MAF >=0.01),])[1]
snp_count$counts[2]  <- dim(variant_5_prime_UTR_premature_start_codon_gain[which(variant_5_prime_UTR_premature_start_codon_gain$MAF >=0.01),])[1]
snp_count$counts[3]  <- dim(variant_5_prime_UTR[which(variant_5_prime_UTR$MAF >=0.01),])[1]
snp_count$counts[4]  <- dim(variant_downstream[which(variant_downstream$MAF >=0.01),])[1]
snp_count$counts[5]  <- dim(variant_initiator_codon[which(variant_initiator_codon$MAF >=0.01),])[1]
snp_count$counts[6]  <- dim(variant_intergenic_region[which(variant_intergenic_region$MAF >=0.01),])[1]
snp_count$counts[7]  <- dim(variant_intragenic[which(variant_intragenic$MAF >=0.01),])[1]
snp_count$counts[8]  <- dim(variant_intron[which(variant_intron$MAF >=0.01),])[1]
snp_count$counts[9]  <- dim(variant_missense[which(variant_missense$MAF >=0.01),])[1]
snp_count$counts[10] <- dim(variant_non_coding_transcript_exon[which(variant_non_coding_transcript_exon$MAF >=0.01),])[1]
snp_count$counts[11] <- dim(variant_splice_acceptor[which(variant_splice_acceptor$MAF >=0.01),])[1]
snp_count$counts[12] <- dim(variant_splice_donor[which(variant_splice_donor$MAF >=0.01),])[1]
snp_count$counts[13] <- dim(variant_splice_region[which(variant_splice_region$MAF >=0.01),])[1]
snp_count$counts[14] <- dim(variant_start_lost[which(variant_start_lost$MAF >=0.01),])[1]
snp_count$counts[15] <- dim(variant_start_retained[which(variant_start_retained$MAF >=0.01),])[1]
snp_count$counts[16] <- dim(variant_stop_gained[which(variant_stop_gained$MAF >=0.01),])[1]
snp_count$counts[17] <- dim(variant_stop_lost[which(variant_stop_lost$MAF >=0.01),])[1]
snp_count$counts[18] <- dim(variant_stop_retained[which(variant_stop_retained$MAF >=0.01),])[1]
snp_count$counts[19] <- dim(variant_synonymous[which(variant_synonymous$MAF >=0.01),])[1]
snp_count$counts[20] <- dim(variant_upstream[which(variant_upstream$MAF >=0.01),])[1]

snp_count$counts[21]  <- dim(variant_3_prime_UTR[which(variant_3_prime_UTR$MAF <0.01),])[1]
snp_count$counts[22]  <- dim(variant_5_prime_UTR_premature_start_codon_gain[which(variant_5_prime_UTR_premature_start_codon_gain$MAF <0.01),])[1]
snp_count$counts[23]  <- dim(variant_5_prime_UTR[which(variant_5_prime_UTR$MAF <0.01),])[1]
snp_count$counts[24]  <- dim(variant_downstream[which(variant_downstream$MAF <0.01),])[1]
snp_count$counts[25]  <- dim(variant_initiator_codon[which(variant_initiator_codon$MAF <0.01),])[1]
snp_count$counts[26]  <- dim(variant_intergenic_region[which(variant_intergenic_region$MAF <0.01),])[1]
snp_count$counts[27]  <- dim(variant_intragenic[which(variant_intragenic$MAF <0.01),])[1]
snp_count$counts[28]  <- dim(variant_intron[which(variant_intron$MAF <0.01),])[1]
snp_count$counts[29]  <- dim(variant_missense[which(variant_missense$MAF <0.01),])[1]
snp_count$counts[30] <- dim(variant_non_coding_transcript_exon[which(variant_non_coding_transcript_exon$MAF <0.01),])[1]
snp_count$counts[31] <- dim(variant_splice_acceptor[which(variant_splice_acceptor$MAF <0.01),])[1]
snp_count$counts[32] <- dim(variant_splice_donor[which(variant_splice_donor$MAF <0.01),])[1]
snp_count$counts[33] <- dim(variant_splice_region[which(variant_splice_region$MAF <0.01),])[1]
snp_count$counts[34] <- dim(variant_start_lost[which(variant_start_lost$MAF <0.01),])[1]
snp_count$counts[35] <- dim(variant_start_retained[which(variant_start_retained$MAF <0.01),])[1]
snp_count$counts[36] <- dim(variant_stop_gained[which(variant_stop_gained$MAF <0.01),])[1]
snp_count$counts[37] <- dim(variant_stop_lost[which(variant_stop_lost$MAF <0.01),])[1]
snp_count$counts[38] <- dim(variant_stop_retained[which(variant_stop_retained$MAF <0.01),])[1]
snp_count$counts[39] <- dim(variant_synonymous[which(variant_synonymous$MAF <0.01),])[1]
snp_count$counts[40] <- dim(variant_upstream[which(variant_upstream$MAF <0.01),])[1]



snp_count$counts_high_r2[1]  <- dim(good_variant_3_prime_UTR[which(good_variant_3_prime_UTR$MAF >=0.01),])[1]
snp_count$counts_high_r2[2]  <- dim(good_variant_5_prime_UTR_premature_start_codon_gain[which(good_variant_5_prime_UTR_premature_start_codon_gain$MAF >=0.01),])[1]
snp_count$counts_high_r2[3]  <- dim(good_variant_5_prime_UTR[which(good_variant_5_prime_UTR$MAF >=0.01),])[1]
snp_count$counts_high_r2[4]  <- dim(good_variant_downstream[which(good_variant_downstream$MAF >=0.01),])[1]
snp_count$counts_high_r2[5]  <- dim(good_variant_initiator_codon[which(good_variant_initiator_codon$MAF >=0.01),])[1]
snp_count$counts_high_r2[6]  <- dim(good_variant_intergenic_region[which(good_variant_intergenic_region$MAF >=0.01),])[1]
snp_count$counts_high_r2[7]  <- dim(good_variant_intragenic[which(good_variant_intragenic$MAF >=0.01),])[1]
snp_count$counts_high_r2[8]  <- dim(good_variant_intron[which(good_variant_intron$MAF >=0.01),])[1]
snp_count$counts_high_r2[9]  <- dim(good_variant_missense[which(good_variant_missense$MAF >=0.01),])[1]
snp_count$counts_high_r2[10] <- dim(good_variant_non_coding_transcript_exon[which(good_variant_non_coding_transcript_exon$MAF >=0.01),])[1]
snp_count$counts_high_r2[11] <- dim(good_variant_splice_acceptor[which(good_variant_splice_acceptor$MAF >=0.01),])[1]
snp_count$counts_high_r2[12] <- dim(good_variant_splice_donor[which(good_variant_splice_donor$MAF >=0.01),])[1]
snp_count$counts_high_r2[13] <- dim(good_variant_splice_region[which(good_variant_splice_region$MAF >=0.01),])[1]
snp_count$counts_high_r2[14] <- dim(good_variant_start_lost[which(good_variant_start_lost$MAF >=0.01),])[1]
snp_count$counts_high_r2[15] <- dim(good_variant_start_retained[which(good_variant_start_retained$MAF >=0.01),])[1]
snp_count$counts_high_r2[16] <- dim(good_variant_stop_gained[which(good_variant_stop_gained$MAF >=0.01),])[1]
snp_count$counts_high_r2[17] <- dim(good_variant_stop_lost[which(good_variant_stop_lost$MAF >=0.01),])[1]
snp_count$counts_high_r2[18] <- dim(good_variant_stop_retained[which(good_variant_stop_retained$MAF >=0.01),])[1]
snp_count$counts_high_r2[19] <- dim(good_variant_synonymous[which(good_variant_synonymous$MAF >=0.01),])[1]
snp_count$counts_high_r2[20] <- dim(good_variant_upstream[which(good_variant_upstream$MAF >=0.01),])[1]

snp_count$counts_high_r2[21]  <- dim(good_variant_3_prime_UTR[which(good_variant_3_prime_UTR$MAF <0.01),])[1]
snp_count$counts_high_r2[22]  <- dim(good_variant_5_prime_UTR_premature_start_codon_gain[which(good_variant_5_prime_UTR_premature_start_codon_gain$MAF <0.01),])[1]
snp_count$counts_high_r2[23]  <- dim(good_variant_5_prime_UTR[which(good_variant_5_prime_UTR$MAF <0.01),])[1]
snp_count$counts_high_r2[24]  <- dim(good_variant_downstream[which(good_variant_downstream$MAF <0.01),])[1]
snp_count$counts_high_r2[25]  <- dim(good_variant_initiator_codon[which(good_variant_initiator_codon$MAF <0.01),])[1]
snp_count$counts_high_r2[26]  <- dim(good_variant_intergenic_region[which(good_variant_intergenic_region$MAF <0.01),])[1]
snp_count$counts_high_r2[27]  <- dim(good_variant_intragenic[which(good_variant_intragenic$MAF <0.01),])[1]
snp_count$counts_high_r2[28]  <- dim(good_variant_intron[which(good_variant_intron$MAF <0.01),])[1]
snp_count$counts_high_r2[29]  <- dim(good_variant_missense[which(good_variant_missense$MAF <0.01),])[1]
snp_count$counts_high_r2[30] <- dim(good_variant_non_coding_transcript_exon[which(good_variant_non_coding_transcript_exon$MAF <0.01),])[1]
snp_count$counts_high_r2[31] <- dim(good_variant_splice_acceptor[which(good_variant_splice_acceptor$MAF <0.01),])[1]
snp_count$counts_high_r2[32] <- dim(good_variant_splice_donor[which(good_variant_splice_donor$MAF <0.01),])[1]
snp_count$counts_high_r2[33] <- dim(good_variant_splice_region[which(good_variant_splice_region$MAF <0.01),])[1]
snp_count$counts_high_r2[34] <- dim(good_variant_start_lost[which(good_variant_start_lost$MAF <0.01),])[1]
snp_count$counts_high_r2[35] <- dim(good_variant_start_retained[which(good_variant_start_retained$MAF <0.01),])[1]
snp_count$counts_high_r2[36] <- dim(good_variant_stop_gained[which(good_variant_stop_gained$MAF <0.01),])[1]
snp_count$counts_high_r2[37] <- dim(good_variant_stop_lost[which(good_variant_stop_lost$MAF <0.01),])[1]
snp_count$counts_high_r2[38] <- dim(good_variant_stop_retained[which(good_variant_stop_retained$MAF <0.01),])[1]
snp_count$counts_high_r2[39] <- dim(good_variant_synonymous[which(good_variant_synonymous$MAF <0.01),])[1]
snp_count$counts_high_r2[40] <- dim(good_variant_upstream[which(good_variant_upstream$MAF <0.01),])[1]

dim(variant_upstream)
dim(variant_upstream[which(variant_upstream$MAF >=0.01),]) 
dim(variant_upstream[which(variant_upstream$MAF <0.01),]) 


dim(good_variant_upstream)
dim(good_variant_upstream[which(good_variant_upstream$MAF >=0.01),]) 
dim(good_variant_upstream[which(good_variant_upstream$MAF <0.01),]) 


snp_count$prop_high_r2 <- snp_count$counts_high_r2 / snp_count$counts * 100

write.table(snp_count,"snp_type_counts.txt", quote=F, row.names=F)


###########        Plot the number of SNPs by SNP types  ###############
all_meta_imp <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_topmed_1000G/snp_type_counts.txt", header=T) 
all_topmed <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/snp_type_counts.txt", header=T)
all_1000g <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/snp_type_counts.txt", header=T) 
all_1000g_38 <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_38/snp_type_counts.txt", header=T) 
all_1000g_37 <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_37/snp_type_counts.txt", header=T) 
all_ga <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_GA/snp_type_counts.txt", header=T) 

dim(all_topmed) # [1] 40  6
dim(all_1000g) # [1] 40  6
dim(all_ga) # [1] 40  6
dim(all_1000g_37) # [1] 40  6
dim(all_1000g_38) # [1] 40  6
dim(all_meta_imp) # [1] 40  6

all_meta_imp$ref_panel <- "TOPMed_1000G"
all_topmed$ref_panel <- "TOPMed"
all_1000g$ref_panel <- "Expanded_1000G"
all_1000g_38$ref_panel <- "1000G_GRCh38"
all_1000g_37$ref_panel <- "1000G_GRCh37_SAS"
all_ga$ref_panel <- "GenomeAsia_Pilot"

snp_count <- rbind(all_meta_imp, all_topmed, all_1000g, all_1000g_38, all_1000g_37, all_ga) # [1] 240   6

snp_count$snp_type <- gsub('_variant', '',snp_count$snp_type)
snp_count$snp_type <- gsub('_gene', '',snp_count$snp_type)
snp_count$snp_type <- gsub('_region', '',snp_count$snp_type)

snp_count_all <- snp_count[,c(1,2,3,4)]
snp_count_well <- snp_count[,c(1,2,3,5)]

colnames(snp_count_well)[colnames(snp_count_well)=="counts_high_r2"] <- "counts"

snp_count_all$snp_cat <- "all"
snp_count_well$snp_cat <- "well_imp"

snp_count_all$ref_panel <- factor(snp_count_all$ref_panel,  levels=c("TOPMed_1000G","TOPMed","Expanded_1000G","1000G_GRCh38", "1000G_GRCh37_SAS", "GenomeAsia_Pilot"))
snp_count_all$maf_cat <- factor(snp_count_all$maf_cat,  levels=c("Common","Rare"))

snp_count_well$ref_panel <- factor(snp_count_well$ref_panel,  levels=c("TOPMed_1000G","TOPMed","Expanded_1000G","1000G_GRCh38", "1000G_GRCh37_SAS", "GenomeAsia_Pilot"))
snp_count_well$maf_cat <- factor(snp_count_well$maf_cat,  levels=c("Common","Rare"))

snp_count_all_high <- snp_count_all[which(snp_count_all$snp_type %in% c("downstream","intergenic","intron","upstream")),]
snp_count_well_high <- snp_count_well[which(snp_count_well$snp_type %in% c("downstream","intergenic","intron","upstream")),]
snp_count_all_high$snp_type <- factor(snp_count_all_high$snp_type,  levels=c("intron","intergenic","downstream","upstream"))
snp_count_well_high$snp_type <- factor(snp_count_well_high$snp_type,  levels=c("intron","intergenic","downstream","upstream"))

snp_count_all_mid <- snp_count_all[which(snp_count_all$snp_type %in% c("3_prime_UTR","5_prime_UTR","missense","","synonymous")),]
snp_count_well_mid <- snp_count_well[which(snp_count_well$snp_type %in% c("3_prime_UTR","5_prime_UTR","missense","non_coding_transcript_exon","synonymous")),]
snp_count_all_mid$snp_type <- factor(snp_count_all_mid$snp_type,  levels=c("non_coding_transcript_exon","3_prime_UTR","missense","synonymous","5_prime_UTR"))
snp_count_well_mid$snp_type <- factor(snp_count_well_mid$snp_type,  levels=c("non_coding_transcript_exon","3_prime_UTR","missense","synonymous","5_prime_UTR"))

snp_count_all_low <- snp_count_all[which(snp_count_all$snp_type %in% c("5_prime_UTR_premature_start_codon_gain","intragenic","splice_acceptor","splice_donor","splice","stop_gained")),]
snp_count_well_low <- snp_count_well[which(snp_count_well$snp_type %in% c("5_prime_UTR_premature_start_codon_gain","intragenic","splice_acceptor","splice_donor","splice","stop_gained")),]
snp_count_all_low$snp_type <- factor(snp_count_all_low$snp_type,  levels=c("splice","5_prime_UTR_premature_start_codon_gain","intragenic","splice_donor","splice_acceptor","stop_gained"))
snp_count_well_low$snp_type <- factor(snp_count_well_low$snp_type,  levels=c("splice","5_prime_UTR_premature_start_codon_gain","intragenic","splice_donor","splice_acceptor","stop_gained"))

snp_count_all_ultra_low <- snp_count_all[which(snp_count_all$snp_type %in% c("initiator_codon","start_lost","start_retained","stop_lost","stop_retained")),]
snp_count_well_ultra_low <- snp_count_well[which(snp_count_well$snp_type %in% c("initiator_codon","start_lost","start_retained","stop_lost","stop_retained")),]
snp_count_all_ultra_low$snp_type <- factor(snp_count_all_ultra_low$snp_type,  levels=c("start_lost","stop_lost","stop_retained","initiator_codon","start_retained"))
snp_count_well_ultra_low$snp_type <- factor(snp_count_well_ultra_low$snp_type,  levels=c("start_lost","stop_lost","stop_retained","initiator_codon","start_retained"))


library(ggplot2)
pdf('snp_type_all_snp_count.pdf',height=10,width=20)
p <- ggplot(data=snp_count_all, aes(x=interaction(maf_cat,snp_type), y=counts, fill=ref_panel)) + xlab("Variant Type") + ylab("Count") + geom_bar(stat="identity", position=position_dodge()) + theme_classic(base_size=20) + scale_x_discrete(guide = guide_axis(angle = 90)) + scale_y_continuous(labels=scales::label_number_si())+theme(axis.text = element_text(size = 20,colour="black"))
p + scale_fill_manual(values=c("#661100", "#D55E00", "#E69F00","#0072B2","#332288","#999999"),name="Reference Panel",labels=c("TOPMed_1000G","TOPMed","Expanded_1000G","1000G_GRCh38", "1000G_GRCh37_SAS", "GenomeAsia_Pilot"))+ theme(plot.title = element_text(hjust = 0.5)) + labs(title="All SNPs")
dev.off()
# +  geom_text(aes(label=scales::label_number_si(accuracy=0.1)(counts)),size=3, position=position_dodge(width=0.9), vjust=-0.25) 

pdf('snp_type_well_imputed_snp_count.pdf',height=10,width=20)
p <- ggplot(data=snp_count_well, aes(x=interaction(maf_cat,snp_type), y=counts, fill=ref_panel)) + xlab("Variant Type") + ylab("Count") + geom_bar(stat="identity", position=position_dodge()) + theme_classic(base_size=20) + scale_x_discrete(guide = guide_axis(angle = 90)) + scale_y_continuous(labels=scales::label_number_si())+theme(axis.text = element_text(size = 20,colour="black"))
p + scale_fill_manual(values=c("#661100", "#D55E00", "#E69F00","#0072B2","#332288","#999999"),name="Reference Panel",labels=c("TOPMed_1000G","TOPMed","Expanded_1000G","1000G_GRCh38", "1000G_GRCh37_SAS", "GenomeAsia_Pilot"))+ theme(legend.position="none",plot.title = element_text(hjust = 0.5)) + labs(title="Well-imputed SNPs")
dev.off()
# +  geom_text(aes(label=scales::label_number_si(accuracy=0.1)(counts)),size=3, position=position_dodge(width=0.9), vjust=-0.25)

# high counts
pdf('snp_type_high_all_snp_count.pdf',height=10,width=20)
p <- ggplot(data=snp_count_all_high, aes(x=interaction(maf_cat,snp_type), y=counts, fill=ref_panel)) + xlab("Variant Type") + ylab("Count") + geom_bar(stat="identity", position=position_dodge()) +  geom_text(aes(label=scales::label_number_si(accuracy=1)(counts)),size=3, position=position_dodge(width=0.9), vjust=-0.25) + theme_classic(base_size=30) + scale_x_discrete(guide = guide_axis(angle = 90)) + scale_y_continuous(labels=scales::label_number_si())+theme(axis.text = element_text(size = 30,colour="black"))
p + scale_fill_manual(values=c("#661100", "#D55E00", "#E69F00","#0072B2","#332288","#999999"),name="Reference Panel",labels=c("TOPMed_1000G","TOPMed","Expanded_1000G","1000G_GRCh38", "1000G_GRCh37_SAS", "GenomeAsia_Pilot"))+ theme(legend.position="none",plot.title = element_text(hjust = 0.5)) + labs(title="All SNPs")
dev.off()

pdf('snp_type_high_well_imputed_snp_count.pdf',height=10,width=20)
p <- ggplot(data=snp_count_well_high, aes(x=interaction(maf_cat,snp_type), y=counts, fill=ref_panel)) + xlab("Variant Type") + ylab("Count") + geom_bar(stat="identity", position=position_dodge()) +  geom_text(aes(label=scales::label_number_si(accuracy=0.1)(counts)),size=2.5, position=position_dodge(width=0.9), vjust=-0.25) + theme_classic(base_size=30) + scale_x_discrete(guide = guide_axis(angle = 90)) + scale_y_continuous(labels=scales::label_number_si())+theme(axis.text = element_text(size = 30,colour="black"))
p + scale_fill_manual(values=c("#661100", "#D55E00", "#E69F00","#0072B2","#332288","#999999"),name="Reference Panel",labels=c("TOPMed_1000G","TOPMed","Expanded_1000G","1000G_GRCh38", "1000G_GRCh37_SAS", "GenomeAsia_Pilot"))+ theme(legend.position="none",plot.title = element_text(hjust = 0.5)) + labs(title="Well-imputed SNPs")
dev.off()

# middle counts
pdf('snp_type_mid_all_snp_count.pdf',height=10,width=20)
p <- ggplot(data=snp_count_all_mid, aes(x=interaction(maf_cat,snp_type), y=counts, fill=ref_panel)) + xlab("Variant Type") + ylab("Count") + geom_bar(stat="identity", position=position_dodge()) +  geom_text(aes(label=scales::label_number_si(accuracy=1)(counts)),size=3, position=position_dodge(width=0.9), vjust=-0.25) + theme_classic(base_size=30) + scale_x_discrete(guide = guide_axis(angle = 90)) + scale_y_continuous(labels=scales::label_number_si())+theme(axis.text = element_text(size = 30,colour="black"))
p + scale_fill_manual(values=c("#661100", "#D55E00", "#E69F00","#0072B2","#332288","#999999"),name="Reference Panel",labels=c("TOPMed_1000G","TOPMed","Expanded_1000G","1000G_GRCh38", "1000G_GRCh37_SAS", "GenomeAsia_Pilot"))+ theme(legend.position="none",plot.title = element_text(hjust = 0.5)) + labs(title="All SNPs")
dev.off()

pdf('snp_type_mid_well_imputed_snp_count.pdf',height=20,width=40)
p <- ggplot(data=snp_count_well_mid, aes(x=interaction(maf_cat,snp_type), y=counts, fill=ref_panel)) + xlab("Variant Type") + ylab("Count") + geom_bar(stat="identity", position=position_dodge()) +  geom_text(aes(label=scales::label_number_si(accuracy=1)(counts)),size=2.5, position=position_dodge(width=0.9), vjust=-0.25) + theme_classic(base_size=30) + scale_x_discrete(guide = guide_axis(angle = 90)) + scale_y_continuous(labels=scales::label_number_si())+theme(axis.text = element_text(size = 30,colour="black"))
p + scale_fill_manual(values=c("#661100", "#D55E00", "#E69F00","#0072B2","#332288","#999999"),name="Reference Panel",labels=c("TOPMed_1000G","TOPMed","Expanded_1000G","1000G_GRCh38", "1000G_GRCh37_SAS", "GenomeAsia_Pilot"))+ theme(legend.position="none",plot.title = element_text(hjust = 0.5)) + labs(title="Well-imputed SNPs")
dev.off()

# low counts
pdf('snp_type_low_all_snp_count.pdf',height=10,width=20)
p <- ggplot(data=snp_count_all_low, aes(x=interaction(maf_cat,snp_type), y=counts, fill=ref_panel)) + xlab("Variant Type") + ylab("Count") + geom_bar(stat="identity", position=position_dodge()) +  geom_text(aes(label=scales::label_number_si(accuracy=1)(counts)),size=2.5, position=position_dodge(width=0.9), vjust=-0.25) + theme_classic(base_size=20) + scale_x_discrete(guide = guide_axis(angle = 90)) + scale_y_continuous(labels=scales::label_number_si())+theme(axis.text = element_text(size = 20,colour="black"))
p + scale_fill_manual(values=c("#661100", "#D55E00", "#E69F00","#0072B2","#332288","#999999"),name="Reference Panel",labels=c("TOPMed_1000G","TOPMed","Expanded_1000G","1000G_GRCh38", "1000G_GRCh37_SAS", "GenomeAsia_Pilot"))+ theme(legend.position="none",plot.title = element_text(hjust = 0.5)) + labs(title="All SNPs")
dev.off()

pdf('snp_type_low_well_imputed_snp_count.pdf',height=10,width=20)
p <- ggplot(data=snp_count_well_low, aes(x=interaction(maf_cat,snp_type), y=counts, fill=ref_panel)) + xlab("Variant Type") + ylab("Count") + geom_bar(stat="identity", position=position_dodge()) +  geom_text(aes(label=scales::label_number_si(accuracy=1)(counts)),size=2.5, position=position_dodge(width=0.9), vjust=-0.25) + theme_classic(base_size=20) + scale_x_discrete(guide = guide_axis(angle = 90)) + scale_y_continuous(labels=scales::label_number_si())+theme(axis.text = element_text(size = 20,colour="black"))
p + scale_fill_manual(values=c("#661100", "#D55E00", "#E69F00","#0072B2","#332288","#999999"),name="Reference Panel",labels=c("TOPMed_1000G","TOPMed","Expanded_1000G","1000G_GRCh38", "1000G_GRCh37_SAS", "GenomeAsia_Pilot"))+ theme(legend.position="none",plot.title = element_text(hjust = 0.5)) + labs(title="Well-imputed SNPs")
dev.off()

# ultra-low counts
pdf('snp_type_ultra_low_all_snp_count.pdf',height=10,width=20)
p <- ggplot(data=snp_count_all_ultra_low, aes(x=interaction(maf_cat,snp_type), y=counts, fill=ref_panel)) + xlab("Variant Type") + ylab("Count") + geom_bar(stat="identity", position=position_dodge()) +  geom_text(aes(label=scales::label_number_si(accuracy=1)(counts)),size=2.5, position=position_dodge(width=0.9), vjust=-0.25) + theme_classic(base_size=30) + scale_x_discrete(guide = guide_axis(angle = 90)) + scale_y_continuous(labels=scales::label_number_si())+theme(axis.text = element_text(size = 30,colour="black"))
p + scale_fill_manual(values=c("#661100", "#D55E00", "#E69F00","#0072B2","#332288","#999999"),name="Reference Panel",labels=c("TOPMed_1000G","TOPMed","Expanded_1000G","1000G_GRCh38", "1000G_GRCh37_SAS", "GenomeAsia_Pilot"))+ theme(legend.position="none",plot.title = element_text(hjust = 0.5)) + labs(title="All SNPs")
dev.off()

pdf('snp_type_ultra_low_well_imputed_snp_count.pdf',height=10,width=20)
p <- ggplot(data=snp_count_well_ultra_low, aes(x=interaction(maf_cat,snp_type), y=counts, fill=ref_panel)) + xlab("Variant Type") + ylab("Count") + geom_bar(stat="identity", position=position_dodge()) +  geom_text(aes(label=scales::label_number_si(accuracy=1)(counts)),size=3, position=position_dodge(width=0.9), vjust=-0.25) + theme_classic(base_size=30) + scale_x_discrete(guide = guide_axis(angle = 90)) + scale_y_continuous(labels=scales::label_number_si())+theme(axis.text = element_text(size = 30,colour="black"))
p + scale_fill_manual(values=c("#661100", "#D55E00", "#E69F00","#0072B2","#332288","#999999"),name="Reference Panel",labels=c("TOPMed_1000G","TOPMed","Expanded_1000G","1000G_GRCh38", "1000G_GRCh37_SAS", "GenomeAsia_Pilot"))+ theme(legend.position="none",plot.title = element_text(hjust = 0.5)) + labs(title="Well-imputed SNPs")
dev.off()



#######################################################################
####             Aggregated Rsq (ultra rare, rare, common)         ####
#######################################################################
# need to change to vcf v4.2 to use vcftools later
gunzip -c  Pakistani_2200_tarSeq.vcf.gz | sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/'  > Pakistani_2200_tarSeq.v4.2.vcf 
# from Dongjing
# 49263 - 28 headers = 49235 SNPs

# change the sample IDs in the tarseq data
ml bcftools
bcftools reheader -s /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/genetic_match_between_TarSeq_and_GSA_2200.tsv  Pakistani_2200_tarSeq.v4.2.vcf  >  Pakistani_2200_tarSeq.v4.2.updateID.vcf

## only keep the 1814 individuals for tarseq data
ml vcftools
vcftools --vcf Pakistani_2200_tarSeq.v4.2.updateID.vcf --recode --recode-INFO-all --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_2433_unphased/pakistan_gsa_included_1814_samples.tsv --out Pakistani_1748_tarSeq.v4.2.updateID
#  49235 variants,  1748 out of 1814 individuals are kept

ml vcftools
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.recode.vcf --chr 1 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.chr1 # 5765 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.recode.vcf --chr 2 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.chr2 # 3308 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.recode.vcf --chr 3 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.chr3 # 3817 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.recode.vcf --chr 4 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.chr4 # 1760 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.recode.vcf --chr 5 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.chr5 # 1420 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.recode.vcf --chr 6 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.chr6 # 2576 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.recode.vcf --chr 7 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.chr7 # 2528 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.recode.vcf --chr 8 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.chr8 # 606 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.recode.vcf --chr 9 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.chr9 # 2379 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.recode.vcf --chr 10 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.chr10 # 1353 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.recode.vcf --chr 11 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.chr11 # 3282 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.recode.vcf --chr 12 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.chr12 # 1265 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.recode.vcf --chr 13 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.chr13 # 1002 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.recode.vcf --chr 14 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.chr14 # 1680 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.recode.vcf --chr 15 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.chr15 # 2745 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.recode.vcf --chr 16 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.chr16 # 2530 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.recode.vcf --chr 17 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.chr17 # 4608 SNPs (GRCh38 version has 4607 SNPs)
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.recode.vcf --chr 18 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.chr18 # 888 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.recode.vcf --chr 19 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.chr19 # 4033 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.recode.vcf --chr 20 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.chr20 # 675 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.recode.vcf --chr 21 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.chr21 # 275 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.recode.vcf --chr 22 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.chr22 # 740 SNPs

grep -v '^#'  Pakistani_1748_tarSeq.chr1.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.chr2.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.chr3.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.chr4.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.chr5.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.chr6.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.chr7.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.chr8.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.chr9.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.chr10.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.chr11.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.chr12.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.chr13.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.chr14.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.chr15.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.chr16.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.chr17.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.chr18.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.chr19.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.chr20.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.chr21.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.chr22.recode.vcf|wc -l


# liftover the vcf from grch37 to grch38 for TOPMed, expanded 1000G and 1000G GRCh38
# https://crossmap.sourceforge.net

# check the manual
CrossMap.py vcf -h

# template script: CrossMap.py vcf input.chain input.vcf refgenome.fa out_vcf
# convert tarseq data from GRCh37 to GRCh38
ml samtools
CrossMap.py vcf /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/hg19ToHg38.over.chain.gz Pakistani_2200_tarSeq.v4.2.updateID.vcf /sc/arion/projects/psychgen/resources/annovar/humandb/hg38_seq/hg38.fa Pakistani_2200_tarSeq.v4.2.updateID.GRCh38.vcf
# 2023-06-07 02:42:03 [INFO]  Total entries: 49235
# 2023-06-07 02:42:03 [INFO]  Failed to map: 1

# grch38 version vcf: 464 headers
wc -l Pakistani_2200_tarSeq.v4.2.updateID.GRCh38.vcf
# 49698 - 464 headers = 49234 SNPs

# after conversion (grch38):  1	6106251	1:6166311:A:G	G	A	.	.	PR	GT	0/0
# before conversion (grch37): 1	6166311	1:6166311:A:G	G	A	.	.	PR	GT	0/0
# position is changed, but not SNP IDs
# rs140903833: 
# 1:6106251 (GRCh38)
# 1:6166311 (GRCh37)

## add chr in the chromosome column
awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' Pakistani_2200_tarSeq.v4.2.updateID.GRCh38.vcf > Pakistani_2200_tarSeq.v4.2.updateID.GRCh38.chr.vcf
# before: 1  	6106251	1:6166311:A:G	G	A	.	.	PR	GT	0/0
# after: chr1	6106251	1:6166311:A:G	G	A	.	.	PR	GT	0/0

## only keep the 1814 individuals for tarseq data
vcftools --vcf Pakistani_2200_tarSeq.v4.2.updateID.GRCh38.chr.vcf --recode --recode-INFO-all --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_2433_unphased/pakistan_gsa_included_1814_samples.tsv --out Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr
# 49234 SNPs,  1748 out of 1814 individuals are kept

# save the 1748 samples 
awk '(NR==464){print $0}' Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.recode.vcf | sed 's/\t/\n/g' > keep_1748_samples
# remove the first 9 columns manually 

ml vcftools
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.recode.vcf --chr chr1 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.GRCh38.chr1 # 5765 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.recode.vcf --chr chr2 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.GRCh38.chr2 # 3308 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.recode.vcf --chr chr3 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.GRCh38.chr3 # 3817 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.recode.vcf --chr chr4 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.GRCh38.chr4 # 1760 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.recode.vcf --chr chr5 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.GRCh38.chr5 # 1420 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.recode.vcf --chr chr6 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.GRCh38.chr6 # 2576 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.recode.vcf --chr chr7 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.GRCh38.chr7 # 2528 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.recode.vcf --chr chr8 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.GRCh38.chr8 # 606 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.recode.vcf --chr chr9 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.GRCh38.chr9 # 2379 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.recode.vcf --chr chr10 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.GRCh38.chr10 # 1353 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.recode.vcf --chr chr11 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.GRCh38.chr11 # 3282 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.recode.vcf --chr chr12 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.GRCh38.chr12 # 1265 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.recode.vcf --chr chr13 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.GRCh38.chr13 # 1002 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.recode.vcf --chr chr14 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.GRCh38.chr14 # 1680 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.recode.vcf --chr chr15 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.GRCh38.chr15 # 2745 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.recode.vcf --chr chr16 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.GRCh38.chr16 # 2530 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.recode.vcf --chr chr17 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.GRCh38.chr17 # 4607 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.recode.vcf --chr chr18 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.GRCh38.chr18 # 888 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.recode.vcf --chr chr19 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.GRCh38.chr19 # 4033 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.recode.vcf --chr chr20 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.GRCh38.chr20 # 675 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.recode.vcf --chr chr21 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.GRCh38.chr21 # 275 SNPs
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.recode.vcf --chr chr22 --recode --recode-INFO-all --out Pakistani_1748_tarSeq.GRCh38.chr22 # 740 SNPs

grep -v '^#'  Pakistani_1748_tarSeq.GRCh38.chr1.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.GRCh38.chr2.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.GRCh38.chr3.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.GRCh38.chr4.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.GRCh38.chr5.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.GRCh38.chr6.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.GRCh38.chr7.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.GRCh38.chr8.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.GRCh38.chr9.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.GRCh38.chr10.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.GRCh38.chr11.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.GRCh38.chr12.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.GRCh38.chr13.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.GRCh38.chr14.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.GRCh38.chr15.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.GRCh38.chr16.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.GRCh38.chr17.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.GRCh38.chr18.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.GRCh38.chr19.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.GRCh38.chr20.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.GRCh38.chr21.recode.vcf|wc -l
grep -v '^#'  Pakistani_1748_tarSeq.GRCh38.chr22.recode.vcf|wc -l
### the ready-to-use tarseq data file are (1) for GRCh37: Pakistani_1748_tarSeq.v4.2.updateID.recode.vcf , (2) for GRCh38: Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.recode.vcf 

# GRCh37 tarseq SNP list per chromosome
sed '/^#/d' Pakistani_1748_tarSeq.chr1.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.chr1.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.chr1.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.chr1.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.chr1.recode.vcf.forward.snplist Pakistani_1748_tarSeq.chr1.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.chr1.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.chr2.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.chr2.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.chr2.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.chr2.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.chr2.recode.vcf.forward.snplist Pakistani_1748_tarSeq.chr2.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.chr2.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.chr3.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.chr3.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.chr3.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.chr3.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.chr3.recode.vcf.forward.snplist Pakistani_1748_tarSeq.chr3.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.chr3.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.chr4.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.chr4.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.chr4.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.chr4.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.chr4.recode.vcf.forward.snplist Pakistani_1748_tarSeq.chr4.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.chr4.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.chr5.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.chr5.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.chr5.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.chr5.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.chr5.recode.vcf.forward.snplist Pakistani_1748_tarSeq.chr5.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.chr5.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.chr6.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.chr6.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.chr6.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.chr6.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.chr6.recode.vcf.forward.snplist Pakistani_1748_tarSeq.chr6.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.chr6.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.chr7.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.chr7.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.chr7.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.chr7.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.chr7.recode.vcf.forward.snplist Pakistani_1748_tarSeq.chr7.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.chr7.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.chr8.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.chr8.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.chr8.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.chr8.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.chr8.recode.vcf.forward.snplist Pakistani_1748_tarSeq.chr8.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.chr8.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.chr9.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.chr9.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.chr9.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.chr9.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.chr9.recode.vcf.forward.snplist Pakistani_1748_tarSeq.chr9.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.chr9.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.chr10.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.chr10.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.chr10.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.chr10.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.chr10.recode.vcf.forward.snplist Pakistani_1748_tarSeq.chr10.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.chr10.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.chr11.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.chr11.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.chr11.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.chr11.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.chr11.recode.vcf.forward.snplist Pakistani_1748_tarSeq.chr11.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.chr11.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.chr12.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.chr12.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.chr12.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.chr12.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.chr12.recode.vcf.forward.snplist Pakistani_1748_tarSeq.chr12.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.chr12.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.chr13.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.chr13.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.chr13.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.chr13.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.chr13.recode.vcf.forward.snplist Pakistani_1748_tarSeq.chr13.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.chr13.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.chr14.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.chr14.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.chr14.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.chr14.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.chr14.recode.vcf.forward.snplist Pakistani_1748_tarSeq.chr14.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.chr14.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.chr15.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.chr15.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.chr15.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.chr15.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.chr15.recode.vcf.forward.snplist Pakistani_1748_tarSeq.chr15.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.chr15.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.chr16.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.chr16.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.chr16.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.chr16.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.chr16.recode.vcf.forward.snplist Pakistani_1748_tarSeq.chr16.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.chr16.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.chr17.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.chr17.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.chr17.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.chr17.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.chr17.recode.vcf.forward.snplist Pakistani_1748_tarSeq.chr17.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.chr17.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.chr18.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.chr18.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.chr18.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.chr18.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.chr18.recode.vcf.forward.snplist Pakistani_1748_tarSeq.chr18.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.chr18.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.chr19.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.chr19.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.chr19.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.chr19.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.chr19.recode.vcf.forward.snplist Pakistani_1748_tarSeq.chr19.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.chr19.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.chr20.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.chr20.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.chr20.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.chr20.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.chr20.recode.vcf.forward.snplist Pakistani_1748_tarSeq.chr20.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.chr20.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.chr21.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.chr21.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.chr21.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.chr21.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.chr21.recode.vcf.forward.snplist Pakistani_1748_tarSeq.chr21.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.chr21.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.chr22.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.chr22.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.chr22.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.chr22.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.chr22.recode.vcf.forward.snplist Pakistani_1748_tarSeq.chr22.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.chr22.recode.vcf.snplist

# GRCh38 tarseq SNP list per chromosome
sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr1.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.GRCh38.chr1.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr1.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.GRCh38.chr1.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.GRCh38.chr1.recode.vcf.forward.snplist Pakistani_1748_tarSeq.GRCh38.chr1.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.GRCh38.chr1.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr2.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.GRCh38.chr2.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr2.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.GRCh38.chr2.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.GRCh38.chr2.recode.vcf.forward.snplist Pakistani_1748_tarSeq.GRCh38.chr2.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.GRCh38.chr2.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr3.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.GRCh38.chr3.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr3.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.GRCh38.chr3.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.GRCh38.chr3.recode.vcf.forward.snplist Pakistani_1748_tarSeq.GRCh38.chr3.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.GRCh38.chr3.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr4.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.GRCh38.chr4.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr4.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.GRCh38.chr4.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.GRCh38.chr4.recode.vcf.forward.snplist Pakistani_1748_tarSeq.GRCh38.chr4.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.GRCh38.chr4.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr5.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.GRCh38.chr5.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr5.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.GRCh38.chr5.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.GRCh38.chr5.recode.vcf.forward.snplist Pakistani_1748_tarSeq.GRCh38.chr5.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.GRCh38.chr5.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr6.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.GRCh38.chr6.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr6.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.GRCh38.chr6.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.GRCh38.chr6.recode.vcf.forward.snplist Pakistani_1748_tarSeq.GRCh38.chr6.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.GRCh38.chr6.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr7.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.GRCh38.chr7.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr7.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.GRCh38.chr7.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.GRCh38.chr7.recode.vcf.forward.snplist Pakistani_1748_tarSeq.GRCh38.chr7.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.GRCh38.chr7.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr8.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.GRCh38.chr8.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr8.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.GRCh38.chr8.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.GRCh38.chr8.recode.vcf.forward.snplist Pakistani_1748_tarSeq.GRCh38.chr8.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.GRCh38.chr8.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr9.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.GRCh38.chr9.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr9.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.GRCh38.chr9.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.GRCh38.chr9.recode.vcf.forward.snplist Pakistani_1748_tarSeq.GRCh38.chr9.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.GRCh38.chr9.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr10.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.GRCh38.chr10.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr10.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.GRCh38.chr10.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.GRCh38.chr10.recode.vcf.forward.snplist Pakistani_1748_tarSeq.GRCh38.chr10.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.GRCh38.chr10.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr11.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.GRCh38.chr11.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr11.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.GRCh38.chr11.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.GRCh38.chr11.recode.vcf.forward.snplist Pakistani_1748_tarSeq.GRCh38.chr11.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.GRCh38.chr11.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr12.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.GRCh38.chr12.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr12.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.GRCh38.chr12.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.GRCh38.chr12.recode.vcf.forward.snplist Pakistani_1748_tarSeq.GRCh38.chr12.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.GRCh38.chr12.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr13.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.GRCh38.chr13.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr13.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.GRCh38.chr13.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.GRCh38.chr13.recode.vcf.forward.snplist Pakistani_1748_tarSeq.GRCh38.chr13.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.GRCh38.chr13.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr14.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.GRCh38.chr14.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr14.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.GRCh38.chr14.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.GRCh38.chr14.recode.vcf.forward.snplist Pakistani_1748_tarSeq.GRCh38.chr14.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.GRCh38.chr14.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr15.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.GRCh38.chr15.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr15.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.GRCh38.chr15.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.GRCh38.chr15.recode.vcf.forward.snplist Pakistani_1748_tarSeq.GRCh38.chr15.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.GRCh38.chr15.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr16.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.GRCh38.chr16.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr16.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.GRCh38.chr16.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.GRCh38.chr16.recode.vcf.forward.snplist Pakistani_1748_tarSeq.GRCh38.chr16.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.GRCh38.chr16.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr17.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.GRCh38.chr17.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr17.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.GRCh38.chr17.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.GRCh38.chr17.recode.vcf.forward.snplist Pakistani_1748_tarSeq.GRCh38.chr17.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.GRCh38.chr17.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr18.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.GRCh38.chr18.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr18.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.GRCh38.chr18.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.GRCh38.chr18.recode.vcf.forward.snplist Pakistani_1748_tarSeq.GRCh38.chr18.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.GRCh38.chr18.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr19.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.GRCh38.chr19.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr19.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.GRCh38.chr19.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.GRCh38.chr19.recode.vcf.forward.snplist Pakistani_1748_tarSeq.GRCh38.chr19.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.GRCh38.chr19.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr20.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.GRCh38.chr20.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr20.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.GRCh38.chr20.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.GRCh38.chr20.recode.vcf.forward.snplist Pakistani_1748_tarSeq.GRCh38.chr20.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.GRCh38.chr20.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr21.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.GRCh38.chr21.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr21.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.GRCh38.chr21.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.GRCh38.chr21.recode.vcf.forward.snplist Pakistani_1748_tarSeq.GRCh38.chr21.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.GRCh38.chr21.recode.vcf.snplist

sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr22.recode.vcf |awk '{print $1":"$2":"$4":"$5}' > Pakistani_1748_tarSeq.GRCh38.chr22.recode.vcf.forward.snplist
sed '/^#/d' Pakistani_1748_tarSeq.GRCh38.chr22.recode.vcf |awk '{print $1":"$2":"$5":"$4}' > Pakistani_1748_tarSeq.GRCh38.chr22.recode.vcf.reverse.snplist
cat Pakistani_1748_tarSeq.GRCh38.chr22.recode.vcf.forward.snplist Pakistani_1748_tarSeq.GRCh38.chr22.recode.vcf.reverse.snplist > Pakistani_1748_tarSeq.GRCh38.chr22.recode.vcf.snplist

######################################################################################
######       Extract SNPs from vcf to those that overlap with tarseq SNPs      #######
######################################################################################


####    Meta-imputation  ####

#!/bin/bash
#BSUB -J pick_snp_overlap_with_tarseq # Job name
#BSUB -n 1
#BSUB -P acc_psychgen # allocation account
#BSUB -q premium # queue
#BSUB -W 48:00 # walltime in HH:MM
#BSUB -R himem
#BSUB -R rusage[mem=96G] # 120 GB of memory 
#BSUB -R span[hosts=1] # all cores from one node
#BSUB -o %J.stdout # output log (%J : JobID)
#BSUB -eo %J.stderr # error log
#BSUB -L /bin/bash # Initialize the execution environment


# limit to 1748 individuals
ml vcftools

vcftools --gzvcf chr1_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr1.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr1_topmed_1000g.meta.tarseq.forward.1748 
# 5751 rows - 16 headers = 5735 SNPs
vcftools --gzvcf chr1_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr1.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr1_topmed_1000g.meta.tarseq.reverse.1748 
# 22 rows - 16 header = 6 SNPs matched with the reverse order of tarseq ()

vcftools --gzvcf chr2_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr2.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr2_topmed_1000g.meta.tarseq.forward.1748 
vcftools --gzvcf chr2_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr2.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr2_topmed_1000g.meta.tarseq.reverse.1748 

vcftools --gzvcf chr3_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr3.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr3_topmed_1000g.meta.tarseq.forward.1748 
vcftools --gzvcf chr3_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr3.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr3_topmed_1000g.meta.tarseq.reverse.1748 

vcftools --gzvcf chr4_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr4.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr4_topmed_1000g.meta.tarseq.forward.1748 
vcftools --gzvcf chr4_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr4.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr4_topmed_1000g.meta.tarseq.reverse.1748 

vcftools --gzvcf chr5_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr5.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr5_topmed_1000g.meta.tarseq.forward.1748 
vcftools --gzvcf chr5_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr5.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr5_topmed_1000g.meta.tarseq.reverse.1748 

vcftools --gzvcf chr6_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr6.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr6_topmed_1000g.meta.tarseq.forward.1748 
vcftools --gzvcf chr6_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr6.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr6_topmed_1000g.meta.tarseq.reverse.1748 

vcftools --gzvcf chr7_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr7.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr7_topmed_1000g.meta.tarseq.forward.1748 
vcftools --gzvcf chr7_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr7.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr7_topmed_1000g.meta.tarseq.reverse.1748 

vcftools --gzvcf chr8_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr8.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr8_topmed_1000g.meta.tarseq.forward.1748 
vcftools --gzvcf chr8_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr8.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr8_topmed_1000g.meta.tarseq.reverse.1748 

vcftools --gzvcf chr9_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr9.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr9_topmed_1000g.meta.tarseq.forward.1748 
vcftools --gzvcf chr9_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr9.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr9_topmed_1000g.meta.tarseq.reverse.1748 

vcftools --gzvcf chr10_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr10.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr10_topmed_1000g.meta.tarseq.forward.1748 
vcftools --gzvcf chr10_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr10.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr10_topmed_1000g.meta.tarseq.reverse.1748 

vcftools --gzvcf chr11_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr11.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr11_topmed_1000g.meta.tarseq.forward.1748 
vcftools --gzvcf chr11_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr11.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr11_topmed_1000g.meta.tarseq.reverse.1748 

vcftools --gzvcf chr12_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr12.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr12_topmed_1000g.meta.tarseq.forward.1748 
vcftools --gzvcf chr12_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr12.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr12_topmed_1000g.meta.tarseq.reverse.1748 

vcftools --gzvcf chr13_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr13.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr13_topmed_1000g.meta.tarseq.forward.1748 
vcftools --gzvcf chr13_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr13.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr13_topmed_1000g.meta.tarseq.reverse.1748 

vcftools --gzvcf chr14_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr14.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr14_topmed_1000g.meta.tarseq.forward.1748 
vcftools --gzvcf chr14_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr14.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr14_topmed_1000g.meta.tarseq.reverse.1748 

vcftools --gzvcf chr15_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr15.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr15_topmed_1000g.meta.tarseq.forward.1748 
vcftools --gzvcf chr15_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr15.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr15_topmed_1000g.meta.tarseq.reverse.1748 

vcftools --gzvcf chr16_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr16.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr16_topmed_1000g.meta.tarseq.forward.1748 
vcftools --gzvcf chr16_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr16.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr16_topmed_1000g.meta.tarseq.reverse.1748 

vcftools --gzvcf chr17_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr17.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr17_topmed_1000g.meta.tarseq.forward.1748 
vcftools --gzvcf chr17_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr17.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr17_topmed_1000g.meta.tarseq.reverse.1748 

vcftools --gzvcf chr18_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr18.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr18_topmed_1000g.meta.tarseq.forward.1748 
vcftools --gzvcf chr18_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr18.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr18_topmed_1000g.meta.tarseq.reverse.1748 

vcftools --gzvcf chr19_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr19.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr19_topmed_1000g.meta.tarseq.forward.1748 
vcftools --gzvcf chr19_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr19.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr19_topmed_1000g.meta.tarseq.reverse.1748 

vcftools --gzvcf chr20_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr20.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr20_topmed_1000g.meta.tarseq.forward.1748 
vcftools --gzvcf chr20_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr20.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr20_topmed_1000g.meta.tarseq.reverse.1748 

vcftools --gzvcf chr21_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr21.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr21_topmed_1000g.meta.tarseq.forward.1748 
vcftools --gzvcf chr21_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr21.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr21_topmed_1000g.meta.tarseq.reverse.1748 

vcftools --gzvcf chr22_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr22.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr22_topmed_1000g.meta.tarseq.forward.1748 
vcftools --gzvcf chr22_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr22.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr22_topmed_1000g.meta.tarseq.reverse.1748 

bsub < pick_snp_overlap_with_tarseq.sh

vcftools --gzvcf chr1_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr1.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr1_topmed_1000g.meta.tarseq.reverse.1748 
vcftools --gzvcf chr2_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr2.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr2_topmed_1000g.meta.tarseq.reverse.1748 
vcftools --gzvcf chr3_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr3.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr3_topmed_1000g.meta.tarseq.reverse.1748 
vcftools --gzvcf chr4_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr4.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr4_topmed_1000g.meta.tarseq.reverse.1748 
vcftools --gzvcf chr5_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr5.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr5_topmed_1000g.meta.tarseq.reverse.1748 
vcftools --gzvcf chr6_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr6.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr6_topmed_1000g.meta.tarseq.reverse.1748 
vcftools --gzvcf chr7_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr7.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr7_topmed_1000g.meta.tarseq.reverse.1748 
vcftools --gzvcf chr8_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr8.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr8_topmed_1000g.meta.tarseq.reverse.1748 
vcftools --gzvcf chr9_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr9.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr9_topmed_1000g.meta.tarseq.reverse.1748 
vcftools --gzvcf chr10_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr10.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr10_topmed_1000g.meta.tarseq.reverse.1748 
vcftools --gzvcf chr11_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr11.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr11_topmed_1000g.meta.tarseq.reverse.1748 
vcftools --gzvcf chr12_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr12.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr12_topmed_1000g.meta.tarseq.reverse.1748 
vcftools --gzvcf chr13_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr13.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr13_topmed_1000g.meta.tarseq.reverse.1748 
vcftools --gzvcf chr14_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr14.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr14_topmed_1000g.meta.tarseq.reverse.1748 
vcftools --gzvcf chr15_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr15.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr15_topmed_1000g.meta.tarseq.reverse.1748 
vcftools --gzvcf chr16_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr16.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr16_topmed_1000g.meta.tarseq.reverse.1748 
vcftools --gzvcf chr17_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr17.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr17_topmed_1000g.meta.tarseq.reverse.1748 
vcftools --gzvcf chr18_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr18.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr18_topmed_1000g.meta.tarseq.reverse.1748 
vcftools --gzvcf chr19_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr19.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr19_topmed_1000g.meta.tarseq.reverse.1748 
vcftools --gzvcf chr20_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr20.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr20_topmed_1000g.meta.tarseq.reverse.1748 
vcftools --gzvcf chr21_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr21.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr21_topmed_1000g.meta.tarseq.reverse.1748 
vcftools --gzvcf chr22_topmed_1000g.meta.metaDose.vcf.gz  --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr22.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr22_topmed_1000g.meta.tarseq.reverse.1748 


####   TOPMed, Expanded 1000G, 1000G GRCh38  ####

#!/bin/bash
#BSUB -J pick_snp_overlap_with_tarseq # Job name
#BSUB -n 1
#BSUB -P acc_psychgen # allocation account
#BSUB -q premium # queue
#BSUB -W 48:00 # walltime in HH:MM
#BSUB -R himem
#BSUB -R rusage[mem=96G] # 120 GB of memory 
#BSUB -R span[hosts=1] # all cores from one node
#BSUB -o %J.stdout # output log (%J : JobID)
#BSUB -eo %J.stderr # error log
#BSUB -L /bin/bash # Initialize the execution environment


# limit to 1748 individuals
ml vcftools

vcftools --gzvcf chr1.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr1.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr1.tarseq.forward.1748 
vcftools --gzvcf chr1.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr1.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr1.tarseq.reverse.1748 

vcftools --gzvcf chr2.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr2.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr2.tarseq.forward.1748 
vcftools --gzvcf chr2.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr2.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr2.tarseq.reverse.1748 

vcftools --gzvcf chr3.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr3.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr3.tarseq.forward.1748 
vcftools --gzvcf chr3.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr3.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr3.tarseq.reverse.1748 

vcftools --gzvcf chr4.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr4.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr4.tarseq.forward.1748 
vcftools --gzvcf chr4.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr4.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr4.tarseq.reverse.1748 

vcftools --gzvcf chr5.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr5.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr5.tarseq.forward.1748 
vcftools --gzvcf chr5.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr5.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr5.tarseq.reverse.1748 

vcftools --gzvcf chr6.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr6.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr6.tarseq.forward.1748 
vcftools --gzvcf chr6.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr6.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr6.tarseq.reverse.1748 

vcftools --gzvcf chr7.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr7.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr7.tarseq.forward.1748 
vcftools --gzvcf chr7.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr7.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr7.tarseq.reverse.1748 

vcftools --gzvcf chr8.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr8.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr8.tarseq.forward.1748 
vcftools --gzvcf chr8.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr8.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr8.tarseq.reverse.1748 

vcftools --gzvcf chr9.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr9.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr9.tarseq.forward.1748 
vcftools --gzvcf chr9.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr9.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr9.tarseq.reverse.1748 

vcftools --gzvcf chr10.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr10.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr10.tarseq.forward.1748 
vcftools --gzvcf chr10.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr10.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr10.tarseq.reverse.1748 

vcftools --gzvcf chr11.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr11.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr11.tarseq.forward.1748 
vcftools --gzvcf chr11.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr11.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr11.tarseq.reverse.1748 

vcftools --gzvcf chr12.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr12.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr12.tarseq.forward.1748 
vcftools --gzvcf chr12.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr12.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr12.tarseq.reverse.1748 

vcftools --gzvcf chr13.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr13.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr13.tarseq.forward.1748 
vcftools --gzvcf chr13.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr13.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr13.tarseq.reverse.1748 

vcftools --gzvcf chr14.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr14.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr14.tarseq.forward.1748 
vcftools --gzvcf chr14.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr14.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr14.tarseq.reverse.1748 

vcftools --gzvcf chr15.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr15.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr15.tarseq.forward.1748 
vcftools --gzvcf chr15.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr15.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr15.tarseq.reverse.1748 

vcftools --gzvcf chr16.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr16.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr16.tarseq.forward.1748 
vcftools --gzvcf chr16.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr16.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr16.tarseq.reverse.1748 

vcftools --gzvcf chr17.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr17.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr17.tarseq.forward.1748 
vcftools --gzvcf chr17.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr17.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr17.tarseq.reverse.1748 

vcftools --gzvcf chr18.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr18.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr18.tarseq.forward.1748 
vcftools --gzvcf chr18.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr18.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr18.tarseq.reverse.1748 

vcftools --gzvcf chr19.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr19.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr19.tarseq.forward.1748 
vcftools --gzvcf chr19.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr19.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr19.tarseq.reverse.1748 

vcftools --gzvcf chr20.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr20.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr20.tarseq.forward.1748 
vcftools --gzvcf chr20.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr20.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr20.tarseq.reverse.1748 

vcftools --gzvcf chr21.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr21.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr21.tarseq.forward.1748 
vcftools --gzvcf chr21.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr21.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr21.tarseq.reverse.1748 

vcftools --gzvcf chr22.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr22.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr22.tarseq.forward.1748 
vcftools --gzvcf chr22.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr22.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr22.tarseq.reverse.1748 

bsub < pick_snp_overlap_with_tarseq.sh

vcftools --gzvcf chr1.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr1.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr1.tarseq.reverse.1748 
vcftools --gzvcf chr2.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr2.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr2.tarseq.reverse.1748 
vcftools --gzvcf chr3.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr3.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr3.tarseq.reverse.1748 
vcftools --gzvcf chr4.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr4.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr4.tarseq.reverse.1748 
vcftools --gzvcf chr5.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr5.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr5.tarseq.reverse.1748 
vcftools --gzvcf chr6.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr6.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr6.tarseq.reverse.1748 
vcftools --gzvcf chr7.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr7.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr7.tarseq.reverse.1748 
vcftools --gzvcf chr8.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr8.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr8.tarseq.reverse.1748 
vcftools --gzvcf chr9.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr9.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr9.tarseq.reverse.1748 
vcftools --gzvcf chr10.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr10.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr10.tarseq.reverse.1748 
vcftools --gzvcf chr11.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr11.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr11.tarseq.reverse.1748 
vcftools --gzvcf chr12.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr12.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr12.tarseq.reverse.1748 
vcftools --gzvcf chr13.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr13.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr13.tarseq.reverse.1748 
vcftools --gzvcf chr14.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr14.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr14.tarseq.reverse.1748 
vcftools --gzvcf chr15.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr15.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr15.tarseq.reverse.1748 
vcftools --gzvcf chr16.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr16.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr16.tarseq.reverse.1748 
vcftools --gzvcf chr17.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr17.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr17.tarseq.reverse.1748 
vcftools --gzvcf chr18.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr18.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr18.tarseq.reverse.1748 
vcftools --gzvcf chr19.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr19.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr19.tarseq.reverse.1748 
vcftools --gzvcf chr20.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr20.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr20.tarseq.reverse.1748 
vcftools --gzvcf chr21.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr21.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr21.tarseq.reverse.1748 
vcftools --gzvcf chr22.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr22.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr22.tarseq.reverse.1748 


####   1000G GRCh37 SAS, Genome Asia  ####

#!/bin/bash
#BSUB -J pick_snp_overlap_with_tarseq # Job name
#BSUB -n 1
#BSUB -P acc_psychgen # allocation account
#BSUB -q premium # queue
#BSUB -W 48:00 # walltime in HH:MM
#BSUB -R himem
#BSUB -R rusage[mem=96G] # 120 GB of memory 
#BSUB -R span[hosts=1] # all cores from one node
#BSUB -o %J.stdout # output log (%J : JobID)
#BSUB -eo %J.stderr # error log
#BSUB -L /bin/bash # Initialize the execution environment


# limit to 1748 individuals
ml vcftools

vcftools --gzvcf chr1.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr1.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr1.tarseq.forward.1748 
vcftools --gzvcf chr1.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr1.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr1.tarseq.reverse.1748 

vcftools --gzvcf chr2.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr2.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr2.tarseq.forward.1748 
vcftools --gzvcf chr2.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr2.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr2.tarseq.reverse.1748 

vcftools --gzvcf chr3.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr3.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr3.tarseq.forward.1748 
vcftools --gzvcf chr3.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr3.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr3.tarseq.reverse.1748 

vcftools --gzvcf chr4.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr4.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr4.tarseq.forward.1748 
vcftools --gzvcf chr4.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr4.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr4.tarseq.reverse.1748 

vcftools --gzvcf chr5.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr5.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr5.tarseq.forward.1748 
vcftools --gzvcf chr5.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr5.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr5.tarseq.reverse.1748 

vcftools --gzvcf chr6.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr6.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr6.tarseq.forward.1748 
vcftools --gzvcf chr6.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr6.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr6.tarseq.reverse.1748 

vcftools --gzvcf chr7.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr7.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr7.tarseq.forward.1748 
vcftools --gzvcf chr7.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr7.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr7.tarseq.reverse.1748 

vcftools --gzvcf chr8.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr8.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr8.tarseq.forward.1748 
vcftools --gzvcf chr8.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr8.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr8.tarseq.reverse.1748 

vcftools --gzvcf chr9.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr9.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr9.tarseq.forward.1748 
vcftools --gzvcf chr9.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr9.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr9.tarseq.reverse.1748 

vcftools --gzvcf chr10.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr10.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr10.tarseq.forward.1748 
vcftools --gzvcf chr10.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr10.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr10.tarseq.reverse.1748 

vcftools --gzvcf chr11.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr11.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr11.tarseq.forward.1748 
vcftools --gzvcf chr11.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr11.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr11.tarseq.reverse.1748 

vcftools --gzvcf chr12.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr12.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr12.tarseq.forward.1748 
vcftools --gzvcf chr12.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr12.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr12.tarseq.reverse.1748 

vcftools --gzvcf chr13.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr13.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr13.tarseq.forward.1748 
vcftools --gzvcf chr13.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr13.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr13.tarseq.reverse.1748 

vcftools --gzvcf chr14.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr14.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr14.tarseq.forward.1748 
vcftools --gzvcf chr14.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr14.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr14.tarseq.reverse.1748 

vcftools --gzvcf chr15.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr15.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr15.tarseq.forward.1748 
vcftools --gzvcf chr15.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr15.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr15.tarseq.reverse.1748 

vcftools --gzvcf chr16.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr16.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr16.tarseq.forward.1748 
vcftools --gzvcf chr16.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr16.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr16.tarseq.reverse.1748 

vcftools --gzvcf chr17.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr17.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr17.tarseq.forward.1748 
vcftools --gzvcf chr17.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr17.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr17.tarseq.reverse.1748 

vcftools --gzvcf chr18.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr18.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr18.tarseq.forward.1748 
vcftools --gzvcf chr18.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr18.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr18.tarseq.reverse.1748 

vcftools --gzvcf chr19.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr19.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr19.tarseq.forward.1748 
vcftools --gzvcf chr19.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr19.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr19.tarseq.reverse.1748 

vcftools --gzvcf chr20.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr20.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr20.tarseq.forward.1748 
vcftools --gzvcf chr20.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr20.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr20.tarseq.reverse.1748 

vcftools --gzvcf chr21.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr21.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr21.tarseq.forward.1748 
vcftools --gzvcf chr21.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr21.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr21.tarseq.reverse.1748 

vcftools --gzvcf chr22.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr22.recode.vcf.forward.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr22.tarseq.forward.1748 
vcftools --gzvcf chr22.dose.vcf.gz --snps /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr22.recode.vcf.reverse.snplist --keep /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples --recode --recode-INFO-all --out chr22.tarseq.reverse.1748 

bsub < pick_snp_overlap_with_tarseq.sh


# switch allele order while printing all the header staring with #
# the output should be tab separated (OFS: output field separator)
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr1_topmed_1000g.meta.tarseq.reverse.1748.recode.vcf > chr1_topmed_1000g.meta.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr2_topmed_1000g.meta.tarseq.reverse.1748.recode.vcf > chr2_topmed_1000g.meta.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr3_topmed_1000g.meta.tarseq.reverse.1748.recode.vcf > chr3_topmed_1000g.meta.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr4_topmed_1000g.meta.tarseq.reverse.1748.recode.vcf > chr4_topmed_1000g.meta.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr5_topmed_1000g.meta.tarseq.reverse.1748.recode.vcf > chr5_topmed_1000g.meta.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr6_topmed_1000g.meta.tarseq.reverse.1748.recode.vcf > chr6_topmed_1000g.meta.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr7_topmed_1000g.meta.tarseq.reverse.1748.recode.vcf > chr7_topmed_1000g.meta.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr8_topmed_1000g.meta.tarseq.reverse.1748.recode.vcf > chr8_topmed_1000g.meta.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr9_topmed_1000g.meta.tarseq.reverse.1748.recode.vcf > chr9_topmed_1000g.meta.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr10_topmed_1000g.meta.tarseq.reverse.1748.recode.vcf > chr10_topmed_1000g.meta.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr11_topmed_1000g.meta.tarseq.reverse.1748.recode.vcf > chr11_topmed_1000g.meta.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr12_topmed_1000g.meta.tarseq.reverse.1748.recode.vcf > chr12_topmed_1000g.meta.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr13_topmed_1000g.meta.tarseq.reverse.1748.recode.vcf > chr13_topmed_1000g.meta.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr14_topmed_1000g.meta.tarseq.reverse.1748.recode.vcf > chr14_topmed_1000g.meta.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr15_topmed_1000g.meta.tarseq.reverse.1748.recode.vcf > chr15_topmed_1000g.meta.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr16_topmed_1000g.meta.tarseq.reverse.1748.recode.vcf > chr16_topmed_1000g.meta.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr17_topmed_1000g.meta.tarseq.reverse.1748.recode.vcf > chr17_topmed_1000g.meta.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr18_topmed_1000g.meta.tarseq.reverse.1748.recode.vcf > chr18_topmed_1000g.meta.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr19_topmed_1000g.meta.tarseq.reverse.1748.recode.vcf > chr19_topmed_1000g.meta.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr20_topmed_1000g.meta.tarseq.reverse.1748.recode.vcf > chr20_topmed_1000g.meta.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr21_topmed_1000g.meta.tarseq.reverse.1748.recode.vcf > chr21_topmed_1000g.meta.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr22_topmed_1000g.meta.tarseq.reverse.1748.recode.vcf > chr22_topmed_1000g.meta.tarseq.reverse.1748.flip.recode.vcf

awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr1.tarseq.reverse.1748.recode.vcf > chr1.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr2.tarseq.reverse.1748.recode.vcf > chr2.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr3.tarseq.reverse.1748.recode.vcf > chr3.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr4.tarseq.reverse.1748.recode.vcf > chr4.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr5.tarseq.reverse.1748.recode.vcf > chr5.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr6.tarseq.reverse.1748.recode.vcf > chr6.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr7.tarseq.reverse.1748.recode.vcf > chr7.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr8.tarseq.reverse.1748.recode.vcf > chr8.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr9.tarseq.reverse.1748.recode.vcf > chr9.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr10.tarseq.reverse.1748.recode.vcf > chr10.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr11.tarseq.reverse.1748.recode.vcf > chr11.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr12.tarseq.reverse.1748.recode.vcf > chr12.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr13.tarseq.reverse.1748.recode.vcf > chr13.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr14.tarseq.reverse.1748.recode.vcf > chr14.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr15.tarseq.reverse.1748.recode.vcf > chr15.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr16.tarseq.reverse.1748.recode.vcf > chr16.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr17.tarseq.reverse.1748.recode.vcf > chr17.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr18.tarseq.reverse.1748.recode.vcf > chr18.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr19.tarseq.reverse.1748.recode.vcf > chr19.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr20.tarseq.reverse.1748.recode.vcf > chr20.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr21.tarseq.reverse.1748.recode.vcf > chr21.tarseq.reverse.1748.flip.recode.vcf
awk '/^#/ {print $0; next} BEGIN { OFS="\t" }{ t = $4; $4 = $5; $5 = t; print; } ' chr22.tarseq.reverse.1748.recode.vcf > chr22.tarseq.reverse.1748.flip.recode.vcf

# no overlap reverse SNP in GenomeAsia

##########################################################################################################################
#############        merge chr1-22 for the extracted SNPs (these SNPs are the ones present in Tarseq)     ################
##########################################################################################################################

############################
####     GenomeAsia     ####
############################

ls chr*.tarseq.forward.1748.recode.vcf|wc -l 
ls chr*.tarseq.reverse.1748.recode.vcf|wc -l
wc -l chr*.tarseq.forward.1748.recode.vcf
wc -l chr*.tarseq.reverse.1748.recode.vcf

ml bcftools
bcftools concat -o all_snp_overlap_tarseq.vcf  chr*.tarseq.forward.1748.recode.vcf
# merged file has 5,394 SNPs for GenomeAsia (44 header lines)
# sum of chr1.tarseq.forward.1748.recode.vcf to chr22.tarseq.forward.1748.recode.vcf is also 5394
# in Genome Asia, there is 0 SNPs for each chromosome in the reverse files (chr*.tarseq.reverse.1748.recode.vcf)

grep -v '^#' all_snp_overlap_tarseq.vcf|wc -l


# merge manually (this does not work because it seems that the vcf header has to follow certain order to be recognized as vcf, but by sorting, it changes the vcf header order)
# grep -h '^##' c*.recode.vcf | sort | uniq > all_snp_overlap_tarseq.vcf
#  -i, --ignore-case         ignore case distinctions
# -h: hide the filename
# grep -m1 '^#CHR' chr9.tarseq.forward.1748.recode.vcf >> all_snp_overlap_tarseq.vcf ## Get the chr header line 
#  --max-count=NUM       stop after NUM matches
# grep -v -h '^#'  c*.recode.vcf   >> all_snp_overlap_tarseq.vcf
# number of overlap SNPs = 11555-41 = 11,514 SNPs

# reorder the vcf to the Participant ID order in keep_1748_samples
ml bcftools
bcftools view -S /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples all_snp_overlap_tarseq.vcf -o all_snp_overlap_tarseq_ordered.vcf

# reorder the vcf so the SNP order starts at chromosome 1
ml vcftools 
vcf-sort  all_snp_overlap_tarseq_ordered.vcf >  all_snp_overlap_tarseq_ordered.sort.vcf
# SNP order: 1:6166382:C:T, 1:6166603:C:T
grep -v '^#'  all_snp_overlap_tarseq_ordered.sort.vcf|wc -l # 5,394 SNPs for GenomeAsia

# aggRSquare calculation
ml git python
ml cmake
/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/aggRSquare/release-build/aggRSquare -v /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.v4.2.updateID.sort.recode.vcf \
-i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_GA/all_snp_overlap_tarseq_ordered.sort.vcf \
-o tarseq_GA_grch37 \
--bins /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_compare/aggRSquare/default_MAF_bins \
--detail

## what about SNPs with Rsq > 0.8
bcftools filter -e 'INFO/R2<0.8' all_snp_overlap_tarseq_ordered.sort.vcf  > all_snp_overlap_tarseq_ordered.highR2.sort.vcf
# -e: exclude
grep -v '^#'  all_snp_overlap_tarseq_ordered.sort.vcf|wc -l  # 5394 SNPs
grep -v '^#'  all_snp_overlap_tarseq_ordered.highR2.sort.vcf|wc -l  # 1445 SNPs

# MAF < 0.01
bcftools filter -i 'INFO/MAF<0.01' all_snp_overlap_tarseq_ordered.sort.vcf |grep -v '^#' |wc -l  # 3654 SNPs with MAF < 0.01
# MAF >= 0.01 
bcftools filter -i 'INFO/MAF>=0.01' all_snp_overlap_tarseq_ordered.sort.vcf |grep -v '^#' |wc -l # 1740 SNPs with MAF >= 0.01

# MAF < 0.01
bcftools filter -i 'INFO/MAF<0.01' all_snp_overlap_tarseq_ordered.highR2.sort.vcf|grep -v '^#' |wc -l  # 277 SNPs with MAF < 0.01
# MAF >= 0.01 
bcftools filter -i 'INFO/MAF>=0.01' all_snp_overlap_tarseq_ordered.highR2.sort.vcf|grep -v '^#' |wc -l # 1168 SNPs with MAF >= 0.01


/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/aggRSquare/release-build/aggRSquare -v /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.v4.2.updateID.sort.recode.vcf \
-i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_GA/all_snp_overlap_tarseq_ordered.highR2.sort.vcf \
-o tarseq_GA_grch37_highR2 \
--bins /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_compare/aggRSquare/default_MAF_bins \
--detail


###### True R2 vs. Estimated R2   ################

awk '(NR>1){print $1}' tarseq_GA_grch37.RSquare  > snp_to_compare_forward # 2126
awk '(NR>1){print $1}' tarseq_GA_grch37.RSquare |sed 's/\:/ /g'|awk '{print $1":"$2":"$4":"$3}'  > snp_to_compare_reverse # 2126
wc -l snp_to_compare_forward
wc -l snp_to_compare_reverse
cat snp_to_compare_forward snp_to_compare_reverse > snp_to_compare # 4252
wc -l snp_to_compare

ml bcftools
bcftools annotate -x ^INFO/R2 all_snp_overlap_tarseq_ordered.sort.vcf > all_snp_overlap_tarseq_ordered.R2.sort.vcf
awk '{print $3"\t"$8}' all_snp_overlap_tarseq_ordered.R2.sort.vcf | grep -f snp_to_compare|sed 's/R2=//' > snp_estimate_r2
wc -l snp_estimate_r2 # 2126 SNPs

# in R
est_r2 <- read.table("snp_estimate_r2",head=F)
true_r2 <- read.table("tarseq_GA_grch37.RSquare",head=F)
reverse_snp <- read.table("snp_to_compare_reverse",head=F)

dim(est_r2) # 2126
dim(true_r2)  # 2126
dim(reverse_snp)  # 2126

colnames(reverse_snp) <- "reverse_SNP"
true_r2 <- cbind(true_r2,reverse_snp)
dim(true_r2)

dup_id <- true_r2[duplicated(true_r2$V1),"V1"] # no duplicated SNP
length(dup_id) # 0 dup SNPs
true_r2 <- true_r2[which(!true_r2$V1 %in% dup_id),] 

dim(est_r2) # 2126
dim(true_r2) # 2126

colnames(est_r2) <- c("SNP","est_R2")
true_r2_forward <- true_r2[,c("V1","V4","V5")]
true_r2_reverse <- true_r2[,c("reverse_SNP","V4","V5")]

head(est_r2)
head(true_r2_forward)
head(true_r2_reverse)

colnames(true_r2_forward) <-  c("SNP","true_R2","Validation_AF")
colnames(true_r2_reverse) <-  c("SNP","true_R2","Validation_AF")

forward_r2 <- merge(true_r2_forward,est_r2,by="SNP") # 8462 
reverse_r2 <- merge(true_r2_reverse,est_r2,by="SNP") # 9 
# sum to 8471 rows
head(forward_r2)
head(reverse_r2)
r2 <- rbind(forward_r2,reverse_r2) # 8471 SNPs

# all SNPs
r2$est_R2 <-as.numeric(r2$est_R2)
r2$true_R2 <-as.numeric(r2$true_R2)

cor(r2$est_R2, r2$true_R2, use = "pairwise.complete.obs")

plot(r2$est_R2, r2$true_R2 , xlab = "Estimated R Square", ylab = "True R Square", pch=20, frame = FALSE)
# Add regression line
abline(lm(true_R2 ~ est_R2 , data = r2), col = "red", lwd=2)

# rare SNPs
rare_r2 <- r2[which(r2$Validation_AF<0.01),] # 6580 SNPs
cor(rare_r2$est_R2, rare_r2$true_R2, use = "pairwise.complete.obs")

# plot 
plot(rare_r2$est_R2, rare_r2$true_R2 , xlab = "Estimated R Square", ylab = "True R Square", pch=20, frame = FALSE)
abline(lm(true_R2 ~ est_R2 , data = rare_r2), col = "red", lwd=2)


########################
####     TopMed     ####
########################

## SNPS in reverse file won't match with Tarseq SNPs and therefore aggRSquare cannot be run

ls chr*.tarseq.forward.1748.recode.vcf|wc -l 
ls chr*.tarseq.reverse.1748.recode.vcf|wc -l
ls chr*.tarseq.reverse.1748.flip.recode.vcf|wc -l
wc -l chr*.tarseq.forward.1748.recode.vcf # 48405 SNPs in total
wc -l chr*.tarseq.reverse.1748.recode.vcf # 50 SNPs in total 
# for example, chr9		138052235		chr9:138052235:CGT:C		CGT		C
wc -l chr*.tarseq.reverse.1748.flip.recode.vcf # 50 SNPs in total 
# for example, chr9		138052235		chr9:138052235:CGT:C		C		CGT

ml bcftools
ml vcftools 
bcftools concat -Oz -o all_forward_snp_overlap_tarseq.vcf.gz  chr*.tarseq.forward.1748.recode.vcf  # 48405 SNPs
zcat all_forward_snp_overlap_tarseq.vcf.gz|grep -v '^#' |wc -l
# use the flip reverse file, instead of the reverse file below
bcftools concat -Oz -o all_reverse_snp_overlap_tarseq.vcf.gz  chr*.tarseq.reverse.1748.flip.recode.vcf  # 50 SNPs
zcat all_reverse_snp_overlap_tarseq.vcf.gz|grep -v '^#' |wc -l
# this contains the flipped alleles

# print out the 50 reverse SNP names
zcat all_reverse_snp_overlap_tarseq.vcf.gz|grep -v '^#' |awk '{print $3}' > all_reverse_snp_overlap_tarseq.snplist # 50 SNPs

tabix -s1 -b2 -e2 all_forward_snp_overlap_tarseq.vcf.gz
tabix -s1 -b2 -e2 all_reverse_snp_overlap_tarseq.vcf.gz

bcftools concat --allow-overlaps all_forward_snp_overlap_tarseq.vcf.gz all_reverse_snp_overlap_tarseq.vcf.gz -o all_snp_overlap_tarseq.vcf # 48455 SNPs
grep -v '^#' all_snp_overlap_tarseq.vcf|wc -l 
# this contains the flip allele orders for the reverse SNPs

# reorder the vcf to the Participant ID order in keep_1748_samples
bcftools view -S /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples all_snp_overlap_tarseq.vcf -o all_snp_overlap_tarseq_ordered.vcf

# reorder the vcf so the SNP order starts at chromosome 1
vcf-sort  all_snp_overlap_tarseq_ordered.vcf >  all_snp_overlap_tarseq_ordered.sort.vcf

# SNP order: 1:6166382:C:T, 1:6166603:C:T
grep -v '^#'  all_snp_overlap_tarseq_ordered.sort.vcf|wc -l # 48,455 SNPs for TOPMed

# aggRSquare calculation
ml git python
ml cmake
/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/aggRSquare/release-build/aggRSquare -v /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.sort.recode.vcf  \
-i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/all_snp_overlap_tarseq_ordered.sort.vcf  \
-o tarseq_topmed_grch38 \
--bins /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_compare/aggRSquare/default_MAF_bins \
--detail

## what about SNPs with Rsq > 0.8
bcftools filter -e 'INFO/R2<0.8' all_snp_overlap_tarseq_ordered.sort.vcf  > all_snp_overlap_tarseq_ordered.highR2.sort.vcf
# -e: exclude
grep -v '^#'  all_snp_overlap_tarseq_ordered.sort.vcf|wc -l  # 48455 SNPs
grep -v '^#'  all_snp_overlap_tarseq_ordered.highR2.sort.vcf|wc -l  # 2735 SNPs

# MAF < 0.01
bcftools filter -i 'INFO/MAF<0.01' all_snp_overlap_tarseq_ordered.sort.vcf |grep -v '^#' |wc -l  # 46566 SNPs with MAF < 0.01
# MAF >= 0.01 
bcftools filter -i 'INFO/MAF>=0.01' all_snp_overlap_tarseq_ordered.sort.vcf |grep -v '^#' |wc -l # 1889 SNPs with MAF >= 0.01

# MAF < 0.01
bcftools filter -i 'INFO/MAF<0.01' all_snp_overlap_tarseq_ordered.highR2.sort.vcf|grep -v '^#' |wc -l  # 1218 SNPs with MAF < 0.01
# MAF >= 0.01 
bcftools filter -i 'INFO/MAF>=0.01' all_snp_overlap_tarseq_ordered.highR2.sort.vcf|grep -v '^#' |wc -l # 1517 SNPs with MAF >= 0.01


/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/aggRSquare/release-build/aggRSquare -v /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.sort.recode.vcf  \
-i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/all_snp_overlap_tarseq_ordered.highR2.sort.vcf  \
-o tarseq_topmed_grch38_highR2 \
--bins /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_compare/aggRSquare/default_MAF_bins \
--detail


# test how the aggRSquare calculate allele frequency in imputation file and validation file
# in validation file: only genotype exists, and therefore the AF is calculated as number of alleles/(2 * sample size) (For example, chr 17, bp 35120032 has 1 allele in the tarseq GT, but 0 in the imputation GT)
# in the imputation file, the AF is calculated as sum(HDS)/(2 * sample size) (for example, the chr17:35120032 AF based on this calculation = 0.000127, very close to 0.000128 in tarseq_topmed_grch38_reverse.RSquare)

cp all_reverse_snp_overlap_tarseq_ordered.sort.vcf  test.vcf

# only keep the genotype format (GT)
bcftools annotate -x ^FORMAT/GT test.vcf > test2.vcf
grep 35120032 test2.vcf

# what about HDS in GT:DS:HDS:GP
bcftools annotate -x ^FORMAT/HDS test.vcf > test3.vcf


###### True R2 vs. Estimated R2   ################

awk '(NR>1){print $1}' tarseq_topmed_grch38.RSquare  > snp_to_compare_forward # 8267
awk '(NR>1){print $1}' tarseq_topmed_grch38.RSquare |sed 's/\:/ /g'|awk '{print $1":"$2":"$4":"$3}'  > snp_to_compare_reverse # 8267
wc -l snp_to_compare_forward
wc -l snp_to_compare_reverse
cat snp_to_compare_forward snp_to_compare_reverse > snp_to_compare # 16534
wc -l snp_to_compare

ml bcftools
bcftools annotate -x ^INFO/R2 all_snp_overlap_tarseq_ordered.sort.vcf > all_snp_overlap_tarseq_ordered.R2.sort.vcf
awk '{print $3"\t"$8}' all_snp_overlap_tarseq_ordered.R2.sort.vcf | grep -f snp_to_compare|sed 's/R2=//' > snp_estimate_r2
wc -l snp_estimate_r2 # 8267 SNPs

# in R
est_r2 <- read.table("snp_estimate_r2",head=F)
true_r2 <- read.table("tarseq_topmed_grch38.RSquare",head=F)
reverse_snp <- read.table("snp_to_compare_reverse",head=F)

dim(est_r2) # 8267
dim(true_r2)  # 8267
dim(reverse_snp)  # 8267

colnames(reverse_snp) <- "reverse_SNP"
true_r2 <- cbind(true_r2,reverse_snp)
dim(true_r2)

dup_id <- true_r2[duplicated(true_r2$V1),"V1"] 
length(dup_id) # 5 dup SNPs
dim(true_r2[which(true_r2$V1 %in% dup_id),]) # 10 rows
true_r2 <- true_r2[which(!true_r2$V1 %in% dup_id),] # 3467 rows 

dim(est_r2) # 8267
dim(true_r2) # 8257

colnames(est_r2) <- c("SNP","est_R2")
true_r2_forward <- true_r2[,c("V1","V4","V5")]
true_r2_reverse <- true_r2[,c("reverse_SNP","V4","V5")]

head(est_r2)
head(true_r2_forward)
head(true_r2_reverse)

colnames(true_r2_forward) <-  c("SNP","true_R2","Validation_AF")
colnames(true_r2_reverse) <-  c("SNP","true_R2","Validation_AF")

forward_r2 <- merge(true_r2_forward,est_r2,by="SNP") # 8462 
reverse_r2 <- merge(true_r2_reverse,est_r2,by="SNP") # 9 
# sum to 8471 rows
head(forward_r2)
head(reverse_r2)
r2 <- rbind(forward_r2,reverse_r2) # 8471 SNPs

# all SNPs
r2$est_R2 <-as.numeric(r2$est_R2)
r2$true_R2 <-as.numeric(r2$true_R2)

cor(r2$est_R2, r2$true_R2, use = "pairwise.complete.obs")

plot(r2$est_R2, r2$true_R2 , xlab = "Estimated R Square", ylab = "True R Square", pch=20, frame = FALSE)
# Add regression line
abline(lm(true_R2 ~ est_R2 , data = r2), col = "red", lwd=2)

# rare SNPs
rare_r2 <- r2[which(r2$Validation_AF<0.01),] # 6580 SNPs
cor(rare_r2$est_R2, rare_r2$true_R2, use = "pairwise.complete.obs")

# plot 
plot(rare_r2$est_R2, rare_r2$true_R2 , xlab = "Estimated R Square", ylab = "True R Square", pch=20, frame = FALSE)
abline(lm(true_R2 ~ est_R2 , data = rare_r2), col = "red", lwd=2)



#######################################
####     1000G GRCh38 Extended     ####
#######################################

## SNPS in reverse file won't match with Tarseq SNPs and therefore aggRSquare cannot be run

ls chr*.tarseq.forward.1748.recode.vcf|wc -l 
ls chr*.tarseq.reverse.1748.flip.recode.vcf|wc -l
wc -l chr*.tarseq.forward.1748.recode.vcf # 11550 SNPs in total
wc -l chr*.tarseq.reverse.1748.flip.recode.vcf # 12 SNPs in total 
wc -l chr*.tarseq.reverse.1748.recode.vcf # 12 SNPs in total 

ml bcftools
ml vcftools 
bcftools concat -Oz -o all_forward_snp_overlap_tarseq.vcf.gz  chr*.tarseq.forward.1748.recode.vcf  # 11550 SNPs
zcat all_forward_snp_overlap_tarseq.vcf.gz|grep -v '^#' |wc -l
# use the flip reverse file, instead of the reverse file below
bcftools concat -Oz -o all_reverse_snp_overlap_tarseq.vcf.gz  chr*.tarseq.reverse.1748.flip.recode.vcf  # 12 SNPs
zcat all_reverse_snp_overlap_tarseq.vcf.gz|grep -v '^#' |wc -l

# print out the 50 reverse SNP names
zcat all_reverse_snp_overlap_tarseq.vcf.gz|grep -v '^#' |awk '{print $3}' > all_reverse_snp_overlap_tarseq.snplist # 12 SNPs


tabix -s1 -b2 -e2 all_forward_snp_overlap_tarseq.vcf.gz
tabix -s1 -b2 -e2 all_reverse_snp_overlap_tarseq.vcf.gz

bcftools concat --allow-overlaps all_forward_snp_overlap_tarseq.vcf.gz all_reverse_snp_overlap_tarseq.vcf.gz -o all_snp_overlap_tarseq.vcf # 11562 SNPs
grep -v '^#' all_snp_overlap_tarseq.vcf|wc -l 

# reorder the vcf to the Participant ID order in keep_1748_samples
bcftools view -S /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples all_snp_overlap_tarseq.vcf -o all_snp_overlap_tarseq_ordered.vcf

# reorder the vcf so the SNP order starts at chromosome 1
vcf-sort  all_snp_overlap_tarseq_ordered.vcf >  all_snp_overlap_tarseq_ordered.sort.vcf
# SNP order: 1:6166382:C:T, 1:6166603:C:T
grep -v '^#'  all_snp_overlap_tarseq_ordered.sort.vcf|wc -l # 48,455 SNPs for TOPMed

# aggRSquare calculation
ml git python
ml cmake
/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/aggRSquare/release-build/aggRSquare -v /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.sort.recode.vcf \
-i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/all_snp_overlap_tarseq_ordered.sort.vcf \
-o tarseq_ex1000G_grch38 \
--bins /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_compare/aggRSquare/default_MAF_bins \
--detail

## what about SNPs with Rsq > 0.8
bcftools filter -e 'INFO/R2<0.8' all_snp_overlap_tarseq_ordered.sort.vcf  > all_snp_overlap_tarseq_ordered.highR2.sort.vcf
# -e: exclude
grep -v '^#'  all_snp_overlap_tarseq_ordered.sort.vcf|wc -l  # 11562 SNPs
grep -v '^#'  all_snp_overlap_tarseq_ordered.highR2.sort.vcf|wc -l  # 2444 SNPs

# MAF < 0.01
bcftools filter -i 'INFO/MAF<0.01' all_snp_overlap_tarseq_ordered.sort.vcf |grep -v '^#' |wc -l  # 9680 SNPs with MAF < 0.01
# MAF >= 0.01 
bcftools filter -i 'INFO/MAF>=0.01' all_snp_overlap_tarseq_ordered.sort.vcf |grep -v '^#' |wc -l # 1882 SNPs with MAF >= 0.01

# MAF < 0.01
bcftools filter -i 'INFO/MAF<0.01' all_snp_overlap_tarseq_ordered.highR2.sort.vcf|grep -v '^#' |wc -l  # 842 SNPs with MAF < 0.01
# MAF >= 0.01 
bcftools filter -i 'INFO/MAF>=0.01' all_snp_overlap_tarseq_ordered.highR2.sort.vcf|grep -v '^#' |wc -l # 1602 SNPs with MAF >= 0.01

/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/aggRSquare/release-build/aggRSquare -v /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.sort.recode.vcf \
-i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/all_snp_overlap_tarseq_ordered.highR2.sort.vcf \
-o tarseq_ex1000G_grch38_highR2 \
--bins /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_compare/aggRSquare/default_MAF_bins \
--detail

###### True R2 vs. Estimated R2   ################

awk '(NR>1){print $1}' tarseq_ex1000G_grch38.RSquare  > snp_to_compare_forward # 3469
awk '(NR>1){print $1}' tarseq_ex1000G_grch38.RSquare |sed 's/\:/ /g'|awk '{print $1":"$2":"$4":"$3}'  > snp_to_compare_reverse # 3469
wc -l snp_to_compare_forward
wc -l snp_to_compare_reverse
cat snp_to_compare_forward snp_to_compare_reverse > snp_to_compare # 6938
wc -l snp_to_compare

ml bcftools
bcftools annotate -x ^INFO/R2 all_snp_overlap_tarseq_ordered.sort.vcf > all_snp_overlap_tarseq_ordered.R2.sort.vcf
awk '{print $3"\t"$8}' all_snp_overlap_tarseq_ordered.R2.sort.vcf | grep -f snp_to_compare|sed 's/R2=//' > snp_estimate_r2
wc -l snp_estimate_r2 # 3469 SNPs

# in R
est_r2 <- read.table("snp_estimate_r2",head=F)
true_r2 <- read.table("tarseq_ex1000G_grch38.RSquare",head=F)
reverse_snp <- read.table("snp_to_compare_reverse",head=F)

dim(est_r2) # 3469
dim(true_r2)  # 3469
dim(reverse_snp)  # 3469

colnames(reverse_snp) <- "reverse_SNP"
true_r2 <- cbind(true_r2,reverse_snp)
dim(true_r2)

dup_id <- true_r2[duplicated(true_r2$V1),"V1"] # no duplicated SNP
length(dup_id) # 1 
dim(true_r2[which(true_r2$V1 %in% dup_id),])
true_r2 <- true_r2[which(!true_r2$V1 %in% dup_id),] # 3467 rows 

dim(est_r2) # 3469
dim(true_r2) # 3467

colnames(est_r2) <- c("SNP","est_R2")
true_r2_forward <- true_r2[,c("V1","V4","V5")]
true_r2_reverse <- true_r2[,c("reverse_SNP","V4","V5")]

head(est_r2)
head(true_r2_forward)
head(true_r2_reverse)

colnames(true_r2_forward) <-  c("SNP","true_R2","Validation_AF")
colnames(true_r2_reverse) <-  c("SNP","true_R2","Validation_AF")

forward_r2 <- merge(true_r2_forward,est_r2,by="SNP") # 8462 
reverse_r2 <- merge(true_r2_reverse,est_r2,by="SNP") # 9 
# sum to 8471 rows
head(forward_r2)
head(reverse_r2)
r2 <- rbind(forward_r2,reverse_r2) # 8471 SNPs

# all SNPs
r2$est_R2 <-as.numeric(r2$est_R2)
r2$true_R2 <-as.numeric(r2$true_R2)

cor(r2$est_R2, r2$true_R2, use = "pairwise.complete.obs")

plot(r2$est_R2, r2$true_R2 , xlab = "Estimated R Square", ylab = "True R Square", pch=20, frame = FALSE)
# Add regression line
abline(lm(true_R2 ~ est_R2 , data = r2), col = "red", lwd=2)

# rare SNPs
rare_r2 <- r2[which(r2$Validation_AF<0.01),] # 6580 SNPs
cor(rare_r2$est_R2, rare_r2$true_R2, use = "pairwise.complete.obs")

# plot 
plot(rare_r2$est_R2, rare_r2$true_R2 , xlab = "Estimated R Square", ylab = "True R Square", pch=20, frame = FALSE)
abline(lm(true_R2 ~ est_R2 , data = rare_r2), col = "red", lwd=2)




##############################
####     1000G GRCh38     ####
##############################

## SNPS in reverse file won't match with Tarseq SNPs and therefore aggRSquare cannot be run

ls chr*.tarseq.forward.1748.recode.vcf|wc -l 
ls chr*.tarseq.reverse.1748.flip.recode.vcf|wc -l
wc -l chr*.tarseq.forward.1748.recode.vcf # 11160 SNPs in total
wc -l chr*.tarseq.reverse.1748.flip.recode.vcf # 12 SNPs in total 
wc -l chr*.tarseq.reverse.1748.recode.vcf # 12 SNPs in total 

ml bcftools
ml vcftools 
bcftools concat -Oz -o all_forward_snp_overlap_tarseq.vcf.gz  chr*.tarseq.forward.1748.recode.vcf  # 11160 SNPs
zcat all_forward_snp_overlap_tarseq.vcf.gz|grep -v '^#' |wc -l
# use the flip reverse file, instead of the reverse file below
bcftools concat -Oz -o all_reverse_snp_overlap_tarseq.vcf.gz  chr*.tarseq.reverse.1748.flip.recode.vcf  # 12 SNPs
zcat all_reverse_snp_overlap_tarseq.vcf.gz|grep -v '^#' |wc -l

# print out the 50 reverse SNP names
zcat all_reverse_snp_overlap_tarseq.vcf.gz|grep -v '^#' |awk '{print $3}' > all_reverse_snp_overlap_tarseq.snplist # 12 SNPs


tabix -s1 -b2 -e2 all_forward_snp_overlap_tarseq.vcf.gz
tabix -s1 -b2 -e2 all_reverse_snp_overlap_tarseq.vcf.gz

bcftools concat --allow-overlaps all_forward_snp_overlap_tarseq.vcf.gz all_reverse_snp_overlap_tarseq.vcf.gz -o all_snp_overlap_tarseq.vcf # 11172 SNPs
grep -v '^#' all_snp_overlap_tarseq.vcf|wc -l 

# reorder the vcf to the Participant ID order in keep_1748_samples
bcftools view -S /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples all_snp_overlap_tarseq.vcf -o all_snp_overlap_tarseq_ordered.vcf

# reorder the vcf so the SNP order starts at chromosome 1
vcf-sort  all_snp_overlap_tarseq_ordered.vcf >  all_snp_overlap_tarseq_ordered.sort.vcf
# SNP order: 1:6166382:C:T, 1:6166603:C:T
grep -v '^#'  all_snp_overlap_tarseq_ordered.sort.vcf|wc -l # 11,172 SNPs

# aggRSquare calculation
ml git python
ml cmake
/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/aggRSquare/release-build/aggRSquare -v /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.sort.recode.vcf \
-i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_38/all_snp_overlap_tarseq_ordered.sort.vcf \
-o tarseq_1000G_OG_grch38 \
--bins /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_compare/aggRSquare/default_MAF_bins \
--detail

## what about SNPs with Rsq > 0.8
bcftools filter -e 'INFO/R2<0.8' all_snp_overlap_tarseq_ordered.sort.vcf  > all_snp_overlap_tarseq_ordered.highR2.sort.vcf
# -e: exclude
grep -v '^#'  all_snp_overlap_tarseq_ordered.sort.vcf|wc -l  # 11172 SNPs
grep -v '^#'  all_snp_overlap_tarseq_ordered.highR2.sort.vcf|wc -l  # 2325 SNPs

# MAF < 0.01
bcftools filter -i 'INFO/MAF<0.01' all_snp_overlap_tarseq_ordered.sort.vcf |grep -v '^#' |wc -l  # 9325 SNPs with MAF < 0.01
# MAF >= 0.01 
bcftools filter -i 'INFO/MAF>=0.01' all_snp_overlap_tarseq_ordered.sort.vcf |grep -v '^#' |wc -l # 1847 SNPs with MAF >= 0.01

# MAF < 0.01
bcftools filter -i 'INFO/MAF<0.01' all_snp_overlap_tarseq_ordered.highR2.sort.vcf|grep -v '^#' |wc -l  # 776 SNPs with MAF < 0.01
# MAF >= 0.01 
bcftools filter -i 'INFO/MAF>=0.01' all_snp_overlap_tarseq_ordered.highR2.sort.vcf|grep -v '^#' |wc -l # 1549 SNPs with MAF >= 0.01

/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/aggRSquare/release-build/aggRSquare -v /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.sort.recode.vcf \
-i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_38/all_snp_overlap_tarseq_ordered.highR2.sort.vcf \
-o tarseq_1000G_OG_grch38_highR2 \
--bins /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_compare/aggRSquare/default_MAF_bins \
--detail

###### True R2 vs. Estimated R2   ################

awk '(NR>1){print $1}' tarseq_1000G_OG_grch38.RSquare  > snp_to_compare_forward # 3361
awk '(NR>1){print $1}' tarseq_1000G_OG_grch38.RSquare |sed 's/\:/ /g'|awk '{print $1":"$2":"$4":"$3}'  > snp_to_compare_reverse # 3361
wc -l snp_to_compare_forward
wc -l snp_to_compare_reverse
cat snp_to_compare_forward snp_to_compare_reverse > snp_to_compare # 6722
wc -l snp_to_compare

ml bcftools
bcftools annotate -x ^INFO/R2 all_snp_overlap_tarseq_ordered.sort.vcf > all_snp_overlap_tarseq_ordered.R2.sort.vcf
awk '{print $3"\t"$8}' all_snp_overlap_tarseq_ordered.R2.sort.vcf | grep -f snp_to_compare|sed 's/R2=//' > snp_estimate_r2
wc -l snp_estimate_r2 # 3361 SNPs

# in R
est_r2 <- read.table("snp_estimate_r2",head=F)
true_r2 <- read.table("tarseq_1000G_OG_grch38.RSquare",head=F)
reverse_snp <- read.table("snp_to_compare_reverse",head=F)

dim(est_r2) # 3361
dim(true_r2)  # 3361
dim(reverse_snp)  # 3361

colnames(reverse_snp) <- "reverse_SNP"
true_r2 <- cbind(true_r2,reverse_snp)
dim(true_r2)

dup_id <- true_r2[duplicated(true_r2$V1),"V1"] # no duplicated SNP
length(dup_id)
true_r2 <- true_r2[which(!true_r2$V1 %in% dup_id),] # 3467 rows 

dim(est_r2) # 3361
dim(true_r2) # 3361

colnames(est_r2) <- c("SNP","est_R2")
true_r2_forward <- true_r2[,c("V1","V4","V5")]
true_r2_reverse <- true_r2[,c("reverse_SNP","V4","V5")]

head(est_r2)
head(true_r2_forward)
head(true_r2_reverse)

colnames(true_r2_forward) <-  c("SNP","true_R2","Validation_AF")
colnames(true_r2_reverse) <-  c("SNP","true_R2","Validation_AF")

forward_r2 <- merge(true_r2_forward,est_r2,by="SNP") # 8462 
reverse_r2 <- merge(true_r2_reverse,est_r2,by="SNP") # 9 
# sum to 8471 rows
head(forward_r2)
head(reverse_r2)
r2 <- rbind(forward_r2,reverse_r2) # 8471 SNPs

# all SNPs
r2$est_R2 <-as.numeric(r2$est_R2)
r2$true_R2 <-as.numeric(r2$true_R2)

cor(r2$est_R2, r2$true_R2, use = "pairwise.complete.obs")

plot(r2$est_R2, r2$true_R2 , xlab = "Estimated R Square", ylab = "True R Square", pch=20, frame = FALSE)
# Add regression line
abline(lm(true_R2 ~ est_R2 , data = r2), col = "red", lwd=2)

# rare SNPs
rare_r2 <- r2[which(r2$Validation_AF<0.01),] # 6580 SNPs
cor(rare_r2$est_R2, rare_r2$true_R2, use = "pairwise.complete.obs")

# plot 
plot(rare_r2$est_R2, rare_r2$true_R2 , xlab = "Estimated R Square", ylab = "True R Square", pch=20, frame = FALSE)
abline(lm(true_R2 ~ est_R2 , data = rare_r2), col = "red", lwd=2)



##############################
####     1000G GRCh37     ####
##############################

## SNPS in reverse file won't match with Tarseq SNPs and therefore aggRSquare cannot be run

ls chr*.tarseq.forward.1748.recode.vcf|wc -l 
ls chr*.tarseq.reverse.1748.flip.recode.vcf|wc -l
wc -l chr*.tarseq.forward.1748.recode.vcf # 11501 SNPs in total
wc -l chr*.tarseq.reverse.1748.flip.recode.vcf # 13 SNPs in total 
wc -l chr*.tarseq.reverse.1748.flip.recode.vcf # 13 SNPs in total 

ml bcftools
ml vcftools 
bcftools concat -Oz -o all_forward_snp_overlap_tarseq.vcf.gz  chr*.tarseq.forward.1748.recode.vcf  # 11501 SNPs
zcat all_forward_snp_overlap_tarseq.vcf.gz|grep -v '^#' |wc -l
# use the flip reverse file, instead of the reverse file below
bcftools concat -Oz -o all_reverse_snp_overlap_tarseq.vcf.gz  chr*.tarseq.reverse.1748.flip.recode.vcf  # 13 SNPs
zcat all_reverse_snp_overlap_tarseq.vcf.gz|grep -v '^#' |wc -l

# print out the 50 reverse SNP names
zcat all_reverse_snp_overlap_tarseq.vcf.gz|grep -v '^#' |awk '{print $3}' > all_reverse_snp_overlap_tarseq.snplist # 13 SNPs


tabix -s1 -b2 -e2 all_forward_snp_overlap_tarseq.vcf.gz
tabix -s1 -b2 -e2 all_reverse_snp_overlap_tarseq.vcf.gz

bcftools concat --allow-overlaps all_forward_snp_overlap_tarseq.vcf.gz all_reverse_snp_overlap_tarseq.vcf.gz -o all_snp_overlap_tarseq.vcf # 11514 SNPs
grep -v '^#' all_snp_overlap_tarseq.vcf|wc -l 

# reorder the vcf to the Participant ID order in keep_1748_samples
bcftools view -S /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples all_snp_overlap_tarseq.vcf -o all_snp_overlap_tarseq_ordered.vcf

# reorder the vcf so the SNP order starts at chromosome 1
vcf-sort  all_snp_overlap_tarseq_ordered.vcf >  all_snp_overlap_tarseq_ordered.sort.vcf
# SNP order: 1:6166382:C:T, 1:6166603:C:T
grep -v '^#'  all_snp_overlap_tarseq_ordered.sort.vcf|wc -l # 48,455 SNPs for TOPMed

# aggRSquare calculation
ml git python
ml cmake
/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/aggRSquare/release-build/aggRSquare -v /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.v4.2.updateID.sort.recode.vcf  \
-i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_37/all_snp_overlap_tarseq_ordered.sort.vcf \
-o tarseq_1000G_OG_grch37 \
--bins /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_compare/aggRSquare/default_MAF_bins \
--detail

## what about SNPs with Rsq > 0.8
bcftools filter -e 'INFO/R2<0.8' all_snp_overlap_tarseq_ordered.sort.vcf  > all_snp_overlap_tarseq_ordered.highR2.sort.vcf
# -e: exclude
grep -v '^#'  all_snp_overlap_tarseq_ordered.sort.vcf|wc -l  # 11514 SNPs
grep -v '^#'  all_snp_overlap_tarseq_ordered.highR2.sort.vcf|wc -l  # 1843 SNPs

# MAF < 0.01
bcftools filter -i 'INFO/MAF<0.01' all_snp_overlap_tarseq_ordered.sort.vcf |grep -v '^#' |wc -l  # 9563 SNPs with MAF < 0.01
# MAF >= 0.01 
bcftools filter -i 'INFO/MAF>=0.01' all_snp_overlap_tarseq_ordered.sort.vcf |grep -v '^#' |wc -l # 1951 SNPs with MAF >= 0.01

# MAF < 0.01
bcftools filter -i 'INFO/MAF<0.01' all_snp_overlap_tarseq_ordered.highR2.sort.vcf|grep -v '^#' |wc -l  # 474 SNPs with MAF < 0.01
# MAF >= 0.01 
bcftools filter -i 'INFO/MAF>=0.01' all_snp_overlap_tarseq_ordered.highR2.sort.vcf|grep -v '^#' |wc -l # 1369 SNPs with MAF >= 0.01

/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/aggRSquare/release-build/aggRSquare -v /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.v4.2.updateID.sort.recode.vcf  \
-i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_37/all_snp_overlap_tarseq_ordered.highR2.sort.vcf \
-o tarseq_1000G_OG_grch37_highR2 \
--bins /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_compare/aggRSquare/default_MAF_bins \
--detail

###### True R2 vs. Estimated R2   ################

awk '(NR>1){print $1}' tarseq_1000G_OG_grch37.RSquare  > snp_to_compare_forward # 3467
awk '(NR>1){print $1}' tarseq_1000G_OG_grch37.RSquare |sed 's/\:/ /g'|awk '{print $1":"$2":"$4":"$3}'  > snp_to_compare_reverse # 3467
cat snp_to_compare_forward snp_to_compare_reverse > snp_to_compare
ml bcftools
bcftools annotate -x ^INFO/R2 all_snp_overlap_tarseq_ordered.sort.vcf > all_snp_overlap_tarseq_ordered.R2.sort.vcf
awk '{print $3"\t"$8}' all_snp_overlap_tarseq_ordered.R2.sort.vcf | grep -f snp_to_compare|sed 's/R2=//' > snp_estimate_r2
wc -l snp_estimate_r2 # 3468 (with 1 more SNPs included from vcf)

# in R
est_r2 <- read.table("snp_estimate_r2",head=F)
true_r2 <- read.table("tarseq_1000G_OG_grch37.RSquare",head=F)
reverse_snp <- read.table("snp_to_compare_reverse",head=F)

dim(est_r2) # 3468
dim(true_r2)  # 3467
dim(reverse_snp)  # 3467

colnames(reverse_snp) <- "reverse_SNP"
true_r2 <- cbind(true_r2,reverse_snp)

dup_id <- true_r2[duplicated(true_r2$V1),"V1"] # no duplicated SNP
dim(true_r2[which(true_r2$V1 %in% dup_id),])
true_r2 <- true_r2[which(!true_r2$V1 %in% dup_id),] # 3467 rows 

dim(est_r2) # 3468
dim(true_r2) # 3467

colnames(est_r2) <- c("SNP","est_R2")
true_r2_forward <- true_r2[,c("V1","V4","V5")]
true_r2_reverse <- true_r2[,c("reverse_SNP","V4","V5")]

head(est_r2)
head(true_r2_forward)
head(true_r2_reverse)

colnames(true_r2_forward) <-  c("SNP","true_R2","Validation_AF")
colnames(true_r2_reverse) <-  c("SNP","true_R2","Validation_AF")

forward_r2 <- merge(true_r2_forward,est_r2,by="SNP") # 8462 
reverse_r2 <- merge(true_r2_reverse,est_r2,by="SNP") # 9 
# sum to 8471 rows
head(forward_r2)
head(reverse_r2)
r2 <- rbind(forward_r2,reverse_r2) # 8471 SNPs

# all SNPs
r2$est_R2 <-as.numeric(r2$est_R2)
r2$true_R2 <-as.numeric(r2$true_R2)

cor(r2$est_R2, r2$true_R2, use = "pairwise.complete.obs")

plot(r2$est_R2, r2$true_R2 , xlab = "Estimated R Square", ylab = "True R Square", pch=20, frame = FALSE)
# Add regression line
abline(lm(true_R2 ~ est_R2 , data = r2), col = "red", lwd=2)

# rare SNPs
rare_r2 <- r2[which(r2$Validation_AF<0.01),] # 6580 SNPs
cor(rare_r2$est_R2, rare_r2$true_R2, use = "pairwise.complete.obs")

# plot 
plot(rare_r2$est_R2, rare_r2$true_R2 , xlab = "Estimated R Square", ylab = "True R Square", pch=20, frame = FALSE)
abline(lm(true_R2 ~ est_R2 , data = rare_r2), col = "red", lwd=2)


##############################
###     meta-imputation   ####
##############################

ls *_topmed_1000g.meta.tarseq.forward.1748.recode.vcf | wc -l
ls *_topmed_1000g.meta.tarseq.reverse.1748.recode.vcf | wc -l
ls *_topmed_1000g.meta.tarseq.reverse.1748.flip.recode.vcf | wc -l
wc -l *_topmed_1000g.meta.tarseq.forward.1748.recode.vcf # 48,854 SNPs
wc -l *_topmed_1000g.meta.tarseq.reverse.1748.recode.vcf  # 53 SNPs
wc -l *_topmed_1000g.meta.tarseq.reverse.1748.flip.recode.vcf # 53 SNPs

ml bcftools
ml vcftools 
bcftools concat -Oz -o all_forward_snp_overlap_tarseq.vcf.gz  *_topmed_1000g.meta.tarseq.forward.1748.recode.vcf   # 48854 SNPs
zcat all_forward_snp_overlap_tarseq.vcf.gz|grep -v '^#' |wc -l
# use the flip reverse file, instead of the reverse file below
bcftools concat -Oz -o all_reverse_snp_overlap_tarseq.vcf.gz  *_topmed_1000g.meta.tarseq.reverse.1748.flip.recode.vcf   # 53 SNPs
zcat all_reverse_snp_overlap_tarseq.vcf.gz|grep -v '^#' |wc -l

# print out the 50 reverse SNP names
zcat all_reverse_snp_overlap_tarseq.vcf.gz|grep -v '^#' |awk '{print $3}' > all_reverse_snp_overlap_tarseq.snplist # 53 SNPs


tabix -s1 -b2 -e2 all_forward_snp_overlap_tarseq.vcf.gz
tabix -s1 -b2 -e2 all_reverse_snp_overlap_tarseq.vcf.gz

bcftools concat --allow-overlaps all_forward_snp_overlap_tarseq.vcf.gz all_reverse_snp_overlap_tarseq.vcf.gz -o all_snp_overlap_tarseq.vcf # 48,907 SNPs
grep -v '^#' all_snp_overlap_tarseq.vcf|wc -l 

# reorder the vcf to the Participant ID order in keep_1748_samples
bcftools view -S /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples all_snp_overlap_tarseq.vcf -o all_snp_overlap_tarseq_ordered.vcf

# reorder the vcf so the SNP order starts at chromosome 1
vcf-sort  all_snp_overlap_tarseq_ordered.vcf >  all_snp_overlap_tarseq_ordered.sort.vcf
# SNP order: 1:6166382:C:T, 1:6166603:C:T
grep -v '^#'  all_snp_overlap_tarseq_ordered.sort.vcf|wc -l # 48,907 SNPs for TOPMed

# aggRSquare calculation
ml git python
ml cmake
/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/aggRSquare/release-build/aggRSquare -v /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.sort.recode.vcf \
-i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_topmed_1000G/all_snp_overlap_tarseq_ordered.sort.vcf \
-o tarseq_meta_imputation \
--bins /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_compare/aggRSquare/default_MAF_bins \
--detail

## what about SNPs with Rsq > 0.8
bcftools filter -e 'INFO/R2<0.8' all_snp_overlap_tarseq_ordered.sort.vcf  > all_snp_overlap_tarseq_ordered.highR2.sort.vcf
# -e: exclude
grep -v '^#'  all_snp_overlap_tarseq_ordered.sort.vcf|wc -l  # 48907 SNPs
grep -v '^#'  all_snp_overlap_tarseq_ordered.highR2.sort.vcf|wc -l  # 2968 SNPs

# MAF < 0.01
bcftools filter -i 'INFO/MAF<0.01' all_snp_overlap_tarseq_ordered.sort.vcf |grep -v '^#' |wc -l  # 46996 SNPs with MAF < 0.01
# MAF >= 0.01 
bcftools filter -i 'INFO/MAF>=0.01' all_snp_overlap_tarseq_ordered.sort.vcf |grep -v '^#' |wc -l # 1911 SNPs with MAF >= 0.01

# MAF < 0.01
bcftools filter -i 'INFO/MAF<0.01' all_snp_overlap_tarseq_ordered.highR2.sort.vcf|grep -v '^#' |wc -l  # 1357 SNPs with MAF < 0.01
# MAF >= 0.01 
bcftools filter -i 'INFO/MAF>=0.01' all_snp_overlap_tarseq_ordered.highR2.sort.vcf|grep -v '^#' |wc -l # 1611 SNPs with MAF >= 0.01


/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/aggRSquare/release-build/aggRSquare -v /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.sort.recode.vcf \
-i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_topmed_1000G/all_snp_overlap_tarseq_ordered.highR2.sort.vcf \
-o tarseq_meta_imputation_highR2 \
--bins /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_compare/aggRSquare/default_MAF_bins \
--detail


#################################################################
######          True R2 vs. Estimated R2         ################
#################################################################

awk '(NR>1){print $1}' tarseq_meta_imputation.RSquare  > snp_to_compare_forward # 8483
awk '(NR>1){print $1}' tarseq_meta_imputation.RSquare |sed 's/\:/ /g'|awk '{print $1":"$2":"$4":"$3}'  > snp_to_compare_reverse # 8483
cat snp_to_compare_forward snp_to_compare_reverse > snp_to_compare
bcftools annotate -x ^INFO/R2 all_snp_overlap_tarseq_ordered.sort.vcf > all_snp_overlap_tarseq_ordered.R2.sort.vcf
awk '{print $3"\t"$8}' all_snp_overlap_tarseq_ordered.R2.sort.vcf | grep -f snp_to_compare|sed 's/R2=//' > snp_estimate_r2
wc -l snp_estimate_r2 # 8483

# in R
est_r2 <- read.table("snp_estimate_r2",head=F)
true_r2 <- read.table("tarseq_meta_imputation.RSquare",head=F)
#SNP.ID	Allele.Frequency	No.Samples	Imputation.R2	Validation.AF	Imputation.AF
reverse_snp <- read.table("snp_to_compare_reverse",head=F)

colnames(reverse_snp) <- "reverse_SNP"
true_r2 <- cbind(true_r2,reverse_snp)

dup_id <- true_r2[duplicated(true_r2$V1),"V1"]
dim(true_r2[which(true_r2$V1 %in% dup_id),])
true_r2 <- true_r2[which(!true_r2$V1 %in% dup_id),] # 8471 rows = 8483-12


dim(est_r2) # 8483
dim(true_r2) # 8471

colnames(est_r2) <- c("SNP","est_R2")
true_r2_forward <- true_r2[,c("V1","V4","V5")]
true_r2_reverse <- true_r2[,c("reverse_SNP","V4","V5")]

head(est_r2)
head(true_r2_forward)
head(true_r2_reverse)

colnames(true_r2_forward) <-  c("SNP","true_R2","Validation_AF")
colnames(true_r2_reverse) <-  c("SNP","true_R2","Validation_AF")

forward_r2 <- merge(true_r2_forward,est_r2,by="SNP") # 8462 
reverse_r2 <- merge(true_r2_reverse,est_r2,by="SNP") # 9 
# sum to 8471 rows
head(forward_r2)
head(reverse_r2)
r2 <- rbind(forward_r2,reverse_r2) # 8471 SNPs

# all SNPs
cor(r2$est_R2, r2$true_R2 )

plot(r2$est_R2, r2$true_R2 , xlab = "Estimated R Square", ylab = "True R Square", pch=20, frame = FALSE)
# Add regression line
abline(lm(true_R2 ~ est_R2 , data = r2), col = "red", lwd=2)

# rare SNPs
rare_r2 <- r2[which(r2$Validation_AF<0.01),] # 6580 SNPs

cor(rare_r2$est_R2, rare_r2$true_R2 )

# plot 
plot(rare_r2$est_R2, rare_r2$true_R2 , xlab = "Estimated R Square", ylab = "True R Square", pch=20, frame = FALSE)
abline(lm(true_R2 ~ est_R2 , data = rare_r2), col = "red", lwd=2)



# meta-imputation: 33 flipped SNPs
# TOPMed: 32 flipped SNPs 
# Expanded 1000G: 6 flipped SNPs (9,132 SNPs matched)
# 1000G GRCh38: 8 flipped SNPs (8,837 SNPs matched)
# 1000G GRCh37 SAS: 13 flipped SNPs (11,501 SNPs matched)
# Genome Asia: 0 flipped SNPs (5,394 SNPs matched)

##########################################################
#######      R2 and EUR MAF/SAS MAF Difference      ######
##########################################################

/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/compare_across_all/MAF_SAS_EUR_1000G_unique_chr_pos
# 84,617,679 lines

####################################
########    add GRCh38 IDs    ######
####################################

# use this liftover website: https://genome.ucsc.edu/cgi-bin/hgLiftOver

# need to prepare a list of SNPs in this format: chr1 6166310 6166311

awk '{print $1}' MAF_SAS_EUR_1000G_unique_chr_pos | sed 's/\:/ /' |awk '(NR>1){print $1,$2-1,$2,$1":"$2}'| sed 's/^/chr/'  > liftover_1000g_snp_grch37_unique.bed # 84,617,678 lines
# different number of lifted SNPs when using $2-1, $2 vs. $2, $2+1
# fourth column: to keep the original grch37 location for easier comparison with grch38

# execute liftover
ml liftover
liftOver liftover_1000g_snp_grch37_unique.bed /sc/arion/projects/psychgen/software/liftover/hg19ToHg38.over.chain.gz lifted_1000g_snp_grch38_unique.bed unlifted.bed
# 36,518 unlifted SNPs
# 84,599,419 lifted SNPs

# checking using one SNP: rs1344706
# grch37: 2:185778428; grch38: 2:184913701

# merge the grch38 position with the MAF_SAS_EUR_1000G_v2 file
awk '{print $1":"$3,$4}' lifted_1000g_snp_grch38_unique.bed >  lifted_1000g_snp_grch38_chr_pos #  84,599,419  lines
# chr1:10177 1:10177
# chr1:10235 1:10235

# check duplicated SNPs
awk '{print $1}' lifted_1000g_snp_grch38_chr_pos|sort | uniq -d > duplicate_1000g_grch38_pos # 1682 duplicated positions
# remove duplicated SNPs
awk 'NR == FNR {a[$1]; next} !($1 in a)' duplicate_1000g_grch38_pos lifted_1000g_snp_grch38_chr_pos > lifted_1000g_snp_grch38_unique_chr_pos
# 84,596,048 lines (some SNPs may have >2 duplicated rows)

# merge the grch38 position with grch37 position, along with SAS MAF and EUR MAF
awk 'NR==FNR {h[$2]=$0; next} {print $0,h[$1]}'  lifted_1000g_snp_grch38_unique_chr_pos  MAF_SAS_EUR_1000G_unique_chr_pos >  tmp # 84617679 lines

awk '(NF==6){print $0}' tmp | awk '{print $1"\t"$5"\t"$2"\t"$3"\t"$4}' |sed '1i SNP_37\tSNP_38\tSAS_AF\tEUR_AF\tDIFF_AF ' > MAF_SAS_EUR_1000G_GRCh37_GRCh38_unique
# 84,596,049 lines

#### Meta-imputation  ###

sed 's/\:/\t/g' tarseq_meta_imputation.RSquare |awk '(NR>1){print $1":"$2"\t"$1":"$2":"$3":"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}'  > tarseq_meta_imputation.edit.RSquare # 8483 lines
# check duplicated SNPs
awk '{print $1}' tarseq_meta_imputation.edit.RSquare|sort | uniq -d > duplicate_snp # 32 duplicated snps
# remove duplicated SNPs
grep -vf duplicate_snp  tarseq_meta_imputation.edit.RSquare >  tarseq_meta_imputation.unique.RSquare
# 8418 lines

awk 'NR==FNR {h[$2]=$0; next} {print $0,h[$1]}' /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/compare_across_all/MAF_SAS_EUR_1000G_GRCh37_GRCh38_unique tarseq_meta_imputation.unique.RSquare |sed '1i CHR_BP\tSNP_ID\tAllele_Frequency\tNo_Samples\tImputation_R2\tValidation_AF\tImputation.AF\tSNP_37\tSNP_38\tSAS_AF\tEUR_AF\tDIFF_AF' |  awk '(NF==12){print $0}'  > tarseq_meta_imputation.unique.SAS_EUR_MAF.RSquare 
# 5914 lines (out of 8418 SNPs, ~70%)
# not all SNPs that are imputed and sequenced can be matched with 1000G SNPs 





#### TOPMed  ###

sed 's/\:/\t/g' tarseq_topmed_grch38.RSquare |awk '(NR>1){print $1":"$2"\t"$1":"$2":"$3":"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}'  > tarseq_topmed_grch38.edit.RSquare # 8267 lines
# check duplicated SNPs
awk '{print $1}' tarseq_topmed_grch38.edit.RSquare|sort | uniq -d > duplicate_snp # 30 duplicated snps
# remove duplicated SNPs
grep -vf duplicate_snp  tarseq_topmed_grch38.edit.RSquare >  tarseq_topmed_grch38.unique.RSquare
# 8206 lines

awk 'NR==FNR {h[$2]=$0; next} {print $0,h[$1]}' /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/compare_across_all/MAF_SAS_EUR_1000G_GRCh37_GRCh38_unique tarseq_topmed_grch38.unique.RSquare |sed '1i CHR_BP\tSNP_ID\tAllele_Frequency\tNo_Samples\tImputation_R2\tValidation_AF\tImputation.AF\tSNP_37\tSNP_38\tSAS_AF\tEUR_AF\tDIFF_AF' |  awk '(NF==12){print $0}'  > tarseq_topmed_grch38.unique.SAS_EUR_MAF.RSquare 
# 5702 lines 

#### Expanded 1000G  ###
sed 's/\:/\t/g' tarseq_ex1000G_grch38.RSquare|awk '(NR>1){print $1":"$2"\t"$1":"$2":"$3":"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}'  > tarseq_ex1000G_grch38.edit.RSquare # 3469 lines
# check duplicated SNPs
awk '{print $1}' tarseq_ex1000G_grch38.edit.RSquare|sort | uniq -d > duplicate_snp # 4 duplicated snps
# remove duplicated SNPs
grep -vf duplicate_snp  tarseq_ex1000G_grch38.edit.RSquare >  tarseq_ex1000G_grch38.unique.RSquare
# 3461 lines

awk 'NR==FNR {h[$2]=$0; next} {print $0,h[$1]}' /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/compare_across_all/MAF_SAS_EUR_1000G_GRCh37_GRCh38_unique tarseq_ex1000G_grch38.unique.RSquare |sed '1i CHR_BP\tSNP_ID\tAllele_Frequency\tNo_Samples\tImputation_R2\tValidation_AF\tImputation.AF\tSNP_37\tSNP_38\tSAS_AF\tEUR_AF\tDIFF_AF' |  awk '(NF==12){print $0}'  > tarseq_ex1000G_grch38.unique.SAS_EUR_MAF.RSquare 
# 3426 lines 

#### 1000G GRCh38 ###
sed 's/\:/\t/g' tarseq_1000G_OG_grch38.RSquare|awk '(NR>1){print $1":"$2"\t"$1":"$2":"$3":"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}'  > tarseq_1000G_OG_grch38.edit.RSquare # 3361 lines
# check duplicated SNPs
awk '{print $1}' tarseq_1000G_OG_grch38.edit.RSquare|sort | uniq -d > duplicate_snp # 0 duplicated snps
# remove duplicated SNPs
grep -vf duplicate_snp  tarseq_1000G_OG_grch38.edit.RSquare >  tarseq_1000G_OG_grch38.unique.RSquare
# 3361 lines

awk 'NR==FNR {h[$2]=$0; next} {print $0,h[$1]}' /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/compare_across_all/MAF_SAS_EUR_1000G_GRCh37_GRCh38_unique tarseq_1000G_OG_grch38.unique.RSquare |sed '1i CHR_BP\tSNP_ID\tAllele_Frequency\tNo_Samples\tImputation_R2\tValidation_AF\tImputation.AF\tSNP_37\tSNP_38\tSAS_AF\tEUR_AF\tDIFF_AF' |  awk '(NF==12){print $0}'  > tarseq_1000G_OG_grch38.unique.SAS_EUR_MAF.RSquare 
# 3334 lines 

#### 1000G GRCh37 ###
sed 's/\:/\t/g' tarseq_1000G_OG_grch37.RSquare|awk '(NR>1){print $1":"$2"\t"$1":"$2":"$3":"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}'  > tarseq_1000G_OG_grch37.edit.RSquare # 3467 lines
# check duplicated SNPs
awk '{print $1}' tarseq_1000G_OG_grch37.edit.RSquare|sort | uniq -d > duplicate_snp # 7 duplicated snps
# remove duplicated SNPs
grep -vf duplicate_snp  tarseq_1000G_OG_grch37.edit.RSquare >  tarseq_1000G_OG_grch37.unique.RSquare
# 3453 lines

awk 'NR==FNR {h[$1]=$0; next} {print $0,h[$1]}' /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/compare_across_all/MAF_SAS_EUR_1000G_GRCh37_GRCh38_unique tarseq_1000G_OG_grch37.unique.RSquare |sed '1i CHR_BP\tSNP_ID\tAllele_Frequency\tNo_Samples\tImputation_R2\tValidation_AF\tImputation.AF\tSNP_37\tSNP_38\tSAS_AF\tEUR_AF\tDIFF_AF' |  awk '(NF==12){print $0}'  > tarseq_1000G_OG_grch37.unique.SAS_EUR_MAF.RSquare 
# 3451 lines 

#### GenomeAsia ###
sed 's/\:/\t/g'  tarseq_GA_grch37.RSquare|awk '(NR>1){print $1":"$2"\t"$1":"$2":"$3":"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}'  > tarseq_GA_grch37.edit.RSquare # 2126 lines
# check duplicated SNPs
awk '{print $1}' tarseq_GA_grch37.edit.RSquare|sort | uniq -d > duplicate_snp # 0 duplicated snps
# remove duplicated SNPs
grep -vf duplicate_snp  tarseq_GA_grch37.edit.RSquare >  tarseq_GA_grch37.unique.RSquare
# 2126 lines

awk 'NR==FNR {h[$1]=$0; next} {print $0,h[$1]}' /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/compare_across_all/MAF_SAS_EUR_1000G_GRCh37_GRCh38_unique tarseq_GA_grch37.unique.RSquare |sed '1i CHR_BP\tSNP_ID\tAllele_Frequency\tNo_Samples\tImputation_R2\tValidation_AF\tImputation.AF\tSNP_37\tSNP_38\tSAS_AF\tEUR_AF\tDIFF_AF' |  awk '(NF==12){print $0}'  > tarseq_GA_grch37.unique.SAS_EUR_MAF.RSquare 
# 2067 lines 

#######      plot the true R2 vs deviation from EUR MAF      ##########
# in R
rsq_meta <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_topmed_1000G/tarseq_meta_imputation.unique.SAS_EUR_MAF.RSquare",head=T)
rsq_topmed  <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/tarseq_topmed_grch38.unique.SAS_EUR_MAF.RSquare",head=T)
rsq_ex1kg  <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/tarseq_ex1000G_grch38.unique.SAS_EUR_MAF.RSquare",head=T)
rsq_1kg38  <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_38/tarseq_1000G_OG_grch38.unique.SAS_EUR_MAF.RSquare",head=T)
rsq_1kg37  <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_37/tarseq_1000G_OG_grch37.unique.SAS_EUR_MAF.RSquare",head=T)
rsq_ga  <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_GA/tarseq_GA_grch37.unique.SAS_EUR_MAF.RSquare",head=T)

reverse_snp_meta <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_topmed_1000G/all_reverse_snp_overlap_tarseq.snplist",head=F)
reverse_snp_topmed  <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/all_reverse_snp_overlap_tarseq.snplist",head=F)
reverse_snp_ex1kg  <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/all_reverse_snp_overlap_tarseq.snplist",head=F)
reverse_snp_1kg38  <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_38/all_reverse_snp_overlap_tarseq.snplist",head=F)
reverse_snp_1kg37  <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_37/all_reverse_snp_overlap_tarseq.snplist",head=F)

dim(reverse_snp_meta) # 53 
dim(reverse_snp_topmed) # 50
dim(reverse_snp_ex1kg) # 12
dim(reverse_snp_1kg38) # 12
dim(reverse_snp_1kg37) # 13 

dim(rsq_meta) # 5913 SNPs 
dim(rsq_topmed) # 5701 SNPs
dim(rsq_ex1kg) # 3425 SNPs
dim(rsq_1kg38) # 3333 SNPs
dim(rsq_1kg37) # 3450 SNPs
dim(rsq_ga) # 2066 SNPs

library(stringr)
reverse_snp_meta[c('chr', 'pos', 'a1', 'a2')] <- str_split_fixed(reverse_snp_meta$V1, ':', 4)
reverse_snp_topmed[c('chr', 'pos', 'a1', 'a2')] <- str_split_fixed(reverse_snp_topmed$V1, ':', 4)
reverse_snp_ex1kg[c('chr', 'pos', 'a1', 'a2')] <- str_split_fixed(reverse_snp_ex1kg$V1, ':', 4)
reverse_snp_1kg38[c('chr', 'pos', 'a1', 'a2')] <- str_split_fixed(reverse_snp_1kg38$V1, ':', 4)
reverse_snp_1kg37[c('chr', 'pos', 'a1', 'a2')] <- str_split_fixed(reverse_snp_1kg37$V1, ':', 4)

reverse_snp_meta$chr_pos <- paste(reverse_snp_meta$chr, reverse_snp_meta$pos, sep=":")
reverse_snp_topmed$chr_pos <- paste(reverse_snp_topmed$chr, reverse_snp_topmed$pos, sep=":")
reverse_snp_ex1kg$chr_pos <- paste(reverse_snp_ex1kg$chr, reverse_snp_ex1kg$pos, sep=":")
reverse_snp_1kg38$chr_pos <- paste(reverse_snp_1kg38$chr, reverse_snp_1kg38$pos, sep=":")
reverse_snp_1kg37$chr_pos <- paste(reverse_snp_1kg37$chr, reverse_snp_1kg37$pos, sep=":")

# check those < -0.2
big_diff_meta <- rsq_meta[which(rsq_meta$Validation_AF < 0.5 & rsq_meta$"Imputation.AF" > 0.5),2] # 14 SNPs (14 -5 = 9 SNPs with big AF difference in sequenced AF and imputed AF))
big_diff_topmed <- rsq_topmed[which(rsq_topmed$Validation_AF < 0.5 & rsq_topmed$"Imputation.AF" > 0.5),2] # 14 SNPs  (14 -5 = 9 SNPs with big AF difference in sequenced AF and imputed AF))
big_diff_ex1kg <- rsq_ex1kg[which(rsq_ex1kg$Validation_AF < 0.5 & rsq_ex1kg$"Imputation.AF" > 0.5),2] # 11 SNPs  (11-5 = 6 SNPs with big AF difference in sequenced AF and imputed AF)
big_diff_1kg38 <- rsq_1kg38[which(rsq_1kg38$Validation_AF < 0.5 & rsq_1kg38$"Imputation.AF" > 0.5),2] # 12 SNPs  (12-5 = 7 SNPs with big AF difference in sequenced AF and imputed AF)
big_diff_1kg37 <- rsq_1kg37[which(rsq_1kg37$Validation_AF < 0.5 & rsq_1kg37$"Imputation.AF" > 0.5),2] # 14 SNPs  (14 -5 = 9 SNPs with big AF difference in sequenced AF and imputed AF)
big_diff_ga <- rsq_ga[which(rsq_ga$Validation_AF < 0.5 & rsq_ga$"Imputation.AF" > 0.5),2] # 5 SNPs  (5-5 = 0 SNPs)

length(big_diff_meta)
length(big_diff_topmed)
length(big_diff_ex1kg)
length(big_diff_1kg38)
length(big_diff_1kg37)
length(big_diff_ga)

rsq_meta$AF_sum <- as.numeric(rsq_meta$Validation_AF) + as.numeric(rsq_meta$"Imputation.AF")
rsq_topmed$AF_sum <- as.numeric(rsq_topmed$Validation_AF) + as.numeric(rsq_topmed$"Imputation.AF")
rsq_ex1kg$AF_sum <- as.numeric(rsq_ex1kg$Validation_AF) + as.numeric(rsq_ex1kg$"Imputation.AF")
rsq_1kg38$AF_sum <- as.numeric(rsq_1kg38$Validation_AF) + as.numeric(rsq_1kg38$"Imputation.AF")
rsq_1kg37$AF_sum <- as.numeric(rsq_1kg37$Validation_AF) + as.numeric(rsq_1kg37$"Imputation.AF")
rsq_ga$AF_sum <- as.numeric(rsq_ga$Validation_AF) + as.numeric(rsq_ga$"Imputation.AF")

rsq_meta[which(rsq_meta$Validation_AF < 0.5 & rsq_meta$"Imputation.AF" > 0.5),]
rsq_topmed[which(rsq_topmed$Validation_AF < 0.5 & rsq_topmed$"Imputation.AF" > 0.5),]
rsq_ex1kg[which(rsq_ex1kg$Validation_AF < 0.5 & rsq_ex1kg$"Imputation.AF" > 0.5),]
rsq_1kg38[which(rsq_1kg38$Validation_AF < 0.5 & rsq_1kg38$"Imputation.AF" > 0.5),]
rsq_1kg37[which(rsq_1kg37$Validation_AF < 0.5 & rsq_1kg37$"Imputation.AF" > 0.5),]
rsq_ga[which(rsq_ga$Validation_AF < 0.5 & rsq_ga$"Imputation.AF" > 0.5),]

list_big_sequence_imputation_DIFF_37 <- c("10:33018351:G:C","11:65389890:G:C","15:90764359:T:G","15:90771750:T:G","17:48656088:G:A")
list_big_sequence_imputation_DIFF_38 <- c("chr10:32729423:G:C","chr11:65622419:G:C","chr15:90221127:T:G","chr15:90228518:T:G","chr17:50578727:G:A")

sum(big_diff_meta %in% list_big_sequence_imputation_DIFF_38) # 5 SNPs
sum(big_diff_topmed %in% list_big_sequence_imputation_DIFF_38) # 5 SNPs
sum(big_diff_ex1kg %in% list_big_sequence_imputation_DIFF_38) # 5 SNPs
sum(big_diff_1kg38 %in% list_big_sequence_imputation_DIFF_38) # 5 SNPs
sum(big_diff_1kg37 %in% list_big_sequence_imputation_DIFF_37) # 5 SNPs
sum(big_diff_ga %in% list_big_sequence_imputation_DIFF_37) # 5 SNPs
# these 5 SNPs are imputed badly for all imputation panels (so perhaps it's called badly by sequencing)


big_diff_meta_chr_pos <- rsq_meta[which(rsq_meta$Validation_AF < 0.5 & rsq_meta$"Imputation.AF" > 0.5),1]  # 14 SNPs
big_diff_topmed_chr_pos  <- rsq_topmed[which(rsq_topmed$Validation_AF < 0.5 & rsq_topmed$"Imputation.AF" > 0.5),1]  # 14 SNPs
big_diff_ex1kg_chr_pos  <- rsq_ex1kg[which(rsq_ex1kg$Validation_AF < 0.5 & rsq_ex1kg$"Imputation.AF" > 0.5),1]   # 11 SNPs
big_diff_1kg38_chr_pos  <- rsq_1kg38[which(rsq_1kg38$Validation_AF < 0.5 & rsq_1kg38$"Imputation.AF" > 0.5),1]  # 12 SNPs
big_diff_1kg37_chr_pos  <- rsq_1kg37[which(rsq_1kg37$Validation_AF < 0.5 & rsq_1kg37$"Imputation.AF" > 0.5),1] # 14 SNPs
big_diff_ga_chr_pos  <- rsq_ga[which(rsq_ga$Validation_AF < 0.5 & rsq_ga$"Imputation.AF" > 0.5),1] # 5 SNPs
 

sum(big_diff_meta_chr_pos %in% reverse_snp_meta$chr_pos) # 0 SNPs with big AF difference in sequencing and imputation (one has AF < 0.5 and another has AF > 0.5) is due to the allele flip 
sum(big_diff_topmed_chr_pos %in% reverse_snp_topmed$chr_pos) # 0 SNPs with big AF difference in sequencing and imputation (one has AF < 0.5 and another has AF > 0.5) is due to the allele flip 
sum(big_diff_ex1kg_chr_pos %in% reverse_snp_ex1kg$chr_pos) # 0 SNPs with big AF difference in sequencing and imputation (one has AF < 0.5 and another has AF > 0.5) is due to the allele flip 
sum(big_diff_1kg38_chr_pos %in% reverse_snp_1kg38$chr_pos) # 0 SNPs with big AF difference in sequencing and imputation (one has AF < 0.5 and another has AF > 0.5) is due to the allele flip 
sum(big_diff_1kg37_chr_pos %in% reverse_snp_1kg37$chr_pos) # 0 SNPs with big AF difference in sequencing and imputation (one has AF < 0.5 and another has AF > 0.5) is due to the allele flip 

# check overlap across different imputation panels
sum(big_diff_meta %in% big_diff_topmed) # 14 out of 14 overlapped
sum(big_diff_ex1kg %in% big_diff_topmed) # 11 out of 11 overlapped
sum(big_diff_1kg38 %in% big_diff_topmed) # 10 out of 12 overlapped
sum(big_diff_1kg37_grch38 %in% big_diff_topmed) # 8 out of 12 overlapped
sum(big_diff_1kg37_grch38 %in% big_diff_1kg38) # 9 out of 12 overlapped

# 14 SNPs in 1kg37 --> convert to grch38 version 
big_diff_1kg37_grch38 <- c("chr10:32729423:G:C", "chr11:65622419:G:C", "chr15:90221127:T:G", "chr15:90228518:T:G", "chr16:10182601:TC:T", "chr17:50578727:G:A", "chr17:50590652:G:C", "chr2:165315807:A:G", "chr2:166035938:A:T", "chr2:166456819:T:G", "chr2:166456827:T:G", "chr2:195716571:T:TGTCTAG")
# 10:21177156:GA:G, chr17:48668021:AT:A" cannot find grch38 location
# check imputation AF and sequencing AF

# check if any of these are in the reverse SNP files


rsq_meta$DIFF_AF_new <- as.numeric(rsq_meta$Validation_AF) - as.numeric(rsq_meta$EUR_AF)
rsq_topmed$DIFF_AF_new <- as.numeric(rsq_topmed$Validation_AF) - as.numeric(rsq_topmed$EUR_AF)
rsq_ex1kg$DIFF_AF_new <- as.numeric(rsq_ex1kg$Validation_AF) - as.numeric(rsq_ex1kg$EUR_AF)
rsq_1kg38$DIFF_AF_new <- as.numeric(rsq_1kg38$Validation_AF) - as.numeric(rsq_1kg38$EUR_AF)
rsq_1kg37$DIFF_AF_new <- as.numeric(rsq_1kg37$Validation_AF) - as.numeric(rsq_1kg37$EUR_AF)
rsq_ga$DIFF_AF_new <- as.numeric(rsq_ga$Validation_AF) - as.numeric(rsq_ga$EUR_AF)

par(mfrow=c(3,2))
plot(rsq_meta$DIFF_AF_new ,rsq_meta$Imputation_R2,   pch=20, frame = FALSE)
plot(rsq_topmed$DIFF_AF_new ,rsq_topmed$Imputation_R2,   pch=20, frame = FALSE)
plot(rsq_ex1kg$DIFF_AF_new ,rsq_ex1kg$Imputation_R2,   pch=20, frame = FALSE)
plot(rsq_1kg38$DIFF_AF_new ,rsq_1kg38$Imputation_R2,   pch=20, frame = FALSE)
plot(rsq_1kg37$DIFF_AF_new ,rsq_1kg37$Imputation_R2,   pch=20, frame = FALSE)
plot(rsq_ga$DIFF_AF_new ,rsq_ga$Imputation_R2,   pch=20, frame = FALSE)

dim(rsq_meta[which(rsq_meta$DIFF_AF_new < -0.4),]) # 7 SNPs 
dim(rsq_topmed[which(rsq_topmed$DIFF_AF_new < -0.4),]) # 7 SNPs
dim(rsq_ex1kg[which(rsq_ex1kg$DIFF_AF_new < -0.4),]) # 6 SNPs
dim(rsq_1kg38[which(rsq_1kg38$DIFF_AF_new < -0.4),]) # 8 SNPs
dim(rsq_1kg37[which(rsq_1kg37$DIFF_AF_new < -0.4),]) # 10 SNPs
dim(rsq_ga[which(rsq_ga$DIFF_AF_new < -0.4),]) # 5 SNPs
# all the imputation AF for these SNPs are similar to those AF in 1000G
# the mean difference exist in between imputation AF and sequencing AF

big_diff_meta_low_R2 <- rsq_meta[which(rsq_meta$DIFF_AF_new < -0.4),2] # 7 SNPs out of 5913 SNPs (0.12 %)
big_diff_topmed_low_R2 <- rsq_topmed[which(rsq_topmed$DIFF_AF_new < -0.4),2] # 7 SNPs out of 5701 SNPs (0.12 %)
big_diff_ex1kg_low_R2 <- rsq_ex1kg[which(rsq_ex1kg$DIFF_AF_new < -0.4),2] # 6 SNPs out of 3425 SNPs (0.18 %)
big_diff_1kg38_low_R2 <- rsq_1kg38[which(rsq_1kg38$DIFF_AF_new < -0.4),2] # 8 SNPs out of 3333 SNPs (0.24 %)
big_diff_1kg37_low_R2 <- rsq_1kg37[which(rsq_1kg37$DIFF_AF_new < -0.4),2] # 10 SNPs out of 3450 SNPs (0.29 %)
big_diff_ga_low_R2 <- rsq_ga[which(rsq_ga$DIFF_AF_new < -0.4),2] # 5 SNPs out of 2066 SNPs (0.24 %)

rsq_meta[which(rsq_meta$DIFF_AF_new < -0.4),] # 7 SNPs out of 5913 SNPs (0.12 %)
rsq_topmed[which(rsq_topmed$DIFF_AF_new < -0.4),] # 7 SNPs out of 5701 SNPs (0.12 %)
rsq_ex1kg[which(rsq_ex1kg$DIFF_AF_new < -0.4),] # 6 SNPs out of 3425 SNPs (0.18 %)
rsq_1kg38[which(rsq_1kg38$DIFF_AF_new < -0.4),] # 8 SNPs out of 3333 SNPs (0.24 %)
rsq_1kg37[which(rsq_1kg37$DIFF_AF_new < -0.4),] # 10 SNPs out of 3450 SNPs (0.29 %)
rsq_ga[which(rsq_ga$DIFF_AF_new < -0.4),] # 5 SNPs out of 2066 SNPs (0.24 %)


sum(big_diff_meta_low_R2 %in% big_diff_topmed_low_R2) # 7 out of 7 overlapped
sum(big_diff_ex1kg_low_R2 %in% big_diff_topmed_low_R2) # 6 out of 6 overlapped
sum(big_diff_topmed_low_R2 %in% big_diff_1kg38_low_R2) # 6 out of 7 overlapped
sum(big_diff_ga_low_R2 %in% big_diff_1kg37_low_R2) # 5 out of 5 overlapped
sum(big_diff_1kg38_low_R2_grch37 %in% big_diff_1kg37_low_R2) # 8 out of 8 overlapped

# conversion of big_diff_1kg38_low_R2 to grch37 
big_diff_1kg38_low_R2_grch37 <- c("10:33018351:G:C" , "10:33123774:T:C" , "10:33123824:G:C", "15:90771750:T:G" , "16:10276458:TC:T", "17:48656088:G:A", "2:166172317:A:G",  "2:167313337:T:G" )


par(mfrow=c(3,2))
plot(abs(rsq_meta$DIFF_AF_new ),rsq_meta$Imputation_R2,   pch=20, frame = FALSE)
plot(abs(rsq_topmed$DIFF_AF_new) ,rsq_topmed$Imputation_R2,   pch=20, frame = FALSE)
plot(abs(rsq_ex1kg$DIFF_AF_new ),rsq_ex1kg$Imputation_R2,   pch=20, frame = FALSE)
plot(abs(rsq_1kg38$DIFF_AF_new) ,rsq_1kg38$Imputation_R2,   pch=20, frame = FALSE)
plot(abs(rsq_1kg37$DIFF_AF_new ),rsq_1kg37$Imputation_R2,   pch=20, frame = FALSE)
plot(abs(rsq_ga$DIFF_AF_new) ,rsq_ga$Imputation_R2,   pch=20, frame = FALSE)


rsq_meta$outlier <- ifelse(abs(rsq_meta$DIFF_AF_new) > 0.4, "yes", "no")
rsq_topmed$outlier <- ifelse(abs(rsq_topmed$DIFF_AF_new) > 0.4, "yes", "no")
rsq_ex1kg$outlier <- ifelse(abs(rsq_ex1kg$DIFF_AF_new) > 0.4, "yes", "no")
rsq_1kg38$outlier <- ifelse(abs(rsq_1kg38$DIFF_AF_new) > 0.4, "yes", "no")
rsq_1kg37$outlier <- ifelse(abs(rsq_1kg37$DIFF_AF_new) > 0.4, "yes", "no")
rsq_ga$outlier <- ifelse(abs(rsq_ga$DIFF_AF_new) > 0.4, "yes", "no")


## plot for the manuscript (outlier dots as red)
library(ggplot2)
ggplot(rsq_meta, aes(x = abs(DIFF_AF_new), y = Imputation_R2 ,group=outlier, color=outlier))  + geom_point(shape=19,size=3)+ theme_classic(base_size=30) + coord_cartesian(ylim = c(0, 1)) +  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + theme(text = element_text(size = 30,colour="black")) + scale_color_manual(values=c("#999999","#D55E00")) + theme(legend.position = "none") + labs(x="Absolute difference from EUR AF", y="True R2") 
ggplot(rsq_topmed, aes(x = abs(DIFF_AF_new), y = Imputation_R2 ,group=outlier, color=outlier))  + geom_point(shape=19,size=3)+ theme_classic(base_size=30) + coord_cartesian(ylim = c(0, 1)) +  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + theme(text = element_text(size = 30,colour="black")) + scale_color_manual(values=c("#999999","#D55E00")) + theme(legend.position = "none") + labs(x="Absolute difference from EUR AF", y="True R2") 
ggplot(rsq_ex1kg, aes(x = abs(DIFF_AF_new), y = Imputation_R2 ,group=outlier, color=outlier))  + geom_point(shape=19,size=3)+ theme_classic(base_size=30) + coord_cartesian(ylim = c(0, 1)) +  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + theme(text = element_text(size = 30,colour="black")) + scale_color_manual(values=c("#999999","#D55E00")) + theme(legend.position = "none") + labs(x="Absolute difference from EUR AF", y="True R2") 
ggplot(rsq_1kg38, aes(x = abs(DIFF_AF_new), y = Imputation_R2 ,group=outlier, color=outlier))  + geom_point(shape=19,size=3)+ theme_classic(base_size=30) + coord_cartesian(ylim = c(0, 1)) +  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + theme(text = element_text(size = 30,colour="black")) + scale_color_manual(values=c("#999999","#D55E00")) + theme(legend.position = "none") + labs(x="Absolute difference from EUR AF", y="True R2") 
ggplot(rsq_1kg37, aes(x = abs(DIFF_AF_new), y = Imputation_R2 ,group=outlier, color=outlier))  + geom_point(shape=19,size=3)+ theme_classic(base_size=30) + coord_cartesian(ylim = c(0, 1)) +  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + theme(text = element_text(size = 30,colour="black")) + scale_color_manual(values=c("#999999","#D55E00")) + theme(legend.position = "none") + labs(x="Absolute difference from EUR AF", y="True R2") 
ggplot(rsq_ga, aes(x = abs(DIFF_AF_new), y = Imputation_R2 ,group=outlier, color=outlier))  + geom_point(shape=19,size=3)+ theme_classic(base_size=30) + coord_cartesian(ylim = c(0, 1)) +  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + theme(text = element_text(size = 30,colour="black")) + scale_color_manual(values=c("#999999","#D55E00")) + theme(legend.position = "none") + labs(x="Absolute difference from EUR AF", y="True R2") 
 
 
summary(lm(Imputation_R2 ~ abs(DIFF_AF_new) , data = rsq_meta)) # p <2e-16
summary(lm(Imputation_R2 ~ abs(DIFF_AF_new) , data = rsq_topmed)) # p <2e-16
summary(lm(Imputation_R2 ~ abs(DIFF_AF_new) , data = rsq_ex1kg)) # p <2e-16
summary(lm(Imputation_R2 ~ abs(DIFF_AF_new) , data = rsq_1kg38)) # p <2e-16
summary(lm(Imputation_R2 ~ abs(DIFF_AF_new) , data = rsq_1kg37)) # p <2e-16
summary(lm(Imputation_R2 ~ abs(DIFF_AF_new) , data = rsq_ga)) # p <2e-16


# the greater the DIFF_AF, the better the imputation_R2 because it is common variants

####### what about common alleles?
rsq_meta_common <- rsq_meta[which(rsq_meta$Validation_AF >= 0.01 ),] # 1877 SNPs
rsq_topmed_common <- rsq_topmed[which(rsq_topmed$Validation_AF >= 0.01 ),] # 1875 SNPs
rsq_ex1kg_common <- rsq_ex1kg[which(rsq_ex1kg$Validation_AF >= 0.01 ),]  # 1288 SNPs
rsq_1kg38_common <- rsq_1kg38[which(rsq_1kg38$Validation_AF >= 0.01 ),]  # 1256 SNPs
rsq_1kg37_common <- rsq_1kg37 [which(rsq_1kg37$Validation_AF >= 0.01 ),] # 1300 SNPs
rsq_ga_common <- rsq_ga[which(rsq_ga$Validation_AF >= 0.01 ),] # 1086 SNPs

dim(rsq_meta_common)
dim(rsq_topmed_common)
dim(rsq_ex1kg_common)
dim(rsq_1kg38_common)
dim(rsq_1kg37_common)
dim(rsq_ga_common)

par(mfrow=c(3,2))
plot(rsq_meta_common$DIFF_AF_new ,rsq_meta_common$Imputation_R2,   pch=20, frame = FALSE)
plot(rsq_topmed_common$DIFF_AF_new ,rsq_topmed_common$Imputation_R2,   pch=20, frame = FALSE)
plot(rsq_ga_common$DIFF_AF_new ,rsq_ga_common$Imputation_R2,   pch=20, frame = FALSE)
plot(rsq_ex1kg_common$DIFF_AF_new ,rsq_ex1kg_common$Imputation_R2,   pch=20, frame = FALSE)
plot(rsq_1kg38_common$DIFF_AF_new ,rsq_1kg38_common$Imputation_R2,   pch=20, frame = FALSE)
plot(rsq_1kg37_common$DIFF_AF_new ,rsq_1kg37_common$Imputation_R2,   pch=20, frame = FALSE) 

# the EUR MAF does not seem to matter as long as the imputation panel has SAS participants

rsq_meta_common$outlier <- ifelse(abs(rsq_meta_common$DIFF_AF_new) > 0.4, "yes", "no")
rsq_topmed_common$outlier <- ifelse(abs(rsq_topmed_common$DIFF_AF_new) > 0.4, "yes", "no")
rsq_ex1kg_common$outlier <- ifelse(abs(rsq_ex1kg_common$DIFF_AF_new) > 0.4, "yes", "no")
rsq_1kg38_common$outlier <- ifelse(abs(rsq_1kg38_common$DIFF_AF_new) > 0.4, "yes", "no")
rsq_1kg37_common$outlier <- ifelse(abs(rsq_1kg37_common$DIFF_AF_new) > 0.4, "yes", "no")
rsq_ga_common$outlier <- ifelse(abs(rsq_ga_common$DIFF_AF_new) > 0.4, "yes", "no")

library(ggplot2)
ggplot(rsq_meta_common, aes(x = abs(DIFF_AF_new), y = Imputation_R2 ,group=outlier, color=outlier))  + geom_point(shape=19,size=3)+ theme_classic(base_size=30) + coord_cartesian(ylim = c(0, 1)) +  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + theme(text = element_text(size = 30,colour="black")) + scale_color_manual(values=c("#999999","#D55E00")) + theme(legend.position = "none") + labs(x="Absolute difference from EUR AF", y="True R2") 
ggplot(rsq_topmed_common, aes(x = abs(DIFF_AF_new), y = Imputation_R2 ,group=outlier, color=outlier))  + geom_point(shape=19,size=3)+ theme_classic(base_size=30) + coord_cartesian(ylim = c(0, 1)) +  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + theme(text = element_text(size = 30,colour="black")) + scale_color_manual(values=c("#999999","#D55E00")) + theme(legend.position = "none") + labs(x="Absolute difference from EUR AF", y="True R2") 
ggplot(rsq_ex1kg_common, aes(x = abs(DIFF_AF_new), y = Imputation_R2 ,group=outlier, color=outlier))  + geom_point(shape=19,size=3)+ theme_classic(base_size=30) + coord_cartesian(ylim = c(0, 1)) +  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + theme(text = element_text(size = 30,colour="black")) + scale_color_manual(values=c("#999999","#D55E00")) + theme(legend.position = "none") + labs(x="Absolute difference from EUR AF", y="True R2") 
ggplot(rsq_1kg38_common, aes(x = abs(DIFF_AF_new), y = Imputation_R2 ,group=outlier, color=outlier))  + geom_point(shape=19,size=3)+ theme_classic(base_size=30) + coord_cartesian(ylim = c(0, 1)) +  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + theme(text = element_text(size = 30,colour="black")) + scale_color_manual(values=c("#999999","#D55E00")) + theme(legend.position = "none") + labs(x="Absolute difference from EUR AF", y="True R2") 
ggplot(rsq_1kg37_common, aes(x = abs(DIFF_AF_new), y = Imputation_R2 ,group=outlier, color=outlier))  + geom_point(shape=19,size=3)+ theme_classic(base_size=30) + coord_cartesian(ylim = c(0, 1)) +  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + theme(text = element_text(size = 30,colour="black")) + scale_color_manual(values=c("#999999","#D55E00")) + theme(legend.position = "none") + labs(x="Absolute difference from EUR AF", y="True R2") 
ggplot(rsq_ga_common, aes(x = abs(DIFF_AF_new), y = Imputation_R2 ,group=outlier, color=outlier))  + geom_point(shape=19,size=3)+ theme_classic(base_size=30) + coord_cartesian(ylim = c(0, 1)) +  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + theme(text = element_text(size = 30,colour="black")) + scale_color_manual(values=c("#999999","#D55E00")) + theme(legend.position = "none") + labs(x="Absolute difference from EUR AF", y="True R2") 
 
summary(lm(Imputation_R2 ~ abs(DIFF_AF_new) , data = rsq_meta_common)) # p =  0.324 
summary(lm(Imputation_R2 ~ abs(DIFF_AF_new) , data = rsq_topmed_common)) # p  = 0.685
summary(lm(Imputation_R2 ~ abs(DIFF_AF_new) , data = rsq_ex1kg_common)) # p =   0.873  
summary(lm(Imputation_R2 ~ abs(DIFF_AF_new), data = rsq_1kg38_common)) # p =   0.146  
summary(lm(Imputation_R2 ~ abs(DIFF_AF_new) , data = rsq_1kg37_common)) # p = 0.191  
summary(lm(Imputation_R2 ~ abs(DIFF_AF_new) , data = rsq_ga_common)) # p =   0.0289

big_diff_meta_low_R2_common <- rsq_meta_common[which(rsq_meta_common$DIFF_AF_new < -0.4),2] # 6 SNPs 
big_diff_topmed_low_R2_common <- rsq_topmed_common[which(rsq_topmed_common$DIFF_AF_new < -0.4),2] # 6 SNPs
big_diff_ex1kg_low_R2_common <- rsq_ex1kg_common[which(rsq_ex1kg_common$DIFF_AF_new < -0.4),2] # 5 SNPs 
big_diff_1kg38_low_R2_common <- rsq_1kg38_common[which(rsq_1kg38_common$DIFF_AF_new < -0.4),2] # 7 SNPs
big_diff_1kg37_low_R2_common<- rsq_1kg37_common[which(rsq_1kg37_common$DIFF_AF_new < -0.4),2] # 7 SNPs
big_diff_ga_low_R2_common <- rsq_ga_common[which(rsq_ga_common$DIFF_AF_new < -0.4),2] # 4 SNPs 

length(big_diff_meta_low_R2_common)
length(big_diff_topmed_low_R2_common)
length(big_diff_ex1kg_low_R2_common)
length(big_diff_1kg38_low_R2_common)
length(big_diff_1kg37_low_R2_common)
length(big_diff_ga_low_R2_common)

tmp <- unique(c(big_diff_meta_low_R2_common,big_diff_topmed_low_R2_common,big_diff_ex1kg_low_R2_common,big_diff_1kg38_low_R2_common,big_diff_1kg37_low_R2_common,big_diff_ga_low_R2_common))


dim(rsq_meta_common[which(rsq_meta_common$DIFF_AF_new < -0.4),]) # 6 SNPs  (all SNPs n = 7, among these, 1 SNPs is rare SNP) 
dim(rsq_topmed_common[which(rsq_topmed_common$DIFF_AF_new < -0.4),]) # 6 SNPs (all SNPs n = 7, among these, 1 SNPs is rare SNP) 
dim(rsq_ex1kg_common[which(rsq_ex1kg_common$DIFF_AF_new < -0.4),]) # 5 SNPs  (all SNPs n = 6, among these, 1 SNPs is rare SNP) 
dim(rsq_1kg38_common[which(rsq_1kg38_common$DIFF_AF_new < -0.4),]) # 7 SNPs (all SNPs n = 8, among these, 1 SNPs is rare SNP) 
dim(rsq_1kg37_common[which(rsq_1kg37_common$DIFF_AF_new < -0.4),]) # 7 SNPs (all SNPs n = 10, among these, 3 SNPs is rare SNP) 
dim(rsq_ga_common[which(rsq_ga_common$DIFF_AF_new < -0.4),]) # 4 SNPs  (all SNPs n = 5, among these, 1 SNPs is rare SNP) 



# what about rare alleles? 
rsq_meta_rare <- rsq_meta[which(rsq_meta$Validation_AF < 0.01 ),] # 4036 SNPs
rsq_topmed_rare <- rsq_topmed[which(rsq_topmed$Validation_AF < 0.01 ),] # 3826 SNPs
rsq_ex1kg_rare <- rsq_ex1kg[which(rsq_ex1kg$Validation_AF < 0.01 ),]  # 2137 SNPs
rsq_1kg38_rare <- rsq_1kg38[which(rsq_1kg38$Validation_AF < 0.01 ),]  # 2077 SNPs
rsq_1kg37_rare <- rsq_1kg37 [which(rsq_1kg37$Validation_AF < 0.01 ),] # 2150 SNPs
rsq_ga_rare <- rsq_ga[which(rsq_ga$Validation_AF < 0.01 ),] # 980 SNPs

dim(rsq_meta_rare)
dim(rsq_topmed_rare)
dim(rsq_ex1kg_rare)
dim(rsq_1kg38_rare)
dim(rsq_1kg37_rare)
dim(rsq_ga_rare)

par(mfrow=c(3,2))
plot(rsq_meta_rare$DIFF_AF_new ,rsq_meta_rare$Imputation_R2,   pch=20, frame = FALSE)
plot(rsq_topmed_rare$DIFF_AF_new ,rsq_topmed_rare$Imputation_R2,   pch=20, frame = FALSE)
plot(rsq_ga_rare$DIFF_AF_new ,rsq_ga_rare$Imputation_R2,   pch=20, frame = FALSE)
plot(rsq_ex1kg_rare$DIFF_AF_new ,rsq_ex1kg_rare$Imputation_R2,   pch=20, frame = FALSE)
plot(rsq_1kg38_rare$DIFF_AF_new ,rsq_1kg38_rare$Imputation_R2,   pch=20, frame = FALSE)
plot(rsq_1kg37_rare$DIFF_AF_new ,rsq_1kg37_rare$Imputation_R2,   pch=20, frame = FALSE) 

# the EUR MAF does not seem to matter as long as the imputation panel has SAS participants

rsq_meta_rare$outlier <- ifelse(abs(rsq_meta_rare$DIFF_AF_new) > 0.1, "yes", "no")
rsq_topmed_rare$outlier <- ifelse(abs(rsq_topmed_rare$DIFF_AF_new) > 0.1, "yes", "no")
rsq_ex1kg_rare$outlier <- ifelse(abs(rsq_ex1kg_rare$DIFF_AF_new) > 0.1, "yes", "no")
rsq_1kg38_rare$outlier <- ifelse(abs(rsq_1kg38_rare$DIFF_AF_new) > 0.1, "yes", "no")
rsq_1kg37_rare$outlier <- ifelse(abs(rsq_1kg37_rare$DIFF_AF_new) > 0.1, "yes", "no")
rsq_ga_rare$outlier <- ifelse(abs(rsq_ga_rare$DIFF_AF_new) > 0.1, "yes", "no")

library(ggplot2)
ggplot(rsq_meta_rare, aes(x = abs(DIFF_AF_new), y = Imputation_R2 ,group=outlier, color=outlier))  + geom_point(shape=19,size=3)+ theme_classic(base_size=30) + coord_cartesian(ylim = c(0, 1)) +  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4)) + theme(text = element_text(size = 30,colour="black")) + scale_color_manual(values=c("#999999","#D55E00")) + theme(legend.position = "none") + labs(x="Absolute difference from EUR AF", y="True R2") 
ggplot(rsq_topmed_rare, aes(x = abs(DIFF_AF_new), y = Imputation_R2 ,group=outlier, color=outlier))  + geom_point(shape=19,size=3)+ theme_classic(base_size=30) + coord_cartesian(ylim = c(0, 1)) +  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4)) + theme(text = element_text(size = 30,colour="black")) + scale_color_manual(values=c("#999999","#D55E00")) + theme(legend.position = "none") + labs(x="Absolute difference from EUR AF", y="True R2") 
ggplot(rsq_ex1kg_rare, aes(x = abs(DIFF_AF_new), y = Imputation_R2 ,group=outlier, color=outlier))  + geom_point(shape=19,size=3)+ theme_classic(base_size=30) + coord_cartesian(ylim = c(0, 1)) +  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4)) + theme(text = element_text(size = 30,colour="black")) + scale_color_manual(values=c("#999999","#D55E00")) + theme(legend.position = "none") + labs(x="Absolute difference from EUR AF", y="True R2") 
ggplot(rsq_1kg38_rare, aes(x = abs(DIFF_AF_new), y = Imputation_R2 ,group=outlier, color=outlier))  + geom_point(shape=19,size=3)+ theme_classic(base_size=30) + coord_cartesian(ylim = c(0, 1)) +  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4)) + theme(text = element_text(size = 30,colour="black")) + scale_color_manual(values=c("#999999","#D55E00")) + theme(legend.position = "none") + labs(x="Absolute difference from EUR AF", y="True R2") 
ggplot(rsq_1kg37_rare, aes(x = abs(DIFF_AF_new), y = Imputation_R2 ,group=outlier, color=outlier))  + geom_point(shape=19,size=3)+ theme_classic(base_size=30) + coord_cartesian(ylim = c(0, 1)) +  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4,0.5,0.6,0.7,0.8,0.9,1.0)) + theme(text = element_text(size = 30,colour="black")) + scale_color_manual(values=c("#999999","#D55E00")) + theme(legend.position = "none") + labs(x="Absolute difference from EUR AF", y="True R2") 
ggplot(rsq_ga_rare, aes(x = abs(DIFF_AF_new), y = Imputation_R2 ,group=outlier, color=outlier))  + geom_point(shape=19,size=3)+ theme_classic(base_size=30) + coord_cartesian(ylim = c(0, 1)) +  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4)) + theme(text = element_text(size = 30,colour="black")) + scale_color_manual(values=c("#999999","#D55E00")) + theme(legend.position = "none") + labs(x="Absolute difference from EUR AF", y="True R2") 
 
summary(lm(Imputation_R2 ~ abs(DIFF_AF_new) , data = rsq_meta_rare)) # p   2.4e-07 
summary(lm(Imputation_R2 ~ abs(DIFF_AF_new) , data = rsq_topmed_rare)) # p  9.12e-07
summary(lm(Imputation_R2 ~ abs(DIFF_AF_new) , data = rsq_ex1kg_rare)) # p =   0.00337 ** 
summary(lm(Imputation_R2 ~ abs(DIFF_AF_new), data = rsq_1kg38_rare)) # p =   0.0197 *  
summary(lm(Imputation_R2 ~ abs(DIFF_AF_new) , data = rsq_1kg37_rare)) # p = 0.851 
summary(lm(Imputation_R2 ~ abs(DIFF_AF_new) , data = rsq_ga_rare)) # p =   0.324  

big_diff_meta_low_R2_rare <- rsq_meta_rare[which(rsq_meta_rare$DIFF_AF_new < -0.1),2] # 3 SNPs 
big_diff_topmed_low_R2_rare <- rsq_topmed_rare[which(rsq_topmed_rare$DIFF_AF_new < -0.1),2] # 3 SNPs
big_diff_ex1kg_low_R2_rare <- rsq_ex1kg_rare[which(rsq_ex1kg_rare$DIFF_AF_new < -0.1),2] # 2 SNPs 
big_diff_1kg38_low_R2_rare <- rsq_1kg38_rare[which(rsq_1kg38_rare$DIFF_AF_new < -0.1),2] # 3 SNPs
big_diff_1kg37_low_R2_rare <- rsq_1kg37_rare[which(rsq_1kg37_rare$DIFF_AF_new < -0.1),2] # 5 SNPs
big_diff_ga_low_R2_rare <- rsq_ga_rare[which(rsq_ga_rare$DIFF_AF_new < -0.1),2] # 1 SNPs 

rsq_meta_rare[which(rsq_meta_rare$DIFF_AF_new < -0.1),] # 3 SNPs 
rsq_topmed_rare[which(rsq_topmed_rare$DIFF_AF_new < -0.1),] # 3 SNPs
rsq_ex1kg_rare[which(rsq_ex1kg_rare$DIFF_AF_new < -0.1),] # 2 SNPs 
rsq_1kg38_rare[which(rsq_1kg38_rare$DIFF_AF_new < -0.1),] # 3 SNPs
rsq_1kg37_rare[which(rsq_1kg37_rare$DIFF_AF_new < -0.1),] # 5 SNPs
rsq_ga_rare[which(rsq_ga_rare$DIFF_AF_new < -0.1),] # 1 SNPs 


length(big_diff_meta_low_R2_rare)
length(big_diff_topmed_low_R2_rare)
length(big_diff_ex1kg_low_R2_rare)
length(big_diff_1kg38_low_R2_rare)
length(big_diff_1kg37_low_R2_rare)
length(big_diff_ga_low_R2_rare)

rsq_meta_rare[which(rsq_meta_rare$DIFF_AF_new < -0.1),] # 3 SNPs (in between -0.1 and -0.4: there are 2 SNPs out of 5913 SNPs + 7 = 9 = 0.15%) 
rsq_topmed_rare[which(rsq_topmed_rare$DIFF_AF_new < -0.1),] # 3 SNPs (in between -0.1 and -0.4: there are 2 SNPs out of 5701 SNPs + 7 = 9 = 0.16%)
rsq_ex1kg_rare[which(rsq_ex1kg_rare$DIFF_AF_new < -0.1),] # 2 SNPs (in between -0.1 and -0.4: there are 1 SNPs out of 3425 SNPs + 6 = 7 = 0.20%)
rsq_1kg38_rare[which(rsq_1kg38_rare$DIFF_AF_new < -0.1),] # 3 SNPs (in between -0.1 and -0.4: there are 2 SNPs out of 3333 SNPs + 8 = 10 = 0.30%)
rsq_1kg37_rare[which(rsq_1kg37_rare$DIFF_AF_new < -0.1),] # 5 SNPs (in between -0.1 and -0.4: there are 2 SNPs out of 3450 SNPs + 10 = 12 = 0.35%)
rsq_ga_rare[which(rsq_ga_rare$DIFF_AF_new < -0.1),] # 1 SNPs (in between -0.1 and -0.4: there are 0 SNPs out of 2066 SNPs + 5 = 5 = 0.24%)

# the imputation AF of these SNPs across different panels are similar to the 1000G SAS AF 

# the 5 SNPS are imputed by all imputation panels, so likely the sequencing is wrong
# but for some SNPs, if it's only imputed by certain panels but not others (e.g. chr15:43518718), then it's likely that the imputation panel has it wrong, and the sequencing could be right

# how many < -0.4
dim(rsq_meta_rare[which(rsq_meta_rare$DIFF_AF_new < -0.4),]) # 1 SNPs 
dim(rsq_topmed_rare[which(rsq_topmed_rare$DIFF_AF_new < -0.4),]) # 1 SNPs 
dim(rsq_ex1kg_rare[which(rsq_ex1kg_rare$DIFF_AF_new < -0.4),]) # 1 SNPs 
dim(rsq_1kg38_rare[which(rsq_1kg38_rare$DIFF_AF_new < -0.4),]) # 1 SNPs 
dim(rsq_1kg37_rare[which(rsq_1kg37_rare$DIFF_AF_new < -0.4),]) # 3 SNPs 
dim(rsq_ga_rare[which(rsq_ga_rare$DIFF_AF_new < -0.4),]) # 1 SNPs 

rsq_meta[which(rsq_meta$CHR_BP =="chr15:43518718"),]
rsq_topmed[which(rsq_topmed$CHR_BP =="chr15:43518718"),]
rsq_ex1kg[which(rsq_ex1kg$CHR_BP =="chr15:43518718"),]
rsq_1kg38[which(rsq_1kg38$CHR_BP =="chr15:43518718"),]
rsq_1kg37[which(rsq_1kg37$CHR_BP =="15:43810916"),]
rsq_ga[which(rsq_ga$CHR_BP =="15:43810916"),]
# this SNP only included in 1k38 and 1k37
# checked in ex1000G and found that this SNP is not imputed by expanded 1000G, topmed, meta-imputation or GenomeAsia

par(mfrow=c(3,2))
plot(abs(rsq_meta_rare$DIFF_AF_new) ,rsq_meta_rare$Imputation_R2,   pch=20, frame = FALSE)
plot(abs(rsq_topmed_rare$DIFF_AF_new)  ,rsq_topmed_rare$Imputation_R2,   pch=20, frame = FALSE)
plot(abs(rsq_ga_rare$DIFF_AF_new)  ,rsq_ga_rare$Imputation_R2,   pch=20, frame = FALSE)
plot(abs(rsq_ex1kg_rare$DIFF_AF_new)  ,rsq_ex1kg_rare$Imputation_R2,   pch=20, frame = FALSE)
plot(abs(rsq_1kg38_rare$DIFF_AF_new)  ,rsq_1kg38_rare$Imputation_R2,   pch=20, frame = FALSE)
plot(abs(rsq_1kg37_rare$DIFF_AF_new)  ,rsq_1kg37_rare$Imputation_R2,   pch=20, frame = FALSE) 



# what about ultra rare alleles? 
rsq_meta_ultra_rare <- rsq_meta[which(rsq_meta$Validation_AF < 0.001 ),] # 2196 SNPs
rsq_topmed_ultra_rare <- rsq_topmed[which(rsq_topmed$Validation_AF < 0.001 ),] # 2062 SNPs
rsq_ex1kg_ultra_rare <- rsq_ex1kg[which(rsq_ex1kg$Validation_AF < 0.001 ),]  # 1014 SNPs
rsq_1kg38_ultra_rare <- rsq_1kg38[which(rsq_1kg38$Validation_AF < 0.001 ),]  # 993 SNPs
rsq_1kg37_ultra_rare <- rsq_1kg37 [which(rsq_1kg37$Validation_AF < 0.001 ),] # 1019 SNPs
rsq_ga_ultra_rare <- rsq_ga[which(rsq_ga$Validation_AF < 0.001 ),] # 374 SNPs

par(mfrow=c(3,2))
plot(rsq_meta_ultra_rare$DIFF_AF_new ,rsq_meta_ultra_rare$Imputation_R2,   pch=20, frame = FALSE)
plot(rsq_topmed_ultra_rare$DIFF_AF_new ,rsq_topmed_ultra_rare$Imputation_R2,   pch=20, frame = FALSE)
plot(rsq_ex1kg_ultra_rare$DIFF_AF_new ,rsq_ex1kg_ultra_rare$Imputation_R2,   pch=20, frame = FALSE)
plot(rsq_1kg38_ultra_rare$DIFF_AF_new ,rsq_1kg38_ultra_rare$Imputation_R2,   pch=20, frame = FALSE)
plot(rsq_1kg37_ultra_rare$DIFF_AF_new ,rsq_1kg37_ultra_rare$Imputation_R2,   pch=20, frame = FALSE) 
plot(rsq_ga_ultra_rare$DIFF_AF_new ,rsq_ga_ultra_rare$Imputation_R2,   pch=20, frame = FALSE)

# the EUR MAF does not seem to matter as long as the imputation panel has SAS participants
rsq_meta_ultra_rare$outlier <- ifelse(abs(rsq_meta_ultra_rare$DIFF_AF_new) > 0.05, "yes", "no")
rsq_topmed_ultra_rare$outlier <- ifelse(abs(rsq_topmed_ultra_rare$DIFF_AF_new) > 0.05, "yes", "no")
rsq_ex1kg_ultra_rare$outlier <- ifelse(abs(rsq_ex1kg_ultra_rare$DIFF_AF_new) > 0.05, "yes", "no")
rsq_1kg38_ultra_rare$outlier <- ifelse(abs(rsq_1kg38_ultra_rare$DIFF_AF_new) > 0.05, "yes", "no")
rsq_1kg37_ultra_rare$outlier <- ifelse(abs(rsq_1kg37_ultra_rare$DIFF_AF_new) > 0.05, "yes", "no")
rsq_ga_ultra_rare$outlier <- ifelse(abs(rsq_ga_ultra_rare$DIFF_AF_new) > 0.05, "yes", "no")

table(rsq_meta_ultra_rare$outlier)
table(rsq_topmed_ultra_rare$outlier)
table(rsq_ex1kg_ultra_rare$outlier)
table(rsq_1kg38_ultra_rare$outlier)
table(rsq_1kg37_ultra_rare$outlier)
table(rsq_ga_ultra_rare$outlier)

library(ggplot2)
ggplot(rsq_meta_ultra_rare, aes(x = abs(DIFF_AF_new), y = Imputation_R2 ,group=outlier, color=outlier))  + geom_point(shape=19,size=3)+ theme_classic(base_size=30) + coord_cartesian(ylim = c(0, 1)) +  scale_x_continuous(breaks = c(0, 0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08)) + theme(text = element_text(size = 30,colour="black")) + scale_color_manual(values=c("#999999","#D55E00")) + theme(legend.position = "none") + labs(x="Absolute difference from EUR AF", y="True R2") 
ggplot(rsq_topmed_ultra_rare, aes(x = abs(DIFF_AF_new), y = Imputation_R2 ,group=outlier, color=outlier))  + geom_point(shape=19,size=3)+ theme_classic(base_size=30) + coord_cartesian(ylim = c(0, 1)) +  scale_x_continuous(breaks = c(0, 0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08)) + theme(text = element_text(size = 30,colour="black")) + scale_color_manual(values=c("#999999","#D55E00")) + theme(legend.position = "none") + labs(x="Absolute difference from EUR AF", y="True R2") 
ggplot(rsq_ex1kg_ultra_rare, aes(x = abs(DIFF_AF_new), y = Imputation_R2 ,group=outlier, color=outlier))  + geom_point(shape=19,size=3)+ theme_classic(base_size=30) + coord_cartesian(ylim = c(0, 1)) +  scale_x_continuous(breaks = c(0, 0.01,0.02)) + theme(text = element_text(size = 30,colour="black")) + scale_color_manual(values=c("#999999","#D55E00")) + theme(legend.position = "none") + labs(x="Absolute difference from EUR AF", y="True R2") 
ggplot(rsq_1kg38_ultra_rare, aes(x = abs(DIFF_AF_new), y = Imputation_R2 ,group=outlier, color=outlier))  + geom_point(shape=19,size=3)+ theme_classic(base_size=30) + coord_cartesian(ylim = c(0, 1)) +  scale_x_continuous(breaks = c(0, 0.01,0.02)) + theme(text = element_text(size = 30,colour="black")) + scale_color_manual(values=c("#999999","#D55E00")) + theme(legend.position = "none") + labs(x="Absolute difference from EUR AF", y="True R2") 
ggplot(rsq_1kg37_ultra_rare, aes(x = abs(DIFF_AF_new), y = Imputation_R2 ,group=outlier, color=outlier))  + geom_point(shape=19,size=3)+ theme_classic(base_size=30) + coord_cartesian(ylim = c(0, 1)) +  scale_x_continuous(breaks = c(0, 0.01,0.02,0.03,0.04,0.05)) + theme(text = element_text(size = 30,colour="black")) + scale_color_manual(values=c("#999999","#D55E00")) + theme(legend.position = "none") + labs(x="Absolute difference from EUR AF", y="True R2") 
ggplot(rsq_ga_ultra_rare, aes(x = abs(DIFF_AF_new), y = Imputation_R2 ,group=outlier, color=outlier))  + geom_point(shape=19,size=3)+ theme_classic(base_size=30) + coord_cartesian(ylim = c(0, 1)) +  scale_x_continuous(breaks = c(0, 0.01,0.02)) + theme(text = element_text(size = 30,colour="black")) + scale_color_manual(values=c("#999999","#D55E00")) + theme(legend.position = "none") + labs(x="Absolute difference from EUR AF", y="True R2") 
 
summary(lm(Imputation_R2 ~ abs(DIFF_AF_new) , data = rsq_meta_ultra_rare)) # p = 0.0146 *  
summary(lm(Imputation_R2 ~ abs(DIFF_AF_new) , data = rsq_topmed_ultra_rare)) # p =  0.0187 *  
summary(lm(Imputation_R2 ~ abs(DIFF_AF_new) , data = rsq_ex1kg_ultra_rare)) # p =   0.13
summary(lm(Imputation_R2 ~ abs(DIFF_AF_new), data = rsq_1kg38_ultra_rare)) # p =   0.091 
summary(lm(Imputation_R2 ~ abs(DIFF_AF_new) , data = rsq_1kg37_ultra_rare)) # p = 0.169
summary(lm(Imputation_R2 ~ abs(DIFF_AF_new) , data = rsq_ga_ultra_rare)) # p =  0.318 

rsq_meta_ultra_rare[which(abs(rsq_meta_ultra_rare$DIFF_AF_new) > 0.05),] # 2 SNPs 
rsq_topmed_ultra_rare[which(abs(rsq_topmed_ultra_rare$DIFF_AF_new) > 0.05),] # 2 SNPs
rsq_ex1kg_ultra_rare[which(abs(rsq_ex1kg_ultra_rare$DIFF_AF_new) > 0.05),] # 0 SNPs 
rsq_1kg38_ultra_rare[which(abs(rsq_1kg38_ultra_rare$DIFF_AF_new) > 0.05),] # 0 SNPs
rsq_1kg37_ultra_rare[which(abs(rsq_1kg37_ultra_rare$DIFF_AF_new) > 0.05),] # 1 SNPs
rsq_ga_ultra_rare[which(abs(rsq_ga_ultra_rare$DIFF_AF_new) > 0.05),] # 0 SNPs 


par(mfrow=c(3,2))
plot(abs(rsq_meta_ultra_rare$DIFF_AF_new) ,rsq_meta_ultra_rare$Imputation_R2,   pch=20, frame = FALSE)
plot(abs(rsq_topmed_ultra_rare$DIFF_AF_new)  ,rsq_topmed_ultra_rare$Imputation_R2,   pch=20, frame = FALSE)
plot(abs(rsq_ex1kg_ultra_rare$DIFF_AF_new)  ,rsq_ex1kg_ultra_rare$Imputation_R2,   pch=20, frame = FALSE)
plot(abs(rsq_1kg38_ultra_rare$DIFF_AF_new)  ,rsq_1kg38_ultra_rare$Imputation_R2,   pch=20, frame = FALSE)
plot(abs(rsq_1kg37_ultra_rare$DIFF_AF_new)  ,rsq_1kg37_ultra_rare$Imputation_R2,   pch=20, frame = FALSE) 
plot(abs(rsq_ga_ultra_rare$DIFF_AF_new)  ,rsq_ga_ultra_rare$Imputation_R2,   pch=20, frame = FALSE)

# check the 14 SNPs that have large AF difference between validation AF (sequencing) and imputation AF 

# topmed, expanded 1000G, 1000G GRCh38
grep -f ../snp_check_sequencing_imputation_grch38 all.chr.info

# meta-imputation
grep -f ../snp_check_sequencing_imputation_grch38  *brief.vcf

# 1000G GRCH37, GenomeAsia
grep -f ../snp_check_sequencing_imputation_grch37 all.chr.info


########################################################################
##########       Calculate MAF for the Tarseq data        ##############
########################################################################


# the original vcf from Oongjing is based on grch37
# therefore use Pakistani_1748_tarSeq.v4.2.updateID.sort.recode.vcf to calculate MAF
ml vcftools
vcftools --vcf Pakistani_1748_tarSeq.v4.2.updateID.sort.recode.vcf  --freq --out Pakistani_1748_tarSeq_grch37

sed 's/\:/\t/g' Pakistani_1748_tarSeq_grch37.frq |sort -g -k6|head
sed 's/\:/\t/g' Pakistani_1748_tarSeq_grch37.frq |sort -g -k6|tail

sed 's/\:/\t/g' Pakistani_1748_tarSeq_grch37.frq |sort -g -k8|head
sed 's/\:/\t/g' Pakistani_1748_tarSeq_grch37.frq |sort -g -k8|tail

# get the number of SNPs with MAF < 0.01
sed 's/\:/\t/g' Pakistani_1748_tarSeq_grch37.frq |awk '$6 < 0.01 || $8 < 0.01' | wc -l # 47338 SNPs with MAF < 0.01
sed 's/\:/\t/g' Pakistani_1748_tarSeq_grch37.frq |awk '$6 >= 0.01 && $8 >= 0.01' |wc -l # 1897 SNPs with MAF >= 0.01

# total number of SNPs
grep -v '^#' Pakistani_1748_tarSeq.v4.2.updateID.sort.recode.vcf |wc -l # 49235 SNPs

# try adding MAF to vcf (GRCH37)
ml bcftools
bcftools +fill-tags Pakistani_1748_tarSeq.v4.2.updateID.sort.recode.vcf --output Pakistani_1748_tarSeq.v4.2.updateID.sort.recode.maf.vcf -- -t AF 

grep -v '^#' Pakistani_1748_tarSeq.v4.2.updateID.sort.recode.vcf |wc -l # 49235 SNPs
grep -v '^#' Pakistani_1748_tarSeq.v4.2.updateID.sort.recode.maf.vcf |wc -l # 49235 SNPs

# check AF of 1 SNP: 9:141016456:A:G (AF of the second allele should be 0.000286205, and indeed it is correct)
grep 141016456  Pakistani_1748_tarSeq.v4.2.updateID.sort.recode.maf.vcf 


bcftools filter -i 'INFO/AF>=0.01 && INFO/AF<=0.99' Pakistani_1748_tarSeq.v4.2.updateID.sort.recode.maf.vcf  |grep -v '^#' | wc -l # 1897 SNPs
bcftools filter -i 'INFO/AF<0.01 || INFO/AF>0.99' Pakistani_1748_tarSeq.v4.2.updateID.sort.recode.maf.vcf  |grep -v '^#' | wc -l # 47338 SNPs

bcftools filter -i 'INFO/AF>=0.01 && INFO/AF<=0.99' Pakistani_1748_tarSeq.v4.2.updateID.sort.recode.maf.vcf  |grep -v '^#' |  awk '{print $3}' > tarseq_common_snp  # 1897 SNPs (no duplicated IDs)
bcftools filter -i 'INFO/AF<0.01 || INFO/AF>0.99' Pakistani_1748_tarSeq.v4.2.updateID.sort.recode.maf.vcf  |grep -v '^#' |  awk '{print $3}' > tarseq_rare_snp  # 47338 SNPs (no duplicated IDs)

sed 's/\:/\t/g' tarseq_common_snp |awk '{print $1":"$2":"$4":"$3}' > tarseq_common_snp_flip # 1897 SNPs
sed 's/\:/\t/g' tarseq_rare_snp |awk '{print $1":"$2":"$4":"$3}' > tarseq_rare_snp_flip # 47338 SNPs

cat tarseq_common_snp tarseq_common_snp_flip |sort  > tarseq_common_snp_double # 3794 SNPs
cat tarseq_rare_snp tarseq_rare_snp_flip |sort  > tarseq_rare_snp_double # 94676 SNPs


# try adding MAF to vcf (GRCH38)
bgzip -c Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.sort.recode.vcf > Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.sort.recode.vcf.gz 
tabix -s1 -b2 -e2 Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.sort.recode.vcf.gz 

ml bcftools
bcftools +fill-tags Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.sort.recode.vcf.gz --output Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.sort.recode.maf.vcf -- -t AF 

grep -v '^#' Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.sort.recode.vcf |wc -l # 49234 SNPs
grep -v '^#' Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.sort.recode.maf.vcf |wc -l # 49234 SNPs

bcftools filter -i 'INFO/AF>=0.01 && INFO/AF<=0.99' Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.sort.recode.maf.vcf |grep -v '^#' | wc -l # 1896 SNPs
bcftools filter -i 'INFO/AF<0.01 || INFO/AF>0.99' Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.sort.recode.maf.vcf |grep -v '^#' | wc -l # 47338 SNPs

bcftools filter -i 'INFO/AF>=0.01 && INFO/AF<=0.99' Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.sort.recode.maf.vcf  |grep -v '^#' |  awk '{print $1":"$2":"$4":"$5}' > tarseq_common_snp_grch38  # 1896 SNPs (no duplicated IDs)
bcftools filter -i 'INFO/AF<0.01 || INFO/AF>0.99' Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.sort.recode.maf.vcf  |grep -v '^#' |  awk '{print $1":"$2":"$4":"$5}' > tarseq_rare_snp_grch38  # 47338 SNPs (no duplicated IDs)

sed 's/\:/\t/g' tarseq_common_snp_grch38 |awk '{print $1":"$2":"$4":"$3}' > tarseq_common_snp_grch38_flip # 1896 SNPs
sed 's/\:/\t/g' tarseq_rare_snp_grch38 |awk '{print $1":"$2":"$4":"$3}' > tarseq_rare_snp_grch38_flip # 47338 SNPs

cat tarseq_common_snp_grch38 tarseq_common_snp_grch38_flip |sort  > tarseq_common_snp_grch38_double # 3792 SNPs
cat tarseq_rare_snp_grch38 tarseq_rare_snp_grch38_flip |sort  > tarseq_rare_snp_grch38_double # 94676 SNPs


## GenomeAsia ##

# find SNP overlap among these vcf files
all_snp_overlap_tarseq_ordered.sort.vcf 
all_snp_overlap_tarseq_ordered.highR2.sort.vcf

head /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/tarseq_common_snp
head /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/tarseq_rare_snp

grep -v '^#' all_snp_overlap_tarseq_ordered.sort.vcf  |awk '{print $3}' |sort >  all_snp_overlap_tarseq_ordered.sort.vcf.snplist  # 5394 SNPS
grep -v '^#' all_snp_overlap_tarseq_ordered.highR2.sort.vcf |awk '{print $3}' |sort >  all_snp_overlap_tarseq_ordered.highR2.sort.vcf.snplist  # 1445 SNPS

comm -12 <(sort all_snp_overlap_tarseq_ordered.sort.vcf.snplist) <(sort /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/tarseq_common_snp_double) |wc -l # 1724 common SNP overlap (in ALL SNPs)
comm -12 <(sort all_snp_overlap_tarseq_ordered.sort.vcf.snplist) <(sort /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/tarseq_rare_snp_double) |wc -l  # 3670 rare SNP overlap (in ALL SNPs)

comm -12 <(sort all_snp_overlap_tarseq_ordered.highR2.sort.vcf.snplist) <(sort /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/tarseq_common_snp_double) |wc -l  # 1162 common SNP overlap (in high quality SNPs)
comm -12 <(sort all_snp_overlap_tarseq_ordered.highR2.sort.vcf.snplist) <(sort /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/tarseq_rare_snp_double) |wc -l  # 283 rare SNP overlap (in high quality SNPs)

## 1000G GRCh 37 ####
grep -v '^#' all_snp_overlap_tarseq_ordered.sort.vcf  |awk '{print $3}' |sort >  all_snp_overlap_tarseq_ordered.sort.vcf.snplist  # 11514 SNPS
grep -v '^#' all_snp_overlap_tarseq_ordered.highR2.sort.vcf |awk '{print $3}' |sort >  all_snp_overlap_tarseq_ordered.highR2.sort.vcf.snplist  # 1843 SNPS

comm -12 <(sort all_snp_overlap_tarseq_ordered.sort.vcf.snplist) <(sort /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/tarseq_common_snp_double) |wc -l # 1885 common SNP overlap (in ALL SNPs)
comm -12 <(sort all_snp_overlap_tarseq_ordered.sort.vcf.snplist) <(sort /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/tarseq_rare_snp_double) |wc -l  # 9629 rare SNP overlap (in ALL SNPs)

comm -12 <(sort all_snp_overlap_tarseq_ordered.highR2.sort.vcf.snplist) <(sort /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/tarseq_common_snp_double) |wc -l  # 1344 common SNP overlap (in high quality SNPs)
comm -12 <(sort all_snp_overlap_tarseq_ordered.highR2.sort.vcf.snplist) <(sort /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/tarseq_rare_snp_double) |wc -l  # 499 rare SNP overlap (in high quality SNPs)

## 1000G GRCh 38 ####
grep -v '^#' all_snp_overlap_tarseq_ordered.sort.vcf  |awk '{print $3}' |sort >  all_snp_overlap_tarseq_ordered.sort.vcf.snplist  # 11172 SNPS
grep -v '^#' all_snp_overlap_tarseq_ordered.highR2.sort.vcf |awk '{print $3}' |sort >  all_snp_overlap_tarseq_ordered.highR2.sort.vcf.snplist  # 2325 SNPS

comm -12 <(sort all_snp_overlap_tarseq_ordered.sort.vcf.snplist) <(sort /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/tarseq_common_snp_grch38_double) |wc -l # 1815 common SNP overlap (in ALL SNPs)
comm -12 <(sort all_snp_overlap_tarseq_ordered.sort.vcf.snplist) <(sort /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/tarseq_rare_snp_grch38_double) |wc -l  # 9357 rare SNP overlap (in ALL SNPs)

comm -12 <(sort all_snp_overlap_tarseq_ordered.highR2.sort.vcf.snplist) <(sort /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/tarseq_common_snp_grch38_double) |wc -l  # 1520 common SNP overlap (in high quality SNPs)
comm -12 <(sort all_snp_overlap_tarseq_ordered.highR2.sort.vcf.snplist) <(sort /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/tarseq_rare_snp_grch38_double) |wc -l  # 805 rare SNP overlap (in high quality SNPs)

## Expanded 1000G  ####
grep -v '^#' all_snp_overlap_tarseq_ordered.sort.vcf  |awk '{print $3}' |sort >  all_snp_overlap_tarseq_ordered.sort.vcf.snplist  # 11562 SNPS
grep -v '^#' all_snp_overlap_tarseq_ordered.highR2.sort.vcf |awk '{print $3}' |sort >  all_snp_overlap_tarseq_ordered.highR2.sort.vcf.snplist  # 2444 SNPS

comm -12 <(sort all_snp_overlap_tarseq_ordered.sort.vcf.snplist) <(sort /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/tarseq_common_snp_grch38_double) |wc -l # 1856 common SNP overlap (in ALL SNPs)
comm -12 <(sort all_snp_overlap_tarseq_ordered.sort.vcf.snplist) <(sort /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/tarseq_rare_snp_grch38_double) |wc -l  # 9706 rare SNP overlap (in ALL SNPs)

comm -12 <(sort all_snp_overlap_tarseq_ordered.highR2.sort.vcf.snplist) <(sort /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/tarseq_common_snp_grch38_double) |wc -l  # 1578 common SNP overlap (in high quality SNPs)
comm -12 <(sort all_snp_overlap_tarseq_ordered.highR2.sort.vcf.snplist) <(sort /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/tarseq_rare_snp_grch38_double) |wc -l  # 866 rare SNP overlap (in high quality SNPs)

## TOPMed  ####
grep -v '^#' all_snp_overlap_tarseq_ordered.sort.vcf  |awk '{print $3}' |sort >  all_snp_overlap_tarseq_ordered.sort.vcf.snplist  # 48455 SNPS
grep -v '^#' all_snp_overlap_tarseq_ordered.highR2.sort.vcf |awk '{print $3}' |sort >  all_snp_overlap_tarseq_ordered.highR2.sort.vcf.snplist  # 2735 SNPS

comm -12 <(sort all_snp_overlap_tarseq_ordered.sort.vcf.snplist) <(sort /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/tarseq_common_snp_grch38_double) |wc -l # 1871 common SNP overlap (in ALL SNPs)
comm -12 <(sort all_snp_overlap_tarseq_ordered.sort.vcf.snplist) <(sort /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/tarseq_rare_snp_grch38_double) |wc -l  # 46584 rare SNP overlap (in ALL SNPs)

comm -12 <(sort all_snp_overlap_tarseq_ordered.highR2.sort.vcf.snplist) <(sort /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/tarseq_common_snp_grch38_double) |wc -l  # 1507 common SNP overlap (in high quality SNPs)
comm -12 <(sort all_snp_overlap_tarseq_ordered.highR2.sort.vcf.snplist) <(sort /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/tarseq_rare_snp_grch38_double) |wc -l  # 1228 rare SNP overlap (in high quality SNPs)

## Meta-imputation  ####
grep -v '^#' all_snp_overlap_tarseq_ordered.sort.vcf  |awk '{print $3}' |sort >  all_snp_overlap_tarseq_ordered.sort.vcf.snplist  # 48907 SNPS
grep -v '^#' all_snp_overlap_tarseq_ordered.highR2.sort.vcf |awk '{print $3}' |sort >  all_snp_overlap_tarseq_ordered.highR2.sort.vcf.snplist  # 2968 SNPS

comm -12 <(sort all_snp_overlap_tarseq_ordered.sort.vcf.snplist) <(sort /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/tarseq_common_snp_grch38_double) |wc -l # 1873 common SNP overlap (in ALL SNPs)
comm -12 <(sort all_snp_overlap_tarseq_ordered.sort.vcf.snplist) <(sort /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/tarseq_rare_snp_grch38_double) |wc -l  # 47034 rare SNP overlap (in ALL SNPs)

comm -12 <(sort all_snp_overlap_tarseq_ordered.highR2.sort.vcf.snplist) <(sort /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/tarseq_common_snp_grch38_double) |wc -l  # 1589 common SNP overlap (in high quality SNPs)
comm -12 <(sort all_snp_overlap_tarseq_ordered.highR2.sort.vcf.snplist) <(sort /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/tarseq_rare_snp_grch38_double) |wc -l  # 1379 rare SNP overlap (in high quality SNPs)





##################################################
########      plot the true Rsq by MAF   #########
##################################################

#####      ALL SNPs   ##### 

# manually remove the # in the header line
# sed '/^#/d' tarseq_metaimp_chr1_test.aggRSquare > tarseq_metaimp_chr1_test.aggRSquare.noheader

# in the path: /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_compare/aggRSquare
# sentence starting with # will not be read by R
metaimp <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_topmed_1000G/tarseq_meta_imputation.aggRSquare",header=F)
ex1000G <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/tarseq_ex1000G_grch38.aggRSquare",header=F)
topmed <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/tarseq_topmed_grch38.aggRSquare",header=F)
GA <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_GA/tarseq_GA_grch37.aggRSquare",header=F)
all_1000G_38 <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_38/tarseq_1000G_OG_grch38.aggRSquare",header=F)
sas_1000G_37 <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_37/tarseq_1000G_OG_grch37.aggRSquare",header=F)

metaimp
topmed
ex1000G
all_1000G_38
sas_1000G_37
GA
# update the SNP ID to the same as tarseq ID (check forward and reverse)

metaimp$ref_panel <- "TOPMed_1000G"
topmed$ref_panel <- "TOPMed"
ex1000G$ref_panel <- "Expanded_1000G"
all_1000G_38$ref_panel <- "1000G_GRCh38"
sas_1000G_37$ref_panel <- "1000G_GRCh37_SAS"
GA$ref_panel <- "GenomeAsia_Pilot"

merged <- rbind(metaimp,topmed,ex1000G,all_1000G_38,sas_1000G_37,GA) #  84 rows, 7 columns
merged$ref_panel <- factor(merged$ref_panel,  levels=c("TOPMed_1000G","TOPMed","Expanded_1000G","1000G_GRCh38", "1000G_GRCh37_SAS", "GenomeAsia_Pilot"))
colnames(merged) <- c("Bin_Aggregated_by_MAF",	"Average_MAF",	"Variant_count"	,"Imputation_R2",	"Gold_MAF","Imputed_MAF"	, "ref_panel")

library(dplyr)
merged$MAF_bin <- case_when(merged$Bin_Aggregated_by_MAF =="(0.000000,0.000500]" ~ "(0,0.0005]", merged$Bin_Aggregated_by_MAF =="(0.000500,0.001000]" ~ "(0.0005,0.001]", merged$Bin_Aggregated_by_MAF =="(0.001000,0.002000]" ~ "(0.001,0.002]", merged$Bin_Aggregated_by_MAF =="(0.002000,0.005000]" ~ "(0.002,0.005]",merged$Bin_Aggregated_by_MAF =="(0.005000,0.010000]" ~ "(0.005,0.01]",merged$Bin_Aggregated_by_MAF =="(0.010000,0.015000]" ~ "(0.01,0.015]",merged$Bin_Aggregated_by_MAF =="(0.015000,0.020000]" ~ "(0.015,0.02]",merged$Bin_Aggregated_by_MAF =="(0.020000,0.035000]" ~ "(0.02,0.035]",merged$Bin_Aggregated_by_MAF =="(0.035000,0.050000]" ~ "(0.035,0.05]",merged$Bin_Aggregated_by_MAF =="(0.050000,0.100000]" ~ "(0.05,0.1]",merged$Bin_Aggregated_by_MAF =="(0.100000,0.200000]" ~ "(0.1,0.2]",merged$Bin_Aggregated_by_MAF =="(0.200000,0.300000]" ~ "(0.2,0.3]",merged$Bin_Aggregated_by_MAF =="(0.300000,0.400000]" ~ "(0.3,0.4]",merged$Bin_Aggregated_by_MAF =="(0.400000,0.500000]" ~ "(0.4,0.5]")

# color blind friendly
colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
scales::show_col(colorBlindGrey8)

library(ggplot2)
library(dplyr)
p <- ggplot(merged, aes(x = MAF_bin, y = Imputation_R2,group=ref_panel, color=ref_panel))  + geom_line() + geom_point(shape=1)+ theme_classic(base_size=20) + coord_cartesian(ylim = c(0, 1)) + scale_x_discrete(guide = guide_axis(angle = 45))  + theme(text = element_text(size = 15,colour="black"))  + geom_hline(yintercept=0.8,  color = "red") + geom_hline(yintercept=0.6, linetype="dashed", color = "red")  + geom_vline(xintercept="(0.01,0.015]",  color = "black") + geom_vline(xintercept="(0.05,0.1]",  color = "black", linetype="dashed") #all font size
p + scale_color_manual(values=c( "#0072B2","#E68F00", "#56B4E9","#CC79A7","#999999","#F0E442"),name="Reference Panel") + labs(x="Minor Allele Frequency", y="Average True R2")


#####      High Quality SNPS (R2 >= 0.8)   ##### 
metaimp <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_topmed_1000G/tarseq_meta_imputation_highR2.aggRSquare",header=F)
ex1000G <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/tarseq_ex1000G_grch38_highR2.aggRSquare",header=F)
topmed <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/tarseq_topmed_grch38_highR2.aggRSquare",header=F)
GA <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_GA/tarseq_GA_grch37_highR2.aggRSquare",header=F)
all_1000G_38 <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_38/tarseq_1000G_OG_grch38_highR2.aggRSquare",header=F)
sas_1000G_37 <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_37/tarseq_1000G_OG_grch37_highR2.aggRSquare",header=F)

# update the SNP ID to the same as tarseq ID (check forward and reverse)

metaimp$ref_panel <- "TOPMed_1000G"
topmed$ref_panel <- "TOPMed"
ex1000G$ref_panel <- "Expanded_1000G"
all_1000G_38$ref_panel <- "1000G_GRCh38"
sas_1000G_37$ref_panel <- "1000G_GRCh37_SAS"
GA$ref_panel <- "GenomeAsia_Pilot"

merged <- rbind(metaimp,topmed,ex1000G,all_1000G_38,sas_1000G_37,GA) #  84 rows, 7 columns
merged$ref_panel <- factor(merged$ref_panel,  levels=c("TOPMed_1000G","TOPMed","Expanded_1000G","1000G_GRCh38", "1000G_GRCh37_SAS", "GenomeAsia_Pilot"))
colnames(merged) <- c("Bin_Aggregated_by_MAF",	"Average_MAF",	"Variant_count"	,"Imputation_R2",	"Gold_MAF","Imputed_MAF"	, "ref_panel")

library(dplyr)
merged$MAF_bin <- case_when(merged$Bin_Aggregated_by_MAF =="(0.000000,0.000500]" ~ "(0,0.0005]", merged$Bin_Aggregated_by_MAF =="(0.000500,0.001000]" ~ "(0.0005,0.001]", merged$Bin_Aggregated_by_MAF =="(0.001000,0.002000]" ~ "(0.001,0.002]", merged$Bin_Aggregated_by_MAF =="(0.002000,0.005000]" ~ "(0.002,0.005]",merged$Bin_Aggregated_by_MAF =="(0.005000,0.010000]" ~ "(0.005,0.01]",merged$Bin_Aggregated_by_MAF =="(0.010000,0.015000]" ~ "(0.01,0.015]",merged$Bin_Aggregated_by_MAF =="(0.015000,0.020000]" ~ "(0.015,0.02]",merged$Bin_Aggregated_by_MAF =="(0.020000,0.035000]" ~ "(0.02,0.035]",merged$Bin_Aggregated_by_MAF =="(0.035000,0.050000]" ~ "(0.035,0.05]",merged$Bin_Aggregated_by_MAF =="(0.050000,0.100000]" ~ "(0.05,0.1]",merged$Bin_Aggregated_by_MAF =="(0.100000,0.200000]" ~ "(0.1,0.2]",merged$Bin_Aggregated_by_MAF =="(0.200000,0.300000]" ~ "(0.2,0.3]",merged$Bin_Aggregated_by_MAF =="(0.300000,0.400000]" ~ "(0.3,0.4]",merged$Bin_Aggregated_by_MAF =="(0.400000,0.500000]" ~ "(0.4,0.5]")

# color blind friendly
colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
scales::show_col(colorBlindGrey8)

library(ggplot2)
library(dplyr)
p <- ggplot(merged, aes(x = MAF_bin, y = Imputation_R2,group=ref_panel, color=ref_panel))  + geom_line(size=1) + geom_point(shape=1)+ theme_classic(base_size=20) + coord_cartesian(ylim = c(0, 1)) + scale_x_discrete(guide = guide_axis(angle = 45))  + theme(text = element_text(size = 15,colour="black"))  + geom_hline(yintercept=0.8,  color = "red") + geom_hline(yintercept=0.6, linetype="dashed", color = "red")  + geom_vline(xintercept="(0.01,0.015]",  color = "black") + geom_vline(xintercept="(0.05,0.1]",  color = "black", linetype="dashed") #all font size
p + scale_color_manual(values=c("#661100", "#D55E00", "#E69F00","#0072B2","#332288","#999999"),name="Reference Panel") + labs(x="Minor Allele Frequency", y="Average True R2")


# what about average of individual true Rsq?

# All SNPs
metaimp <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_topmed_1000G/tarseq_meta_imputation.RSquare",header=F) # 8483 SNPs
topmed <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/tarseq_topmed_grch38.RSquare",header=F) # 8267 SNPs
ex1000G <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/tarseq_ex1000G_grch38.RSquare",header=F) # 3469 SNPs
all_1000G_38 <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_38/tarseq_1000G_OG_grch38.RSquare",header=F) # 3361 SNPs
sas_1000G_37 <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_37/tarseq_1000G_OG_grch37.RSquare",header=F) # 3467 SNPs
GA <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_GA/tarseq_GA_grch37.RSquare",header=F) # 2126 SNPs

metaimp$ref_panel <- "TOPMed_1000G"
topmed$ref_panel <- "TOPMed"
ex1000G$ref_panel <- "Expanded_1000G"
all_1000G_38$ref_panel <- "1000G_GRCh38"
sas_1000G_37$ref_panel <- "1000G_GRCh37_SAS"
GA$ref_panel <- "GenomeAsia_Pilot"

merged <- rbind(metaimp,topmed,ex1000G,all_1000G_38,sas_1000G_37,GA) #  29173 rows 
colnames(merged) <- c("SNP_ID",	"Allele_Frequency",	"No_Samples"	,"Imputation_R2",	"Validation_AF","Imputation_AF"	, "ref_panel")

for (i in 1:dim(merged)[1]) {
 if (merged$Validation_AF[i] <= 0.5) {
    merged$MAF[i] <- merged$Validation_AF[i] 
    }
 else if (merged$Validation_AF[i] > 0.5) {
    merged$MAF[i] <- 1- merged$Validation_AF[i] 
    }
}

merged$MAF <- as.numeric(merged$MAF)
summary(merged$MAF)

library(dplyr)
merged$MAF_bin <- case_when(merged$MAF > 0 & merged$MAF <= 0.0005 ~ "(0,0.0005]", merged$MAF > 0.0005 & merged$MAF <= 0.001 ~ "(0.0005,0.001]", merged$MAF > 0.001 & merged$MAF <= 0.002 ~ "(0.001,0.002]", merged$MAF > 0.002 & merged$MAF <= 0.005 ~ "(0.002,0.005]",merged$MAF > 0.005 & merged$MAF <= 0.01 ~ "(0.005,0.01]",merged$MAF > 0.01 & merged$MAF <= 0.015 ~ "(0.01,0.015]",merged$MAF > 0.015 & merged$MAF <= 0.02 ~ "(0.015,0.02]", merged$MAF > 0.02 & merged$MAF <= 0.035 ~ "(0.02,0.035]", merged$MAF > 0.035 & merged$MAF <= 0.05 ~ "(0.035,0.05]", merged$MAF > 0.05 & merged$MAF <= 0.1 ~ "(0.05,0.1]", merged$MAF > 0.1 & merged$MAF <= 0.2  ~ "(0.1,0.2]", merged$MAF > 0.2 & merged$MAF <= 0.3  ~ "(0.2,0.3]", merged$MAF > 0.3 & merged$MAF <= 0.4  ~ "(0.3,0.4]", merged$MAF > 0.4 & merged$MAF <= 0.5  ~ "(0.4,0.5]")
table(merged$MAF_bin)

# install.packages("plotrix")
# library(plotrix)
true_r2_maf <- aggregate(Imputation_R2~ref_panel+MAF_bin, data=merged, mean)
# FUN = function(x) c(mean = mean(x), se = std.error(x))
true_r2_maf$ref_panel <- factor(true_r2_maf$ref_panel,  levels=c("TOPMed_1000G","TOPMed","Expanded_1000G","1000G_GRCh38", "1000G_GRCh37_SAS", "GenomeAsia_Pilot"))
true_r2_maf$MAF_bin <- factor(true_r2_maf$MAF_bin,  levels=c("(0,0.0005]","(0.0005,0.001]","(0.001,0.002]","(0.002,0.005]","(0.005,0.01]","(0.01,0.015]","(0.015,0.02]","(0.02,0.035]","(0.035,0.05]","(0.05,0.1]","(0.1,0.2]","(0.2,0.3]","(0.3,0.4]","(0.4,0.5]"))

library(ggplot2)
library(dplyr)
pdf('true_r2_maf_all_snp.pdf',height=10,width=15)
p <- ggplot(true_r2_maf, aes(x = MAF_bin, y = Imputation_R2,group=ref_panel, color=ref_panel))  + geom_line(size=1) + geom_point(shape=1)+ theme_classic(base_size=30) + coord_cartesian(ylim = c(0, 1)) + scale_x_discrete(guide = guide_axis(angle = 45))  + theme(text = element_text(size = 30,colour="black"))  + geom_hline(yintercept=0.8,  color = "dark grey") + geom_hline(yintercept=0.6, linetype="dashed", color = "dark grey")  + geom_vline(xintercept="(0.01,0.015]",  color = "dark grey",linetype="dashed") + geom_vline(xintercept="(0.05,0.1]",  color = "dark grey") #all font size
p + scale_color_manual(values=c("#661100", "#D55E00", "#E69F00","#0072B2","#332288","#999999"),name="Reference Panel") + labs(x="Minor Allele Frequency", y="Average True R2") + 
  annotate('text', x = "(0,0.0005]", y = 0.45, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.0005,0.001]", y = 0.4, label='"**"', parse=TRUE, color = "red",size=10) +
  annotate('text', x = "(0.001,0.002]", y = 0.45, label='"**"', parse=TRUE, color = "red",size=10) +
  annotate('text', x = "(0.002,0.005]", y = 0.6, label='"*"', parse=TRUE, color = "red",size=10) +
  annotate('text', x = "(0.005,0.01]", y = 0.7, label='"*"', parse=TRUE, color = "red",size=10) +
  annotate('text', x = "(0.01,0.015]", y = 0.75, label='"*"', parse=TRUE, color = "red",size=10) +
  annotate('text', x = "(0.4,0.5]", y = 0.95, label='"*"', parse=TRUE, color = "red",size=10)
dev.off()
# one way ANOVA
merged_1 <- merged[which(merged$MAF_bin == "(0,0.0005]"),] # 7631 rows
merged_2 <- merged[which(merged$MAF_bin == "(0.0005,0.001]"),] # 4654 rows
merged_3 <- merged[which(merged$MAF_bin == "(0.001,0.002]"),] # 2749 rows
merged_4 <- merged[which(merged$MAF_bin == "(0.002,0.005]"),] # 3494 rows
merged_5 <- merged[which(merged$MAF_bin == "(0.005,0.01]"),] # 2020 rows
merged_6 <- merged[which(merged$MAF_bin == "(0.01,0.015]"),] # 902 rows
merged_7 <- merged[which(merged$MAF_bin == "(0.015,0.02]"),] # 464 rows
merged_8 <- merged[which(merged$MAF_bin == "(0.02,0.035]"),] # 1029 rows
merged_9 <- merged[which(merged$MAF_bin == "(0.035,0.05]"),] # 597 rows
merged_10 <- merged[which(merged$MAF_bin == "(0.05,0.1]"),] # 1271 rows
merged_11 <- merged[which(merged$MAF_bin == "(0.1,0.2]"),] # 1488 rows
merged_12 <- merged[which(merged$MAF_bin == "(0.2,0.3]"),] # 1167 rows
merged_13 <- merged[which(merged$MAF_bin == "(0.3,0.4]"),] # 1007 rows
merged_14 <- merged[which(merged$MAF_bin == "(0.4,0.5]"),] # 700 rows

table(merged$MAF_bin)
# testing if the aggregate function calculate the average true R2 correctly: YES! 
table(merged_7$ref_panel)
test <- merged_7[which(merged_7$ref_panel=="GenomeAsia_Pilot"),]
mean(test$Imputation_R2) # the average R2 for GenomeAsia in the MAF bin of (0.015,0.02] is 0.6271125 (matched with the value in true_r2_maf)

# Non-parametric alternative to one-way ANOVA test
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_1) # P  < 2.2e-16 (**) (probably due to topmed has way many low-quality rare variants imputed, therefore drag down the average R2)
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_2) # P  < 2.2e-16 (**)
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_3) # P  0.0008228 (**)
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_4) # P  0.00391 (*)
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_5) # P   0.02918 (*)
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_6) # P  0.00701 (*)
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_7) # P  0.3296
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_8) # P  0.5176
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_9) # P  0.7334
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_10) # P  0.8644
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_11) # P  0.1805
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_12) # P  0.7145
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_13) # P 0.3076
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_14) # P 0.02616 (*)
# 0.05/14 = 0.003571 (** Bonferroni significant )
# 0.05 (* nominal significant)

# http://www.sthda.com/english/wiki/one-way-anova-test-in-r

# check for pairwise difference
# Pairwise t-tests with no assumption of equal variances
pairwise.t.test(merged_1$Imputation_R2, merged_1$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) # GA > ex1KG ~ 1KG38 ~ 1KG37 > meta ~ TOPMed
# genome_asia better than any other panels
# GenomeAsia significantly the best
# meta-imputation and topmed are significantly the worst

#                 1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     0.9939           -            -              -               
# Expanded_1000G   0.5703           0.5777       -              -               
# GenomeAsia_Pilot 0.0060           0.0062       0.0182         -               
# TOPMed           9.1e-10          1.2e-09      6.0e-12        6.0e-11         
# TOPMed_1000G     1.8e-08          2.3e-08      1.7e-10        3.5e-10         
#                  TOPMed
# 1000G_GRCh38     -     
# Expanded_1000G   -     
# GenomeAsia_Pilot -     
# TOPMed           -     
# TOPMed_1000G     0.3369

pairwise.t.test(merged_2$Imputation_R2, merged_2$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) # ex1KG ~ GA ~ 1KG38 ~ 1KG37 > meta ~ TOPMed
# expanded 1000G better than TOPMed and meta-imputation
# only meta-imputation and topmed are significantly the worst

#                 1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     0.95648          -            -              -               
# Expanded_1000G   0.20201          0.23092      -              -               
# GenomeAsia_Pilot 0.60187          0.63657      0.62802        -               
# TOPMed           0.00072          0.00080      1e-06          0.00277         
# TOPMed_1000G     0.02753          0.02726      0.00018        0.03098         
#                 TOPMed 
# 1000G_GRCh38     -      
# Expanded_1000G   -      
# GenomeAsia_Pilot -      
# TOPMed           -      
# TOPMed_1000G     0.08319

pairwise.t.test(merged_3$Imputation_R2, merged_3$ref_panel,  p.adjust.method = "none", pool.sd = FALSE)  # ex1KG ~ meta ~ 1KG38 > 1KG37 ~ TOPMed ~ GA
# expanded 1000G better than TOPMed and genome_asia
# only meta-imputation and expanded 1000G are significantly better
# GenomeAsia significantly the worst

#                 1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     0.4687           -            -              -               
# Expanded_1000G   0.0507           0.2161       -              -               
# GenomeAsia_Pilot 0.1012           0.0281       0.0018         -               
# TOPMed           0.7813           0.2593       0.0103         0.1121          
# TOPMed_1000G     0.1057           0.4402       0.4935         0.0034          
#                  TOPMed
# 1000G_GRCh38     -     
# Expanded_1000G   -     
# GenomeAsia_Pilot -     
# TOPMed           -     
# TOPMed_1000G     0.0163

pairwise.t.test(merged_4$Imputation_R2, merged_4$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) #  ex1KG ~ meta > 1KG37 ~ 1KG38 ~ GA ~ TOPMed
# expanded 1000G better than TOPMed and 1000G GRCh38 and 1000G SAS GRCh37
# only meta-imputation and expanded 1000G are significantly better

#                  1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     0.97538          -            -              -               
# Expanded_1000G   0.03715          0.03441      -              -               
# GenomeAsia_Pilot 0.84275          0.86326      0.05097        -               
# TOPMed           0.30612          0.32297      0.00063        0.53571         
# TOPMed_1000G     0.08334          0.07730      0.53250        0.10513         
#                  TOPMed 
# 1000G_GRCh38     -      
# Expanded_1000G   -      
# GenomeAsia_Pilot -      
# TOPMed           -      
# TOPMed_1000G     0.00107

pairwise.t.test(merged_5$Imputation_R2, merged_5$ref_panel,  p.adjust.method = "none", pool.sd = FALSE)  # meta-imputation ~ ex1KG > 1KG37 ~ 1KG38 ~ GA ~ TOPMed 
# meta-imputation better than TOPMed, 1000G GRCh38, 1000G SAS GRCh37, GenomeAsia
# only meta-imputation and expanded 1000G are significantly better

#                  1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     0.73637          -            -              -               
# Expanded_1000G   0.27043          0.15357      -              -               
# GenomeAsia_Pilot 0.53192          0.74942      0.10923        -               
# TOPMed           0.42828          0.68187      0.04319        0.98875         
# TOPMed_1000G     0.03351          0.01384      0.35298        0.01241         
#                  TOPMed 
# 1000G_GRCh38     -      
# Expanded_1000G   -      
# GenomeAsia_Pilot -      
# TOPMed           -      
# TOPMed_1000G     0.00091

pairwise.t.test(merged_6$Imputation_R2, merged_6$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) # meta-imputation > TOPMed (mean: 0.6588, p=0.049) ~ ex1KG (mean: 0.6585, p=0.071, possibly fewer SNPs therefore larger P) ~ 1KG38 ~ 1KG37 ~ GA 
# meta-imputation better than TOPMed, 1000G GRCh38, 1000G SAS GRCh37, GenomeAsia
# expanded 1000G as good as meta-imputation

#                 1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     0.6633           -            -              -               
# Expanded_1000G   0.1994           0.4104       -              -               
# GenomeAsia_Pilot 0.9797           0.6612       0.2156         -               
# TOPMed           0.1670           0.3741       0.9913         0.1858          
# TOPMed_1000G     0.0017           0.0088       0.0708         0.0031          
#                  TOPMed
# 1000G_GRCh38     -     
# Expanded_1000G   -     
# GenomeAsia_Pilot -     
# TOPMed           -     
# TOPMed_1000G     0.0489

pairwise.t.test(merged_7$Imputation_R2, merged_7$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) # meta ~ TOPMed ~ Ex1KG ~ 1K38 ~ GA > 1K37

#                  1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     0.905            -            -              -               
# Expanded_1000G   0.659            0.758        -              -               
# GenomeAsia_Pilot 0.974            0.935        0.699          -               
# TOPMed           0.237            0.318        0.502          0.285           
# TOPMed_1000G     0.044            0.073        0.134          0.065           
#                  TOPMed
# 1000G_GRCh38     -     
# Expanded_1000G   -     
# GenomeAsia_Pilot -     
# TOPMed           -     
# TOPMed_1000G     0.312 

pairwise.t.test(merged_8$Imputation_R2, merged_8$ref_panel,  p.adjust.method = "none", pool.sd = FALSE)  # no sig pair, meta ~ ex1KG ~ 1G37 ~ TOPMed ~ 1K38 ~ GA

#                  1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     0.971            -            -              -               
# Expanded_1000G   0.595            0.575        -              -               
# GenomeAsia_Pilot 0.856            0.884        0.493          -               
# TOPMed           0.970            1.000        0.548          0.876           
# TOPMed_1000G     0.131            0.127        0.351          0.106           
#                  TOPMed
# 1000G_GRCh38     -     
# Expanded_1000G   -     
# GenomeAsia_Pilot -     
# TOPMed           -     
# TOPMed_1000G     0.097 

pairwise.t.test(merged_9$Imputation_R2, merged_9$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) # no sig pair, meta ~ TOPMed ~ 1K38 ~ ex1KG ~ GA ~ 1KG37

#                  1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     0.68             -            -              -               
# Expanded_1000G   0.74             0.94         -              -               
# GenomeAsia_Pilot 0.90             0.79         0.84           -               
# TOPMed           0.41             0.71         0.65           0.53            
# TOPMed_1000G     0.12             0.26         0.23           0.20            
#                  TOPMed
# 1000G_GRCh38     -     
# Expanded_1000G   -     
# GenomeAsia_Pilot -     
# TOPMed           -     
# TOPMed_1000G     0.39 

pairwise.t.test(merged_10$Imputation_R2, merged_10$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) # no sig pair, ex1KG ~ meta ~ GA ~ 1K38 ~ 1K37 ~ TOPMed

#                  1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     0.90             -            -              -               
# Expanded_1000G   0.69             0.78         -              -               
# GenomeAsia_Pilot 0.76             0.86         0.92           -               
# TOPMed           0.59             0.52         0.34           0.40            
# TOPMed_1000G     0.75             0.85         0.92           1.00            
#                  TOPMed
# 1000G_GRCh38     -     
# Expanded_1000G   -     
# GenomeAsia_Pilot -     
# TOPMed           -     
# TOPMed_1000G     0.36 

pairwise.t.test(merged_11$Imputation_R2, merged_11$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) # no sig pair, meta ~ ex1KG ~ TOPMed ~ GA ~ 1KG38 ~ 1KG37

#                  1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     0.68             -            -              -               
# Expanded_1000G   0.33             0.58         -              -               
# GenomeAsia_Pilot 0.62             0.94         0.64           -               
# TOPMed           0.42             0.71         0.84           0.77            
# TOPMed_1000G     0.11             0.24         0.53           0.26            
#                  TOPMed
# 1000G_GRCh38     -     
# Expanded_1000G   -     
# GenomeAsia_Pilot -     
# TOPMed           -     
# TOPMed_1000G     0.37  

pairwise.t.test(merged_12$Imputation_R2, merged_12$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) # no sig pair, meta ~ GA ~ ex1KG ~ TOPMed ~ 1KG38 ~ 1KG37

#                  1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     0.79             -            -              -               
# Expanded_1000G   0.39             0.57         -              -               
# GenomeAsia_Pilot 0.24             0.37         0.76           -               
# TOPMed           0.44             0.64         0.85           0.59            
# TOPMed_1000G     0.17             0.28         0.63           0.87            
#                  TOPMed
# 1000G_GRCh38     -     
# Expanded_1000G   -     
# GenomeAsia_Pilot -     
# TOPMed           -     
# TOPMed_1000G     0.44  

pairwise.t.test(merged_13$Imputation_R2, merged_13$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) # no sig pair, meta ~ TOPMed ~ ex1KG ~ 1K37 ~ 1K38 ~ Pilot

#                  1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     0.91             -            -              -               
# Expanded_1000G   0.82             0.74         -              -               
# GenomeAsia_Pilot 0.52             0.60         0.40           -               
# TOPMed           0.61             0.54         0.79           0.25            
# TOPMed_1000G     0.30             0.26         0.43           0.10            
#                  TOPMed
# 1000G_GRCh38     -     
# Expanded_1000G   -     
# GenomeAsia_Pilot -     
# TOPMed           -     
# TOPMed_1000G     0.55  

pairwise.t.test(merged_14$Imputation_R2, merged_14$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) # meta-imputation ~ TOPMed ~ ex1KG ~ 1KG38 ~ 1KG37 > GA
# expanded 1000G as good as meta-imputation

#                  1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     0.955            -            -              -               
# Expanded_1000G   0.768            0.814        -              -               
# GenomeAsia_Pilot 0.694            0.659        0.501          -               
# TOPMed           0.186            0.216        0.313          0.101           
# TOPMed_1000G     0.060            0.074        0.114          0.032           
#                  TOPMed
# 1000G_GRCh38     -     
# Expanded_1000G   -     
# GenomeAsia_Pilot -     
# TOPMed           -     
# TOPMed_1000G     0.491 


# testing anova
aov_1 <- aov(Imputation_R2 ~ ref_panel, data = merged_1)
library(car)
leveneTest(Imputation_R2 ~ ref_panel, data = merged_1) # unequal variance across groups
plot(aov_1, 2) # unlike normal distribution 
# therefore using non parametric test

# check distribution of R2 in (0,0.0005)

par(mfrow=c(3,2))
hist(merged_1$Imputation_R2[which(merged_1$ref_panel=="TOPMed_1000G")])
hist(merged_1$Imputation_R2[which(merged_1$ref_panel=="TOPMed")])
hist(merged_1$Imputation_R2[which(merged_1$ref_panel=="Expanded_1000G")])
hist(merged_1$Imputation_R2[which(merged_1$ref_panel=="1000G_GRCh38")])
hist(merged_1$Imputation_R2[which(merged_1$ref_panel=="1000G_GRCh37_SAS")])
hist(merged_1$Imputation_R2[which(merged_1$ref_panel=="GenomeAsia_Pilot")])


#####      High Quality SNPS (R2 >= 0.8)   ##### 
metaimp <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_topmed_1000G/tarseq_meta_imputation_highR2.RSquare",header=F) # 1684 SNPs (has ~40 indels, some indel could have perfect true R2 value close to 1, e.g., chr19:2111874:AC:A has a MAF of 0.08, with a true R2 of 0.989504)
topmed <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/tarseq_topmed_grch38_highR2.RSquare",header=F) # 1549 SNPs
ex1000G <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/tarseq_ex1000G_grch38_highR2.RSquare",header=F) # 1435 SNPs
all_1000G_38 <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_38/tarseq_1000G_OG_grch38_highR2.RSquare",header=F) # 1328 SNPs
sas_1000G_37 <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_37/tarseq_1000G_OG_grch37_highR2.RSquare",header=F) # 1058 SNPs
GA <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_GA/tarseq_GA_grch37_highR2.RSquare",header=F) # 854 SNPs

dim(metaimp)
dim(topmed)
dim(ex1000G)
dim(all_1000G_38)
dim(sas_1000G_37)
dim(GA)

metaimp$ref_panel <- "TOPMed_1000G"
topmed$ref_panel <- "TOPMed"
ex1000G$ref_panel <- "Expanded_1000G"
all_1000G_38$ref_panel <- "1000G_GRCh38"
sas_1000G_37$ref_panel <- "1000G_GRCh37_SAS"
GA$ref_panel <- "GenomeAsia_Pilot"

merged <- rbind(metaimp,topmed,ex1000G,all_1000G_38,sas_1000G_37,GA) #  7908 rows 
colnames(merged) <- c("SNP_ID",	"Allele_Frequency",	"No_Samples"	,"Imputation_R2",	"Validation_AF","Imputation_AF"	, "ref_panel")

for (i in 1:dim(merged)[1]) {
 if (merged$Validation_AF[i] <= 0.5) {
    merged$MAF[i] <- merged$Validation_AF[i] 
    }
 else if (merged$Validation_AF[i] > 0.5) {
    merged$MAF[i] <- 1- merged$Validation_AF[i] 
    }
}

merged$MAF <- as.numeric(merged$MAF)
summary(merged$MAF)

library(dplyr)
merged$MAF_bin <- case_when(merged$MAF > 0 & merged$MAF <= 0.0005 ~ "(0,0.0005]", merged$MAF > 0.0005 & merged$MAF <= 0.001 ~ "(0.0005,0.001]", merged$MAF > 0.001 & merged$MAF <= 0.002 ~ "(0.001,0.002]", merged$MAF > 0.002 & merged$MAF <= 0.005 ~ "(0.002,0.005]",merged$MAF > 0.005 & merged$MAF <= 0.01 ~ "(0.005,0.01]",merged$MAF > 0.01 & merged$MAF <= 0.015 ~ "(0.01,0.015]",merged$MAF > 0.015 & merged$MAF <= 0.02 ~ "(0.015,0.02]", merged$MAF > 0.02 & merged$MAF <= 0.035 ~ "(0.02,0.035]", merged$MAF > 0.035 & merged$MAF <= 0.05 ~ "(0.035,0.05]", merged$MAF > 0.05 & merged$MAF <= 0.1 ~ "(0.05,0.1]", merged$MAF > 0.1 & merged$MAF <= 0.2  ~ "(0.1,0.2]", merged$MAF > 0.2 & merged$MAF <= 0.3  ~ "(0.2,0.3]", merged$MAF > 0.3 & merged$MAF <= 0.4  ~ "(0.3,0.4]", merged$MAF > 0.4 & merged$MAF <= 0.5  ~ "(0.4,0.5]")
table(merged$MAF_bin) # all MAF_bin add up to 7908 SNPs

# install.packages("plotrix")
# library(plotrix)
true_r2_maf <- aggregate(Imputation_R2~ref_panel+MAF_bin, data=merged, mean)
# FUN = function(x) c(mean = mean(x), se = std.error(x))
true_r2_maf$ref_panel <- factor(true_r2_maf$ref_panel,  levels=c("TOPMed_1000G","TOPMed","Expanded_1000G","1000G_GRCh38", "1000G_GRCh37_SAS", "GenomeAsia_Pilot"))
true_r2_maf$MAF_bin <- factor(true_r2_maf$MAF_bin,  levels=c("(0,0.0005]","(0.0005,0.001]","(0.001,0.002]","(0.002,0.005]","(0.005,0.01]","(0.01,0.015]","(0.015,0.02]","(0.02,0.035]","(0.035,0.05]","(0.05,0.1]","(0.1,0.2]","(0.2,0.3]","(0.3,0.4]","(0.4,0.5]"))

library(ggplot2)
library(dplyr)
pdf('true_r2_maf_well_imputed.pdf',height=10,width=15)
p <- ggplot(true_r2_maf, aes(x = MAF_bin, y = Imputation_R2,group=ref_panel, color=ref_panel))  + geom_line(size=1) + geom_point(shape=1)+ theme_classic(base_size=30) + coord_cartesian(ylim = c(0, 1)) + scale_x_discrete(guide = guide_axis(angle = 45))  + theme(text = element_text(size = 30,colour="black"))  + geom_hline(yintercept=0.8,  color = "dark grey") + geom_hline(yintercept=0.6, linetype="dashed", color = "dark grey")  + geom_vline(xintercept="(0.01,0.015]",  color = "dark grey",linetype="dashed") + geom_vline(xintercept="(0.05,0.1]",  color = "dark grey") #all font size
p + scale_color_manual(values=c("#661100", "#D55E00", "#E69F00","#0072B2","#332288","#999999"),name="Reference Panel") + labs(x="Minor Allele Frequency", y="Average True R2") +
  annotate('text', x = "(0.001,0.002]", y = 0.95, label='"*"', parse=TRUE, color = "red",size=10) +
  annotate('text', x = "(0.002,0.005]", y = 0.95, label='"**"', parse=TRUE, color = "red",size=10) +
  annotate('text', x = "(0.005,0.01]", y = 0.95, label='"**"', parse=TRUE, color = "red",size=10) 
dev.off()

# one way ANOVA
merged_1 <- merged[which(merged$MAF_bin == "(0,0.0005]"),] # 391 rows
merged_2 <- merged[which(merged$MAF_bin == "(0.0005,0.001]"),] # 487 rows
merged_3 <- merged[which(merged$MAF_bin == "(0.001,0.002]"),] # 329 rows
merged_4 <- merged[which(merged$MAF_bin == "(0.002,0.005]"),] # 816 rows
merged_5 <- merged[which(merged$MAF_bin == "(0.005,0.01]"),] # 683 rows
merged_6 <- merged[which(merged$MAF_bin == "(0.01,0.015]"),] # 296 rows
merged_7 <- merged[which(merged$MAF_bin == "(0.015,0.02]"),] # 170 rows
merged_8 <- merged[which(merged$MAF_bin == "(0.02,0.035]"),] # 576 rows
merged_9 <- merged[which(merged$MAF_bin == "(0.035,0.05]"),] # 360 rows
merged_10 <- merged[which(merged$MAF_bin == "(0.05,0.1]"),] # 857 rows
merged_11 <- merged[which(merged$MAF_bin == "(0.1,0.2]"),] # 972 rows
merged_12 <- merged[which(merged$MAF_bin == "(0.2,0.3]"),] # 812 rows
merged_13 <- merged[which(merged$MAF_bin == "(0.3,0.4]"),] # 701 rows
merged_14 <- merged[which(merged$MAF_bin == "(0.4,0.5]"),] # 458 rows

table(merged$MAF_bin)

# Non-parametric alternative to one-way ANOVA test
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_1) # P  0.6524
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_2) # P  0.4132
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_3) # P  0.006679 (*)
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_4) # P  6.878e-08 (**)
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_5) # P  0.0005763 (**)
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_6) # P  0.1958
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_7) # P  0.4084
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_8) # P  0.05995
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_9) # P  0.3902
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_10) # P  0.7652
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_11) # P  0.3168
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_12) # P  0.5183
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_13) # P  0.7797
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_14) # P  0.107
# 0.05/14 = 0.003571 (** Bonferroni significant )
# 0.05 (* nominal significant)

# check for pairwise difference
# Pairwise t-tests with no assumption of equal variances
pairwise.t.test(merged_3$Imputation_R2, merged_3$ref_panel,  p.adjust.method = "none", pool.sd = FALSE)  # TOPMed lower than all other panels

#                 1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     0.570            -            -              -               
# Expanded_1000G   0.452            0.788        -              -               
# GenomeAsia_Pilot 0.798            0.452        0.368          -               
# TOPMed           0.023            0.012        0.026          0.038           
# TOPMed_1000G     0.299            0.497        0.703          0.262           
#                  TOPMed
# 1000G_GRCh38     -     
# Expanded_1000G   -     
# GenomeAsia_Pilot -     
# TOPMed           -     
# TOPMed_1000G     0.041 

pairwise.t.test(merged_4$Imputation_R2, merged_4$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) # GenomeAsia better than all, except 1000G SAS GRCh37, TOPMed worse than all except meta-imputation

#                 1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
#  1000G_GRCh38     0.11509          -            -              -               
#  Expanded_1000G   0.11450          0.92426      -              -               
#  GenomeAsia_Pilot 0.15711          0.00210      0.00163        -               
#  TOPMed           0.00021          0.01389      0.00663        4.8e-07         
#  TOPMed_1000G     0.00688          0.22440      0.15709        2.9e-05         
#                   TOPMed 
#  1000G_GRCh38     -      
#  Expanded_1000G   -      
#  GenomeAsia_Pilot -      
#  TOPMed           -      
#  TOPMed_1000G     0.15124

pairwise.t.test(merged_5$Imputation_R2, merged_5$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) # GenomeAsia better than all other panels

#                  1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
# 1000G_GRCh38     0.04669          -            -              -               
# Expanded_1000G   0.02477          0.81502      -              -               
# GenomeAsia_Pilot 0.03773          3.7e-05      9.9e-06        -               
# TOPMed           0.00224          0.17011      0.23657        1.5e-06         
# TOPMed_1000G     0.17458          0.46766      0.32721        0.00031         
#                  TOPMed 
# 1000G_GRCh38     -      
# Expanded_1000G   -      
# GenomeAsia_Pilot -      
# TOPMed           -      
# TOPMed_1000G     0.04274

######################################################################
########      plot the number of SNPs overlap with Tarseq    #########
######################################################################

dat <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/snp_count_overlap",head=T)
dat$ref_panel <- factor(dat$ref_panel,  levels=c("Targeted_Sequencing","TOPMed_1000G","TOPMed","Expanded_1000G","1000G_GRCh38", "1000G_GRCh37_SAS", "GenomeAsia_Pilot"))

#####      ALL SNPs   ##### 
library(ggplot2)
# rare
# pdf(paste(plottingdir,paste(trait,p_type, 'ccgwas_GTEXv8_Manhattan.pdf', sep='_'), sep='/'),height=10,width=20)
pdf('number_snp_overlap_tarseq_rare_all.pdf',height=15,width=20)
p <- ggplot(data=dat, aes(x=ref_panel, y=snp_count_rare, fill=ref_panel)) + xlab("Reference Panel") + ylab("Count") + geom_bar(stat="identity") + theme_classic(base_size=30) +  geom_text(aes(label=scales::label_number_si(accuracy=0.01)(snp_count_rare)),size=10, position=position_dodge(width=0.9), vjust=-0.25)  + scale_x_discrete(guide = guide_axis(angle = 45))  +theme(axis.text = element_text(size = 30,colour="black"))
p + scale_fill_manual(values=c("#DDCC77","#661100", "#D55E00", "#E69F00","#0072B2","#332288","#999999"))+ theme(legend.position = "none")
dev.off()
# common
pdf('number_snp_overlap_tarseq_common_all.pdf',height=15,width=20)
p <- ggplot(data=dat, aes(x=ref_panel, y=snp_count_common, fill=ref_panel)) + xlab("Reference Panel") + ylab("Count") + geom_bar(stat="identity") + theme_classic(base_size=30) +  geom_text(aes(label=scales::label_number_si(accuracy=0.01)(snp_count_common)),size=10, position=position_dodge(width=0.9), vjust=-0.25)  + scale_x_discrete(guide = guide_axis(angle = 45))  +theme(axis.text = element_text(size = 30,colour="black"))
p + scale_fill_manual(values=c("#DDCC77", "#661100", "#D55E00", "#E69F00","#0072B2","#332288","#999999"))+ theme(legend.position = "none")
dev.off()

#####      High Quality SNPS (R2 >= 0.8)   ##### 
# rare
pdf('number_snp_overlap_tarseq_rare_well_imputed.pdf',height=15,width=20)
p <- ggplot(data=dat, aes(x=ref_panel, y=snp_count_rare_high_r2, fill=ref_panel)) + xlab("Reference Panel") + ylab("Count") + geom_bar(stat="identity") + theme_classic(base_size=30) +  geom_text(aes(label=scales::label_number_si(accuracy=0.01)(snp_count_rare_high_r2)),size=10, position=position_dodge(width=0.9), vjust=-0.25)  + scale_x_discrete(guide = guide_axis(angle = 45))  +theme(axis.text = element_text(size = 30,colour="black"))
p + scale_fill_manual(values=c("#DDCC77", "#661100", "#D55E00", "#E69F00","#0072B2","#332288","#999999"))+ theme(legend.position = "none")
dev.off()
# common
pdf('number_snp_overlap_tarseq_common_well_imputed.pdf',height=15,width=20)
p <- ggplot(data=dat, aes(x=ref_panel, y=snp_count_common_high_r2, fill=ref_panel)) + xlab("Reference Panel") + ylab("Count") + geom_bar(stat="identity") + theme_classic(base_size=30) +  geom_text(aes(label=scales::label_number_si(accuracy=0.01)(snp_count_common_high_r2)),size=10, position=position_dodge(width=0.9), vjust=-0.25)  + scale_x_discrete(guide = guide_axis(angle = 45)) +theme(axis.text = element_text(size = 30,colour="black"))
p + scale_fill_manual(values=c("#DDCC77","#661100", "#D55E00", "#E69F00","#0072B2","#332288","#999999"))+ theme(legend.position = "none")
dev.off()


#################################################################################################
#############        Estimated R2 by MAF for those SNPs with true R2 values      ################
#################################################################################################

metaimp <- read.table("",header=F) # 8483 SNPs
topmed <- read.table("",header=F) # 8267 SNPs
ex1000G <- read.table("",header=F) # 3469 SNPs
all_1000G_38 <- read.table("",header=F) # 3361 SNPs
sas_1000G_37 <- read.table("",header=F) # 3467 SNPs
GA <- read.table("",header=F) # 2126 SNPs

#### # meta-imputation

# get the SNP list
awk '(NR>1){print $1}' /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_topmed_1000G/tarseq_meta_imputation.RSquare  > /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_topmed_1000G/tarseq_meta_imputation.RSquare.snplist # 8483 SNPs (there are 6 duplicated SNPs, so 12 duplicated rows, 0.14%)
# check duplicated SNP names
sort -k1 tarseq_meta_imputation.RSquare.snplist|uniq -d > tarseq_meta_imputation.RSquare.snplist.dup
# remove duplicated SNP names
grep -vf tarseq_meta_imputation.RSquare.snplist.dup tarseq_meta_imputation.RSquare  > tarseq_meta_imputation.nodup.RSquare # 8471 SNPs 
# get the new SNP list 
awk '(NR>1){print $1}'  tarseq_meta_imputation.nodup.RSquare > tarseq_meta_imputation.nodup.RSquare.snplist # 8471 SNPs
# create a reverse SNP list
sed 's/\:/\t/g' tarseq_meta_imputation.nodup.RSquare.snplist | awk '{print $1":"$2":"$4":"$3}' > tarseq_meta_imputation.nodup.RSquare.snplist.reverse
# grep the forward SNPs from INFO file
grep -w -F -f tarseq_meta_imputation.nodup.RSquare.snplist all.chr.sort.info | awk '{print $1, $3}'> all.chr.trueR2.forward.info # 8462 SNPs
# -F: search for fixed strings (plain text) rather than regular expressions; -w: match whole word, -f: takes patterns from file, one per line.
# grep the reverse SNPs from INFO file
grep -w -F -f tarseq_meta_imputation.nodup.RSquare.snplist.reverse all.chr.sort.info | awk '{print $1, $3}' > all.chr.trueR2.reverse.info # 9 SNPs
# change the allele order in the reverse info file
sed 's/\:/\t/g' all.chr.trueR2.reverse.info | awk '{print $1":"$2":"$4":"$3"\t"$5}'  > all.chr.trueR2.reverse.flip.info 
# combine the forward and reverse INFO file
cat all.chr.trueR2.forward.info  all.chr.trueR2.reverse.flip.info >  all.chr.trueR2.info  # 8471 SNPs
# merge with Rsquare file
awk 'NR==FNR {h[$1]=$2; next} {print $0"\t"h[$1]}'   all.chr.trueR2.info  tarseq_meta_imputation.nodup.RSquare | awk '(NR>1){print $0}' |sed '1i SNP_ID\tAllele_Frequency\tNo_Samples\tImputation_R2\tValidation_AF\tImputation_AF\tRsq'> tarseq_meta_imputation.nodup.RSquare.Rsq # 8471 SNPs with 1 header line


# TOPMed
# get the SNP list
awk '(NR>1){print $1}' /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/tarseq_topmed_grch38.RSquare  >  /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/tarseq_topmed_grch38.RSquare.snplist # 8267 SNPs
# check duplicated SNP names
sort -k1 tarseq_topmed_grch38.RSquare.snplist |uniq -d > tarseq_topmed_grch38.RSquare.snplist.dup
# remove duplicated SNP names
grep -vf tarseq_topmed_grch38.RSquare.snplist.dup tarseq_topmed_grch38.RSquare  > tarseq_topmed_grch38.nodup.RSquare # 8257 SNPs 
# get the new SNP list 
awk '(NR>1){print $1}'  tarseq_topmed_grch38.nodup.RSquare > tarseq_topmed_grch38.nodup.RSquare.snplist # 8257 SNPs
# create a reverse SNP list
sed 's/\:/\t/g' tarseq_topmed_grch38.nodup.RSquare.snplist | awk '{print $1":"$2":"$4":"$3}' > tarseq_topmed_grch38.nodup.RSquare.snplist.reverse # 8257 SNPs
# grep the forward SNPs from INFO file
grep -w -F -f tarseq_topmed_grch38.nodup.RSquare.snplist all.chr.sort.info | awk '{print $1, $7}'> all.chr.trueR2.forward.info # 8248 SNPs
# -F: search for fixed strings (plain text) rather than regular expressions; -w: match whole word, -f: takes patterns from file, one per line.
# grep the reverse SNPs from INFO file
grep -w -F -f tarseq_topmed_grch38.nodup.RSquare.snplist.reverse all.chr.sort.info | awk '{print $1, $7}' > all.chr.trueR2.reverse.info # 9 SNPs
# change the allele order in the reverse info file
sed 's/\:/\t/g' all.chr.trueR2.reverse.info | awk '{print $1":"$2":"$4":"$3"\t"$5}'  > all.chr.trueR2.reverse.flip.info 
# combine the forward and reverse INFO file
cat all.chr.trueR2.forward.info  all.chr.trueR2.reverse.flip.info >  all.chr.trueR2.info  # 8257 SNPs
# merge with Rsquare file
awk 'NR==FNR {h[$1]=$2; next} {print $0"\t"h[$1]}'   all.chr.trueR2.info  tarseq_topmed_grch38.nodup.RSquare | awk '(NR>1){print $0}' |sed '1i SNP_ID\tAllele_Frequency\tNo_Samples\tImputation_R2\tValidation_AF\tImputation_AF\tRsq'> tarseq_topmed_grch38.nodup.RSquare.Rsq # 8257 SNPs with 1 header line


# expanded 1000G
# get the SNP list
awk '(NR>1){print $1}' /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/tarseq_ex1000G_grch38.RSquare  > /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/tarseq_ex1000G_grch38.RSquare.snplist # 3469 SNPs
# check duplicated SNP names
sort -k1 tarseq_ex1000G_grch38.RSquare.snplist |uniq -d > tarseq_ex1000G_grch38.RSquare.snplist.dup # 1 duplicated SNP
# remove duplicated SNP names
grep -vf tarseq_ex1000G_grch38.RSquare.snplist.dup tarseq_ex1000G_grch38.RSquare  > tarseq_ex1000G_grch38.nodup.RSquare # 3467 SNPs 
# get the new SNP list 
awk '(NR>1){print $1}'  tarseq_ex1000G_grch38.nodup.RSquare > tarseq_ex1000G_grch38.nodup.RSquare.snplist # 3467 SNPs
# create a reverse SNP list
sed 's/\:/\t/g' tarseq_ex1000G_grch38.nodup.RSquare.snplist  | awk '{print $1":"$2":"$4":"$3}' > tarseq_ex1000G_grch38.nodup.RSquare.snplist.reverse # 3467 SNPs
# grep the forward SNPs from INFO file
grep -w -F -f tarseq_ex1000G_grch38.nodup.RSquare.snplist all.chr.sort.info | awk '{print $1, $7}'> all.chr.trueR2.forward.info # 3465 SNPs
# -F: search for fixed strings (plain text) rather than regular expressions; -w: match whole word, -f: takes patterns from file, one per line.
# grep the reverse SNPs from INFO file
grep -w -F -f tarseq_ex1000G_grch38.nodup.RSquare.snplist.reverse all.chr.sort.info | awk '{print $1, $7}' > all.chr.trueR2.reverse.info # 2 SNPs
# change the allele order in the reverse info file
sed 's/\:/\t/g' all.chr.trueR2.reverse.info | awk '{print $1":"$2":"$4":"$3"\t"$5}'  > all.chr.trueR2.reverse.flip.info 
# combine the forward and reverse INFO file
cat all.chr.trueR2.forward.info  all.chr.trueR2.reverse.flip.info >  all.chr.trueR2.info  # 3467 SNPs
# merge with Rsquare file
awk 'NR==FNR {h[$1]=$2; next} {print $0"\t"h[$1]}'   all.chr.trueR2.info  tarseq_ex1000G_grch38.nodup.RSquare | awk '(NR>1){print $0}' |sed '1i SNP_ID\tAllele_Frequency\tNo_Samples\tImputation_R2\tValidation_AF\tImputation_AF\tRsq'> tarseq_ex1000G_grch38.nodup.RSquare.Rsq # 3467 SNPs with 1 header line


# 1000G GRCh38
# get the SNP list
awk '(NR>1){print $1}' /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_38/tarseq_1000G_OG_grch38.RSquare  > /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_38/tarseq_1000G_OG_grch38.RSquare.snplist # 3361 SNPs
# check duplicated SNP names
sort -k1 tarseq_1000G_OG_grch38.RSquare.snplist |uniq -d # no duplicated SNP
# create a reverse SNP list
sed 's/\:/\t/g' tarseq_1000G_OG_grch38.RSquare.snplist | awk '{print $1":"$2":"$4":"$3}' > tarseq_1000G_OG_grch38.RSquare.snplist.reverse # 3361 SNPs
# grep the forward SNPs from INFO file
grep -w -F -f tarseq_1000G_OG_grch38.RSquare.snplist  all.chr.sort.info | awk '{print $1, $7}'> all.chr.trueR2.forward.info # 3357 SNPs
# -F: search for fixed strings (plain text) rather than regular expressions; -w: match whole word, -f: takes patterns from file, one per line.
# grep the reverse SNPs from INFO file
grep -w -F -f tarseq_1000G_OG_grch38.RSquare.snplist.reverse  all.chr.sort.info | awk '{print $1, $7}' > all.chr.trueR2.reverse.info # 4 SNPs
# change the allele order in the reverse info file
sed 's/\:/\t/g' all.chr.trueR2.reverse.info | awk '{print $1":"$2":"$4":"$3"\t"$5}'  > all.chr.trueR2.reverse.flip.info 
# combine the forward and reverse INFO file
cat all.chr.trueR2.forward.info  all.chr.trueR2.reverse.flip.info >  all.chr.trueR2.info  # 3361 SNPs
# merge with Rsquare file
awk 'NR==FNR {h[$1]=$2; next} {print $0"\t"h[$1]}'   all.chr.trueR2.info  tarseq_1000G_OG_grch38.RSquare | awk '(NR>1){print $0}' |sed '1i SNP_ID\tAllele_Frequency\tNo_Samples\tImputation_R2\tValidation_AF\tImputation_AF\tRsq'> tarseq_1000G_OG_grch38.nodup.RSquare.Rsq # 3361 SNPs with 1 header line


# 1000G GRCh37 SAS
# get the SNP list
awk '(NR>1){print $1}' /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_37/tarseq_1000G_OG_grch37.RSquare > /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_37/tarseq_1000G_OG_grch37.RSquare.snplist # 3467 SNPs
# check duplicated SNP names
sort -k1 tarseq_1000G_OG_grch37.RSquare.snplist |uniq -d # no duplicated SNP
# create a reverse SNP list
sed 's/\:/\t/g' tarseq_1000G_OG_grch37.RSquare.snplist | awk '{print $1":"$2":"$4":"$3}' > tarseq_1000G_OG_grch37.RSquare.snplist.reverse # 3467 SNPs
# grep the forward SNPs from INFO file
grep -w -F -f tarseq_1000G_OG_grch37.RSquare.snplist  all.chr.sort.info | awk '{print $1, $7}'> all.chr.trueR2.forward.info # 3465 SNPs
# -F: search for fixed strings (plain text) rather than regular expressions; -w: match whole word, -f: takes patterns from file, one per line.
# grep the reverse SNPs from INFO file
grep -w -F -f tarseq_1000G_OG_grch37.RSquare.snplist.reverse  all.chr.sort.info | awk '{print $1, $7}' > all.chr.trueR2.reverse.info # 2 SNPs
# change the allele order in the reverse info file
sed 's/\:/\t/g' all.chr.trueR2.reverse.info | awk '{print $1":"$2":"$4":"$3"\t"$5}'  > all.chr.trueR2.reverse.flip.info 
# combine the forward and reverse INFO file
cat all.chr.trueR2.forward.info  all.chr.trueR2.reverse.flip.info >  all.chr.trueR2.info  # 3467 SNPs
# merge with Rsquare file
awk 'NR==FNR {h[$1]=$2; next} {print $0"\t"h[$1]}'   all.chr.trueR2.info  tarseq_1000G_OG_grch37.RSquare | awk '(NR>1){print $0}' |sed '1i SNP_ID\tAllele_Frequency\tNo_Samples\tImputation_R2\tValidation_AF\tImputation_AF\tRsq'> tarseq_1000G_OG_grch37.nodup.RSquare.Rsq # 3467 SNPs with 1 header line


# GenomeAsia
# get the SNP list
awk '(NR>1){print $1}' /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_GA/tarseq_GA_grch37.RSquare > /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_GA/tarseq_GA_grch37.RSquare.snplist # 2126 SNPs
# check duplicated SNP names
sort -k1 tarseq_GA_grch37.RSquare.snplist |uniq -d # no duplicated SNP
# create a reverse SNP list
sed 's/\:/\t/g' tarseq_GA_grch37.RSquare.snplist | awk '{print $1":"$2":"$4":"$3}' > tarseq_GA_grch37.RSquare.snplist.reverse # 2126 SNPs
# grep the forward SNPs from INFO file
grep -w -F -f tarseq_GA_grch37.RSquare.snplist  all.chr.sort.info | awk '{print $1, $7}'> all.chr.trueR2.info # 2126 SNPs
# -F: search for fixed strings (plain text) rather than regular expressions; -w: match whole word, -f: takes patterns from file, one per line.
# grep the reverse SNPs from INFO file
grep -w -F -f tarseq_GA_grch37.RSquare.snplist.reverse  all.chr.sort.info | awk '{print $1, $7}' > all.chr.trueR2.reverse.info # 0 SNPs
# merge with Rsquare file
awk 'NR==FNR {h[$1]=$2; next} {print $0"\t"h[$1]}'   all.chr.trueR2.info  tarseq_GA_grch37.RSquare | awk '(NR>1){print $0}' |sed '1i SNP_ID\tAllele_Frequency\tNo_Samples\tImputation_R2\tValidation_AF\tImputation_AF\tRsq'> tarseq_GA_grch37.nodup.RSquare.Rsq # 2126 SNPs with 1 header line


# in R
metaimp <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_topmed_1000G/tarseq_meta_imputation.nodup.RSquare.Rsq",header=T) # 8483 SNPs
topmed <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/tarseq_topmed_grch38.nodup.RSquare.Rsq",header=T) # 8267 SNPs
ex1000G <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/tarseq_ex1000G_grch38.nodup.RSquare.Rsq",header=T) # 3469 SNPs
all_1000G_38 <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_38/tarseq_1000G_OG_grch38.nodup.RSquare.Rsq",header=T) # 3361 SNPs
sas_1000G_37 <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_37/tarseq_1000G_OG_grch37.nodup.RSquare.Rsq",header=T) # 3467 SNPs
GA <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_GA/tarseq_GA_grch37.nodup.RSquare.Rsq",header=T) # 2126 SNPs

metaimp$ref_panel <- "TOPMed_1000G"
topmed$ref_panel <- "TOPMed"
ex1000G$ref_panel <- "Expanded_1000G"
all_1000G_38$ref_panel <- "1000G_GRCh38"
sas_1000G_37$ref_panel <- "1000G_GRCh37_SAS"
GA$ref_panel <- "GenomeAsia_Pilot"

dim(metaimp) # [1]  8471    8
dim(topmed) # [1] 8257    8
dim(ex1000G) # [1] 3467    8
dim(all_1000G_38) # [1] 3361    8
dim(sas_1000G_37) # [1] 3467    8
dim(GA) # [1] 2126    8

merged <- rbind(metaimp,topmed,ex1000G,all_1000G_38,sas_1000G_37,GA) #  [1] 29149     8

for (i in 1:dim(merged)[1]) {
 if (merged$Validation_AF[i] <= 0.5) {
    merged$MAF[i] <- merged$Validation_AF[i] 
    }
 else if (merged$Validation_AF[i] > 0.5) {
    merged$MAF[i] <- 1- merged$Validation_AF[i] 
    }
}

merged$MAF <- as.numeric(merged$MAF)
summary(merged$MAF)

library(dplyr)
merged$MAF_bin <- case_when(merged$MAF > 0 & merged$MAF <= 0.0005 ~ "(0,0.0005]", merged$MAF > 0.0005 & merged$MAF <= 0.001 ~ "(0.0005,0.001]", merged$MAF > 0.001 & merged$MAF <= 0.002 ~ "(0.001,0.002]", merged$MAF > 0.002 & merged$MAF <= 0.005 ~ "(0.002,0.005]",merged$MAF > 0.005 & merged$MAF <= 0.01 ~ "(0.005,0.01]",merged$MAF > 0.01 & merged$MAF <= 0.015 ~ "(0.01,0.015]",merged$MAF > 0.015 & merged$MAF <= 0.02 ~ "(0.015,0.02]", merged$MAF > 0.02 & merged$MAF <= 0.035 ~ "(0.02,0.035]", merged$MAF > 0.035 & merged$MAF <= 0.05 ~ "(0.035,0.05]", merged$MAF > 0.05 & merged$MAF <= 0.1 ~ "(0.05,0.1]", merged$MAF > 0.1 & merged$MAF <= 0.2  ~ "(0.1,0.2]", merged$MAF > 0.2 & merged$MAF <= 0.3  ~ "(0.2,0.3]", merged$MAF > 0.3 & merged$MAF <= 0.4  ~ "(0.3,0.4]", merged$MAF > 0.4 & merged$MAF <= 0.5  ~ "(0.4,0.5]")
table(merged$MAF_bin) # add up to 29149

merged$Rsq <- as.numeric(merged$Rsq)

# install.packages("plotrix")
# library(plotrix)
est_r2_maf <- aggregate(Rsq~ref_panel+MAF_bin, data=merged, mean)
# FUN = function(x) c(mean = mean(x), se = std.error(x))
est_r2_maf$ref_panel <- factor(est_r2_maf$ref_panel,  levels=c("TOPMed_1000G","TOPMed","Expanded_1000G","1000G_GRCh38", "1000G_GRCh37_SAS", "GenomeAsia_Pilot"))
est_r2_maf$MAF_bin <- factor(est_r2_maf$MAF_bin,  levels=c("(0,0.0005]","(0.0005,0.001]","(0.001,0.002]","(0.002,0.005]","(0.005,0.01]","(0.01,0.015]","(0.015,0.02]","(0.02,0.035]","(0.035,0.05]","(0.05,0.1]","(0.1,0.2]","(0.2,0.3]","(0.3,0.4]","(0.4,0.5]"))

library(ggplot2)
library(dplyr)
pdf('estimated_r2_maf_snps_with_true_r2.pdf',height=10,width=15)
p <- ggplot(est_r2_maf, aes(x = MAF_bin, y = Rsq,group=ref_panel, color=ref_panel))  + geom_line(size=1) + geom_point(shape=1)+ theme_classic(base_size=30) + coord_cartesian(ylim = c(0, 1)) + scale_x_discrete(guide = guide_axis(angle = 45))  + theme(text = element_text(size = 30,colour="black"))  + geom_hline(yintercept=0.8,  color = "dark grey") + geom_hline(yintercept=0.6, linetype="dashed", color = "dark grey")  + geom_vline(xintercept="(0.01,0.015]",  color = "dark grey",linetype="dashed") + geom_vline(xintercept="(0.05,0.1]",  color = "dark grey") #all font size
p + scale_color_manual(values=c("#661100", "#D55E00", "#E69F00","#0072B2","#332288","#999999"),name="Reference Panel") + labs(x="Minor Allele Frequency", y="Average Estimated R2") + theme(plot.title = element_text(hjust = 0.5)) + labs(title="SNPs (imputed and sequenced)") +
  annotate('text', x = "(0,0.0005]", y = 0.4, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.0005,0.001]", y = 0.5, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.001,0.002]", y = 0.6, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.002,0.005]", y = 0.7, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.005,0.01]", y = 0.78, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.01,0.015]", y = 0.85, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.015,0.02]", y = 0.87, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.02,0.035]", y = 0.9, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.035,0.05]", y = 0.95, label='"*"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.05,0.1]", y = 0.96, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.1,0.2]", y = 0.98, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.2,0.3]", y = 0.98, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.3,0.4]", y = 0.98, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.4,0.5]", y = 0.98, label='"**"', parse=TRUE, color = "red",size=10) 
dev.off()

# one way ANOVA
merged_1 <- merged[which(merged$MAF_bin == "(0,0.0005]"),] # 7619 rows
merged_2 <- merged[which(merged$MAF_bin == "(0.0005,0.001]"),] # 4654 rows
merged_3 <- merged[which(merged$MAF_bin == "(0.001,0.002]"),] # 2745 rows
merged_4 <- merged[which(merged$MAF_bin == "(0.002,0.005]"),] # 3494 rows
merged_5 <- merged[which(merged$MAF_bin == "(0.005,0.01]"),] # 2016 rows
merged_6 <- merged[which(merged$MAF_bin == "(0.01,0.015]"),] # 902 rows
merged_7 <- merged[which(merged$MAF_bin == "(0.015,0.02]"),] # 464 rows
merged_8 <- merged[which(merged$MAF_bin == "(0.02,0.035]"),] # 1029 rows
merged_9 <- merged[which(merged$MAF_bin == "(0.035,0.05]"),] # 597 rows
merged_10 <- merged[which(merged$MAF_bin == "(0.05,0.1]"),] # 1271 rows
merged_11 <- merged[which(merged$MAF_bin == "(0.1,0.2]"),] # 1484 rows
merged_12 <- merged[which(merged$MAF_bin == "(0.2,0.3]"),] # 1167 rows
merged_13 <- merged[which(merged$MAF_bin == "(0.3,0.4]"),] # 1007 rows
merged_14 <- merged[which(merged$MAF_bin == "(0.4,0.5]"),] # 700 rows

dim(merged_1)
dim(merged_2)
dim(merged_3)
dim(merged_4)
dim(merged_5)
dim(merged_6)
dim(merged_7)
dim(merged_8)
dim(merged_9)
dim(merged_10)
dim(merged_11)
dim(merged_12)
dim(merged_13)
dim(merged_14)

table(merged$MAF_bin)
# testing if the aggregate function calculate the average true R2 correctly: YES! 

# Non-parametric alternative to one-way ANOVA test
kruskal.test(Rsq ~ ref_panel, data = merged_1) # P  < 2.2e-16 (**) 
kruskal.test(Rsq ~ ref_panel, data = merged_2) # P  < 2.2e-16 (**)
kruskal.test(Rsq ~ ref_panel, data = merged_3) # P  < 2.2e-16 (**)
kruskal.test(Rsq ~ ref_panel, data = merged_4) # P  < 2.2e-16 (**)
kruskal.test(Rsq ~ ref_panel, data = merged_5) # P  < 2.2e-16 (**)
kruskal.test(Rsq ~ ref_panel, data = merged_6) # P  4.626e-13 (**)
kruskal.test(Rsq ~ ref_panel, data = merged_7) # P  5.064e-05 (**)
kruskal.test(Rsq ~ ref_panel, data = merged_8) # P  8.883e-07 (**)
kruskal.test(Rsq ~ ref_panel, data = merged_9) # P  0.003399 (*)
kruskal.test(Rsq ~ ref_panel, data = merged_10) # P  4.018e-08 (**)
kruskal.test(Rsq ~ ref_panel, data = merged_11) # P  < 2.2e-16 (**)
kruskal.test(Rsq ~ ref_panel, data = merged_12) # P  4.307e-16 (**)
kruskal.test(Rsq ~ ref_panel, data = merged_13) # P  < 2.2e-16 (**)
kruskal.test(Rsq ~ ref_panel, data = merged_14) # P  4.015e-13 (**)
# 0.05/15 = 0.003 (** Bonferroni significant )
# 0.05 (* nominal significant)

# http://www.sthda.com/english/wiki/one-way-anova-test-in-r

# check for pairwise difference
# Pairwise t-tests with no assumption of equal variances
pairwise.t.test(merged_1$Rsq, merged_1$ref_panel,  p.adjust.method = "none", pool.sd = FALSE) 
# true R2: GA > ex1KG ~ 1KG38 ~ 1KG37 > meta ~ TOPMed
# est R2: 1KG38 ~ ex1KG > GA ~ 1K37 > meta ~ TOPMed

                 1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
1000G_GRCh38     2.2e-11          -            -              -               
Expanded_1000G   2.6e-11          0.90528      -              -               
GenomeAsia_Pilot 0.44155          6.6e-05      8.4e-05        -               
TOPMed           1.1e-05          < 2e-16      < 2e-16        0.00061         
TOPMed_1000G     0.00612          < 2e-16      < 2e-16        0.01378         
                 TOPMed 
1000G_GRCh38     -      
Expanded_1000G   -      
GenomeAsia_Pilot -      
TOPMed           -      
TOPMed_1000G     0.00810

pairwise.t.test(merged_2$Rsq, merged_2$ref_panel,  p.adjust.method = "none", pool.sd = FALSE)
# true R2: ex1KG ~ GA ~ 1KG38 ~ 1KG37 > meta ~ TOPMed
# est R2: ex1KG ~ 1KG38 > meta > TOPMed ~ 1KG37 ~ GA

                 1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
1000G_GRCh38     2.1e-14          -            -              -               
Expanded_1000G   < 2e-16          0.2997       -              -               
GenomeAsia_Pilot 0.2480           2.0e-12      3.4e-15        -               
TOPMed           0.5131           7.5e-16      < 2e-16        0.0851          
TOPMed_1000G     0.0047           5.8e-10      4.6e-14        0.0011          
                 TOPMed
1000G_GRCh38     -     
Expanded_1000G   -     
GenomeAsia_Pilot -     
TOPMed           -     
TOPMed_1000G     0.0070

pairwise.t.test(merged_3$Rsq, merged_3$ref_panel,  p.adjust.method = "none", pool.sd = FALSE)
# true R2: ex1KG ~ meta ~ 1KG38 > 1KG37 ~ TOPMed ~ GA
# est R2: ex1KG > 1KG38 ~ meta > TOMed  > 1KG37 > GA

                 1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
1000G_GRCh38     1.4e-08          -            -              -               
Expanded_1000G   2.1e-15          0.02934      -              -               
GenomeAsia_Pilot 0.00038          2.6e-15      < 2e-16        -               
TOPMed           0.00145          0.00125      6.8e-09        7.0e-10         
TOPMed_1000G     2.2e-08          0.35247      0.00046        4.9e-15         
                 TOPMed 
1000G_GRCh38     -      
Expanded_1000G   -      
GenomeAsia_Pilot -      
TOPMed           -      
TOPMed_1000G     0.00523

pairwise.t.test(merged_4$Rsq, merged_4$ref_panel,  p.adjust.method = "none", pool.sd = FALSE)
# true R2: ex1KG ~ meta > 1KG37 ~ 1KG38 ~ GA ~ TOPMed
# est R2: ex1KG ~ meta > 1KG38 ~ TOPMed > 1KG37 > GA

                 1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
1000G_GRCh38     1.4e-10          -            -              -               
Expanded_1000G   < 2e-16          0.01698      -              -               
GenomeAsia_Pilot 0.00812          5.7e-15      < 2e-16        -               
TOPMed           1.2e-08          0.24037      0.00018        4.3e-13         
TOPMed_1000G     < 2e-16          0.11100      0.30305        < 2e-16         
                 TOPMed 
1000G_GRCh38     -      
Expanded_1000G   -      
GenomeAsia_Pilot -      
TOPMed           -      
TOPMed_1000G     0.00232

pairwise.t.test(merged_5$Rsq, merged_5$ref_panel,  p.adjust.method = "none", pool.sd = FALSE)
# true R2: meta-imputation ~ ex1KG > 1KG37 ~ 1KG38 ~ GA ~ TOPMed
# est R2: meta ~ ex1KG > 1K38 ~ TOPMed >  1KG37  > GA

                 1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
1000G_GRCh38     2.2e-06          -            -              -               
Expanded_1000G   4.3e-10          0.12602      -              -               
GenomeAsia_Pilot 0.00593          3.4e-11      7.9e-15        -               
TOPMed           2.0e-06          0.71415      0.03894        2.9e-11         
TOPMed_1000G     7.8e-14          0.00973      0.33781        < 2e-16         
                 TOPMed 
1000G_GRCh38     -      
Expanded_1000G   -      
GenomeAsia_Pilot -      
TOPMed           -      
TOPMed_1000G     0.00084

pairwise.t.test(merged_6$Rsq, merged_6$ref_panel,  p.adjust.method = "none", pool.sd = FALSE)
# true R2: meta-imputation > TOPMed (mean: 0.6588, p=0.049) ~ ex1KG (mean: 0.6585, p=0.071, possibly fewer SNPs therefore larger P) ~ 1KG38 ~ 1KG37 ~ GA
# est R2: meta ~ ex1KG ~ TOPMed > 1KG38 > 1K37 >  GA 

                 1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
1000G_GRCh38     0.00151          -            -              -               
Expanded_1000G   0.00013          0.54917      -              -               
GenomeAsia_Pilot 0.01448          3.3e-07      1.6e-08        -               
TOPMed           0.00013          0.66731      0.83140        1.5e-08         
TOPMed_1000G     2.7e-07          0.06811      0.21251        3.7e-11         
                 TOPMed 
1000G_GRCh38     -      
Expanded_1000G   -      
GenomeAsia_Pilot -      
TOPMed           -      
TOPMed_1000G     0.11389

pairwise.t.test(merged_7$Rsq, merged_7$ref_panel,  p.adjust.method = "none", pool.sd = FALSE)
# true R2: meta ~ TOPMed ~ Ex1KG ~ 1K38 ~ GA > 1K37
# est R2: meta ~TOPMed ~ ex1KG ~1K38 > 1KG37 ~ GA 

                 1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
1000G_GRCh38     0.00957          -            -              -               
Expanded_1000G   0.00496          0.79963      -              -               
GenomeAsia_Pilot 0.46843          0.00233      0.00123        -               
TOPMed           0.00077          0.54210      0.75032        0.00023         
TOPMed_1000G     5.7e-05          0.15360      0.25610        2.4e-05         
                 TOPMed 
1000G_GRCh38     -      
Expanded_1000G   -      
GenomeAsia_Pilot -      
TOPMed           -      
TOPMed_1000G     0.31629

pairwise.t.test(merged_8$Rsq, merged_8$ref_panel,  p.adjust.method = "none", pool.sd = FALSE)
# true R2: no sig pair, meta ~ ex1KG ~ 1G37 ~ TOPMed ~ 1K38 ~ GA
# est R2: meta ~ex1KG ~1K38 ~ TOPMed >  1KG37 ~ GA

                 1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
1000G_GRCh38     0.00144          -            -              -               
Expanded_1000G   0.00034          0.66628      -              -               
GenomeAsia_Pilot 0.13340          2.6e-05      6.0e-06        -               
TOPMed           0.00156          0.77043      0.44378        2.7e-05         
TOPMed_1000G     1.2e-05          0.28866      0.55693        2.7e-07         
                 TOPMed 
1000G_GRCh38     -      
Expanded_1000G   -      
GenomeAsia_Pilot -      
TOPMed           -      
TOPMed_1000G     0.12768

pairwise.t.test(merged_9$Rsq, merged_9$ref_panel,  p.adjust.method = "none", pool.sd = FALSE)
# true R2: no sig pair, meta ~ TOPMed ~ 1K38 ~ ex1KG ~ GA ~ 1KG37
# est R2: meta ~ ex1KG ~ 1K38 ~ TOPMED > 1K37 ~ GA

                 1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
1000G_GRCh38     0.00755          -            -              -               
Expanded_1000G   0.00373          0.79831      -              -               
GenomeAsia_Pilot 0.58328          0.00452      0.00252        -               
TOPMed           0.00522          0.98370      0.76156        0.00348         
TOPMed_1000G     0.00027          0.25974      0.38274        0.00033         
                 TOPMed 
1000G_GRCh38     -      
Expanded_1000G   -      
GenomeAsia_Pilot -      
TOPMed           -      
TOPMed_1000G     0.18960

pairwise.t.test(merged_10$Rsq, merged_10$ref_panel,  p.adjust.method = "none", pool.sd = FALSE)
# true R2: no sig pair, ex1KG ~ meta ~ GA ~ 1K38 ~ 1K37 ~ TOPMed
# est R2: ex1KG ~ meta ~ 1K38 ~ TOPMed > 1K37 ~ GA

                 1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
1000G_GRCh38     0.00096          -            -              -               
Expanded_1000G   0.00028          0.75982      -              -               
GenomeAsia_Pilot 0.50489          0.00018      5.1e-05        -               
TOPMed           0.02627          0.13727      0.06270        0.00550         
TOPMed_1000G     0.00031          0.97874      0.75361        5.7e-05         
                 TOPMed 
1000G_GRCh38     -      
Expanded_1000G   -      
GenomeAsia_Pilot -      
TOPMed           -      
TOPMed_1000G     0.08656

pairwise.t.test(merged_11$Rsq, merged_11$ref_panel,  p.adjust.method = "none", pool.sd = FALSE)
# true R2: no sig pair, meta ~ ex1KG ~ TOPMed ~ GA ~ 1KG38 ~ 1KG37
# est R2: ex1KG  ~ meta ~ 1K38 > TOPMed > 1K37 ~ GA

                 1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
1000G_GRCh38     8.1e-07          -            -              -               
Expanded_1000G   1.2e-08          0.47497      -              -               
GenomeAsia_Pilot 0.07580          2.1e-09      4.3e-11        -               
TOPMed           0.00119          0.01905      0.00110        3.6e-06         
TOPMed_1000G     5.4e-09          0.59161      0.79453        2.6e-11         
                 TOPMed 
1000G_GRCh38     -      
Expanded_1000G   -      
GenomeAsia_Pilot -      
TOPMed           -      
TOPMed_1000G     0.00083

pairwise.t.test(merged_12$Rsq, merged_12$ref_panel,  p.adjust.method = "none", pool.sd = FALSE)
# true R2: no sig pair, meta ~ GA ~ ex1KG ~ TOPMed ~ 1KG38 ~ 1KG37
# est R2: ex1KG ~ meta ~ 1K38 > TOPMed > 1K37 ~ GA

                 1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
1000G_GRCh38     0.00020          -            -              -               
Expanded_1000G   1.7e-06          0.40145      -              -               
GenomeAsia_Pilot 0.50756          3.0e-05      2.8e-07        -               
TOPMed           0.01646          0.03976      0.00058        0.00288         
TOPMed_1000G     5.9e-06          0.78188      0.44862        9.8e-07         
                 TOPMed 
1000G_GRCh38     -      
Expanded_1000G   -      
GenomeAsia_Pilot -      
TOPMed           -      
TOPMed_1000G     0.00219

pairwise.t.test(merged_13$Rsq, merged_13$ref_panel,  p.adjust.method = "none", pool.sd = FALSE)
# true R2: no sig pair, meta ~ TOPMed ~ ex1KG ~ 1K37 ~ 1K38 ~ GA
# est R2: ex1KG ~ meta ~ 1K38 > TOPMed > 1K37 ~ GA

                 1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
1000G_GRCh38     0.00017          -            -              -               
Expanded_1000G   1.6e-05          0.65641      -              -               
GenomeAsia_Pilot 0.09220          2.5e-06      2.9e-07        -               
TOPMed           0.03183          0.03467      0.00639        0.00050         
TOPMed_1000G     1.2e-05          0.88429      0.71439        2.7e-07         
                 TOPMed 
1000G_GRCh38     -      
Expanded_1000G   -      
GenomeAsia_Pilot -      
TOPMed           -      
TOPMed_1000G     0.00709

pairwise.t.test(merged_14$Rsq, merged_14$ref_panel,  p.adjust.method = "none", pool.sd = FALSE)
# true R2: meta-imputation ~ TOPMed ~ ex1KG ~ 1KG38 ~ 1KG37 > GA
# est R2: meta ~ ex1KG ~ 1K38 ~ TOPMed > 1K37 ~ GA

                 1000G_GRCh37_SAS 1000G_GRCh38 Expanded_1000G GenomeAsia_Pilot
1000G_GRCh38     0.02550          -            -              -               
Expanded_1000G   0.01076          0.80283      -              -               
GenomeAsia_Pilot 0.05013          0.00025      9.1e-05        -               
TOPMed           0.04528          0.52749      0.33845        0.00037         
TOPMed_1000G     0.00044          0.38107      0.53862        4.4e-06         
                 TOPMed 
1000G_GRCh38     -      
Expanded_1000G   -      
GenomeAsia_Pilot -      
TOPMed           -      
TOPMed_1000G     0.05104





###############################################################
#############        aggRSquare Comparison     ################
###############################################################
ml git python
ml cmake
/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/aggRSquare/release-build/aggRSquare -v /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.v4.2.updateID.sort.recode.vcf -i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_GA/all_snp_overlap_tarseq_ordered.sort.vcf -o tarseq_GA_grch37  --bins /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_compare/aggRSquare/default_MAF_bins --detail


### meta-imputation

#!/bin/bash
#BSUB -J chr1_aggR2_tarseq_metaimp_forward_reverse # Job name
#BSUB -n 1
#BSUB -P acc_psychgen # allocation account
#BSUB -q premium # queue
#BSUB -W 12:00 # walltime in HH:MM
#BSUB -R himem
#BSUB -R rusage[mem=96G] # 120 GB of memory 
#BSUB -R span[hosts=1] # all cores from one node
#BSUB -o %J.stdout # output log (%J : JobID)
#BSUB -eo %J.stderr # error log
#BSUB -L /bin/bash # Initialize the execution environment

ml git python
ml cmake
# meta-imputation
/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/aggRSquare/release-build/aggRSquare -v /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr1.recode.vcf -i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_topmed_1000G/chr1_topmed_1000g.meta.tarseq.forward.1748.recode.vcf.gz -o tarseq_metaimp_chr1_forward_test  --bins /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_compare/aggRSquare/default_MAF_bins --detail
# -- Found 5735 common variants  -- Analysis included 1025 variants.

# because it's the reverse order, so the common variants cannot be found
# after the SNP name changed (with allele order fliped), common variants still cannot be found
# the REF and ALT alleles are flipped, but the individual genotypes are still matched... (so it is the REF ALT itself need to be flipped)
# After flipping the ref and alt alleles manually in the meta-imputation file, it finally matches!! 
/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/aggRSquare/release-build/aggRSquare -v /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr1.recode.vcf -i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_topmed_1000G/chr1_topmed_1000g.meta.tarseq.reverse.1748.recode.vcf.gz -o tarseq_metaimp_chr1_reverse_test  --bins /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_compare/aggRSquare/default_MAF_bins --detail
# aggRSquare don't match SNPs based on $3 (SNP Names), instead it match SNPS based on chr, pos, ref allele, alt allele (so allele order matters)
#  -- Found 6 common variants,    -- Analysis included 2 variants.




bsub < chr1_aggR2_tarseq_metaimp_forward_reverse.sh



## limit the meta-imputation to the 1748 individuals

#!/bin/bash
#BSUB -J vcf_gz_remove_sample # Job name
#BSUB -n 1
#BSUB -P acc_psychgen # allocation account
#BSUB -q premium # queue
#bsub -n 1
#BSUB -W 48:00 # walltime in HH:MM
#BSUB -R himem
#BSUB -R rusage[mem=96G] # 120 GB of memory 
#BSUB -R span[hosts=1] # all cores from one node
#BSUB -o %J.stdout # output log (%J : JobID)
#BSUB -eo %J.stderr # error log
#BSUB -L /bin/bash # Initialize the execution environment

ml bcftools
bcftools view -S /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples -Oz -o chr1_topmed_1000g.meta.1748.metaDose.vcf.gz chr1_topmed_1000g.meta.metaDose.vcf.gz
# reorder the sample order to be the same as the order in keep_1748_samples


#!/bin/bash
#BSUB -J vcf_gz_remove_sample # Job name
#BSUB -n 1
#BSUB -P acc_psychgen # allocation account
#BSUB -q premium # queue
#bsub -n 1
#BSUB -W 48:00 # walltime in HH:MM
#BSUB -R himem
#BSUB -R rusage[mem=96G] # 120 GB of memory 
#BSUB -R span[hosts=1] # all cores from one node
#BSUB -o %J.stdout # output log (%J : JobID)
#BSUB -eo %J.stderr # error log
#BSUB -L /bin/bash # Initialize the execution environment

ml bcftools
bcftools view -S /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples -Oz -o chr1.1748.dose.vcf.gz chr1.dose.vcf.gz

bsub < vcf_gz_remove_sample.sh

# ALL Overlap SNPs with Tarseq across 22 chromosomes in GenomeAsia 
# in gsa_1814_mocha_GA folder

# sort the tarseq file (GRCh37)
ml vcftools 
vcf-sort  Pakistani_1748_tarSeq.v4.2.updateID.recode.vcf  >  Pakistani_1748_tarSeq.v4.2.updateID.sort.recode.vcf

# sort the tarseq file (GRCh38)
ml vcftools 
vcf-sort  Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.recode.vcf   > Pakistani_1748_tarSeq.v4.2.updateID.GRCh38.chr.sort.recode.vcf 
 

#############        aggRSquare Comparison     ################

ml git python
ml cmake
/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/aggRSquare/release-build/aggRSquare -v /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.v4.2.updateID.sort.recode.vcf -i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_GA/all_snp_overlap_tarseq_ordered.sort.vcf -o tarseq_GA_grch37  --bins /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_compare/aggRSquare/default_MAF_bins --detail



# in TOPMed

#!/bin/bash
#BSUB -J vcf_gz_remove_sample # Job name
#BSUB -n 1
#BSUB -P acc_psychgen # allocation account
#BSUB -q premium # queue
#bsub -n 1
#BSUB -W 48:00 # walltime in HH:MM
#BSUB -R himem
#BSUB -R rusage[mem=120G] # 120 GB of memory 
#BSUB -R span[hosts=1] # all cores from one node
#BSUB -o %J.stdout # output log (%J : JobID)
#BSUB -eo %J.stderr # error log
#BSUB -L /bin/bash # Initialize the execution environment

ml bcftools
bcftools view -S /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/keep_1748_samples -Oz -o chr1.1748.dose.vcf.gz chr1.dose.vcf.gz

bsub < vcf_gz_remove_sample.sh

### test comparison: imputation accuracy on chromosome 1 
ml git python
ml cmake
git clone https://github.com/yukt/aggRSquare.git
cd aggRSquare
bash install.sh

### meta-imputation

#!/bin/bash
#BSUB -J chr1_aggR2_tarseq_metaimp # Job name
#BSUB -n 1
#BSUB -P acc_psychgen # allocation account
#BSUB -q premium # queue
#BSUB -W 12:00 # walltime in HH:MM
#BSUB -R himem
#BSUB -R rusage[mem=96G] # 120 GB of memory 
#BSUB -R span[hosts=1] # all cores from one node
#BSUB -o %J.stdout # output log (%J : JobID)
#BSUB -eo %J.stderr # error log
#BSUB -L /bin/bash # Initialize the execution environment

ml git python
ml cmake
# meta-imputation
/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/aggRSquare/release-build/aggRSquare -v /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr1.recode.vcf -i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_topmed_1000G/chr1_topmed_1000g.meta.1748.metaDose.vcf.gz -o tarseq_metaimp_chr1_test  --bins /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_compare/aggRSquare/default_MAF_bins --detail
# --validationFormat --imputationFormat
# -v: validation file

bsub < chr1_aggR2_tarseq_metaimp.sh

# expanded 1000G

#!/bin/bash
#BSUB -J chr1_aggR2_tarseq_1000G # Job name
#BSUB -n 1
#BSUB -P acc_psychgen # allocation account
#BSUB -q premium # queue
#BSUB -W 12:00 # walltime in HH:MM
#BSUB -R himem
#BSUB -R rusage[mem=96G] # 120 GB of memory 
#BSUB -R span[hosts=1] # all cores from one node
#BSUB -o %J.stdout # output log (%J : JobID)
#BSUB -eo %J.stderr # error log
#BSUB -L /bin/bash # Initialize the execution environment

ml git python
ml cmake
/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/aggRSquare/release-build/aggRSquare -v /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr1.recode.vcf -i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/chr1.1748.dose.vcf.gz -o tarseq_1000G_chr1_test  --bins /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_compare/aggRSquare/default_MAF_bins --detail

bsub < chr1_aggR2_tarseq_1000G.sh

## TOPMed 

#!/bin/bash
#BSUB -J chr1_aggR2_tarseq_topmed # Job name
#BSUB -n 1
#BSUB -P acc_psychgen # allocation account
#BSUB -q premium # queue
#BSUB -W 12:00 # walltime in HH:MM
#BSUB -R himem
#BSUB -R rusage[mem=96G] # 120 GB of memory 
#BSUB -R span[hosts=1] # all cores from one node
#BSUB -o %J.stdout # output log (%J : JobID)
#BSUB -eo %J.stderr # error log
#BSUB -L /bin/bash # Initialize the execution environment

ml git python
ml cmake
/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/aggRSquare/release-build/aggRSquare -v /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr1.recode.vcf -i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/chr1.1748.dose.vcf.gz -o tarseq_topmed_chr1_test  --bins /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_compare/aggRSquare/default_MAF_bins --detail

bsub < chr1_aggR2_tarseq_topmed.sh


## GenomeAsia 

#!/bin/bash
#BSUB -J chr1_aggR2_tarseq_GA # Job name
#BSUB -n 1
#BSUB -P acc_psychgen # allocation account
#BSUB -q premium # queue
#BSUB -W 12:00 # walltime in HH:MM
#BSUB -R himem
#BSUB -R rusage[mem=96G] # 120 GB of memory 
#BSUB -R span[hosts=1] # all cores from one node
#BSUB -o %J.stdout # output log (%J : JobID)
#BSUB -eo %J.stderr # error log
#BSUB -L /bin/bash # Initialize the execution environment

ml git python
ml cmake
/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/aggRSquare/release-build/aggRSquare -v /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr1.recode.vcf -i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_GA/chr1.1748.dose.vcf.gz -o tarseq_GA_chr1_test  --bins /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_compare/aggRSquare/default_MAF_bins --detail

bsub < chr1_aggR2_tarseq_GA.sh

## 1000G GRCh38 

#!/bin/bash
#BSUB -J chr1_aggR2_tarseq_1000G_grch38 # Job name
#BSUB -n 1
#BSUB -P acc_psychgen # allocation account
#BSUB -q premium # queue
#BSUB -W 12:00 # walltime in HH:MM
#BSUB -R himem
#BSUB -R rusage[mem=96G] # 120 GB of memory 
#BSUB -R span[hosts=1] # all cores from one node
#BSUB -o %J.stdout # output log (%J : JobID)
#BSUB -eo %J.stderr # error log
#BSUB -L /bin/bash # Initialize the execution environment

ml git python
ml cmake
/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/aggRSquare/release-build/aggRSquare -v /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr1.recode.vcf -i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_38/chr1.1748.dose.vcf.gz -o tarseq_1000G_GRCh38_chr1_test  --bins /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_compare/aggRSquare/default_MAF_bins --detail

bsub < chr1_aggR2_tarseq_1000G_grch38.sh

## 1000G GRCh37 

#!/bin/bash
#BSUB -J chr1_aggR2_tarseq_1000G_grch37 # Job name
#BSUB -n 1
#BSUB -P acc_psychgen # allocation account
#BSUB -q premium # queue
#BSUB -W 12:00 # walltime in HH:MM
#BSUB -R himem
#BSUB -R rusage[mem=96G] # 120 GB of memory 
#BSUB -R span[hosts=1] # all cores from one node
#BSUB -o %J.stdout # output log (%J : JobID)
#BSUB -eo %J.stderr # error log
#BSUB -L /bin/bash # Initialize the execution environment

ml git python
ml cmake
/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/aggRSquare/release-build/aggRSquare -v /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.chr1.recode.vcf -i /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_37/chr1.1748.dose.vcf.gz -o tarseq_1000G_GRCh37_chr1_test  --bins /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_compare/aggRSquare/default_MAF_bins --detail

bsub < chr1_aggR2_tarseq_1000G_grch37.sh



# match by position, not by SNP name
# for example
# chr1	6106575	chr1:6106575:C:T	C	T	.	PASS	AF=0.00055;MAF=0.00055;R2=0.31768;IMPUTED;AC=1;AN=3496	GT:DS:HDS:GP	0|0:0.001:0,0.001:0.999,0.001,0	0|0:0:0,0:1,0,0 (1000G imputed)
# chr1	6106575	1:6166635:C:T	C	T	.	.	PR	GT	0/0	0/0 (tarseq)

# manually check how many overlap between tarseq chr1 and 1000G chr1

# in the targeted_seq_data folder
awk '(NR>464){print $1":"$2":"$4":"$5}' Pakistani_1748_tarSeq.GRCh38.chr1.recode.vcf > Pakistani_1748_tarSeq.GRCh38.chr1.recode.vcf.snplist
# 6229 - 464 = 5765 SNPs
awk '(NR>464){print $1":"$2":"$5":"$4}' Pakistani_1748_tarSeq.GRCh38.chr1.recode.vcf > Pakistani_1748_tarSeq.GRCh38.chr1.recode.vcf.flip.snplist

# in the gsa_1814_mocha_1000G folder
zcat chr1.1748.dose.vcf.gz| awk '(NR>25){print $3}' > chr1.1748.dose.vcf.snplist
# 3,898,075 - 25 = 3,898,050 SNPs

comm -12 <(sort /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr1.recode.vcf.snplist) <(sort /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/chr1.1748.dose.vcf.snplist) > overlap_tarseq_1000G_chr1
# 1434 overlap variants
comm -12 <(sort /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr1.recode.vcf.flip.snplist) <(sort /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/chr1.1748.dose.vcf.snplist) > overlap_tarseq_1000G_chr1_flip
# 1 overlap variants when the alleles are flipped

# check the Rsq manually

grep 6110293 /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr1.recode.vcf > tarseq_chr1_6110293
awk '(NR==464){print $0}' /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr1.recode.vcf > tarseq_header
cat tarseq_header tarseq_chr1_6110293 > tarseq_chr1_6110293_header # manually removed the # sign

zcat /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/chr1.1748.dose.vcf.gz | grep chr1:6110293:C:T   > expanded_1000g_chr1_6110293
zcat /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/chr1.1748.dose.vcf.gz | awk '(NR==25){print $0}' > expanded_1000g_header
cat expanded_1000g_header expanded_1000g_chr1_6110293 > expanded_1000g_chr1_6110293_header

# in R
tarseq <- read.table("tarseq_chr1_6110293_header",header=F)
ex1000G <- read.table("expanded_1000g_chr1_6110293_header",header=F)

tarseq <- as.data.frame(t(as.matrix(tarseq)))
tarseq <- tarseq[-c(1:9),]
table(tarseq$V2)
# ./.  0/0  0/1  1/1 
#   8 1247  392  101 
# 8 missing genotype, 1740 have genotype data

tarseq_nomiss <- tarseq[which(tarseq$V2 != "./."),] # 1740 

library(dplyr)
tarseq_nomiss$geno <- case_when(tarseq_nomiss$V2 =="0/0" ~ 0, tarseq_nomiss$V2 =="0/1" ~ 1,tarseq_nomiss$V2 =="1/1" ~ 2)
table(tarseq_nomiss$geno,tarseq_nomiss$V2)
# AF = (392+101*2)/(1740*2) = 594/3480 = 0.170690
id_keep <- tarseq_nomiss[,1]

ex1000G <- as.data.frame(t(as.matrix(ex1000G)))
ex1000G <- ex1000G[-c(1:9),]

ex1000G_nomiss <- ex1000G[which(ex1000G$V1 %in% id_keep),] # 1740 samples

library(stringr)
# Split name column into firstname and last name
ex1000G_nomiss[c('GT', 'DS', 'HDS', 'GP')] <- str_split_fixed(ex1000G_nomiss$V2, ':', 4)
ex1000G_nomiss <- ex1000G_nomiss[,c("V1","DS")]
sum(as.numeric(ex1000G_nomiss$DS))/(1740*2)
# AF = 0.1590147
merge_tarseq_ex1000G <- merge(tarseq_nomiss,ex1000G_nomiss, by="V1" )


cor(as.numeric(merge_tarseq_ex1000G$geno),as.numeric(merge_tarseq_ex1000G$DS), method = "pearson")
#  pearson =    0.9629336
# squared pearson =  0.9629336 *  0.9629336 = 0.9272411 (matched with the Imputation.R2 value)


### check a second SNP

grep 6106575 /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/targeted_seq_data/Pakistani_1748_tarSeq.GRCh38.chr1.recode.vcf > tarseq_chr1_6106575
cat tarseq_header tarseq_chr1_6106575 > tarseq_chr1_6106575_header # manually removed the # sign

zcat /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G/chr1.1748.dose.vcf.gz | grep chr1:6106575:C:T  > expanded_1000g_chr1_6106575
cat expanded_1000g_header expanded_1000g_chr1_6106575 > expanded_1000g_chr1_6106575_header

# in R
tarseq <- read.table("tarseq_chr1_6106575_header",header=F)
ex1000G <- read.table("expanded_1000g_chr1_6106575_header",header=F)

tarseq <- as.data.frame(t(as.matrix(tarseq)))
tarseq <- tarseq[-c(1:9),]
table(tarseq$V2)
# ./.  0/0  0/1  1/1 
#   2 1745    1 
# 2 missing genotype, 1746 have genotype data

tarseq_nomiss <- tarseq[which(tarseq$V2 != "./."),] # 1746

library(dplyr)
tarseq_nomiss$geno <- case_when(tarseq_nomiss$V2 =="0/0" ~ 0, tarseq_nomiss$V2 =="0/1" ~ 1,tarseq_nomiss$V2 =="1/1" ~ 2)
table(tarseq_nomiss$geno,tarseq_nomiss$V2)
# AF = 1/(1746*2) = 0.000286
id_keep <- tarseq_nomiss[,1]

ex1000G <- as.data.frame(t(as.matrix(ex1000G)))
ex1000G <- ex1000G[-c(1:9),]

ex1000G_nomiss <- ex1000G[which(ex1000G$V1 %in% id_keep),] # 1740 samples

library(stringr)
# Split name column into firstname and last name
ex1000G_nomiss[c('GT', 'DS', 'HDS', 'GP')] <- str_split_fixed(ex1000G_nomiss$V2, ':', 4)
ex1000G_nomiss <- ex1000G_nomiss[,c("V1","DS")]
sum(as.numeric(ex1000G_nomiss$DS))/(1746*2)
# AF = 0.0005506873
merge_tarseq_ex1000G <- merge(tarseq_nomiss,ex1000G_nomiss, by="V1" )


cor(as.numeric(merge_tarseq_ex1000G$geno),as.numeric(merge_tarseq_ex1000G$DS), method = "pearson")
#  pearson =   -0.001358449
# squared pearson = -0.001358449 * -0.001358449 = 1.845384e-06 ~ 0.000002  (matched with the Imputation.R2 value)

# yes the R2 calculation are correct in aggRSquare (some minor deviation from manual calculation)

##################################################
########      plot the true Rsq by MAF   #########
##################################################

# manually remove the # in the header line
sed '/^#/d' tarseq_metaimp_chr1_test.aggRSquare > tarseq_metaimp_chr1_test.aggRSquare.noheader
sed '/^#/d' tarseq_1000G_chr1_test.aggRSquare > tarseq_1000G_chr1_test.aggRSquare.noheader

sed '/^#/d' tarseq_topmed_chr1_test.aggRSquare > tarseq_topmed_chr1_test.aggRSquare.noheader
sed '/^#/d' tarseq_GA_chr1_test.aggRSquare > tarseq_GA_chr1_test.aggRSquare.noheader
sed '/^#/d' tarseq_1000G_GRCh38_chr1_test.aggRSquare > tarseq_1000G_GRCh38_chr1_test.aggRSquare.noheader
sed '/^#/d' tarseq_1000G_GRCh37_chr1_test.aggRSquare > tarseq_1000G_GRCh37_chr1_test.aggRSquare.noheader

# in the path: /sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_compare/aggRSquare
metaimp <- read.table("tarseq_metaimp_chr1_test.aggRSquare.noheader",header=T)
ex1000G <- read.table("tarseq_1000G_chr1_test.aggRSquare.noheader",header=T)

topmed <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_topmed/tarseq_topmed_chr1_test.aggRSquare.noheader",header=T)
GA <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_GA/tarseq_GA_chr1_test.aggRSquare.noheader",header=T)
all_1000G_38 <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_38/tarseq_1000G_GRCh38_chr1_test.aggRSquare.noheader",header=T)
sas_1000G_37 <- read.table("/sc/arion/projects/psychgen/HUCKINS_LAB_DONT_DELETE/jiayi_folders/Pakistani_GSA_impute/gsa_1814_mocha_1000G_OG_37/tarseq_1000G_GRCh37_chr1_test.aggRSquare.noheader",header=T)

# update the SNP ID to the same as tarseq ID (check forward and reverse)


metaimp$ref_panel <- "TOPMed_1000G"
topmed$ref_panel <- "TOPMed"
ex1000G$ref_panel <- "Expanded_1000G"
all_1000G_38$ref_panel <- "1000G_GRCh38"
sas_1000G_37$ref_panel <- "1000G_GRCh37_SAS"
GA$ref_panel <- "GenomeAsia_Pilot"

merged <- rbind(metaimp,topmed,ex1000G,all_1000G_38,sas_1000G_37,GA) #  84 rows, 7 columns
merged$ref_panel <- factor(merged$ref_panel,  levels=c("TOPMed_1000G","TOPMed","Expanded_1000G","1000G_GRCh38", "1000G_GRCh37_SAS", "GenomeAsia_Pilot"))

# color blind friendly
colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
scales::show_col(colorBlindGrey8)

library(ggplot2)
library(dplyr)
p <- ggplot(merged, aes(x = Bin.Aggregated.by.MAF, y = Imputation.R2,group=ref_panel, color=ref_panel))  + geom_line() + geom_point(shape=1)+ theme_classic(base_size=20) + coord_cartesian(ylim = c(0, 1)) + scale_x_discrete(guide = guide_axis(angle = 45))  + theme(text = element_text(size = 15,colour="black"))  + geom_hline(yintercept=0.8, linetype="dashed", color = "red") #all font size
p + scale_color_manual(values=c( "#0072B2","#E68F00", "#56B4E9","#CC79A7","#999999","#F0E442"),name="Reference Panel") + labs(x="Minor Allele Frequency", y="Average True R2")
#009E73

merged_meta_1000g_topmed <- merged[which(merged$ref_panel %in% c("TOPMed_1000G","TOPMed","Expanded_1000G")),]
p <- ggplot(merged_meta_1000g_topmed, aes(x = Bin.Aggregated.by.MAF, y = Imputation.R2,group=ref_panel, color=ref_panel))  + geom_line() + geom_point(shape=1)+ theme_classic(base_size=20) + coord_cartesian(ylim = c(0, 1)) + scale_x_discrete(guide = guide_axis(angle = 45))  + theme(text = element_text(size = 15,colour="black"))  + geom_hline(yintercept=0.8, linetype="dashed", color = "red") #all font size
p + scale_color_manual(values=c( "#0072B2","#E68F00", "#56B4E9"),name="Reference Panel") + labs(x="Minor Allele Frequency", y="Average True R2")

# how the aggregated Rsq is calculated
aggR <- read.table("tarseq_metaimp_chr1_test.aggRSquare.noheader",header=T)
indR <- read.table("tarseq_metaimp_chr1_test.RSquare", header=F)

colnames(indR) <- c("SNP.ID","Allele.Frequency","No.Samples","Imputation.R2","Validation.AF","Imputation.AF")

indR$"Allele.Frequency" <- as.numeric(indR$"Allele.Frequency")

for (i in 1:dim(indR)[1]) {
  if (indR$"Allele.Frequency"[i] <=0.5) {
    indR$MAF[i] <- indR$"Allele.Frequency"[i]
  }
  else {indR$MAF[i] <- 1- indR$"Allele.Frequency"[i] }
}

maf0_05 <- indR[which(indR$MAF <= 0.0005),] # 345 (aggR: 345)
maf0_1 <- indR[which(indR$MAF <= 0.001 & indR$MAF > 0.0005),] # 178  (aggR: 178)
maf0_2 <- indR[which(indR$MAF <= 0.002 & indR$MAF > 0.001),] # 75  (aggR: 75)
maf0_5 <- indR[which(indR$MAF <= 0.005 & indR$MAF> 0.002),] # 119 (aggR: 119) 
maf1 <- indR[which(indR$MAF <= 0.01 & indR$MAF > 0.005 ),] # 63  (aggR: 63)
maf1_5 <- indR[which(indR$MAF <= 0.015 & indR$MAF > 0.01),] # 16  (aggR: 16)
maf2 <- indR[which(indR$MAF <= 0.02 & indR$MAF > 0.015 ),] # 12  (aggR: 12)
maf3_5 <- indR[which(indR$MAF <= 0.035 & indR$MAF > 0.02 ),] # 25  (aggR: 25)
maf5 <- indR[which(indR$MAF <= 0.05 & indR$MAF> 0.035  ),] # 20  (aggR: 20)
maf10 <- indR[which(indR$MAF <= 0.1 &  indR$MAF > 0.05),] # 21  (aggR: 21)
maf20 <- indR[which(indR$MAF <= 0.2 &  indR$MAF > 0.1),] # 57  (aggR: 57) 
maf30 <- indR[which(indR$MAF <= 0.3 & indR$MAF > 0.2 ),] # 43  (aggR: 43)  
maf40 <- indR[which(indR$MAF <= 0.4 &  indR$MAF  > 0.3),] # 30  (aggR: 30)  
maf50 <- indR[which(indR$MAF <= 0.5 &  indR$MAF  > 0.4),] # 21 (aggR: 21)  

dim(maf0_05) 
dim(maf0_1)
dim(maf0_2)
dim(maf0_5)
dim(maf1)
dim(maf1_5)
dim(maf2)
dim(maf3_5)
dim(maf5)
dim(maf10)
dim(maf20)
dim(maf30)
dim(maf40)
dim(maf50)
 
# question: number of SNPs, calculation of Pearson R2, calculation of aggregate R square
 
sum(maf0_05$"Imputation.R2")/dim(maf0_05)[1] # average R2 = 0.1817488 (aggR file R2 =  0.071575)
sum(maf0_1$"Imputation.R2")/dim(maf0_1)[1] # average R2 = 0.2909237 (aggR file R2 =  0.188234)
sum(maf0_2$"Imputation.R2")/dim(maf0_2)[1] # average R2 =  0.3410112 (aggR file R2 =  0.246433)
sum(maf0_5$"Imputation.R2")/dim(maf0_5)[1] # average R2 =  0.5722874 (aggR file R2 =  0.495111)
sum(maf1$"Imputation.R2")/dim(maf1)[1] # average R2 = 0.6593659 (aggR file R2 =  0.657491)
sum(maf1_5$"Imputation.R2")/dim(maf1_5)[1] # average R2 = 0.7415713 (aggR file R2 =  0.730690)
sum(maf2$"Imputation.R2")/dim(maf2)[1] # average R2 = 0.704127 (aggR file R2 =   0.690520)
sum(maf3_5$"Imputation.R2")/dim(maf3_5)[1] # average R2 = 0.8163102 (aggR file R2 =  0.763792)
sum(maf5$"Imputation.R2")/dim(maf5)[1] # average R2 = 0.9659496 (aggR file R2 =   0.962478)
sum(maf10$"Imputation.R2")/dim(maf10)[1] # average R2 = 0.8815348 (aggR file R2 =  0.889672)
sum(maf20$"Imputation.R2")/dim(maf20)[1] # average R2 = 0.906451 (aggR file R2 = 0.901500)
sum(maf30$"Imputation.R2")/dim(maf30)[1] # average R2 = 0.9550922 (aggR file R2 =  0.924139)
sum(maf40$"Imputation.R2")/dim(maf40)[1] # average R2 = 0.9484828 (aggR file R2 =  0.952891)
sum(maf50$"Imputation.R2")/dim(maf50)[1] # average R2 =  0.9338291 (aggR file R2 =  0.884489)


sum(maf0_05$MAF)/dim(maf0_05)[1] # average MAF = 0.0002862754 (aggR average MAF = 0.000286)
sum(maf0_1$MAF)/dim(maf0_1)[1] # average MAF = 0.0006725506 (aggR average MAF = 0.000673)
sum(maf0_2$MAF)/dim(maf0_2)[1] # average MAF = 0.00140276 (aggR average MAF = 0.001403)
sum(maf0_5$MAF)/dim(maf0_5)[1] # average MAF = 0.003143185 (aggR average MAF = 0.003143) 
sum(maf1$MAF)/dim(maf1)[1] # average MAF = 0.006572524  (aggR average MAF = 0.006572) 
sum(maf1_5$MAF)/dim(maf1_5)[1] # average MAF = 0.01282763 (aggR average MAF =  0.012828)
sum(maf2$MAF)/dim(maf2)[1] # average MAF =  0.01725983 (aggR average MAF = 0.017260)
sum(maf3_5$MAF)/dim(maf3_5)[1] # average MAF = 0.02566952 (aggR average MAF =  0.025670)
sum(maf5$MAF)/dim(maf5)[1] # average MAF = 0.04476265 (aggR average MAF = 0.044763)
sum(maf10$MAF)/dim(maf10)[1] # average MAF =  0.07260005 (aggR average MAF = 0.072600)
sum(maf20$MAF)/dim(maf20)[1] # average MAF =   0.1446488 (aggR average MAF = 0.144649)
sum(maf30$MAF)/dim(maf30)[1] # average MAF = 0.2461621 (aggR average MAF =  0.246162) 
sum(maf40$MAF)/dim(maf40)[1] # average MAF = 0.3597733 (aggR average MAF = 0.359773) 
sum(maf50$MAF)/dim(maf50)[1] # average MAF =  0.445615 (aggR average MAF = 0.445615) 

# Bin.Aggregated.by.MAF : MAF bin in the format of (min, max]
# Average.MAF : the average MAF of all the variants in the bin