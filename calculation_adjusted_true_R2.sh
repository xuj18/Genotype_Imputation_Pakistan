#####################################################################################
########              R script: Calculation of adjusted true R2             #########
#####################################################################################


# combine true R2 dataset and estimated R2 dataset for sequenced SNPs
dat <- merge(true_r2_maf_seq,est_r2_maf_seq, by=c("ref_panel","MAF_bin"))

dat <- merge(dat,est_r2_maf_gw,by=c("ref_panel","MAF_bin")) # 14*6 = 84 rows 

dat$ratio <- dat$Imputation_R2/dat$Estimated_R2

write.table(dat,"ratio_true_r2_vs_est_r2_seq_data",sep='\t', row.names=FALSE,quote=FALSE)

# calculate the adjusted true R2 for each imputation panel at each MAF bin

metaimp <- as.data.frame(metaimp)
topmed <- as.data.frame(topmed)
ex1000G <- as.data.frame(ex1000G)
all_1000G_38 <- as.data.frame(all_1000G_38)
sas_1000G_37 <- as.data.frame(sas_1000G_37)
GA <- as.data.frame(GA)
dat <- as.data.frame(dat)

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

dim(metaimp) 
dim(topmed) 
dim(ex1000G) 
dim(all_1000G_38)  
dim(sas_1000G_37) 
dim(GA)  


ratio_metaimp      <- dat[which(dat$MAF_bin=="(0.0005,0.001]" & dat$ref_panel=="TOPMed_1000G"),"ratio"]
ratio_topmed       <- dat[which(dat$MAF_bin=="(0.0005,0.001]" & dat$ref_panel=="TOPMed"),"ratio"]
ratio_ex1000G      <- dat[which(dat$MAF_bin=="(0.0005,0.001]" & dat$ref_panel=="Expanded_1000G"),"ratio"]
ratio_all_1000G_38 <- dat[which(dat$MAF_bin=="(0.0005,0.001]" & dat$ref_panel=="1000G_GRCh38"),"ratio"]
ratio_sas_1000G_37 <- dat[which(dat$MAF_bin=="(0.0005,0.001]" & dat$ref_panel=="1000G_GRCh37_SAS"),"ratio"]
ratio_GA           <- dat[which(dat$MAF_bin=="(0.0005,0.001]" & dat$ref_panel=="GenomeAsia_Pilot"),"ratio"]


metaimp$adjusted_true_r2      <- metaimp$Rsq * ratio_metaimp
topmed$adjusted_true_r2       <- topmed$Rsq * ratio_topmed
ex1000G$adjusted_true_r2      <- ex1000G$Rsq * ratio_ex1000G
all_1000G_38$adjusted_true_r2 <- all_1000G_38$Rsq * ratio_all_1000G_38
sas_1000G_37$adjusted_true_r2 <- sas_1000G_37$Rsq * ratio_sas_1000G_37
GA$adjusted_true_r2           <- GA$Rsq * ratio_GA

merged <- rbind(metaimp,topmed,ex1000G,all_1000G_38,sas_1000G_37,GA)  # [1] 6243917       5

kruskal.test(adjusted_true_r2 ~ ref_panel, data = merged)
pairwise.wilcox.test(merged$adjusted_true_r2, g = merged$ref_panel, p.adjust.method = "bonferroni")

