###############################################################################
####           R script:  Classify SNPs into SNVs and indels               ####
###############################################################################

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