################################################################################################
####           R script:  Classify SNPs as multi-allelic or bi-allelic (rare, common)       ####
################################################################################################


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
snp_type$counts[3] <- dim(imputed[which(imputed$dup==1 & imputed$MAF < 0.01 & imputed$MAF > 0),])[1]
snp_type$counts[4] <- dim(imputed[which(imputed$dup==0 & imputed$MAF < 0.01 & imputed$MAF > 0),])[1]
snp_type$counts[5] <- dim(imputed[which(imputed$dup==1 & imputed$MAF >= 0.01),])[1]
snp_type$counts[6] <- dim(imputed[which(imputed$dup==0 & imputed$MAF >= 0.01),])[1]


snp_type$counts_high_r2[1] <- dim(imputed[which(imputed$dup==1 & imputed$Rsq >= 0.8),])[1]
snp_type$counts_high_r2[2] <- dim(imputed[which(imputed$dup==0 & imputed$Rsq >= 0.8),])[1]
snp_type$counts_high_r2[3] <- dim(imputed[which(imputed$dup==1 & imputed$MAF < 0.01 & imputed$MAF > 0 & imputed$Rsq >= 0.8),])[1]
snp_type$counts_high_r2[4] <- dim(imputed[which(imputed$dup==0 & imputed$MAF < 0.01 & imputed$MAF > 0 & imputed$Rsq >= 0.8),])[1]
snp_type$counts_high_r2[5] <- dim(imputed[which(imputed$dup==1 & imputed$MAF >= 0.01 & imputed$Rsq >= 0.8),])[1]
snp_type$counts_high_r2[6] <- dim(imputed[which(imputed$dup==0 & imputed$MAF >= 0.01 & imputed$Rsq >= 0.8),])[1]

snp_type$prop_high_r2 <- snp_type$counts_high_r2 / snp_type$counts * 100

write.table(snp_type,"snp_type_multi_allele_bi_allele.txt", quote=F, row.names=F)
