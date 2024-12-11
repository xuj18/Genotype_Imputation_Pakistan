##############################################################
########      R script: plot the true Rsq by MAF     #########
##############################################################

metaimp$ref_panel <- "TOPMed_1000G"
topmed$ref_panel <- "TOPMed"
ex1000G$ref_panel <- "Expanded_1000G"
all_1000G_38$ref_panel <- "1000G_GRCh38"
sas_1000G_37$ref_panel <- "1000G_GRCh37_SAS"
GA$ref_panel <- "GenomeAsia_Pilot"

dim(metaimp) # 8483
dim(topmed) # 8267
dim(ex1000G) # 5038
dim(all_1000G_38) # 4905
dim(sas_1000G_37) # 5072
dim(GA) # 3423

library(dplyr)
# remove duplicates
dup_meta <- metaimp[duplicated(metaimp$V1),1]  # 6 duplicated SNPs
length(dup_meta)
metaimp_dedup <- metaimp[which(!metaimp$V1 %in% dup_meta),] # 8483 - 12 = 8471 SNPs

dup_topmed <- topmed[duplicated(topmed$V1),1]  # 5 duplicated SNPs
length(dup_topmed)
topmed_dedup <- topmed[which(!topmed$V1 %in% dup_topmed),] # 8267 - 10 = 8257 SNPs

dup_ex1000G <- ex1000G[duplicated(ex1000G$V1),1]  # 1 duplicated SNP
length(dup_ex1000G)
ex1000G_dedup <- ex1000G[which(!ex1000G$V1 %in% dup_ex1000G),] # 5038 - 2 = 5036 SNPs

dup_1000G_38 <- all_1000G_38[duplicated(all_1000G_38$V1),1]  # 0
length(dup_1000G_38)
all_1000G_38_dedup <- all_1000G_38[which(!all_1000G_38$V1 %in% dup_1000G_38),] # 4905 - 0 = 4905 SNPs

dup_1000G_37 <- sas_1000G_37[duplicated(sas_1000G_37$V1),1]  # 0
length(dup_1000G_37)
sas_1000G_37_dedup <- sas_1000G_37[which(!sas_1000G_37$V1 %in% dup_1000G_37),] # 5072 - 0 = 5072 SNPs

dup_GA <- GA[duplicated(GA$V1),1]  # 0
length(dup_GA)
GA_dedup <- GA[which(!GA$V1 %in% dup_GA),] # 3423 - 0 = 3423 SNPs

# all add up to 35164

# check the MAF of duplicated SNPs (column: v2)

# meta-imputation: 6 SNPs (1 common, 5 rare)
metaimp[which(metaimp$V1 %in% dup_meta),] 

# TOPMed: 5 SNPs (1 common, 4 rare)
topmed[which(topmed$V1 %in% dup_topmed),] 

# expanded 1000G: 1 SNP (1 rare)
ex1000G[which(ex1000G$V1 %in% dup_ex1000G),] 



# merge across all panels
merged <- rbind(metaimp_dedup,topmed_dedup,ex1000G_dedup,all_1000G_38_dedup,sas_1000G_37_dedup,GA_dedup) #  35164 rows 
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
summary(merged$MAF) # MAF ranges from 0.000286 to 0.5

# what is the second smallest MAF: 0.000287
# library(tidyverse)
# test <- unique(merged$MAF)
# test[order(test)]
# Validation_AF: the allele frequency calculated from genotype data in the input validation file

library(dplyr)
merged$MAF_bin <- case_when(merged$MAF > 0 & merged$MAF <= 0.0005 ~ "(0,0.0005]", merged$MAF > 0.0005 & merged$MAF <= 0.001 ~ "(0.0005,0.001]", merged$MAF > 0.001 & merged$MAF <= 0.002 ~ "(0.001,0.002]", merged$MAF > 0.002 & merged$MAF <= 0.005 ~ "(0.002,0.005]",merged$MAF > 0.005 & merged$MAF <= 0.01 ~ "(0.005,0.01]",merged$MAF > 0.01 & merged$MAF <= 0.015 ~ "(0.01,0.015]",merged$MAF > 0.015 & merged$MAF <= 0.02 ~ "(0.015,0.02]", merged$MAF > 0.02 & merged$MAF <= 0.035 ~ "(0.02,0.035]", merged$MAF > 0.035 & merged$MAF <= 0.05 ~ "(0.035,0.05]", merged$MAF > 0.05 & merged$MAF <= 0.1 ~ "(0.05,0.1]", merged$MAF > 0.1 & merged$MAF <= 0.2  ~ "(0.1,0.2]", merged$MAF > 0.2 & merged$MAF <= 0.3  ~ "(0.2,0.3]", merged$MAF > 0.3 & merged$MAF <= 0.4  ~ "(0.3,0.4]", merged$MAF > 0.4 & merged$MAF <= 0.5  ~ "(0.4,0.5]")
table(merged$MAF_bin)
table(merged$MAF_bin,merged$ref_panel)
table(merged$ref_panel)
# install.packages("plotrix")
# library(plotrix)
true_r2_maf <- aggregate(Imputation_R2~ref_panel+MAF_bin, data=merged, mean)
# FUN = function(x) c(mean = mean(x), se = std.error(x))
true_r2_maf$ref_panel <- factor(true_r2_maf$ref_panel,  levels=c("TOPMed_1000G","TOPMed","Expanded_1000G","1000G_GRCh38", "1000G_GRCh37_SAS", "GenomeAsia_Pilot"))
true_r2_maf$MAF_bin <- factor(true_r2_maf$MAF_bin,  levels=c("(0,0.0005]","(0.0005,0.001]","(0.001,0.002]","(0.002,0.005]","(0.005,0.01]","(0.01,0.015]","(0.015,0.02]","(0.02,0.035]","(0.035,0.05]","(0.05,0.1]","(0.1,0.2]","(0.2,0.3]","(0.3,0.4]","(0.4,0.5]"))

write.table(true_r2_maf,"average_true_r2_maf_bin_sequenced_snps",row=FALSE,quote=FALSE)


library(ggplot2)
library(dplyr)
pdf('true_r2_maf_all_snp.pdf',height=10,width=15) # Figure 1a
p <- ggplot(true_r2_maf, aes(x = MAF_bin, y = Imputation_R2,group=ref_panel, color=ref_panel))  + geom_line(size=1) + geom_point(shape=1)+ theme_classic(base_size=30) + coord_cartesian(ylim = c(0, 1)) + scale_x_discrete(guide = guide_axis(angle = 45))  + theme(text = element_text(size = 30,colour="black"))  + geom_hline(yintercept=0.8,  color = "dark grey") + geom_hline(yintercept=0.6, linetype="dashed", color = "dark grey")  + geom_vline(xintercept="(0.01,0.015]",  color = "dark grey",linetype="dashed") + geom_vline(xintercept="(0.05,0.1]",  color = "dark grey") #all font size
p + scale_color_manual(values=c("#661100", "#D55E00", "#E69F00","#0072B2","#332288","#999999"),name="Reference Panel") + labs(x="Minor Allele Frequency", y="Average True R2") + 
  annotate('text', x = "(0,0.0005]", y = 0.45, label='"**"', parse=TRUE, color = "red",size=10) +  
  annotate('text', x = "(0.0005,0.001]", y = 0.4, label='"**"', parse=TRUE, color = "red",size=10) +
  annotate('text', x = "(0.001,0.002]", y = 0.45, label='"**"', parse=TRUE, color = "red",size=10) +
  annotate('text', x = "(0.002,0.005]", y = 0.6, label='"**"', parse=TRUE, color = "red",size=10) +
  annotate('text', x = "(0.005,0.01]", y = 0.7, label='"*"', parse=TRUE, color = "red",size=10) 
dev.off()
