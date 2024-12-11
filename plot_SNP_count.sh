#############################################################################################
####         R script:    plot number of SNPs and well-imputed SNPs (rare, common)       ####
#############################################################################################

# get the numbers first
library(data.table)
dat <- as.data.frame(fread("all.chr.info", select = c("SNP","MAF","Rsq","Genotyped")))
dat <- dat[which(dat$Genotyped!="Typed_Only"),]

rare <- dat[which(dat$MAF < 0.01 & dat $MAF > 0),] # including monomorphic SNPs: 35,213,874 SNPs, excluding monomorphic SNPs: 30,442,265 SNPs
common <- dat[which(dat$MAF >= 0.01 ),]

dim(dat)
dim(rare)
dim(common)
# rare + common add up to total

dim(rare[which(rare$Rsq >=0.8),])
dim(common[which(common$Rsq >=0.8),])


# plot the numbers
snp_count <- read.table("snp_count_wo_monomorphic",header=T)
table(snp_count$ref_panel)
table(snp_count$maf_cat)
table(snp_count$snp_cat)

snp_count$ref_panel <- factor(snp_count$ref_panel,  levels=c("TOPMed_1000G","TOPMed","Expanded_1000G","1000G_GRCh38", "1000G_GRCh37_SAS", "GenomeAsia_Pilot"))
snp_count$maf_cat <- factor(snp_count$maf_cat,  levels=c("common","rare"))

snp_count_all <- snp_count[which(snp_count$snp_cat=="all_snps"),]
snp_count_well <- snp_count[which(snp_count$snp_cat=="well_imputed_snps"),]


library(ggplot2)
pdf('all_snp_count.pdf',height=10,width=15)
p <- ggplot(data=snp_count_all, aes(x=maf_cat, y=count, fill=ref_panel)) + xlab("Variant Type") + ylab("Count") + geom_bar(stat="identity", position=position_dodge()) +  geom_text(aes(label=scales::label_number_si(accuracy=0.01)(count)),size=4, position=position_dodge(width=0.9), vjust=-0.25) + theme_classic(base_size=30)  + scale_x_discrete(labels=c("Common","Rare")) + scale_y_continuous(labels=scales::label_number_si())+theme(axis.text = element_text(size = 30,colour="black")) 
p + scale_fill_manual(values=c("#661100", "#D55E00", "#E69F00","#0072B2","#332288","#999999"),name="Reference Panel",labels=c("TOPMed_1000G","TOPMed","Expanded_1000G","1000G_GRCh38", "1000G_GRCh37_SAS", "GenomeAsia_Pilot")) + theme(legend.position = "none" ,plot.title = element_text(hjust = 0.5)) + labs(title="All SNPs")
dev.off()

pdf('well_imputed_snp_count.pdf',height=10,width=15)
p <- ggplot(data=snp_count_well, aes(x=maf_cat, y=count, fill=ref_panel)) + xlab("Variant Type") + ylab("Count") + geom_bar(stat="identity", position=position_dodge()) +  geom_text(aes(label=scales::label_number_si(accuracy=0.01)(count)),size=4, position=position_dodge(width=0.9), vjust=-0.25) + theme_classic(base_size=30)  + scale_x_discrete(labels=c("Common","Rare")) + scale_y_continuous(labels=scales::label_number_si())+theme(axis.text = element_text(size = 30,colour="black")) 
p + scale_fill_manual(values=c("#661100", "#D55E00", "#E69F00","#0072B2","#332288","#999999"),name="Reference Panel",labels=c("TOPMed_1000G","TOPMed","Expanded_1000G","1000G_GRCh38", "1000G_GRCh37_SAS", "GenomeAsia_Pilot")) + theme(legend.position = "none",plot.title = element_text(hjust = 0.5)) + labs(title="Well-imputed SNPs")
dev.off()