####################################################################
####                      R script: R2 by MAF                   ####
####################################################################


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
   dat <- subset(imputed, MAF > 0 & MAF <= 0.0005)
   mean_r2[i,1] <- "(0,0.0005]"
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
