###########################################################################################
########     R script: compare true Rsq across imputation panels at each MAF      #########
###########################################################################################

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


# Non-parametric alternative to one-way ANOVA test
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_1) # P  < 2.2e-16 (**) (probably due to topmed has way many low-quality rare variants imputed, therefore drag down the average R2)
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_2) # P  < 2.2e-16 (**)
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_3) # P  0.0001138 (**)
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_4) # P  0.00164 (**)
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_5) # P   0.02117 (*)
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_6) # P  0.5666 
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_7) # P  0.877
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_8) # P  0.6122
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_9) # P  0.894
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_10) # P  0.8434
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_11) # P  0.4808
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_12) # P  0.7344
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_13) # P 0.839
kruskal.test(Imputation_R2 ~ ref_panel, data = merged_14) # P 0.7879 
# 0.05/14 = 0.003571 (** Bonferroni significant )
# 0.05 (* nominal significant)

# check for pairwise difference
to_save <- true_r2_maf  %>%  arrange(Imputation_R2) %>%  arrange( MAF_bin)
write.table(to_save,"figure_1a_true_r2_by_maf_by_panel.txt", quote=F, col=T, row=F,sep="\t")


pairwise.wilcox.test(merged_1$Imputation_R2, merged_1$ref_panel, p.adjust.method = "bonferroni")
pairwise.wilcox.test(merged_2$Imputation_R2, merged_2$ref_panel, p.adjust.method = "bonferroni")
pairwise.wilcox.test(merged_3$Imputation_R2, merged_3$ref_panel, p.adjust.method = "bonferroni")
pairwise.wilcox.test(merged_4$Imputation_R2, merged_4$ref_panel, p.adjust.method = "bonferroni")
pairwise.wilcox.test(merged_5$Imputation_R2, merged_5$ref_panel, p.adjust.method = "bonferroni")
pairwise.wilcox.test(merged_6$Imputation_R2, merged_6$ref_panel, p.adjust.method = "bonferroni")
pairwise.wilcox.test(merged_7$Imputation_R2, merged_7$ref_panel, p.adjust.method = "bonferroni")
pairwise.wilcox.test(merged_8$Imputation_R2, merged_8$ref_panel, p.adjust.method = "bonferroni")
pairwise.wilcox.test(merged_9$Imputation_R2, merged_9$ref_panel, p.adjust.method = "bonferroni")
pairwise.wilcox.test(merged_10$Imputation_R2, merged_10$ref_panel, p.adjust.method = "bonferroni")
pairwise.wilcox.test(merged_11$Imputation_R2, merged_11$ref_panel, p.adjust.method = "bonferroni")
pairwise.wilcox.test(merged_12$Imputation_R2, merged_12$ref_panel, p.adjust.method = "bonferroni")
pairwise.wilcox.test(merged_13$Imputation_R2, merged_13$ref_panel, p.adjust.method = "bonferroni")
pairwise.wilcox.test(merged_14$Imputation_R2, merged_14$ref_panel, p.adjust.method = "bonferroni")