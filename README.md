# File Dictionary

1. R2_and_wilcox_test_statistics.xlsx
   
This spreadsheet includes the imputation accuracy represented by different types of R2 (squared correlation) (e.g. true R2, adjusted true R2, estimated R2) for each imputation reference panel at each bin of minor allele frequency (MAF) in the Xu et al. 2024 paper (see below for citation). In addition, this spreadsheet includes pairwise comparison of imputation accuracy (i.e., measured by true R2, adjusted true R2, and/or estimated R2) across different imputation panels at each MAF bin, with significance tested by Wilcoxon rank-sum test (p-values are shown). 

2. R2_by_MAF.sh

This file contains the R script to calculate average R2 for each MAF bin.

3. calculation_adjusted_true_R2.sh

This file contains the R script to calculate adjusted true R2 for each MAF bin for each imputation panel.

4. compare_Rsq_across_imputation_panels.sh

This file contains the R script to compare Rsq across different imputation panels at each MAF bin using Kruskal-Wallis rank sum test for overall significance and Wilcoxon rank-sum test for pairwise comparison.

5. find_SNV_indel.sh

This file contains the R script to classify SNPs into single nucleotide variant versus insertion/deletion variant.

6. find_multi_and_bi_allelic_SNP.sh

This file contains the R script to classify SNPs into biallelic variant versus multiallelic variant.

7. mocha_idat_to_vcf.sh

This file contains the Google Cloud script to covert Illumina idat files to vcf files for the genotype data.

8. plot_Rsq_by_MAF.sh

This file contains the R script to plot the R2 across the spectrum of MAF for all imputation panels.

9. plot_SNP_count.sh

This file contains the R script to count and plot the SNP counts.

10. snp_annotation.sh

This file contains the R script to annotate the SNPs with different variant type (e.g., introns, extrons, missense).


# Citation
1. Xu J, Liu D, Hassan A, Genovese G, Cote AC, Fennessy B, Cheng E, Charney AW, Knowles JA, Ayub M, Peterson RE, Bigdeli TB, Huckins LM. Evaluation of imputation performance of multiple reference panels in a Pakistani population. Human Genetics and Genomics Advances 2024 (in print).
   
