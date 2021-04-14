## Additional analyses for male vs female PD

```
Containing:
1. LDSCORE on full stats...
2. LDSCORE on stats without UKB proxy...
3. Polygenic risk scoring in NEUROX_DBGAP only...
4. Polygenic risk scoring in all datasets...
5. Test differences of beta's...
6. Male/Female specific polygenic risk scoring in UKB...

```


### 1. LDSCORE on full stats...

```
# main purpose of this:
# to check for genome wide correlation between male and female PD summary statistics

cd /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/RESULTS/to_meta_files/

mkdir LDSCORE
cd LDSCORE

module load ldsc/1.0.0-101-g89c13a7
# get reference data
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
tar -jxvf eur_w_ld_chr.tar.bz2

### copy data from folder
scp ../MALE_PD_filtered_sumstats_no_multi_allelics_RSID.txt .
scp ../FEMALE_PD_filtered_sumstats_no_multi_allelics_RSID.txt .

# The columns should be -> "snpid","A1","A2","Zscore", "N", "P.value"
# and need rs-ids...
cut -f 2- MALE_PD_filtered_sumstats_no_multi_allelics_RSID.txt > output.txt

munge_sumstats.py --out MALE --sumstats output.txt \
--merge-alleles eur_w_ld_chr/w_hm3.snplist \
--N 110616 --snp ID --a1 Allele1 --a2 Allele2 --p P-value --frq Freq1

cut -f 2- FEMALE_PD_filtered_sumstats_no_multi_allelics_RSID.txt > output.txt

munge_sumstats.py --out FEMALE --sumstats output.txt \
--merge-alleles eur_w_ld_chr/w_hm3.snplist \
--N 104082 --snp ID --a1 Allele1 --a2 Allele2 --p P-value --frq Freq1

# do final comparison

ldsc.py --rg MALE.sumstats.gz,FEMALE.sumstats.gz --out PD_male_female --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/

*********************************************************************
* LD Score Regression (LDSC)
* Version 1.0.0
* (C) 2014-2015 Brendan Bulik-Sullivan and Hilary Finucane
* Broad Institute of MIT and Harvard / MIT Department of Mathematics
* GNU General Public License v3
*********************************************************************
Call: 
./ldsc.py \
--ref-ld-chr eur_w_ld_chr/ \
--out PD_male_female \
--rg MALE.sumstats.gz,FEMALE.sumstats.gz \
--w-ld-chr eur_w_ld_chr/ 

Beginning analysis at Tue Dec 29 13:49:35 2020
Reading summary statistics from MALE.sumstats.gz ...
Read summary statistics for 1155229 SNPs.
Reading reference panel LD Score from eur_w_ld_chr/[1-22] ...
Read reference panel LD Scores for 1290028 SNPs.
Removing partitioned LD Scores with zero variance.
Reading regression weight LD Score from eur_w_ld_chr/[1-22] ...
Read regression weight LD Scores for 1290028 SNPs.
After merging with reference panel LD, 1149983 SNPs remain.
After merging with regression SNP LD, 1149983 SNPs remain.
Computing rg for phenotype 2/2
Reading summary statistics from FEMALE.sumstats.gz ...
Read summary statistics for 1217311 SNPs.
After merging with summary statistics, 1149983 SNPs remain.
1149058 SNPs with valid alleles.

Heritability of phenotype 1
---------------------------
Total Observed scale h2: 0.0632 (0.0069)
Lambda GC: 1.0588
Mean Chi^2: 1.0934
Intercept: 0.9548 (0.0077)
Ratio < 0 (usually indicates GC correction).

Heritability of phenotype 2/2
-----------------------------
Total Observed scale h2: 0.0538 (0.0078)
Lambda GC: 1.0255
Mean Chi^2: 1.0645
Intercept: 0.9516 (0.0077)
Ratio < 0 (usually indicates GC correction).

Genetic Covariance
------------------
Total Observed scale gencov: 0.0511 (0.0069)
Mean z1*z2: 0.1037
Intercept: -0.0079 (0.0066)

Genetic Correlation
-------------------
Genetic Correlation: 0.877 (0.0699)
Z-score: 12.5448
P: 4.2444e-36


Summary of Genetic Correlation Results
p1                  p2     rg      se        z           p  h2_obs  h2_obs_se  h2_int  h2_int_se  gcov_int  gcov_int_se
MALE.sumstats.gz  FEMALE.sumstats.gz  0.877  0.0699  12.5448  4.2444e-36  0.0538     0.0078  0.9516     0.0077   -0.0079       0.0066

```

### 2. LDSCORE on stats without UKB proxy...

```
# main purpose of this:
# to check for genome wide correlation between male and female PD summary statistics

cd /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/RESULTS/to_meta_files/

mkdir LDSCORE
cd LDSCORE

module load ldsc/1.0.0-101-g89c13a7
# get reference data
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
tar -jxvf eur_w_ld_chr.tar.bz2

### copy data from folder
scp ../MALE_PD_filtered_sumstats_NO_PROXY_MAF_no_multi_allelics_RSID.txt .
scp ../FEMALE_PD_filtered_sumstats_NO_PROXY_no_multi_allelics_RSID.txt .

# The columns should be -> "snpid","A1","A2","Zscore", "N", "P.value"
# and need rs-ids...
cut -f 2- MALE_PD_filtered_sumstats_NO_PROXY_MAF_no_multi_allelics_RSID.txt > output.txt

munge_sumstats.py --out MALE_no_proxy --sumstats output.txt \
--merge-alleles eur_w_ld_chr/w_hm3.snplist \
--N 32356 --snp ID --a1 Allele1 --a2 Allele2 --p P-value --frq Freq1

cut -f 2- FEMALE_PD_filtered_sumstats_NO_PROXY_no_multi_allelics_RSID.txt > output.txt

munge_sumstats.py --out FEMALE_no_proxy --sumstats output.txt \
--merge-alleles eur_w_ld_chr/w_hm3.snplist \
--N 28277 --snp ID --a1 Allele1 --a2 Allele2 --p P-value --frq Freq1

# do final comparison

ldsc.py --rg MALE_no_proxy.sumstats.gz,FEMALE_no_proxy.sumstats.gz --out PD_male_female_no_proxy --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/

*********************************************************************
* LD Score Regression (LDSC)
* Version 1.0.0
* (C) 2014-2015 Brendan Bulik-Sullivan and Hilary Finucane
* Broad Institute of MIT and Harvard / MIT Department of Mathematics
* GNU General Public License v3
*********************************************************************
Call: 
./ldsc.py \
--ref-ld-chr eur_w_ld_chr/ \
--out PD_male_female_no_proxy \
--rg MALE_no_proxy.sumstats.gz,FEMALE_no_proxy.sumstats.gz \
--w-ld-chr eur_w_ld_chr/ 

Beginning analysis at Thu Dec 31 13:40:07 2020
Reading summary statistics from MALE_no_proxy.sumstats.gz ...
Read summary statistics for 1150592 SNPs.
Reading reference panel LD Score from eur_w_ld_chr/[1-22] ...
Read reference panel LD Scores for 1290028 SNPs.
Removing partitioned LD Scores with zero variance.
Reading regression weight LD Score from eur_w_ld_chr/[1-22] ...
Read regression weight LD Scores for 1290028 SNPs.
After merging with reference panel LD, 1145796 SNPs remain.
After merging with regression SNP LD, 1145796 SNPs remain.
Computing rg for phenotype 2/2
Reading summary statistics from FEMALE_no_proxy.sumstats.gz ...
Read summary statistics for 1217311 SNPs.
After merging with summary statistics, 1145796 SNPs remain.
1144782 SNPs with valid alleles.

Heritability of phenotype 1
---------------------------
Total Observed scale h2: 0.2077 (0.0238)
Lambda GC: 1.0802
Mean Chi^2: 1.1171
Intercept: 0.9825 (0.008)
Ratio < 0 (usually indicates GC correction).

Heritability of phenotype 2/2
-----------------------------
Total Observed scale h2: 0.1857 (0.0304)
Lambda GC: 1.0466
Mean Chi^2: 1.081
Intercept: 0.9735 (0.0081)
Ratio < 0 (usually indicates GC correction).

Genetic Covariance
------------------
Total Observed scale gencov: 0.1738 (0.0241)
Mean z1*z2: 0.1021
Intercept: -0.0058 (0.0065)

Genetic Correlation
-------------------
Genetic Correlation: 0.8849 (0.07)
Z-score: 12.649
P: 1.1331e-36


Summary of Genetic Correlation Results
p1                           p2      rg    se       z           p  h2_obs  h2_obs_se  h2_int  h2_int_se  gcov_int  gcov_int_se
MALE_no_proxy.sumstats.gz  FEMALE_no_proxy.sumstats.gz  0.8849  0.07  12.649  1.1331e-36  0.1857     0.0304  0.9735     0.0081   -0.0058       0.0065

```

### 3. Polygenic risk scoring in NEUROX_DBGAP...

```
# main purpose of this:
# to check whether the genetic risk score is different between males and females with PD 

module load plink
module load R
module load flashpca

# create working dir
mkdir GRS_differences
cd /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/GRS_differences/

##### make new covariate from cases only
cd /PATH/TO/CORNELIS_TEMP/PD_FINAL_PLINK_2018/

awk '$10 == "NEUROX_DBGAP"' IPDGC_all_samples_covariates.txt | cut -f 1,2 > NEUROX_DBGAP_ONLY.txt 

plink --bfile UMIMPUTED_PD_september_2018_no_duplicates --maf 0.01 --geno 0.01 --hwe 5e-6 --make-bed --out filter_better \
--exclude range /PATH/TO/GENERAL/exclusion_regions_hg19.txt --memory 230000 --threads 19 --keep NEUROX_DBGAP_ONLY.txt --filter-cases
# 38945 variants and 5405 people pass filters and QC
plink --bfile filter_better --indep-pairwise 1000 10 0.02 --out pruned_better_and_more --memory 230000 --threads 19
plink --bfile filter_better --extract pruned_better_and_more.prune.in --make-bed --out input_pca --memory 230000 --threads 19
# 13024 variants and 5405 people pass filters and QC.
flashpca --bfile input_pca --suffix _NEUROX_DBGAP_ONLY_cases.txt --numthreads 19
mv *_NEUROX_DBGAP_ONLY_cases.txt /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/GRS_differences/

# now merge with OG covariates
cd /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/GRS_differences/

module load R
R
all_cov <- read.table("/PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/covariates/IPDGC_all_samples_covariates_SHORT.txt",header=T)
neuroX <- read.table("pcs_NEUROX_DBGAP_ONLY_cases.txt",header=T)
neuroX$IID <- NULL
MM = merge(all_cov,neuroX,by='FID')
write.table(MM, file="COV_NEUROX_DBGAP_ONLY_cases.txt", quote=FALSE,row.names=F,sep="\t")

##### calculate polygenic risk score
module load plink

cd /PATH/TO/CORNELIS_TEMP/PD_FINAL_PLINK_2018/
plink --bfile HARDCALLS_PD_september_2018_no_cousins \
--keep /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/GRS_differences/COV_NEUROX_DBGAP_ONLY_cases.txt \
--score /PATH/TO/GENERAL/META5_GRS_chr_bp_no_GBA_no_G2019S.txt \
--out /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/GRS_differences/SCORE_FILE.txt
# 5405 people remaining

##### test difference between male and female risk score...
cd /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/GRS_differences/

module load R
R
# Load in the necessary packages 
library(dplyr)
library(ggplot2)

# Load in data from the step above
data1 = read.table("SCORE_FILE.txt.profile",header=T)

# Also load in covariates (this can be done in PLINK or using flashPCA)
data2 = read.table("COV_NEUROX_DBGAP_ONLY_cases.txt",header=T)

# Drop the IID column to prevent double columns
data2$IID <- NULL

# Merge by FID column, can also make shorter in this case -> by='FID'
MM = merge(data2,data1,by.x='FID',by.y='FID')
temp <- anti_join(data1, data2, by='FID')
data <- rbind(temp, MM)

# Convert score values to Z-score
meanGRS <- mean(data$SCORE)
sdGRS <- sd(data$SCORE)
data$SCOREZ <- (data$SCORE - meanGRS)/sdGRS

# Simple plot Generate boxplot case-control data 
pdf("MALE_FEMALE_GRS_difference.pdf",width=6)
boxplot(data$SCOREZ~data$sex,col=c('powderblue', 'pink1'),xlab="1 Male, 2 Female",ylab="GRS Z-score")
grid()
dev.off()

# More complex plot Generate boxplot case-control data 
pdf("MALE_FEMALE_GRS_difference_fancy.pdf",width=6)
p <- ggplot(data, aes(x=as.factor(sex), y=SCOREZ, fill=as.factor(sex))) + 
  geom_violin(trim=FALSE)
# Add boxplot
p2 <- p+geom_boxplot(width=0.4, fill="white" ) + theme_minimal()
p2 + scale_fill_manual(values=c("powderblue", "pink1")) + theme_bw() + 
labs(title="IPDGC PD GRS",x="1 Male PD-case, 2 Female PD-case", y = "GRS Z-score") + theme(legend.position="none")
dev.off()

# Calculate if the GRS is statistical significant between cases and controls using AGE, SEX and first 5 principal components (PC)
thisFormula1 <- formula(paste("SCOREZ ~ SEX_COV + AGE + PC1 + PC2 + PC3 + PC4 + PC5"))
model1 <- glm(thisFormula1, data = data)
print(summary(model1))

Coefficients:
             Estimate Std. Error t value Pr(>|t|)    
(Intercept)  0.382510   0.068485   5.585 2.45e-08 ***
SEX_COV      0.071746   0.028377   2.528   0.0115 *  
AGE         -0.006659   0.001076  -6.189 6.49e-10 ***
PC1          1.173007   0.225946   5.192 2.16e-07 ***
PC2         -0.388056   0.433114  -0.896   0.3703    
PC3          0.468774   0.522476   0.897   0.3696    
PC4          0.878961   0.586377   1.499   0.1339    
PC5          0.292534   0.597333   0.490   0.6243    

```

### 4. Polygenic risk scoring in all datasets...

```
# main purpose of this:
# to check whether the genetic risk score is different between males and females with PD 

module load plink
module load R
module load flashpca

# create working dir (already done)
# mkdir GRS_differences
cd /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/GRS_differences/

##### make new covariate from cases only
cd /PATH/TO/CORNELIS_TEMP/PD_FINAL_PLINK_2018/

awk '$10 == "DUTCH"' IPDGC_all_samples_covariates.txt | awk '$6 == "2"' | cut -f 1,2 > DUTCH_cases.txt 
awk '$10 == "FINLAND"' IPDGC_all_samples_covariates.txt | awk '$6 == "2"' | cut -f 1,2 > FINLAND_cases.txt 
awk '$10 == "GERMANY"' IPDGC_all_samples_covariates.txt | awk '$6 == "2"' | cut -f 1,2 > GERMANY_cases.txt 
awk '$10 == "HBS"' IPDGC_all_samples_covariates.txt | awk '$6 == "2"' | cut -f 1,2 > HBS_cases.txt 
awk '$10 == "MCGILL"' IPDGC_all_samples_covariates.txt | awk '$6 == "2"' | cut -f 1,2 > MCGILL_cases.txt 
awk '$10 == "MF"' IPDGC_all_samples_covariates.txt | awk '$6 == "2"' | cut -f 1,2 > MF_cases.txt 
awk '$10 == "NEUROX_DBGAP"' IPDGC_all_samples_covariates.txt | awk '$6 == "2"' | cut -f 1,2 > NEUROX_DBGAP_cases.txt 
awk '$10 == "NIA"' IPDGC_all_samples_covariates.txt | awk '$6 == "2"' | cut -f 1,2 > NIA_cases.txt 
awk '$10 == "OSLO"' IPDGC_all_samples_covariates.txt | awk '$6 == "2"' | cut -f 1,2 > OSLO_cases.txt 
awk '$10 == "PDBP"' IPDGC_all_samples_covariates.txt | awk '$6 == "2"' | cut -f 1,2 > PDBP_cases.txt 
awk '$10 == "PPMI"' IPDGC_all_samples_covariates.txt | awk '$6 == "2"' | cut -f 1,2 > PPMI_cases.txt 
awk '$10 == "SHULMAN"' IPDGC_all_samples_covariates.txt | awk '$6 == "2"' | cut -f 1,2 > SHULMAN_cases.txt 
awk '$10 == "SPAIN3"' IPDGC_all_samples_covariates.txt | awk '$6 == "2"' | cut -f 1,2 > SPAIN3_cases.txt 
awk '$10 == "SPAIN4"' IPDGC_all_samples_covariates.txt | awk '$6 == "2"' | cut -f 1,2 > SPAIN4_cases.txt 
awk '$10 == "TUBI"' IPDGC_all_samples_covariates.txt | awk '$6 == "2"' | cut -f 1,2 > TUBI_cases.txt 
awk '$10 == "UK_GWAS"' IPDGC_all_samples_covariates.txt | awk '$6 == "2"' | cut -f 1,2 > UK_GWAS_cases.txt 
awk '$10 == "VANCE"' IPDGC_all_samples_covariates.txt | awk '$6 == "2"' | cut -f 1,2 > VANCE_cases.txt 

# make actual PC's
cut -f 10 IPDGC_all_samples_covariates.txt | sort -u | grep -v PROBAND | grep -v PROPARK | grep -v DATASET > NEW_list.txt

cat NEW_list.txt  | while read line
do
	plink --bfile UMIMPUTED_PD_september_2018_no_duplicates --maf 0.01 --geno 0.01 --hwe 5e-6 --make-bed --out filter_better \
	--exclude range /PATH/TO/GENERAL/exclusion_regions_hg19.txt --memory 230000 --threads 19 --keep "$line"_cases.txt --filter-cases
	plink --bfile filter_better --indep-pairwise 1000 10 0.02 --out pruned_better_and_more --memory 230000 --threads 19
	plink --bfile filter_better --extract pruned_better_and_more.prune.in --make-bed --out input_pca --memory 230000 --threads 19
	flashpca --bfile input_pca --suffix _$line.txt --numthreads 19
	mv *_$line.txt /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/GRS_differences/
done

# now merge with OG covariates
cd /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/GRS_differences/

module load R
R
for (dataset in c("DUTCH","FINLAND","GERMANY","HBS","MCGILL","MF","NEUROX_DBGAP","NIA","OSLO","PDBP","PPMI","SHULMAN","SPAIN3","SPAIN4","TUBI","UK_GWAS","VANCE")){
all_cov <- read.table("/PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/covariates/IPDGC_all_samples_covariates_SHORT.txt",header=T)
case <- read.table(paste("pcs_",dataset,".txt",sep=""),header=T)
case$IID <- NULL
CASE2 = merge(all_cov,case,by='FID')
write.table(CASE2, file=paste(dataset,"_cov_case.txt",sep=""), quote=FALSE,row.names=F,sep="\t")
}
q()
n

# merge all together
cat *_cov_case.txt | grep -v PHENO_PLINK > temp.txt
head -1 DUTCH_cov_case.txt > header.txt
cat header.txt temp.txt > ALL_CASE_COVARIATES.txt

# clean-up
mkdir other_files
mv eigenvalues_* other_files/
mv pve_* other_files/
mv pcs_* other_files/
mv eigenvectors_* other_files/

# done...

##### calculate polygenic risk score
module load plink

cd /PATH/TO/CORNELIS_TEMP/PD_FINAL_PLINK_2018/

plink --bfile HARDCALLS_PD_september_2018_no_cousins \
--keep /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/GRS_differences/ALL_CASE_COVARIATES.txt \
--score /PATH/TO/GENERAL/META5_GRS_chr_bp_no_GBA_no_G2019Sv2.txt \
--out /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/GRS_differences/SCORE_FILE_ALL_CASE.txt
# 19438 people remaining

# now run the association tests... 

cd /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/GRS_differences/

module load R
R
all_cov <- read.table("ALL_CASE_COVARIATES.txt",header=T)
all_score <- read.table("SCORE_FILE_ALL_CASE.txt.profile",header=T)
all_score$IID <- NULL
data = merge(all_cov,all_score,by='FID')

# loop in R for datasets with AGE information
for (dataset in c("DUTCH","GERMANY","HBS","MCGILL","NEUROX_DBGAP","NIA","OSLO","PDBP","PPMI","SHULMAN","SPAIN3","SPAIN4","TUBI","UK_GWAS")){
print(dataset)
newdata <- subset(data, DATASET == dataset)
print(dim(newdata))
# Convert score values to Z-score
meanGRS <- mean(newdata$SCORE)
sdGRS <- sd(newdata$SCORE)
newdata$SCOREZ <- (newdata$SCORE - meanGRS)/sdGRS
print(dim(newdata))
thisFormula1 <- formula(paste("SCOREZ ~ SEX_COV + AGE + PC1 + PC2 + PC3 + PC4 + PC5"))
model1 <- glm(thisFormula1, data = newdata)
sink(paste("RESULTS",dataset,".txt",sep=""))
print(summary(model1))
sink()
print(summary(model1))
print("DONE")
}

# loop in R for datasets WITHOUT AGE information
for (dataset in c("MF","VANCE","FINLAND")){
print(dataset)
newdata <- subset(data, DATASET == dataset)
print(dim(newdata))
# Convert score values to Z-score
meanGRS <- mean(newdata$SCORE)
sdGRS <- sd(newdata$SCORE)
newdata$SCOREZ <- (newdata$SCORE - meanGRS)/sdGRS
print(dim(newdata))
thisFormula1 <- formula(paste("SCOREZ ~ SEX_COV + PC1 + PC2 + PC3 + PC4 + PC5"))
model1 <- glm(thisFormula1, data = newdata)
sink(paste("RESULTS",dataset,".txt",sep=""))
print(summary(model1))
sink()
print(summary(model1))
print("DONE")
}

# combine results...

ls | grep RESULTS > ORDER_of_results.txt
echo "dataset" > temp.txt
cat temp.txt ORDER_of_results.txt > col1.txt
cat RESULTSDUTCH.txt RESULTSFINLAND.txt RESULTSGERMANY.txt RESULTSHBS.txt RESULTSMCGILL.txt RESULTSMF.txt RESULTSNEUROX_DBGAP.txt RESULTSNIA.txt RESULTSOSLO.txt RESULTSPDBP.txt RESULTSPPMI.txt RESULTSSHULMAN.txt RESULTSSPAIN3.txt RESULTSSPAIN4.txt RESULTSTUBI.txt RESULTSUK_GWAS.txt RESULTSVANCE.txt | grep SEX_COV > COMBO_resssults.txt
grep Estimate RESULTSDUTCH.txt > header_resssults.txt
cat header_resssults.txt COMBO_resssults.txt > final_resssults.txt
paste col1.txt final_resssults.txt > final_resssultsv2.txt
# fix up file in editor...

# FOREST results... => Thanks Julie!
module load R
R
library(metafor)
data <- read.table("final_resssultsv2.txt", header = T)
labs <- gsub(".*\\.","", data$dataset)
yi   <- data$Estimate
sei  <- data$StdError
resFe  <- rma(yi=yi, sei=sei, method="FE")
resRe  <- rma(yi=yi, sei=sei)
print(summary(resFe))
print(summary(resRe))
pdf(file ="final_plot_random_effects.pdf", width = 8, height = 6)
#tiff(file="test.tiff", res=300, width = 1600, height = 1600)
Pvalue <- formatC(resRe$pval, digits=4)
forest(resRe, xlim=c(-2,2),atransf=exp, xlab=paste("Odds Ratio (95%CI) for SNP",sep=""), slab=labs, mlab="Random Effects", col = "red", border = "red", cex=.9, at=log(c(0.5, 1, 2)))
mtext("GRS ~ male/female status", line=.3, cex=1.2, font=2)
mtext(paste("P=",Pvalue,sep=""), line=-1, cex=1, font=2)
dev.off()

### DONE...
```


### 5. Test differences of beta's...

```

# main purpose of this:
# to check whether any of the beta's are significantly different between males and females...

module load R
cd /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/RESULTS
R
data <- read.table("SUP_table2.txt",header=T)
data$Z = (data$Effect_male-data$Effect_female) / sqrt(data$StdErr_male^2 + data$StdErr_female^2)
data$P_dif = 2*pnorm(-abs(data$Z))
write.table(data,file="SUP_table2_with_comparison_test.txt",sep="\t",row.names=F,quote=FALSE)
q()
n

# Done...
```

### 6. Male/Female specific polygenic risk scoring in UKB...

```
# main purpose of this:
# to check whether prediction of PD status using GRS is difference between males and females...

cd /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/RESULTS/to_meta_files/

# filter for MAF >1%, at least present in 12 datasets out of 17 (~2/3) this makes sure NeuroX isnt a problem, and I2 of 80%
head -1 META_NEW_MALE_NO_UKB_AT_ALL1.tbl > header.txt
awk '{if ($6 > 0.01) print $0;}' META_NEW_MALE_NO_UKB_AT_ALL1.tbl | awk '{if ($6 < 0.99) print $0;}' | awk '{if ($14 > 11) print $0;}' | awk '{if ($12 < 80) print $0;}' > temp
# filter for genome-wide significance
awk '{if ($10 < 0.00000005) print $0;}' temp > temp2
# add header back in
cat header.txt temp2 > META_NEW_MALE_NO_UKB_AT_ALL1_filtered_sumstats.txt

# filter for MAF >1%, at least present in 12 datasets out of 17 (~2/3) this makes sure NeuroX isnt a problem, and I2 of 80%
head -1 META_NEW_FEMALE_NO_UKB_AT_ALL1.tbl > header.txt
awk '{if ($6 > 0.01) print $0;}' META_NEW_FEMALE_NO_UKB_AT_ALL1.tbl | awk '{if ($6 < 0.99) print $0;}' | awk '{if ($14 > 11) print $0;}' | awk '{if ($12 < 80) print $0;}' > temp
# filter for genome-wide significance
awk '{if ($10 < 0.00000005) print $0;}' temp > temp2
# add header back in
cat header.txt temp2 > META_NEW_FEMALE_NO_UKB_AT_ALL1_filtered_sumstats.txt

# create male and female specific score...

### option1: Use Nalls et al 2019 variants and assign sex-specific effect sizes
### option2: Use only top variant from each genomic region passing significance in each sex-specific GWAS


### option1: Use Nalls et al 2019 variants and assign sex-specific effect sizes
module load R
R
require("data.table")
male <- fread("META_NEW_MALE_NO_UKB_AT_ALL1.tbl",header=T)
female <- fread("META_NEW_FEMALE_NO_UKB_AT_ALL1.tbl",header=T)
GRS <- fread("/PATH/TO/GENERAL/META5_GRS_chr_bp_no_GBA_no_G2019Sv2.txt",header=F)
GRS$V2 <- NULL
GRS$V3 <- NULL
MM = merge(GRS,male,by.x='V1',by.y="MarkerName")
FF = merge(GRS,female,by.x='V1',by.y="MarkerName")
MM2 <- data.frame(MM$V1,MM$Allele1,MM$Effect)
FF2 <- data.frame(FF$V1,FF$Allele1,FF$Effect)
# make allele capitalized
MM2$A1 <- toupper(MM2$MM.Allele1)
FF2$A1 <- toupper(FF2$FF.Allele1)
# rename and reorder...
MM2$MM.Allele1 <- NULL
colnames(MM2)[1] <- "SNP"
colnames(MM2)[2] <- "Effect"
colnames(MM2)[3] <- "Allele"
MM3 <- MM2[,c(1,3,2)]
FF2$FF.Allele1 <- NULL
colnames(FF2)[1] <- "SNP"
colnames(FF2)[2] <- "Effect"
colnames(FF2)[3] <- "Allele"
FF3 <- FF2[,c(1,3,2)]
# load in variant IDs
UKB <- fread("/PATH/TO/UKBIOBANK/PD_GRS_ONLY/convert_sheet.txt",header=F)
MM4 = merge(MM3,UKB,by.x='SNP',by.y="V7")
FF4 = merge(FF3,UKB,by.x='SNP',by.y="V7")
MM5 <- MM4[,c(5,2,3)]
FF5 <- FF4[,c(5,2,3)]
colnames(MM5)[1] <- "SNP"
colnames(FF5)[1] <- "SNP"
# save...
write.table(MM5, file="META5_score_file_MALE_weights.txt", quote=FALSE,row.names=F,sep="\t")
write.table(FF5, file="META5_score_file_FEMALE_weights.txt", quote=FALSE,row.names=F,sep="\t")
q()
n

# calculate risk score...
module load plink
cd /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/RESULTS/to_meta_files/

plink --bfile /PATH/TO/UKBIOBANK/PD_GRS_ONLY/PD_GRS_UKBB_EURO_only_temp \
--score META5_score_file_MALE_weights.txt \
--keep /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/UKBB/COV_UKB_PD_MALE_cases_control_over60.txt \
--out META5_score_file_MALE_weights_OUTPUT_MALES

plink --bfile /PATH/TO/UKBIOBANK/PD_GRS_ONLY/PD_GRS_UKBB_EURO_only_temp \
--score META5_score_file_MALE_weights.txt \
--keep /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/UKBB/COV_UKB_PD_FEMALE_cases_control_over60.txt \
--out META5_score_file_MALE_weights_OUTPUT_FEMALES

plink --bfile /PATH/TO/UKBIOBANK/PD_GRS_ONLY/PD_GRS_UKBB_EURO_only_temp \
--score META5_score_file_FEMALE_weights.txt \
--keep /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/UKBB/COV_UKB_PD_MALE_cases_control_over60.txt \
--out META5_score_file_FEMALE_weights_OUTPUT_MALES

plink --bfile /PATH/TO/UKBIOBANK/PD_GRS_ONLY/PD_GRS_UKBB_EURO_only_temp \
--score META5_score_file_FEMALE_weights.txt \
--keep /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/UKBB/COV_UKB_PD_FEMALE_cases_control_over60.txt \
--out META5_score_file_FEMALE_weights_OUTPUT_FEMALES

plink --bfile /PATH/TO/UKBIOBANK/PD_GRS_ONLY/PD_GRS_UKBB_EURO_only_temp \
--score /PATH/TO/GENERAL/META5_no_UKBB_GRS_no_GBA_no_G2019S.txt \
--keep /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/UKBB/COV_UKB_PD_MALE_cases_control_over60.txt \
--out META5_score_file_NORMAL_weights_OUTPUT_MALES

plink --bfile /PATH/TO/UKBIOBANK/PD_GRS_ONLY/PD_GRS_UKBB_EURO_only_temp \
--score /PATH/TO/GENERAL/META5_no_UKBB_GRS_no_GBA_no_G2019S.txt \
--keep /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/UKBB/COV_UKB_PD_FEMALE_cases_control_over60.txt \
--out META5_score_file_NORMAL_weights_OUTPUT_FEMALES

### calculations in R
module load R
R
library("pROC")
male_male = read.table("META5_score_file_MALE_weights_OUTPUT_MALES.profile",header=T)
male_female = read.table("META5_score_file_MALE_weights_OUTPUT_FEMALES.profile",header=T)
female_male = read.table("META5_score_file_FEMALE_weights_OUTPUT_MALES.profile",header=T)
female_female = read.table("META5_score_file_FEMALE_weights_OUTPUT_FEMALES.profile",header=T)
norm_male = read.table("META5_score_file_NORMAL_weights_OUTPUT_MALES.profile",header=T)
norm_female = read.table("META5_score_file_NORMAL_weights_OUTPUT_FEMALES.profile",header=T)
COV_male = read.table("/PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/UKBB/COV_UKB_PD_MALE_cases_control_over60.txt",header=T)
COV_female = read.table("/PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/UKBB/COV_UKB_PD_FEMALE_cases_control_over60.txt",header=T)
COV_male$IID <- NULL
COV_female$IID <- NULL

# Merge by FID column, can also make shorter in this case -> by='FID'
MM1 = merge(male_male,COV_male,by.x='FID',by.y='FID')
MM2 = merge(female_male,COV_male,by.x='FID',by.y='FID')
MM3 = merge(norm_male,COV_male,by.x='FID',by.y='FID')
FF1 = merge(female_female,COV_female,by.x='FID',by.y='FID')
FF2 = merge(male_female,COV_female,by.x='FID',by.y='FID')
FF3 = merge(norm_female,COV_female,by.x='FID',by.y='FID')

Do 6 times....
MM1 => P=4.569776e-25 Beta= 0.361536 SE=0.034959 AUC=0.6045 CI=95% CI: 0.5854-0.6235
MM2 => P=2.120856e-25 Beta= 0.366467 SE=0.035187 AUC=0.6059 CI=95% CI: 0.5871-0.6247
MM3 => P=1.206051e-26 Beta= 0.374805 SE=0.035080 AUC=0.6079 CI=95% CI: 0.589-0.6268 
FF1 => P=2.352693e-20 Beta= 0.41421 SE=0.04480 AUC=0.6203 CI=95% CI: 0.5963-0.6442 
FF2 => P=1.280007e-19 Beta= 0.40444 SE=0.04463 AUC=0.6179 CI=95% CI: 0.5941-0.6417 
FF3 => P=3.557093e-21 Beta= 0.42245 SE=0.04473 AUC=0.6229 CI=95% CI: 0.5989-0.6468

# Convert score values to Z-score
data <- FF3
meanGRS <- mean(data$SCORE)
sdGRS <- sd(data$SCORE)
data$SCOREZ <- (data$SCORE - meanGRS)/sdGRS
# Calculate if the GRS is statistical significant between cases and controls using AGE, SEX and first 5 principal components (PC)
data$PHENO <- data$STATUS-1
thisFormula1 <- formula(paste("PHENO ~ SCOREZ + GENETIC_SEX + AGE_OF_RECRUIT + PC1 + PC2 + PC3 + PC4 + PC5"))
model1 <- glm(thisFormula1, data = data, family=binomial)
print(summary(model1))
summary(model1)$coefficients[,4] 
rocAuc <- roc(data$PHENO, data$SCOREZ)
auc(rocAuc)
ci(rocAuc, of="auc")
coords(rocAuc, "best")


### option2: Use only top variant from each genomic region passing significance in each sex-specific GWAS

Checking metal output using NO UKBIOBANK data...
Then taking top variant from each genomic region passing 5E-8 significant from Figure1

Males: 13 variants
Females: 10 variants

Then extracting these from UKB biobank...



```


![Alt Text](https://media3.giphy.com/media/l0MYw3oeYCUJhj5FC/200.gif)



