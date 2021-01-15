## Adding in UK Biobank


### First sorting out sample inclusion files

```
cd /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/UKBB/

# PD cases
scp /PATH/TO/UKBIOBANK/PHENOTYPE_DATA/disease_groups/PD.txt . 
# proxies
scp /PATH/TO/UKBIOBANK/PHENOTYPE_DATA/disease_groups/PD_parent_no_PD.txt .
# Samples to NOT use for controls are here:
scp /PATH/TO/UKBIOBANK/PHENOTYPE_DATA/disease_groups/PD_case_or_PD_parent.txt .

# Full covariates are here:
scp /PATH/TO/UKBIOBANK/PHENOTYPE_DATA/covariates_phenome_to_use.txt .
to include in GWAS are:
PC1-5 (to be generated), TOWNSERND, AGE of RECRUIT and SEX

# To create PCs from here:
/PATH/TO/UKBIOBANK/raw_genotypes_no_cousins/UKBB_raw_data_no_cousins.*

# Create controls from here:
/PATH/TO/UKBIOBANK/raw_genotypes_no_cousins/UKBB_raw_data_no_cousins.*
# RANDOM sample over 60 years old people from here:
/PATH/TO/UKBIOBANK/ICD10_UKBB/Covariates/all_samples_60andover.txt
# but make sure you exclude the PD cases and proxies from there....
# Ideally case-control ration 1:10

# get origin of proxy from here
scp /PATH/TO/UKBIOBANK/OLD/illness_father_PD_CASES_20107.txt /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/UKBB/
scp /PATH/TO/UKBIOBANK/OLD/illness_mother_PD_CASES_20110.txt /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/UKBB/

wc -l illness_father_PD_CASES_20107.txt # 10831
wc -l illness_mother_PD_CASES_20110.txt # 7533
---- > 18364
cat illness_father_PD_CASES_20107.txt illness_mother_PD_CASES_20110.txt | sed -e 's/"//g' | sort -u | wc -l
---- > 18164
# so 200 overlapping with both parents with PD... 
# removing these from father given the higher N there...
grep -v -f illness_mother_PD_CASES_20110.txt illness_father_PD_CASES_20107.txt > illness_father_PD_CASES_20107_no_duplicates.txt
wc -l illness_father_PD_CASES_20107_no_duplicates.txt # 10631
wc -l illness_mother_PD_CASES_20110.txt # 7533
cat illness_father_PD_CASES_20107_no_duplicates.txt illness_mother_PD_CASES_20110.txt | sed -e 's/"//g' | sort -u | wc -l
---- > 18164
# remove ""
sed -i -e 's/"//g' illness_father_PD_CASES_20107_no_duplicates.txt
sed -i -e 's/"//g' illness_mother_PD_CASES_20110.txt

# remove PD cases...
cut -f 1 PD.txt > first_row_cases_only.txt
grep -v -f first_row_cases_only.txt illness_father_PD_CASES_20107_no_duplicates.txt > illness_father_PD_CASES_20107_no_duplicates_no_PD.txt
grep -v -f first_row_cases_only.txt illness_mother_PD_CASES_20110.txt > illness_mother_PD_CASES_20110_no_PD.txt

# pre QC:
number of samples with father PD => 10531
number of samples with mother PD => 7464

number of males with PD => XX
number of females with PD => XX

```

### Make covariate files

```
# TO make covariates files: N=4
PD vs control males
PD vs control females
PD proxy vs control father PD
PD proxy vs control mother PD

# make subset files...
module load R
cd /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/UKBB/
R
require("data.table")
require("dplyr")
cov <- fread("covariates_phenome_to_use.txt",header=T)
PD <- fread("PD.txt",header=T)
Mother <- fread("illness_mother_PD_CASES_20110_no_PD.txt",header=F)
Father <- fread("illness_father_PD_CASES_20107_no_duplicates_no_PD.txt",header=F)
keep <- fread("/PATH/TO/UKBIOBANK/raw_genotypes_no_cousins/UKBB_raw_data_no_cousins.fam",header=F)
not_control <- fread("/PATH/TO/UKBIOBANK/PHENOTYPE_DATA/disease_groups/PD_case_or_PD_parent.txt",header=T)
# remove columns
PDshort <-  data.frame(PD$eid)
covshort <- data.frame(cov$FID, cov$AGE_OF_RECRUIT, cov$GENETIC_SEX)
keepshort <- data.frame(keep$V1)
not_controlshort <- data.frame(not_control$eid)
# rename
names(PDshort) <- c("FID")
names(Father) <- c("FID")
names(Mother) <- c("FID")
names(covshort) <- c("FID", "AGE", "SEX")
names(keepshort) <- c("FID")
names(not_controlshort) <- c("FID")
# merging
PD2 = merge(PDshort,covshort,by.x='FID',by.y='FID')
Father2 = merge(Father,covshort,by.x='FID',by.y='FID')
Mother2 = merge(Mother,covshort,by.x='FID',by.y='FID')
keep2 = merge(keepshort,covshort,by.x='FID',by.y='FID')
not_control2 = merge(not_controlshort,covshort,by.x='FID',by.y='FID')
# samples to keep
keep2$AGE <- NULL
keep2$SEX <- NULL
# subset males and females in PD2
PD_male <- subset(PD2, SEX==1)
PD_female <- subset(PD2, SEX==0)
# merge with keep to only keep "good" samples
PD_male2 = merge(PD_male,keep2,by.x='FID',by.y='FID')
PD_female2 = merge(PD_female,keep2,by.x='FID',by.y='FID')
Father3 = merge(Father2,keep2,by.x='FID',by.y='FID')
Mother3 = merge(Mother2,keep2,by.x='FID',by.y='FID')
not_control3 = merge(not_control2,keep2,by.x='FID',by.y='FID')
# check sizes
dim(PD_male2) # 967
dim(PD_female2) # 563
dim(Father3) # 7948
dim(Mother3) # 5487
dim(not_control3) # 14960
# get old controls....
keep2 = merge(keepshort,covshort,by.x='FID',by.y='FID')
keep4 <- anti_join(keep2,not_control3, by = c('FID'))
keep5 <- subset(keep4, AGE >= 60)
dim(keep5)# 156208
# add column on status
PD_male2$STATUS <- "PD"
PD_female2$STATUS <- "PD"
Father3$STATUS <- "PROXY"
Mother3$STATUS <- "PROXY"
keep5$STATUS <- "CONTROL"
# get random controls... 
PD3_controls <- sample_n(keep5, 15300,replace = FALSE)
Proxy3_controls <- anti_join(keep5,PD3_controls, by = c('FID'))
dim(PD3_controls) # 15300
dim(Proxy3_controls) # 140908
# sort out PD male-female
PD_male_control <- subset(PD3_controls, SEX==1)
PD_female_control <- subset(PD3_controls, SEX==0)
PD_male_FINAL <- rbind(PD_male2, PD_male_control)
PD_female_FINAL <- rbind(PD_female2, PD_female_control)
# sort out proxy mother-father
dim(Proxy3_controls) # 140908
Proxy3_controls_mother <- sample_n(Proxy3_controls, 70454,replace = FALSE)
Proxy3_controls_father <- anti_join(Proxy3_controls,Proxy3_controls_mother, by = c('FID'))
Mother4 <- rbind(Mother3, Proxy3_controls_mother)
Father4 <- rbind(Father3, Proxy3_controls_father)
# add some things to make plink love these files
PD_male_FINAL$IID <- PD_male_FINAL$FID
PD_female_FINAL$IID <- PD_female_FINAL$FID
PD_male_FINAL <- PD_male_FINAL[c(1,5,3,4,2)]
PD_female_FINAL <- PD_female_FINAL[c(1,5,3,4,2)]
# save files
write.table(PD_male_FINAL, file="UKB_PD_MALE_cases_control_over60.txt", quote=FALSE,row.names=F,sep="\t")
write.table(PD_female_FINAL, file="UKB_PD_FEMALE_cases_control_over60.txt", quote=FALSE,row.names=F,sep="\t")
# add some things to make plink love these files v2
Mother4$IID <- Mother4$FID
Father4$IID <- Father4$FID 
Mother4 <- Mother4[,c(1,5,3,4,2)]
Father4 <- Father4[,c(1,5,3,4,2)]
# save files
write.table(Father4, file="UKB_PROXY_FATHER_cases_control_over60.txt", quote=FALSE,row.names=F,sep="\t")
write.table(Mother4, file="UKB_PROXY_MOTHER_cases_control_over60.txt", quote=FALSE,row.names=F,sep="\t")
q() n

# Total numbers:
UKB_PD_FEMALE_cases_control_over60.txt # 8518
cases: 563
controls: 7954
UKB_PD_MALE_cases_control_over60.txt # 8314
cases: 967
controls: 7346
UKB_PROXY_FATHER_cases_control_over60.txt # 78403
proxies: 7948
controls: 70454
UKB_PROXY_MOTHER_cases_control_over60.txt # 75942
proxies: 5487
controls: 70454

# Done for now...
```

### Calculate PC's

```
cd /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/UKBB/

# calculate PC's
https://github.com/gabraham/flashpca
cd 
module load flashpca
module load plink

# make file for loop...
ls | grep UKB_P > PC_files.txt

# start loop
cat PC_files.txt  | while read line
do 
	plink --bfile /PATH/TO/UKBIOBANK/raw_genotypes_no_cousins/UKBB_raw_data_no_cousins \
	--maf 0.05 --geno 0.01 --hwe 5e-6 --autosome \
	--exclude /PATH/TO/GENERAL/exclusion_regions_hg19.txt --make-bed --out FILENAME_2 \
	--keep $line
	plink --bfile FILENAME_2 --indep-pairwise 1000 10 0.02 --autosome --out pruned_data
	plink --bfile FILENAME_2 --extract pruned_data.prune.in --make-bed --out FILENAME_3 
	flashpca --bfile FILENAME_3 --suffix _$line --numthreads 19
done

# sanity check, note header is present in UKB_PD/PROXY files...:

UKB_PD_FEMALE_cases_control_over60.txt # 8518
pcs_UKB_PD_FEMALE_cases_control_over60.txt # 8518
UKB_PD_MALE_cases_control_over60.txt # 8314
pcs_UKB_PD_MALE_cases_control_over60.txt # 8314 
UKB_PROXY_FATHER_cases_control_over60.txt # 78403
pcs_UKB_PROXY_FATHER_cases_control_over60.txt # 78403 
UKB_PROXY_MOTHER_cases_control_over60.txt # 75942
pcs_UKB_PROXY_MOTHER_cases_control_over60.txt # 75942 


# merge files in R

module load R
R
require(data.table)
cov <- fread("covariates_phenome_to_use.txt",header=T)
pc1 <- fread("pcs_UKB_PD_MALE_cases_control_over60.txt",header=T)
pc2 <- fread("pcs_UKB_PD_FEMALE_cases_control_over60.txt",header=T)
pc3 <- fread("pcs_UKB_PROXY_FATHER_cases_control_over60.txt",header=T)
pc4 <- fread("pcs_UKB_PROXY_MOTHER_cases_control_over60.txt",header=T)
pc1$IID <- NULL
pc2$IID <- NULL
pc3$IID <- NULL
pc4$IID <- NULL
pheno1 <- fread("UKB_PD_MALE_cases_control_over60.txt",header=T)
pheno2 <- fread("UKB_PD_FEMALE_cases_control_over60.txt",header=T)
pheno3 <- fread("UKB_PROXY_FATHER_cases_control_over60.txt",header=T)
pheno4 <- fread("UKB_PROXY_MOTHER_cases_control_over60.txt",header=T)
pheno1$IID <- pheno1$SEX <- pheno1$AGE <- NULL
pheno2$IID <- pheno2$SEX <- pheno2$AGE <- NULL
pheno3$IID <- pheno3$SEX <- pheno3$AGE <- NULL
pheno4$IID <- pheno4$SEX <- pheno4$AGE <- NULL
MM1 = merge(cov,pheno1,by='FID')
MM2 = merge(cov,pheno2,by='FID')
MM3 = merge(cov,pheno3,by='FID')
MM4 = merge(cov,pheno4,by='FID')
MM11 = merge(MM1,pc1,by='FID')
MM22 = merge(MM2,pc2,by='FID')
MM33 = merge(MM3,pc3,by='FID')
MM44 = merge(MM4,pc4,by='FID')
write.table(MM11, file="COV_UKB_PD_MALE_cases_control_over60.txt", quote=FALSE,row.names=F,sep="\t")
write.table(MM22, file="COV_UKB_PD_FEMALE_cases_control_over60.txt", quote=FALSE,row.names=F,sep="\t")
write.table(MM33, file="COV_UKB_PROXY_FATHER_cases_control_over60.txt", quote=FALSE,row.names=F,sep="\t")
write.table(MM44, file="COV_UKB_PROXY_MOTHER_cases_control_over60.txt", quote=FALSE,row.names=F,sep="\t")
q()
n

# sanity check, note header is present in UKB_PD/PROXY files...:
   8314 COV_UKB_PD_MALE_cases_control_over60.txt
   8518 COV_UKB_PD_FEMALE_cases_control_over60.txt
   78403 COV_UKB_PROXY_FATHER_cases_control_over60.txt
   75942 COV_UKB_PROXY_MOTHER_cases_control_over60.txt.txt

# replace PD/CONTROL/PROXY with numbers
sed -i 's/PD/2/g' COV_UKB_PD_MALE_cases_control_over60.txt
sed -i 's/CONTROL/1/g' COV_UKB_PD_MALE_cases_control_over60.txt
sed -i 's/PD/2/g' COV_UKB_PD_FEMALE_cases_control_over60.txt
sed -i 's/CONTROL/1/g' COV_UKB_PD_FEMALE_cases_control_over60.txt
sed -i 's/PROXY/2/g' COV_UKB_PROXY_FATHER_cases_control_over60.txt
sed -i 's/CONTROL/1/g' COV_UKB_PROXY_FATHER_cases_control_over60.txt
sed -i 's/PROXY/2/g' COV_UKB_PROXY_MOTHER_cases_control_over60.txt
sed -i 's/CONTROL/1/g' COV_UKB_PROXY_MOTHER_cases_control_over60.txt
```

### Start GWAS

```
cd /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/UKBB/
mkdir GWAS_output
module load plink/2.0-dev-20191128

# note data prefiltered for:
R2 =< 0_8
--geno 0.1
--hwe 1e-6
--keep EUROPEAN.txt
--maf 0.01
--mind 0.1

# 1) COV_UKB_PD_MALE_cases_control_over60
# 1 binary phenotype loaded (966 cases, 7337 controls)
for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
  do
plink2 --pfile /PATH/TO/UKBIOBANK/FILTER_IMPUTED_DATA/chr$chnum.UKBB.EU.filtered \
--pheno-name STATUS --pheno COV_UKB_PD_MALE_cases_control_over60.txt \
--covar COV_UKB_PD_MALE_cases_control_over60.txt --memory 235000 \
--glm hide-covar firth-fallback cols=+a1freq,+a1freqcc,+a1count,+totallele,+a1countcc,+totallelecc,+err \
--out GWAS_output/UKB_PD_MALE_cases_control_over60_chr$chnum --covar-name AGE_OF_RECRUIT,TOWNSEND,PC1,PC2,PC3,PC4,PC5 --covar-variance-standardize
done

------------------

# 2) COV_UKB_PD_FEMALE_cases_control_over60
# 1 binary phenotype loaded (563 cases, 7941 controls).
for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
  do
plink2 --pfile /PATH/TO/UKBIOBANK/FILTER_IMPUTED_DATA/chr$chnum.UKBB.EU.filtered \
--pheno-name STATUS --pheno COV_UKB_PD_FEMALE_cases_control_over60.txt \
--covar COV_UKB_PD_FEMALE_cases_control_over60.txt --memory 235000 \
--glm hide-covar firth-fallback cols=+a1freq,+a1freqcc,+a1count,+totallele,+a1countcc,+totallelecc,+err \
--out GWAS_output/UKB_PD_FEMALE_cases_control_over60_chr$chnum --covar-name AGE_OF_RECRUIT,TOWNSEND,PC1,PC2,PC3,PC4,PC5 --covar-variance-standardize
done

------------------

# 3) COV_UKB_PROXY_FATHER_cases_control_over60
# 1 binary phenotype loaded (7936 cases, 70324 controls)

for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
  do
plink2 --pfile /PATH/TO/UKBIOBANK/FILTER_IMPUTED_DATA/chr$chnum.UKBB.EU.filtered \
--pheno-name STATUS --pheno COV_UKB_PROXY_FATHER_cases_control_over60.txt \
--covar COV_UKB_PROXY_FATHER_cases_control_over60.txt --memory 235000 \
--glm hide-covar firth-fallback cols=+a1freq,+a1freqcc,+a1count,+totallele,+a1countcc,+totallelecc,+err \
--out GWAS_output/UKB_PROXY_FATHER_cases_control_over60_chr$chnum --covar-name AGE_OF_RECRUIT,TOWNSEND,PC1,PC2,PC3,PC4,PC5 --covar-variance-standardize
done

# script:
#!/bin/bash
# sbatch --cpus-per-task=10 --mem=235g --time=24:00:00 PROXY_FATHER_GWAS.sh 1
CHNUM=$1
module load plink/2.0-dev-20191128
cd /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/UKBB/
plink2 --pfile /PATH/TO/UKBIOBANK/FILTER_IMPUTED_DATA/chr$CHNUM.UKBB.EU.filtered \
--pheno-name STATUS --pheno COV_UKB_PROXY_FATHER_cases_control_over60.txt \
--covar COV_UKB_PROXY_FATHER_cases_control_over60.txt --memory 235000 \
--glm hide-covar firth-fallback cols=+a1freq,+a1freqcc,+a1count,+totallele,+a1countcc,+totallelecc,+err \
--out GWAS_output/UKB_PROXY_FATHER_cases_control_over60_chr$CHNUM --covar-name AGE_OF_RECRUIT,TOWNSEND,PC1,PC2,PC3,PC4,PC5 --covar-variance-standardize

for chnum in {2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
  do
sbatch --cpus-per-task=10 --mem=235g --time=24:00:00 PROXY_FATHER_GWAS.sh $chnum
done

------------------

# 4) COV_UKB_PROXY_MOTHER_cases_control_over60
# 1 binary phenotype loaded (5473 cases, 70332 controls)
for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
  do
plink2 --pfile /PATH/TO/UKBIOBANK/FILTER_IMPUTED_DATA/chr$chnum.UKBB.EU.filtered \
--pheno-name STATUS --pheno COV_UKB_PROXY_MOTHER_cases_control_over60.txt \
--covar COV_UKB_PROXY_MOTHER_cases_control_over60.txt --memory 235000 \
--glm hide-covar firth-fallback cols=+a1freq,+a1freqcc,+a1count,+totallele,+a1countcc,+totallelecc,+err \
--out GWAS_output/UKB_PROXY_MOTHER_cases_control_over60_chr$chnum --covar-name AGE_OF_RECRUIT,TOWNSEND,PC1,PC2,PC3,PC4,PC5 --covar-variance-standardize
done

# script:
#!/bin/bash
# sbatch --cpus-per-task=10 --mem=235g --time=24:00:00 PROXY_MOTHER_GWAS.sh 1
CHNUM=$1
module load plink/2.0-dev-20191128
cd /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/UKBB/
plink2 --pfile /PATH/TO/UKBIOBANK/FILTER_IMPUTED_DATA/chr$CHNUM.UKBB.EU.filtered \
--pheno-name STATUS --pheno COV_UKB_PROXY_MOTHER_cases_control_over60.txt \
--covar COV_UKB_PROXY_MOTHER_cases_control_over60.txt --memory 235000 \
--glm hide-covar firth-fallback cols=+a1freq,+a1freqcc,+a1count,+totallele,+a1countcc,+totallelecc,+err \
--out GWAS_output/UKB_PROXY_MOTHER_cases_control_over60_chr$CHNUM --covar-name AGE_OF_RECRUIT,TOWNSEND,PC1,PC2,PC3,PC4,PC5 --covar-variance-standardize

for chnum in {2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
  do
sbatch --cpus-per-task=10 --mem=235g --time=24:00:00 PROXY_MOTHER_GWAS.sh $chnum
done

------------------

```

### Reformat output prep for meta-analysis

```
cd /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/UKBB/GWAS_output/

# merge data
cat UKB_PD_FEMALE_cases_control_over60_chr*hybrid | grep -v "ERRCODE" > temp.txt
head -1 UKB_PD_FEMALE_cases_control_over60_chr11*hybrid > header.txt
cat header.txt temp.txt > UKB_PD_FEMALE_cases_control_GWAS.txt

cat UKB_PD_MALE_cases_control_over60_chr*hybrid | grep -v "ERRCODE" > temp.txt
head -1 UKB_PD_MALE_cases_control_over60_chr11*hybrid > header.txt
cat header.txt temp.txt > UKB_PD_MALE_cases_control_GWAS.txt

cat UKB_PROXY_FATHER_cases_control_over60_chr*hybrid | grep -v "ERRCODE" > temp.txt
head -1 UKB_PROXY_FATHER_cases_control_over60_chr11*hybrid > header.txt
cat header.txt temp.txt > UKB_PROXY_FATHER_control_GWAS.txt

cat UKB_PROXY_MOTHER_cases_control_over60_chr*hybrid | grep -v "ERRCODE" > temp.txt
head -1 UKB_PROXY_FATHER_cases_control_over60_chr11*hybrid > header.txt
cat header.txt temp.txt > UKB_PROXY_MOTHER_control_GWAS.txt

-----------------------

# reformat data
module load python/3.6
python /PATH/TO/projects/CHR_X/UKBB/RESULTS/reformat_plink2_results.py --infile UKB_PD_FEMALE_cases_control_GWAS.txt \
--outfile toMeta.UKB_PD_FEMALE_cases_control_GWAS.txt --B-or-C B    

python /PATH/TO/projects/CHR_X/UKBB/RESULTS/reformat_plink2_results.py --infile UKB_PD_MALE_cases_control_GWAS.txt \
--outfile toMeta.UKB_PD_MALE_cases_control_GWAS.txt --B-or-C B    

python /PATH/TO/projects/CHR_X/UKBB/RESULTS/reformat_plink2_results.py --infile UKB_PROXY_FATHER_control_GWAS.txt \
--outfile toProxy.UKB_PROXY_FATHER_control_GWAS.txt --B-or-C B    

python /PATH/TO/projects/CHR_X/UKBB/RESULTS/reformat_plink2_results.py --infile UKB_PROXY_MOTHER_control_GWAS.txt \
--outfile toProxy.UKB_PD_PROXY_MOTHER_control_GWAS.txt --B-or-C B    

# convert proxies to "normal"

# make .csv
module load R
R
require("data.table")
data1 <- fread("toProxy.UKB_PROXY_FATHER_control_GWAS.txt",header=T)
data2 <- fread("toProxy.UKB_PD_PROXY_MOTHER_control_GWAS.txt",header=T)
write.table(data1, file="toConvert.UKB_PROXY_FATHER_control_GWAS.csv",quote=F,row.names=F,sep=",")
write.table(data2, file="toConvert.UKB_PD_PROXY_MOTHER_control_GWAS.csv",quote=F,row.names=F,sep=",")
q()
n

python /PATH/TO/projects/CHR_X/UKBB/RESULTS/Proxy_conversion/proxy_gwas_gwaxStyle.py \
--infile toConvert.UKB_PROXY_FATHER_control_GWAS.csv --beta-proxy beta \
--se-proxy se --p-proxy P --outfile toMeta.UKB_PROXY_FATHER_control_GWAS.csv

python /PATH/TO/projects/CHR_X/UKBB/RESULTS/Proxy_conversion/proxy_gwas_gwaxStyle.py \
--infile toConvert.UKB_PD_PROXY_MOTHER_control_GWAS.csv --beta-proxy beta \
--se-proxy se --p-proxy P --outfile toMeta.UKB_PD_PROXY_MOTHER_control_GWAS.csv

# make .txt
module load R
R
require("data.table")
data1 <- fread("toMeta.UKB_PROXY_FATHER_control_GWAS.csv",header=T)
data2 <- fread("toMeta.UKB_PD_PROXY_MOTHER_control_GWAS.csv",header=T)
write.table(data1, file="toMeta.UKB_PROXY_FATHER_control_GWAS.txt",quote=F,row.names=F,sep="\t")
write.table(data2, file="toMeta.UKB_PD_PROXY_MOTHER_control_GWAS.txt",quote=F,row.names=F,sep="\t")
q()
n

# prep the UKB files for forest

## normal PD cases
cut -f 1,2,3 toMeta.UKB_PD_MALE_cases_control_GWAS.txt > part1.txt
cut -f 4,5,6 toMeta.UKB_PD_MALE_cases_control_GWAS.txt > part3.txt
cut -f 7 toMeta.UKB_PD_MALE_cases_control_GWAS.txt > part2.txt
paste part1.txt part2.txt part3.txt > toMeta.SHORT_UKB_PD_MALE_cases_control_GWAS.txt

cut -f 1,2,3 toMeta.UKB_PD_FEMALE_cases_control_GWAS.txt > part1.txt
cut -f 4,5,6 toMeta.UKB_PD_FEMALE_cases_control_GWAS.txt > part3.txt
cut -f 7 toMeta.UKB_PD_FEMALE_cases_control_GWAS.txt > part2.txt
paste part1.txt part2.txt part3.txt > toMeta.SHORT_UKB_PD_FEMALE_cases_control_GWAS.txt

## proxy PDs
cut -f 1,2,3,7,15,16,17 toMeta.UKB_PROXY_FATHER_control_GWAS.txt > toMeta.SHORT_UKB_PROXY_FATHER_control_GWAS.txt
cut -f 1,2,3,7,15,16,17 toMeta.UKB_PD_PROXY_MOTHER_control_GWAS.txt > toMeta.SHORT_UKB_PROXY_MOTHER_control_GWAS.txt

## sanity checks...

grep 4:90626111 toMeta.SHORT*
looks OK!

# update variant names of UKB....

ls | grep toMeta.SHORT | grep txt > list.txt

cat list.txt  | while read line
do 
   	sed -i 's/:A//g' $line
	sed -i 's/:T//g' $line
	sed -i 's/:C//g' $line
	sed -i 's/:G//g' $line
done

### DONE for now...

# copy to metal
cd /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/UKBB/GWAS_output/

scp toMeta.SHORT_UKB_* /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/RESULTS/to_meta_files/

```



![Alt Text](https://media4.giphy.com/media/RkDX47fpp2nHlaZdjY/giphy.gif)



