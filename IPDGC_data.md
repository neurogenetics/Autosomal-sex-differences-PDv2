## Processing of the IPDGC data


```
Working folders:

Results folder => /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/

```


### Datasets to use: N=17 !!!
```
DUTCH
FINLAND
GERMANY
HBS
MCGILL
MF
NEUROX_DBGAP
NIA
OSLO
PDBP
PPMI
SHULMAN
SPAIN3
SPAIN4
TUBI
UK_GWAS
VANCE

# numbers:
male_case=12054
male_control=11999
female_case=7384
female_control=12389
ALL=43826

samples stored in:
females_only.txt
males_only.txt
```


#### Making covariate files

```
module load plink

# make file per dataset and per male/female
# females
awk '$10 == "DUTCH"' IPDGC_all_samples_covariates.txt | awk '$5 == "2"' | cut -f 1,2 > DUTCH_females.txt 
awk '$10 == "FINLAND"' IPDGC_all_samples_covariates.txt | awk '$5 == "2"' | cut -f 1,2 > FINLAND_females.txt 
awk '$10 == "GERMANY"' IPDGC_all_samples_covariates.txt | awk '$5 == "2"' | cut -f 1,2 > GERMANY_females.txt 
awk '$10 == "HBS"' IPDGC_all_samples_covariates.txt | awk '$5 == "2"' | cut -f 1,2 > HBS_females.txt 
awk '$10 == "MCGILL"' IPDGC_all_samples_covariates.txt | awk '$5 == "2"' | cut -f 1,2 > MCGILL_females.txt 
awk '$10 == "MF"' IPDGC_all_samples_covariates.txt | awk '$5 == "2"' | cut -f 1,2 > MF_females.txt 
awk '$10 == "NEUROX_DBGAP"' IPDGC_all_samples_covariates.txt | awk '$5 == "2"' | cut -f 1,2 > NEUROX_DBGAP_females.txt 
awk '$10 == "NIA"' IPDGC_all_samples_covariates.txt | awk '$5 == "2"' | cut -f 1,2 > NIA_females.txt 
awk '$10 == "OSLO"' IPDGC_all_samples_covariates.txt | awk '$5 == "2"' | cut -f 1,2 > OSLO_females.txt 
awk '$10 == "PDBP"' IPDGC_all_samples_covariates.txt | awk '$5 == "2"' | cut -f 1,2 > PDBP_females.txt 
awk '$10 == "PPMI"' IPDGC_all_samples_covariates.txt | awk '$5 == "2"' | cut -f 1,2 > PPMI_females.txt 
awk '$10 == "SHULMAN"' IPDGC_all_samples_covariates.txt | awk '$5 == "2"' | cut -f 1,2 > SHULMAN_females.txt 
awk '$10 == "SPAIN3"' IPDGC_all_samples_covariates.txt | awk '$5 == "2"' | cut -f 1,2 > SPAIN3_females.txt 
awk '$10 == "SPAIN4"' IPDGC_all_samples_covariates.txt | awk '$5 == "2"' | cut -f 1,2 > SPAIN4_females.txt 
awk '$10 == "TUBI"' IPDGC_all_samples_covariates.txt | awk '$5 == "2"' | cut -f 1,2 > TUBI_females.txt 
awk '$10 == "UK_GWAS"' IPDGC_all_samples_covariates.txt | awk '$5 == "2"' | cut -f 1,2 > UK_GWAS_females.txt 
awk '$10 == "VANCE"' IPDGC_all_samples_covariates.txt | awk '$5 == "2"' | cut -f 1,2 > VANCE_females.txt 
# males
awk '$10 == "DUTCH"' IPDGC_all_samples_covariates.txt | awk '$5 == "1"' | cut -f 1,2 > DUTCH_males.txt 
awk '$10 == "FINLAND"' IPDGC_all_samples_covariates.txt | awk '$5 == "1"' | cut -f 1,2 > FINLAND_males.txt 
awk '$10 == "GERMANY"' IPDGC_all_samples_covariates.txt | awk '$5 == "1"' | cut -f 1,2 > GERMANY_males.txt 
awk '$10 == "HBS"' IPDGC_all_samples_covariates.txt | awk '$5 == "1"' | cut -f 1,2 > HBS_males.txt 
awk '$10 == "MCGILL"' IPDGC_all_samples_covariates.txt | awk '$5 == "1"' | cut -f 1,2 > MCGILL_males.txt 
awk '$10 == "MF"' IPDGC_all_samples_covariates.txt | awk '$5 == "1"' | cut -f 1,2 > MF_males.txt 
awk '$10 == "NEUROX_DBGAP"' IPDGC_all_samples_covariates.txt | awk '$5 == "1"' | cut -f 1,2 > NEUROX_DBGAP_males.txt 
awk '$10 == "NIA"' IPDGC_all_samples_covariates.txt | awk '$5 == "1"' | cut -f 1,2 > NIA_males.txt 
awk '$10 == "OSLO"' IPDGC_all_samples_covariates.txt | awk '$5 == "1"' | cut -f 1,2 > OSLO_males.txt 
awk '$10 == "PDBP"' IPDGC_all_samples_covariates.txt | awk '$5 == "1"' | cut -f 1,2 > PDBP_males.txt 
awk '$10 == "PPMI"' IPDGC_all_samples_covariates.txt | awk '$5 == "1"' | cut -f 1,2 > PPMI_males.txt 
awk '$10 == "SHULMAN"' IPDGC_all_samples_covariates.txt | awk '$5 == "1"' | cut -f 1,2 > SHULMAN_males.txt 
awk '$10 == "SPAIN3"' IPDGC_all_samples_covariates.txt | awk '$5 == "1"' | cut -f 1,2 > SPAIN3_males.txt 
awk '$10 == "SPAIN4"' IPDGC_all_samples_covariates.txt | awk '$5 == "1"' | cut -f 1,2 > SPAIN4_males.txt 
awk '$10 == "TUBI"' IPDGC_all_samples_covariates.txt | awk '$5 == "1"' | cut -f 1,2 > TUBI_males.txt 
awk '$10 == "UK_GWAS"' IPDGC_all_samples_covariates.txt | awk '$5 == "1"' | cut -f 1,2 > UK_GWAS_males.txt 
awk '$10 == "VANCE"' IPDGC_all_samples_covariates.txt | awk '$5 == "1"' | cut -f 1,2 > VANCE_males.txt 

# make actual PC's
module load plink
module load flashpca

# prep LOOOOP to speed things up
cd /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/
# make list file
ls | grep "txt" | grep -v "IPDGC" | sed 's/.txt//g' > list.txt
scp list.txt ../PD_FINAL_PLINK_2018/

# LOOOOP
cd /PATH/TO/CORNELIS_TEMP/PD_FINAL_PLINK_2018/

cat list.txt  | while read line
do 
	plink --bfile UMIMPUTED_PD_september_2018_no_duplicates --maf 0.01 --geno 0.01 --hwe 5e-6 --make-bed --out filter_better \
	--exclude range /data/AND/GENERAL/exclusion_regions_hg19.txt --memory 230000 --threads 19 --keep /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/$line.txt
	plink --bfile filter_better --indep-pairwise 1000 10 0.02 --out pruned_better_and_more --memory 230000 --threads 19
	plink --bfile filter_better --extract pruned_better_and_more.prune.in --make-bed --out input_pca --memory 230000 --threads 19
	flashpca --bfile input_pca --suffix _$line.txt --numthreads 19
	mv *_$line.txt /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/covariates/
done

# now merge with OG covariates   
module load R
cd /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/covariates/

R

for (dataset in c("DUTCH","FINLAND","GERMANY","HBS","MCGILL","MF","NEUROX_DBGAP","NIA","OSLO","PDBP","PPMI","SHULMAN","SPAIN3","SPAIN4","TUBI","UK_GWAS","VANCE")){
all_cov <- read.table("IPDGC_all_samples_covariates_SHORT.txt",header=T)
male <- read.table(paste("pcs_",dataset,"_males.txt",sep=""),header=T)
female <- read.table(paste("pcs_",dataset,"_females.txt",sep=""),header=T)
male$IID <- NULL
female$IID <- NULL
Male = merge(all_cov,male,by='FID')
Female = merge(all_cov,female,by='FID')
write.table(Male, file=paste(dataset,"_cov_male.txt",sep=""), quote=FALSE,row.names=F,sep="\t")
write.table(Female, file=paste(dataset,"_cov_female.txt",sep=""), quote=FALSE,row.names=F,sep="\t")
}

# clean-up
mkdir other_files
mv eigenvalues_* other_files/
mv pve_* other_files/
mv pcs_* other_files/
mv eigenvectors_* other_files/

```

### RUN GWAS

```
RVTESTS => Program version: 20190205

# STORED IN BATCH SCRIPTS:
SEX_GWAS_final.sh
SEX_GWAS_final_NO_AGE.sh

#!/bin/bash

# sbatch --cpus-per-task=10 --mem=10g --time=24:00:00 SEX_GWAS_final.sh DUTCH male 22

dataset=$1
sex=$2
chr=$3

cd /PATH/TO/PD/imputed_data/$dataset

/PATH/TO/GENERAL/executable/rvtest --noweb --hide-covar --rangeFile maf001rsq03minimums_chr"$chr".txt \
--out /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/RESULTS/"$dataset"_"$sex"_"$chr" --single wald \
--inVcf chr"$chr".dose.vcf.gz --dosage DS --pheno /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/covariates/"$dataset"_cov_"$sex".txt \
--pheno-name PHENO_PLINK --covar /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/covariates/"$dataset"_cov_"$sex".txt \
--covar-name AGE,PC1,PC2,PC3,PC4,PC5 --imputeCov --numThread 9

'LOOOOOP'
# NO age available of datasets
cat dataset_no_age.txt | while read line
do
	for sex in {"male","female"}; 
	do
		for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
		do
			echo "bash SEX_GWAS_final_NO_AGE.sh $line $sex $chnum" >> gwas_no_age.swarm
done
done
done

# note dataset_no_age.txt => VANCE, MF and FINLAND

# age available of datasets
cat dataset_with_age.txt | while read line
do
	for sex in {"male","female"}; 
	do
		for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
		do
			echo "bash SEX_GWAS_final.sh $line $sex $chnum" >> gwas_with_age.swarm
done
done
done

# start swarm....
swarm -f gwas_no_age.swarm -g 10 -t 10 --time=24:00:00 --devel
swarm -f gwas_with_age.swarm -g 10 -t 10 --time=24:00:00 --devel

#######
'ERRORS'
#######

'error1' with all files forgot to $ sex so all files were called sex in stead of male/female....

'error2' with SPAIN4... because sample names do not match
in .vcf samples are without SPAIN4_
scp SPAIN4_cov_male.txt SPAIN4_cov_male_original.txt
scp SPAIN4_cov_female.txt SPAIN4_cov_female_original.txt

cut -f 1,2 SPAIN4_cov_male_original.txt | sed 's/SPAIN4_//g' > temp
paste temp SPAIN4_cov_male_original.txt | cut -f1,2,5- > SPAIN4_cov_male.txt

cut -f 1,2 SPAIN4_cov_female_original.txt | sed 's/SPAIN4_//g' > temp
paste temp SPAIN4_cov_female_original.txt | cut -f1,2,5- > SPAIN4_cov_female.txt

## rerun:
for sex in {"male","female"}; 
	do
			for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
			do
			echo "bash SEX_GWAS_final.sh SPAIN4 $sex $chnum" >> gwas_with_age_redo_SPAIN4.swarm
done
done

swarm -f gwas_with_age_redo_SPAIN4.swarm -g 10 -t 10 --time=24:00:00

'error3' Finland also needs to be without age....

## rerun:
for sex in {"male","female"}; 
do
	for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
	do
	echo "bash SEX_GWAS_final_NO_AGE.sh FINLAND $sex $chnum" >> gwas_no_age_redo_FINLAND.swarm
done
done

swarm -f gwas_no_age_redo_FINLAND.swarm -g 10 -t 10 --time=12:00:00

cat FINLAND_female*.SingleWald.assoc | grep -v 'N_INFORMATIVE' > allChrs_FINLAND_female.assoc 
cat FINLAND_male*.SingleWald.assoc | grep -v 'N_INFORMATIVE' > allChrs_FINLAND_male.assoc 

# then merge etc.
```


### Start prepping for meta-analyses

```
# merge all GWAS files per chromosome from GWAS per dataset
cd /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/RESULTS/

cat ../all_datasets.txt | while read line
do
for sex in {"male","female"}; 
do
cat "$line"_"$sex"*.SingleWald.assoc | grep -v 'N_INFORMATIVE' > allChrs_"$line"_"$sex".assoc 
done
done

# clean up
mkdir results_files_per_CHR
mv *.log results_files_per_CHR/
mv *.SingleWald.assoc results_files_per_CHR/

# merge all info files per chromosome from imputed data
'NOTE' cannot do separatedly between males and females because data is imputed together...

cd /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/
cat dataset_no_age.txt dataset_with_age.txt > all_datasets.txt

cat all_datasets.txt | while read line
do
	cd /data/AND/PD/imputed_data/$line
	cat maf001rsq03minimums_chr*.info | grep -v 'Rsq' > /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/RESULTS/allChrs_"$line".info
done

# Then merge in R
module load R
run like this:
R < make_to_meta_files.R --no-save  

#!/usr/bin/env Rscript
require("data.table")
# LOOOOP
for (dataset in c("DUTCH","FINLAND","GERMANY","HBS","MCGILL","MF","NEUROX_DBGAP","NIA","OSLO","PDBP","PPMI","SHULMAN","SPAIN3","SPAIN4","TUBI","UK_GWAS","VANCE")){
infos <- fread(paste("allChrs_",dataset,".info",sep=""))
colnames(infos) <- c("SNP","ALT_Frq","Rsq")
#  MALE
assoc <- fread(paste("allChrs_",dataset,"_male.assoc",sep=""))
colnames(assoc) <- c("CHROM","POS","REF","ALT","N_INFORMATIVE","Test","Beta","SE","Pvalue")
data <- merge(infos, assoc, by.x = "SNP", by.y = "Test", all.y = T)
dat <- subset(data, Beta < 5 & Beta > -5 & !is.na(data$Pvalue))
dat$chr <- paste("chr",dat$CHROM, sep = "")
dat$markerID <- paste(dat$chr,dat$POS, sep = ":")
dat$minorAllele <- ifelse(dat$ALT_Frq <= 0.5, as.character(dat$ALT), as.character(dat$REF))
dat$majorAllele <- ifelse(dat$ALT_Frq <= 0.5, as.character(dat$REF), as.character(dat$ALT))
dat$beta <- ifelse(dat$ALT_Frq <= 0.5, dat$Beta, dat$Beta*-1)
dat$se <- dat$SE
dat$maf <- ifelse(dat$ALT_Frq <= 0.5, dat$ALT_Frq, 1 - dat$ALT_Frq)
dat$P <- dat$Pvalue
dat0 <- dat[,c("markerID","minorAllele","majorAllele","beta","se","maf","P")]
write.table(dat0, file=(paste("toMeta.",dataset,"male.tab",sep="")), quote = F, sep = "\t", row.names = F)
# FEMALE
assoc <- fread(paste("allChrs_",dataset,"_female.assoc",sep=""))
colnames(assoc) <- c("CHROM","POS","REF","ALT","N_INFORMATIVE","Test","Beta","SE","Pvalue")
data <- merge(infos, assoc, by.x = "SNP", by.y = "Test", all.y = T)
dat <- subset(data, Beta < 5 & Beta > -5 & !is.na(data$Pvalue))
dat$chr <- paste("chr",dat$CHROM, sep = "")
dat$markerID <- paste(dat$chr,dat$POS, sep = ":")
dat$minorAllele <- ifelse(dat$ALT_Frq <= 0.5, as.character(dat$ALT), as.character(dat$REF))
dat$majorAllele <- ifelse(dat$ALT_Frq <= 0.5, as.character(dat$REF), as.character(dat$ALT))
dat$beta <- ifelse(dat$ALT_Frq <= 0.5, dat$Beta, dat$Beta*-1)
dat$se <- dat$SE
dat$maf <- ifelse(dat$ALT_Frq <= 0.5, dat$ALT_Frq, 1 - dat$ALT_Frq)
dat$P <- dat$Pvalue
dat0 <- dat[,c("markerID","minorAllele","majorAllele","beta","se","maf","P")]
write.table(dat0, file=(paste("toMeta.",dataset,"female.tab",sep="")), quote = F, sep = "\t", row.names = F)
}

# clean-up
mkdir to_meta_files
mv toMeta* to_meta_files/
mkdir merged_files_per_dataset
mv *info merged_files_per_dataset/
mv *assoc merged_files_per_dataset/
```


![Alt Text](https://media0.giphy.com/media/26u4lOMA8JKSnL9Uk/giphy.gif)



