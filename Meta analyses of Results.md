## Meta-analysing GWAS results

### Start preparing files for metal
```
cd /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/RESULTS/to_meta_files/

# fix variant naming issue...
mkdir ORIGINAL_TOMETA_FILES
scp toMeta.* ORIGINAL_TOMETA_FILES
# make list
ls | grep tab > list2.txt
# loop over list
cat list2.txt  | while read line
do 
	sed -i 's/chr//g' $line
done

# remove multi-allelics (if present)
ls | grep toMeta > list_to_remove_multi_allelic.txt

cat list_to_remove_multi_allelic.txt  | while read line
do 
   echo "$line"
   awk '{a[$1]++;b[$1]=$0}END{for(x in a)if(a[x]==1)print b[x]}' "$line" > no_multi_allelics_"$line"
done

# add back in header... Ooops
cat list_to_remove_multi_allelic.txt  | while read line
do 
   echo "$line"
   head -1 "$line" > header.txt
   cat header.txt no_multi_allelics_"$line" > temp.txt
   mv temp.txt no_multi_allelics_"$line"
done
# note in the end this is not needed due to step later on, when filtering straight from HRC for multi allelics... 
# thanks for trying though

```

### Run metal

```
# start meta-analyzing...
module load metal

# 19 files each.... 17x IPDGC + 2x UKB...

metal metal_females.txt
metal metal_males.txt

metal metal_females_no_multi_allelics.txt
metal metal_males_no_multi_allelics.txt

# also create sumstats withouth UKB proxies

metal metal_females_no_UKB_proxy.txt
metal metal_males_no_UKB_proxy.txt

# also create sumstats withouth UKB data at all...

metal metal_females_no_UKB.txt
metal metal_males_no_UKB.txt

# note likely still some multi allelics in there...
# so some more filtering:
cd /PATH/TO/GENERAL/
cut -f 1,2 HRC.r1-1.GRCh37.wgs.mac5.sites.tab | sed -e 's/\t/:/g' | sort | uniq -u > no_multi_allelics_HRC.txt
# keep for later filtering....
```


### Post metal filtering and sorting

```
# to sort meta-analyzed file
sort -gk 10 META_NEW_MALE1.tbl | head
# SNCA and MAPT top-hits => Yippieeee

sort -gk 10 META_NEW_FEMALE1.tbl | head
# SNCA and MAPT top-hits => Yippieeee

# filter for MAF >1%, at least present in 13 datasets out of 19 (~2/3) this makes sure NeuroX isnt a problem, and I2 of 80%
head -1 META_NEW_MALE1.tbl > header.txt
awk '{if ($6 > 0.01) print $0;}' META_NEW_MALE1.tbl | awk '{if ($6 < 0.99) print $0;}' | awk '{if ($14 > 12) print $0;}' | awk '{if ($12 < 80) print $0;}' > temp
cat header.txt temp > MALE_PD_filtered_sumstats.txt

head -1 META_NEW_FEMALE1.tbl > header.txt
awk '{if ($6 > 0.01) print $0;}' META_NEW_FEMALE1.tbl | awk '{if ($6 < 0.99) print $0;}' | awk '{if ($14 > 12) print $0;}' | awk '{if ($12 < 80) print $0;}' > temp
cat header.txt temp > FEMALE_PD_filtered_sumstats.txt

# MAF 5% just to check lambda
head -1 META_NEW_MALE1.tbl > header.txt
awk '{if ($6 > 0.05) print $0;}' META_NEW_MALE1.tbl | awk '{if ($6 < 0.95) print $0;}' | awk '{if ($14 > 12) print $0;}' | awk '{if ($12 < 80) print $0;}' > temp
cat header.txt temp > MALE_PD_filtered_sumstats_5percent_MAF.txt

# MAF 5% just to check lambda
head -1 META_NEW_FEMALE1.tbl > header.txt
awk '{if ($6 > 0.05) print $0;}' META_NEW_FEMALE1.tbl | awk '{if ($6 < 0.95) print $0;}' | awk '{if ($14 > 12) print $0;}' | awk '{if ($12 < 80) print $0;}' > temp
cat header.txt temp > FEMALE_PD_filtered_sumstats_5percent_MAF.txt

# filter no proxy ones...
head -1 META_NEW_MALE_NO_PROXY1.tbl > header.txt
awk '{if ($6 > 0.01) print $0;}' META_NEW_MALE_NO_PROXY1.tbl | awk '{if ($6 < 0.99) print $0;}' | awk '{if ($14 > 12) print $0;}' | awk '{if ($12 < 80) print $0;}' > temp
cat header.txt temp > MALE_PD_filtered_sumstats_NO_PROXY.txt

head -1 META_NEW_FEMALE_NO_PROXY1.tbl > header.txt
awk '{if ($6 > 0.01) print $0;}' META_NEW_FEMALE_NO_PROXY1.tbl | awk '{if ($6 < 0.99) print $0;}' | awk '{if ($14 > 12) print $0;}' | awk '{if ($12 < 80) print $0;}' > temp
cat header.txt temp > FEMALE_PD_filtered_sumstats_NO_PROXY.txt

# filter no UKB data at all...
head -1 META_NEW_MALE_NO_UKB_AT_ALL1.tbl > header.txt
awk '{if ($6 > 0.01) print $0;}' META_NEW_MALE_NO_UKB_AT_ALL1.tbl | awk '{if ($6 < 0.99) print $0;}' | awk '{if ($14 > 11) print $0;}' | awk '{if ($12 < 80) print $0;}' > temp
cat header.txt temp > MALE_PD_filtered_sumstats_NO_UKB_AT_ALL.txt

head -1 META_NEW_FEMALE_NO_UKB_AT_ALL1.tbl > header.txt
awk '{if ($6 > 0.01) print $0;}' META_NEW_FEMALE_NO_UKB_AT_ALL1.tbl | awk '{if ($6 < 0.99) print $0;}' | awk '{if ($14 > 11) print $0;}' | awk '{if ($12 < 80) print $0;}' > temp
cat header.txt temp > FEMALE_PD_filtered_sumstats_NO_UKB_AT_ALL.txt


# another quick filter for multiallelics...
module load R
R
require("data.table")
# load in data
male <- fread("MALE_PD_filtered_sumstats.txt",header=T)
male2 <- fread("MALE_PD_filtered_sumstats_5percent_MAF.txt",header=T)
male3 <- fread("MALE_PD_filtered_sumstats_NO_PROXY.txt",header=T)
male4 <- fread("MALE_PD_filtered_sumstats_NO_UKB_AT_ALL.txt",header=T)
female <- fread("FEMALE_PD_filtered_sumstats.txt",header=T)
female2 <- fread("FEMALE_PD_filtered_sumstats_5percent_MAF.txt",header=T)
female3 <- fread("FEMALE_PD_filtered_sumstats_NO_PROXY.txt",header=T)
female4 <- fread("FEMALE_PD_filtered_sumstats_NO_UKB_AT_ALL.txt",header=T)
HRC <- fread("/PATH/TO/GENERAL/no_multi_allelics_HRC.txt",header=F)
# merge files
MM1 <- merge(HRC,male, by.x="V1",by.y="MarkerName")
MM2 <- merge(HRC,male2, by.x="V1",by.y="MarkerName")
MM3 <- merge(HRC,male3, by.x="V1",by.y="MarkerName")
MM4 <- merge(HRC,male4, by.x="V1",by.y="MarkerName")
FF1 <- merge(HRC,female, by.x="V1",by.y="MarkerName")
FF2 <- merge(HRC,female2, by.x="V1",by.y="MarkerName")
FF3 <- merge(HRC,female3, by.x="V1",by.y="MarkerName")
FF4 <- merge(HRC,female4, by.x="V1",by.y="MarkerName")
# rename column name
colnames(MM1)[1] <- "MarkerName"
colnames(MM2)[1] <- "MarkerName"
colnames(MM3)[1] <- "MarkerName"
colnames(MM4)[1] <- "MarkerName"
colnames(FF1)[1] <- "MarkerName"
colnames(FF2)[1] <- "MarkerName"
colnames(FF3)[1] <- "MarkerName"
colnames(FF4)[1] <- "MarkerName"
# save files
write.table(MM1, file="MALE_PD_filtered_sumstats_no_multi_allelics.txt", quote = F, sep = "\t", row.names = F)
write.table(MM2, file="MALE_PD_filtered_sumstats_5percent_MAF_no_multi_allelics.txt", quote = F, sep = "\t", row.names = F)
write.table(MM3, file="MALE_PD_filtered_sumstats_NO_PROXY_MAF_no_multi_allelics.txt", quote = F, sep = "\t", row.names = F)
write.table(MM4, file="MALE_PD_filtered_sumstats_NO_UKB_AT_ALL_no_multi_allelics.txt", quote = F, sep = "\t", row.names = F)
write.table(FF1, file="FEMALE_PD_filtered_sumstats_no_multi_allelics.txt", quote = F, sep = "\t", row.names = F)
write.table(FF2, file="FEMALE_PD_filtered_sumstats_5percent_MAF_no_multi_allelics.txt", quote = F, sep = "\t", row.names = F)
write.table(FF3, file="FEMALE_PD_filtered_sumstats_NO_PROXY_no_multi_allelics.txt", quote = F, sep = "\t", row.names = F)
write.table(FF4, file="FEMALE_PD_filtered_sumstats_NO_UKB_AT_ALL_no_multi_allelics.txt", quote = F, sep = "\t", row.names = F)
q()
n

wc -l MALE_PD_filtered_sumstats_no_multi_allelics.txt  # 7153507
wc -l MALE_PD_filtered_sumstats.txt  # 7172164
wc -l MALE_PD_filtered_sumstats_5percent_MAF_no_multi_allelics.txt  # 5268882
wc -l MALE_PD_filtered_sumstats_5percent_MAF.txt # 5282609
wc -l MALE_PD_filtered_sumstats_NO_PROXY_MAF_no_multi_allelics.txt # 7016524
wc -l MALE_PD_filtered_sumstats_NO_UKB_AT_ALL_no_multi_allelics.txt # 7234142

wc -l FEMALE_PD_filtered_sumstats_no_multi_allelics.txt  # 7141404
wc -l FEMALE_PD_filtered_sumstats.txt  # 7141404
wc -l FEMALE_PD_filtered_sumstats_5percent_MAF_no_multi_allelics.txt  # 5269024
wc -l FEMALE_PD_filtered_sumstats_5percent_MAF.txt  # 5269024
wc -l FEMALE_PD_filtered_sumstats_NO_PROXY_no_multi_allelics.txt # 7076952
wc -l FEMALE_PD_filtered_sumstats_NO_UKB_AT_ALL_no_multi_allelics.txt # 7309890


# move forward with:
MALE_PD_filtered_sumstats_no_multi_allelics.txt
MALE_PD_filtered_sumstats_5percent_MAF_no_multi_allelics.txt
MALE_PD_filtered_sumstats_NO_PROXY_MAF_no_multi_allelics.txt
MALE_PD_filtered_sumstats_NO_UKB_AT_ALL_no_multi_allelics.txt
FEMALE_PD_filtered_sumstats_no_multi_allelics.txt
FEMALE_PD_filtered_sumstats_5percent_MAF_no_multi_allelics.txt
FEMALE_PD_filtered_sumstats_NO_PROXY_no_multi_allelics.txt
FEMALE_PD_filtered_sumstats_NO_UKB_AT_ALL_no_multi_allelics.txt

```



### Additional QC on meta-analyses results

```
# lambda in R

cd /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/RESULTS/to_meta_files/

module load R
R
require("data.table")
data = fread("MALE_PD_filtered_sumstats_no_multi_allelics.txt",header=T)
#data = fread("FEMALE_PD_filtered_sumstats_no_multi_allelics.txt",header=T)
#data = fread("MALE_PD_filtered_sumstats_5percent_MAF_no_multi_allelics.txt",header=T)
#data = fread("FEMALE_PD_filtered_sumstats_5percent_MAF_no_multi_allelics.txt",header=T)
#data = fread("MALE_PD_filtered_sumstats_NO_PROXY_MAF_no_multi_allelics.txt",header=T)
#data = fread("FEMALE_PD_filtered_sumstats_NO_PROXY_no_multi_allelics.txt",header=T)
#data = fread("MALE_PD_filtered_sumstats_NO_UKB_AT_ALL_no_multi_allelics.txt",header=T)
#data = fread("FEMALE_PD_filtered_sumstats_NO_UKB_AT_ALL_no_multi_allelics.txt",header=T)


p <- data$"P-value"
n <- length(data$P)
x2obs <- qchisq(as.numeric(as.character(p)), 1, lower.tail = FALSE)
x2exp <- qchisq((1:n)/n, 1, lower.tail = FALSE)
lambda <- median(x2obs)/median(x2exp)
lambda # the closer to 1 to better
# include UKB proxy
#lambda1000_males <- 1 + ( lambda -1 ) * (1/20956 + 1/89660)/(1/1000 + 1/1000)
#lambda1000_females <- 1 + ( lambda -1 ) * (1/13420 + 1/90662)/(1/1000 + 1/1000)
# NOT include UKB proxy
#lambda1000_males <- 1 + ( lambda -1 ) * (1/13020 + 1/19336)/(1/1000 + 1/1000)
#lambda1000_females <- 1 + ( lambda -1 ) * (1/7947 + 1/20330)/(1/1000 + 1/1000)
# NOT include UKB at all
#lambda1000_males <- 1 + ( lambda -1 ) * (1/12054 + 1/11999)/(1/1000 + 1/1000)
#lambda1000_females <- 1 + ( lambda -1 ) * (1/7384 + 1/12389)/(1/1000 + 1/1000)

lambda1000

Lambda MAF 1% normal
MALES: 1.04312
FEMALES: 1.019263

Lambda 1000 => MAF 1% normal
MALES: 1.001269
FEMALES: 1.000824

Lambda 5% normal
MALES: 1.066897
FEMALES: 1.031141

Lambda 1000 => MAF 5% normal
MALES: 1.001969
FEMALES: 1.001332

## very reasonable Lambda values (!!)

No proxy beta's
Lambda MAF 1% normal
MALES: 1.065429
FEMALES: 1.038796

Lambda 1000 => MAF 1% normal
MALES: 1.004205
FEMALES: 1.003395

No UKB data at all
Lambda MAF 1% normal
MALES: 1.071307
FEMALES: 1.039275

Lambda 1000 => MAF 1% normal
MALES: 1.005929
FEMALES: 1.004245

```

### QQ plot in R

```
cd /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/RESULTS/to_meta_files/

module load R
R
require("data.table")
data = fread("MALE_PD_filtered_sumstats_no_multi_allelics.txt",header=T)
#data = fread("FEMALE_PD_filtered_sumstats_no_multi_allelics.txt",header=T)
p <- data$"P-value"
observed <- sort(p)
lobs <- -(log10(observed))
expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))
# depending on p-values update the dimensions of the plot which is in the example below -> 10
png("qqplot_male.png")
plot(c(0,40), c(0,40), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,50), ylim=c(0,50), las=1, xaxs="i", yaxs="i", bty="l")
points(lexp, lobs, pch=23, cex=.4, bg="black") 
dev.off()

png("qqplot_female.png")
plot(c(0,40), c(0,40), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,35), ylim=c(0,35), las=1, xaxs="i", yaxs="i", bty="l")
points(lexp, lobs, pch=23, cex=.4, bg="black") 
dev.off()
```



### add RS-IDS for locuszoom and for better usability of data in future...

```
cd /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/RESULTS/to_meta_files/

module load R
R
require("data.table")
male = fread("MALE_PD_filtered_sumstats_no_multi_allelics.txt",header=T)
male2 = fread("MALE_PD_filtered_sumstats_NO_PROXY_MAF_no_multi_allelics.txt",header=T)
male3 = fread("MALE_PD_filtered_sumstats_NO_UKB_AT_ALL_no_multi_allelics.txt",header=T)
female = fread("FEMALE_PD_filtered_sumstats_no_multi_allelics.txt",header=T)
female2 = fread("FEMALE_PD_filtered_sumstats_NO_PROXY_no_multi_allelics.txt",header=T)
female3 = fread("FEMALE_PD_filtered_sumstats_NO_UKB_AT_ALL_no_multi_allelics.txt",header=T)
HRC <- fread("/PATH/TO/GENERAL/no_multi_allelics_HRC.txt",header=F)
RS <- fread("/PATH/TO/GENERAL/HRC_RS_conversion_final.txt",header=T)

# merge
MM <- merge(male, RS, by.x = "MarkerName", by.y = "POS",all.x=TRUE)
MM2 <- merge(male2, RS, by.x = "MarkerName", by.y = "POS",all.x=TRUE)
MM3 <- merge(male3, RS, by.x = "MarkerName", by.y = "POS",all.x=TRUE)
FF <- merge(female, RS, by.x = "MarkerName", by.y = "POS",all.x=TRUE)
FF2 <- merge(female2, RS, by.x = "MarkerName", by.y = "POS",all.x=TRUE)
FF3 <- merge(female3, RS, by.x = "MarkerName", by.y = "POS",all.x=TRUE)

write.table(MM, file="MALE_PD_filtered_sumstats_no_multi_allelics_RSID.txt", quote = F, sep = "\t", row.names = F)
write.table(MM2, file="MALE_PD_filtered_sumstats_NO_PROXY_MAF_no_multi_allelics_RSID.txt", quote = F, sep = "\t", row.names = F)
write.table(MM3, file="MALE_PD_filtered_sumstats_NO_UKB_AT_ALL_no_multi_allelics_RSID.txt", quote = F, sep = "\t", row.names = F)
write.table(FF, file="FEMALE_PD_filtered_sumstats_no_multi_allelics_RSID.txt", quote = F, sep = "\t", row.names = F)
write.table(FF2, file="FEMALE_PD_filtered_sumstats_NO_PROXY_no_multi_allelics_RSID.txt", quote = F, sep = "\t", row.names = F)
write.table(FF3, file="FEMALE_PD_filtered_sumstats_NO_UKB_AT_ALL_no_multi_allelics_RSID.txt", quote = F, sep = "\t", row.names = F)


q()
n

```


```
# store data for long-term use

scp MALE_PD_filtered_sumstats_no_multi_allelics_RSID.txt /PATH/TO/PD/summary_stats/sex_specific/MALE/
scp MALE_PD_filtered_sumstats_NO_PROXY_MAF_no_multi_allelics_RSID.txt /PATH/TO/PD/summary_stats/sex_specific/MALE/
scp MALE_PD_filtered_sumstats_NO_UKB_AT_ALL_no_multi_allelics_RSID.txt /PATH/TO/PD/summary_stats/sex_specific/MALE/
scp FEMALE_PD_filtered_sumstats_no_multi_allelics_RSID.txt /PATH/TO/PD/summary_stats/sex_specific/FEMALE/
scp FEMALE_PD_filtered_sumstats_NO_PROXY_no_multi_allelics_RSID.txt /PATH/TO/PD/summary_stats/sex_specific/FEMALE/
scp FEMALE_PD_filtered_sumstats_NO_UKB_AT_ALL_no_multi_allelics_RSID.txt /PATH/TO/PD/summary_stats/sex_specific/FEMALE/

scp toMeta.*male.tab /PATH/TO/PD/summary_stats/sex_specific/MALE/cohort_level/
scp toMeta.SHORT_UKB_PD_MALE_cases_control_GWAS.txt /PATH/TO/PD/summary_stats/sex_specific/MALE/cohort_level/
scp toMeta.SHORT_UKB_PROXY_FATHER_control_GWAS.txt /PATH/TO/PD/summary_stats/sex_specific/MALE/cohort_level/

scp toMeta.*female.tab /PATH/TO/PD/summary_stats/sex_specific/FEMALE/cohort_level/
scp toMeta.SHORT_UKB_PD_FEMALE_cases_control_GWAS.txt /PATH/TO/PD/summary_stats/sex_specific/FEMALE/cohort_level/
scp toMeta.SHORT_UKB_PROXY_MOTHER_control_GWAS.txt /PATH/TO/PD/summary_stats/sex_specific/FEMALE/cohort_level/

scp META_NEW_MALE1.tbl.info /PATH/TO/PD/summary_stats/sex_specific/MALE/
scp META_NEW_MALE_NO_PROXY1.tbl.info /PATH/TO/PD/summary_stats/sex_specific/MALE/
scp META_NEW_MALE_NO_UKB_AT_ALL1.tbl.info /PATH/TO/PD/summary_stats/sex_specific/MALE/

scp META_NEW_FEMALE1.tbl.info /PATH/TO/PD/summary_stats/sex_specific/FEMALE/
scp META_NEW_FEMALE_NO_PROXY1.tbl.info /PATH/TO/PD/summary_stats/sex_specific/FEMALE/
scp META_NEW_FEMALE_NO_UKB_AT_ALL1.tbl.info /PATH/TO/PD/summary_stats/sex_specific/FEMALE/

```


### META5 variants

```
cd /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/RESULTS/to_meta_files/

module load R
R
require("data.table")
male = fread("MALE_PD_filtered_sumstats_no_multi_allelics_RSID.txt",header=T)
female = fread("FEMALE_PD_filtered_sumstats_no_multi_allelics_RSID.txt",header=T)
data2 = fread("/PATH/TO/GENERAL/META5_GRS_chr_bp.txt",header=F)
MM <- merge(male, data2, by.x = "MarkerName", by.y = "V1")
FF <- merge(female, data2, by.x = "MarkerName", by.y = "V1")
write.table(MM, file="MALE_META5_BETA.txt", quote = F, sep = "\t", row.names = F)
write.table(FF, file="FEMALE_META5_BETA.txt", quote = F, sep = "\t", row.names = F)
q()
n

-----------

cd /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/RESULTS/to_meta_files/

module load R
R
require("data.table")
male = fread("MALE_META5_BETA.txt",header=T)
female = fread("FEMALE_META5_BETA.txt",header=T)
colnames(male)[8] <- "Effect_male"
colnames(male)[9] <- "StdErr_male"
colnames(male)[10] <- "P_male"
colnames(female)[8] <- "Effect_female"
colnames(female)[9] <- "StdErr_female"
colnames(female)[10] <- "P_female"
colnames(female)[20] <- "Effect_META5"

MM <- merge(male, female, by.x = "MarkerName", by.y = "MarkerName")

cor.test(MM$Effect_male,MM$Effect_female)

pdf("male_PD_vs_female_PD.pdf")
plot(MM$Effect_male,MM$Effect_female, pch=16,main="Correlation = 0.951",xlab="Male PD GWAS",ylab="Female PD GWAS")
grid()
dev.off()

cor.test(MM$Effect_male,MM$Effect_female) # 0.9507056
cor.test(MM$Effect_male,MM$Effect_META5) # 0.983635
cor.test(MM$Effect_female,MM$Effect_META5) # 0.9645291 

# more fancy plot

library("ggplot2")
library("data.table")
library("ggrepel")

plotBs <- ggplot(MM, aes(x=Effect_male, y=Effect_female)) +
geom_point() + 
theme_bw() + geom_smooth(se = T, method = lm) + 
xlim(-0.3,0.25) + ylim(-0.4,0.25) + geom_vline(xintercept = 0, linetype = 2) + 
geom_hline(yintercept = 0, linetype = 2)
pdf("BETA_BETA_PLOT_META5_ONLY.pdf", width=10, height=8)
print(plotBs + xlab("Male PD GWAS") + ylab("Female PD GWAS")) 
graphics.off()

# save data for table
write.table(MM, file="DATA_FOR_INPUT_SUP_TABLE2_BETA_BETA.txt", quote = F, sep = "\t", row.names = F)


```

### Mirrored Manhattan plot

```
cd /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/RESULTS/to_meta_files/

# mirror Manhattan plot

# prep data:
# males
cut -f 1 MALE_PD_filtered_sumstats_no_multi_allelics.txt > VARIANT.txt
cut -f 1 MALE_PD_filtered_sumstats_no_multi_allelics.txt | sed -e 's/:/\t/g' > VARIANT2.txt
cut -f 10 MALE_PD_filtered_sumstats_no_multi_allelics.txt > P.txt
paste VARIANT.txt VARIANT2.txt P.txt | grep -v "MarkerName" > input_males_mirror_plot.txt
# females
cut -f 1 FEMALE_PD_filtered_sumstats_no_multi_allelics.txt > VARIANT.txt
cut -f 1 FEMALE_PD_filtered_sumstats_no_multi_allelics.txt | sed -e 's/:/\t/g' > VARIANT2.txt
cut -f 10 FEMALE_PD_filtered_sumstats_no_multi_allelics.txt > P.txt
paste VARIANT.txt VARIANT2.txt P.txt | grep -v "MarkerName" > input_females_mirror_plot.txt

# run mirror plot code

module load R
R
devtools::install_github('anastasia-lucas/hudson')
library(hudson)
require("data.table")
data.t = fread("input_males_mirror_plot.txt",header=F)
data.b = fread("input_females_mirror_plot.txt",header=F)
colnames(data.t) <- c("SNP", "CHR", "POS", "pvalue")
colnames(data.b) <- c("SNP", "CHR", "POS", "pvalue")
# format:
#SNP CHR  POS     pvalue
#1 rs1   1    1 0.06396809
#2 rs2   1  201 0.90213420
gmirror(top=data.t, bottom=data.b, tline=5E-8, bline=5E-8, 
toptitle="PD Male GWAS", bottomtitle = "PD Female GWAS", 
highlight_p = c(5E-8,5E-8), highlighter="green")
q()
n

# files saved as => gmirror.png

# add in manually the hits...

```

### Searching for potential differences between males and females...

```
cd /PATH/TO/CORNELIS_TEMP/MALE_FEMALE_PD/RESULTS/to_meta_files/

# prep data....
awk '{if ($10 < 0.00000005) print $0;}' MALE_PD_filtered_sumstats_no_multi_allelics.txt | cut -f 1,2,3,4,8,9,10,11 > temp1.txt 
awk '{if ($10 < 0.00000005) print $0;}' FEMALE_PD_filtered_sumstats_no_multi_allelics.txt | cut -f 1,2,3,4,8,9,10,11 > temp2.txt 
head -1 FEMALE_PD_filtered_sumstats_no_multi_allelics.txt | cut -f 1,2,3,4,8,9,10,11 > header.txt
cat header.txt temp1.txt > MALE_PD_genome_wide_signit.txt
cat header.txt temp2.txt > FEMALE_PD_genome_wide_signit.txt

# wc -l MALE_PD_genome_wide_signit.txt # 3384
# wc -l FEMALE_PD_genome_wide_signit.txt # 3232

# got to R and merge
R
require("data.table")
male_gwas = fread('MALE_PD_filtered_sumstats_no_multi_allelics.txt',header=T)
females_gwas = fread('FEMALE_PD_filtered_sumstats_no_multi_allelics.txt',header=T)
male_signi = fread('MALE_PD_genome_wide_signit.txt',header=T)
female_signi = fread('FEMALE_PD_genome_wide_signit.txt',header=T)
MM <- merge(male_gwas, female_signi, by.x = "MarkerName", by.y = "MarkerName")
MM2 <- merge(females_gwas, male_signi, by.x = "MarkerName", by.y = "MarkerName")

# MM
beta1 = MM$Effect.x
beta2 = MM$Effect.y
SE1 = MM$StdErr.x
SE2 = MM$StdErr.y
Z = (beta1-beta2) / sqrt(SE1^2 + SE2^2)
MM$P = 2*pnorm(-abs(Z))

# MM2
beta1 = MM2$Effect.x
beta2 = MM2$Effect.y
SE1 = MM2$StdErr.x
SE2 = MM2$StdErr.y
Z = (beta1-beta2) / sqrt(SE1^2 + SE2^2)
MM2$P = 2*pnorm(-abs(Z))

write.table(MM, file="male_gwas_with_female_signi_hits.txt", quote = F, sep = "\t", row.names = F)
write.table(MM2, file="female_gwas_with_male_signi_hits.txt", quote = F, sep = "\t", row.names = F)

q()
n
# close and restart

R
require("data.table")

male_gwas = fread('male_gwas_with_female_signi_hits.txt',header=T)
female_gwas = fread('female_gwas_with_male_signi_hits.txt',header=T)

pdf("compare_SIGNI_HITS.pdf")

plot(male_gwas$Effect.x,male_gwas$Effect.y,pch=16,main="Female PD GWAS 5E-8 variants in male PD GWAS",xlab="Beta Male GWAS", ylab="Beta Female GWAS")
grid()
plot(female_gwas$Effect.x,female_gwas$Effect.y,pch=16,main="Male PD GWAS 5E-8 variants in female PD GWAS",xlab="Beta Female GWAS", ylab="Beta male GWAS")
grid()
dev.off()

# fancy plots.....

library("ggplot2")
library("data.table")
library("ggrepel")

plotBs <- ggplot(male_gwas, aes(x=Effect.x, y=Effect.y)) +
geom_point() + 
theme_bw() + geom_smooth(se = T, method = lm) + 
xlim(-0.65,0.55) + ylim(-1.2,1) + geom_vline(xintercept = 0, linetype = 2) + 
geom_hline(yintercept = 0, linetype = 2)
pdf("male_gwas_with_female_signi_hits.pdf", width=10, height=8)
print(plotBs + xlab("Beta Male GWAS") + ylab("Beta Female GWAS") + ggtitle("Female PD GWAS 5E-8 variants in male PD GWAS")) 
graphics.off()


plotBs <- ggplot(female_gwas, aes(x=Effect.x, y=Effect.y)) +
geom_point() + 
theme_bw() + geom_smooth(se = T, method = lm) + 
xlim(-0.65,0.55) + ylim(-1.2,1) + geom_vline(xintercept = 0, linetype = 2) + 
geom_hline(yintercept = 0, linetype = 2)
pdf("female_gwas_with_male_signi_hits.pdf", width=10, height=8)
print(plotBs + xlab("Beta Female GWAS") + ylab("Beta Male GWAS") + ggtitle("Male PD GWAS 5E-8 variants in female PD GWAS")) 
graphics.off()

```


![Alt Text](https://media1.tenor.com/images/ee479a8148cf18ba7582bc16e891eec1/tenor.gif?itemid=16045121)


