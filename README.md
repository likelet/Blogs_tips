# Blogs_tips

## Table of Contents

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->

- [Subset bamfile with chromosome names and convert into paired fastq](#Subset-bamfile-with-chromosome-names-and-convert-into-paired-fastq)
- [Cluster management](#Cluster-management)
- [R code for ploting nomograph from competing risk survival analysis model](#R-code-for-ploting-nomograph-from-competing-risk-survival-analysis-model)
- [Setting docker download mirror site](#Setting-docker-download-mirror-site)
- [Install bioconductor R package using VPS](#Install-bioconductor-R-package-using-VPS)
- [install bioconductor R package using mirror at UTSC](#install-bioconductor-R-package-using-mirror-at-UTSC)



## Subset bamfile with chromosome names and convert into paired fastq  
* software required: **[sambamba](https://github.com/lomereiter/sambamba)** and **[bam2fastx](https://github.com/infphilo/tophat)** from tophat binary distribution.<br>

 > sambamba usages should refer to https://github.com/lomereiter/sambamba/wiki/%5Bsambamba-view%5D-Filter-expression-syntax#basic-conditions-for-fields

```shell 
#using star output bamfile as example 
#!/bin/sh
bamin=$1
#extract reads aligned to chr2
sambamba view -F "ref_id==1" -f bam $bamin -o ${bamin%%Aligned.sortedByCoord.out.bam}_chr2.bam
#sort reads by names if not presorted by software
sambamba sort -n ${bamin%%Aligned.sortedByCoord.out.bam}_chr2.bam -o ${bamin%%Aligned.sortedByCoord.out.bam}_chr2.sort.bam
#bam2fastq
bam2fastx -PANQ -o ${bamin%%Aligned.sortedByCoord.out.bam}_chr2.fq.gz ${bamin%%Aligned.sortedByCoord.out.bam}_chr2.sort.bam

```
**PS**: the numbers specified in `ref_id` means the ref order list in header from bamfle, which can be checked by 
`samtools view -H your.bam` if samtools was installed. 


## Cluster management 
* 1. shudown system 
Shut down computational node 
```shell
#!/bin/sh
for i in `seq 1 3`
do
 ssh cu0$i "hostname;init 0"
done
```
umount storage 
```shell
umount /home
```
shutdown login node 
```shell
poweroff
```
## R code for ploting nomograph from competing risk survival analysis model 
```R
library(cmprsk)
library(rms)
### add path 
setwd("C:\\Users\\hh\\Desktop\\nomo")
rt<-read.csv("Stomach.csv")
rt
View(rt)
attach(rt) 
#change variable names

cov<-cbind(sexC, Age, AJCC_T,AJCC_N,AJCC_M,Surgery)
for (i in 1:6)
{
  cov[,i]<-factor(cov[,i])
}
status<-factor(status)
z <- crr(time,status,cov)
z.p <- predict(z,cov)
n=60#suppose I want to predict the probability of event at time 60(an order)
df<-data.frame(y=z.p[n,-1],cov)
ddist <- datadist(df)  
options(datadist='ddist') 
lmod<-ols(y~(sexC)+(Age)+(AJCC_T)+(AJCC_N)+(AJCC_M)+(Surgery),data=df)#
nom<-nomogram(lmod)
plot(nom,lplabel=paste("prob. of incidence T",round(z.p[n,1],2),sep="="))
```
## Setting docker download mirror site 
Sometimes you may find that it's extrimely painfull to pull docker image from docker.io in china. So this tip can help you to set a mirror site locally in your docker pull command.  
* 1. First, find the file `/etc/docker/daemon.json` and modify it with root authority.
```{javascript}
{
  "registry-mirrors": ["https://registry.docker-cn.com"]
}
```
* 2. Secondly, restart your docker service. 
## Install bioconductor R package using VPS.   

    proxychains4 Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("BSgenome")'

## install bioconductor R package using mirror at UTSC. 

    source("http://bioconductor.org/biocLite.R")
    options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
    biocLite("your package")

## GVC output format information. 
* This tips record the GVC software output format 
1. `cnv.simp`

| Name | Description| Example |
|-----------|--------------|-------------|
|Chromosome1|染色体号|1,2,…,X,Y|
|Position1|起始坐标|1350000|
|Chromosome2|染色体号|1,2,…,X,Y|
|Position2|终止坐标|5345511|
|Value|区域的CNV信号值：LRR=log2(肿瘤样品覆盖度/对照样品覆盖度）；|LOSS: Value < (-0.5);GAIN: Value > 0.5;NORMAL: Value >= (-0.5) & Value <= 0.5|
|Type|突变类型|LOSS：缺失；GAIN：重复；NORMAL：正常；|
|SampleID|样品号|TCGA-21-1081-01A|
|Num_Probes|探测的CNV区域内的单位（如bed区域）个数，若无此值，可设为“-”；|121|
2. `cnv.CBS`

| Name | Description| Example |
|-----------|--------------|-------------|
|SampleID|样品号|TCGA-21-1081-01A|
|Chromosome|染色体号|1,2,…,X,Y|
|start|起始坐标|39175351|
|end|终止坐标|39957876|
|Num_Probes|探测的CNV区域内的单位（如bed区域）个数，若无此值，可设为“-”；|172|
|seg.mean|区域的CNV信号平均值：LRR=log2(肿瘤样品覆盖度/对照样品覆盖度）；|-0.8416|
|seg.sd|区域的CNV信号值的标准差|0.3835|
|seg.median|区域的CNV信号中值：LRR=log2(肿瘤样品覆盖度/对照样品覆盖度）；|-0.8635|
|seg.mad|区域的CNV信号值的平均绝对偏差|0.3432|
|bstat|二元分割的最大统计量|11.07|
|pval|区域分割的pvalue|2.45E-26|
|lcl|α/ 2置信下限（基因组位置）|39946512|
|ucl|α/ 2置信上限（基因组位置）|39981424|
3. somatic `snv.simp` formt 

| Name | Description| Example |
|-----------|--------------|-------------|
|Chromosome|染色体号|1,2,…,X,Y|
|Position|坐标位置|1350000|
|Ref|参考碱基|G|
|Mut|突变碱基|A|
|SampleID|样品号|TCGA-21-1081-01A|
|Depth|肿瘤样品里此点的测序深度|100|
|Alt_count|肿瘤样品里此点突变碱基的测序深度|30|
|MutAF_tumor|肿瘤样品里此点突变碱基的频率|0.3|
|Alleles_tumor|肿瘤样品里此点的碱基及每个碱基的测序深度|G:70,A:30,A:0,G:0|
|Depth_normal|对照样品里此点的测序深度|174|
|Alt_count_normal|对照样品里此点突变碱基的测序深度|0|
|MutAF_normal|对照样品里此点突变碱基的频率|0|
|Alleles_normal|对照样品里此点的碱基及每个碱基的测序深度|C:174,A:0,G:0,T:0|
|consensus_base_tumor|肿瘤样品里此点的一致性碱基|Y|
|consensus_quality_tumor|肿瘤样品里此点的一致性质量|228|
|snp_quality_tumor|肿瘤样品里此点的SNP质量|182|
|Function_region|突变功能区域|exonic,splicing,intergenic和intronic等|
|Mutation_type|突变类型|nonsynonymous,synonymous和stopgain等|
|Gene_name|突变所在的基因的hugo symbol|TP53|
|Gene_ID|突变所在的基因的Ensemble ID|ENSG00000141510|
|Codon_change|Codon改变|AGA->AGT|
|AA_change|氨基酸改变|R280S|
|CDNA_change|cDNA改变|A840T|
|ExonicFunc|Exonic区域突变的详细信息|ENSG00000141510:ENST00000359597:exon7:c.A840T#AGA->AGT#R->S#:p.R280S|
|GeneDetail|||
|Primary_site|COSMIC库注释||
|Cancer_type|||
|COSMIC_ID|||
|1KG_MAF_ALL|千人基因组注释||
|1KG_MAF_AFR|||
|1KG_MAF_AMR|||
|1KG_MAF_EAS|||
|1KG_MAF_EUR|||
|1KG_MAF_SAS|||
|GeneSymbol|clivar库注释||
|ClinicalSignificance|||
|PhenotypeIDs|||
|SIFT_pred|突变的有害性注释||
|Polyphen2_HDIV_pred|||
|Polyphen2_HVAR_pred|||
|LRT_pred|||
|MutationTaster_pred|||
|MutationAssessor_pred|||
|FATHMM_pred|||
|PROVEAN_pred|||
|fathmm-MKL_coding_pred|||
|MetaSVM_pred|||
|MetaLR_pred|||
|GERP++_RS|||
|rsID|dbsnp库注释||
