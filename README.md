# Blogs_tips
## Subset bamfile with chromosome names and convert into paired fastq  
* software required:[sambamba](https://github.com/lomereiter/sambamba) and [bam2fastx](https://github.com/infphilo/tophat) from tophat binary distribution.<br>
sambamba usages should refer to https://github.com/lomereiter/sambamba/wiki/%5Bsambamba-view%5D-Filter-expression-syntax#basic-conditions-for-fields

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
PS: the numbers specified in `ref_id` means the ref order list in header from bamfle, which can be checked by `samtools view -H your.bam` if samtools installed 
