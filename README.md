```python

```

## 1. Data
<font style="font-family: Verdana; font-size:1.2em;">
 
### Sample info

| Type | value |
| --- | --- |
| Sample type | DNA from FFPE blocks |
| No. of samples | 40 |
| Capture type | Whole exome |
| Sequencing platform | Illumina HiSeq |
| Library type | Paired end (150*2) |

<br>

    
=======================================================================
<br>**Library preparation**<br>
SureselectXT Target Enrichment system kit


**Adapters used in this study are:**<br>
**Illumina Universal Adapters**<br>
`5’-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT-3’` <br><br>
**Index Adapter:** <br>
`5’-GATCGGAAGAGCACACGTCTGAACTCCAGTCAC [INDEX] ATCTCGTATGCCGTCTTCTGCTTG-3’` <br>

=======================================================================
    
</font>

### Reference Genome
<font style="font-family: Verdana; font-size:1.2em;">
    
**Reference Genome :** hg19 <br>
<font color=brown>**1. Download reference genome data using ftp**</font>

&emsp;&emsp;&emsp;&emsp;location: ftp.broadinstitute.org/bundle <br>
&emsp;&emsp;&emsp;&emsp;username: gsapubftp-anonymous <br>
&emsp;&emsp;&emsp;&emsp;password:<br>
    
<font color=brown>**2. Access data via google bucket** </font><br>
**hg19**<br>https://storage.cloud.google.com/broad-references/hg19/v0/Homo_sapiens_assembly19.fasta<br>

**hg38**<br>
https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/<br>
    

<font color=brown>**3. Web browser** </font><br>
Broad Institute <br>
ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/ <br>
<br>UCSC Genome<br>
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/ <br>
<br>Igenome Illumina<br>
http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Homo_sapiens/Ensembl/GRCh37/Homo_sapiens_Ensembl_GRCh37.tar.gz  <br>

=======================================================================
    
</font>


## 2. Softwares

<font style="font-family: Verdana; font-size:1.2em;">
    
<br>    
<font color=brown>**FASTQC** </font><br>
Home Page: <br>
https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ <br>
    
Download Location: <br>
https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip <br>
    
Manual Page: <br>
https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/ <br>

=======================================================================

<font color=brown>**Trim Galore** </font><br>
Download Location: <br>
https://github.com/FelixKrueger/TrimGalore <br>

Manual Page: <br>
https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md <br>

<font color=brown>**Cut Adapt** </font><br>
Download Location: <br>
https://github.com/marcelm/cutadapt/ <br>
run in the shell `pip install cutadapt`

Manual Page: <br>
https://cutadapt.readthedocs.io/en/stable/

=======================================================================

<font color=brown>**GATK** </font><br>
Download Location:<br>
https://github.com/broadinstitute/gatk/releases/download/4.1.4.1/gatk-4.1.4.1.zip <br>

Manual Page:<br>
https://software.broadinstitute.org/gatk/documentation/ 

=======================================================================

<font color=brown>**BWA** </font><br>
Download Location:<br>
https://liquidtelecom.dl.sourceforge.net/project/bio-bwa/bwa-0.7.17.tar.bz2<br>

Manual page:<br>
http://bio-bwa.sourceforge.net/<br>
http://bio-bwa.sourceforge.net/bwa.shtml<br>
    
<font color=brown>**Bowtie2**</font> <br>
Home Page: <br>
http://bowtie-bio.sourceforge.net/bowtie2/index.shtml <br>

Download Page: <br>
https://sourceforge.net/projects/bowtie-bio/files/latest/download <br>
https://liquidtelecom.dl.sourceforge.net/project/bowtie-bio/bowtie/1.2.3/bowtie-src-x86_64.zip <br>

=======================================================================
   
<font color=brown>**PICARD** </font><br>
Download Location:<br>
https://github.com/broadinstitute/picard/releases/download/2.21.4/picard.jar <br>

Manual Page:<br>
https://broadinstitute.github.io/picard/command-line-overview.html <br>

Home Page:<br>
https://broadinstitute.github.io/picard/ 

=======================================================================

<font color=brown>**SAMTools** </font> <br>
Download Location:<br>
https://liquidtelecom.dl.sourceforge.net/project/samtools/samtools/1.10/samtools-1.10.tar.bz2 <br>
     
Manual Page:<br>
http://www.htslib.org/doc/samtools.1.html <br>
http://quinlanlab.org/tutorials/samtools/samtools.html<br>

=======================================================================

<font color=brown>**SAMBAMBA** </font><br>
Download Location:<br>
https://github.com/biod/sambamba<br>
**installation**<br>
`git clone --recursive https://github.com/lomereiter/sambamba.git` <br>
 `cd sambamba` <br>
 `make sambamba-ldmd2-64` <br>

Manual Page:<br>   
https://lomereiter.github.io/sambamba/docs/sambamba-view.html<br>

=======================================================================

<font color=brown>**Bedtool** </font> <br>
Download Location:<br>
https://github.com/arq5x/bedtools2/releases  <br>
**installation**<br>
`wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz` <br>
`tar -zxvf bedtools-2.29.1.tar.gz` <br>
`cd bedtools2` <br>
`make` <br>

Manual Page: <br> 
https://bedtools.readthedocs.io/en/latest/<br>

=======================================================================

<font color=brown>**SNPeff** </font><br>
Download Location:<br>
https://github.com/pcingola/SnpEff <br>
https://excellmedia.dl.sourceforge.net/project/snpeff/snpEff_latest_core.zip <br>
http://snpeff.sourceforge.net/download.html <br>

Manual Page: <br>
http://snpeff.sourceforge.net/SnpEff_manual.html


<font color=brown>**Annovar** </font><br>
http://annovar.openbioinformatics.org/en/latest/ <br> 
**VarAFT** <br>
https://varaft.eu/PackageForInstall/varaft-2.16.deb<br>

=======================================================================
</font>



# Pipeline for variant calling

<font style="font-family: Verdana; font-size:1.2em;">

### 1. Indexing the reference genome

<br>**Using bwa**<br>
`bwa index -a bwtsw hg19.fa`

<br>**Using samtools**<br>
`samtools faidx reference.fa`
    
<br>**Create sequnce dictionary**<br>

`java -jar picard.jar CreateSequenceDictionary`
      `R=hg19.fa`
      `O=hg19.dict`

----------------------------------
`R` : Reference Genome <br>
`O` : Output dictionary

</font>

<font style="font-family: Verdana; font-size:1.2em;">


### 2. Checking the Quality of sequence
`./fastqc sample1_R1.fastq` #forward strand <br>
`./fastqc sample1_R2.fastq` #reverse strand
    
</font>

<font style="font-family: Verdana; font-size:1.2em;">

### 3. Removing the adaptors using trim_galore
**command to run**

`trim_galore \
-–paired \
-–phred33 \
-–quality 20 \
–a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT \
–stringency 5 –e 0.1 \
–t \
–r1 35 \
–r2 35 \
sample1_R1.fastq sample1_R2.fastq`

    
--------------------------------------------------------------------

**parameters**

`--paired` : For paired end reads <br>

--------------------------------------------------------------------
`-–phred33` : Instructs Cutadapt to use ASCII+33 quality scores as Phred scores (Sanger/Illumina 1.9+ encoding) for quality trimming.<br> The other option is `--phred64` : Instructs Cutadapt to use ASCII+64 quality scores as Phred scores (Illumina 1.5 encoding) for quality trimming. <br>

--------------------------------------------------------------------
`-q` or `--quality` : Trim low-quality ends from reads in addition to adapter removal. For RRBS samples, quality trimming will be performed first, and adapter trimming is carried in a second round. Other files are quality and adapter trimmed in a single pass. The algorithm is the same as the one used by BWA (Subtract INT from all qualities; compute partial sums from all indices to the end of the sequence; cut sequence at the index at which the sum is minimal).

--------------------------------------------------------------------
`-a` or `--adapter` <br>
`-a2` or `adapter2` <br>

`-a` <br>
Adapter sequence to be trimmed. If not specified explicitly, Trim Galore will try to auto-detect whether the Illumina universal, Nextera transposase or Illumina small RNA adapter sequence was used

`-a2` <br>
Optional adapter sequence to be trimmed off read 2 of paired-end files. This option requires --paired to be specified as well. If the libraries to be trimmed are smallRNA then a2 will be set to the Illumina small RNA 5' adapter automatically (`GATCGTCGGACT`). A single base may also be given as e.g. `-a2` A{10}, to be expanded to `-a2` `AAAAAAAAAA`

--------------------------------------------------------------------

`-t` or `--trim1` : Trims 1 bp off every read from its 3' end.
This may be needed for FastQ files that are to be aligned as paired-end data with Bowtie 1. This is because Bowtie (1) regards alignments like this as invalid (whenever a start/end coordinate is contained within the other read):

--------------------------------------------------------------------

`-r1` or `--length_1` : Unpaired single-end read length cutoff needed for read 1 to be written to .unpaired_1.fq output file. These reads may be mapped in single-end mode.
Default: 35 bp

`-r2` or `--length_2` : Unpaired single-end read length cutoff needed for read 2 to be written to .unpaired_2.fq output file. These reads may be mapped in single-end mode.
Default: 35 bp 

--------------------------------------------------------------------

`-s` or `--stringency` : Overlap with adapter sequence required to trim a sequence.
Defaults to a very stringent setting of 1, i.e. even a single base pair of overlapping sequence will be trimmed of the 3' end of any read.

    
</font>

<font style="font-family: Verdana; font-size:1.2em;">

### 4. Re-checking the Quality of sequence

`./fastqc sample1_R1_trimmed.fastq` #forward strand <br>
`./fastqc sample1_R2_trimmed.fastq` #reverse strand
    
</font>


```python

```

<font style="font-family: Verdana; font-size:1.2em;">

### 5. Map to Reference

`bwa mem \
-M \
-t 8 \
hg19.fa \
sample1_R1_trimmed.fastq sample1_R2_trimmed.fastq > aligned_paired_reads.sam` <br>

**Converting SAM to BAM** <br>
`samtools view -Sb  aligned_paired_reads.sam  >  aligned_paired_reads.bam` <br>

<font color=brown>**Both the above steps can be achieved in following single step**</font>
    
`bwa mem \
-M \
-t 8 \
hg19.fa \
sample1_R1_trimmed.fastq sample1_R2_trimmed.fastq | samtools view -1 -bS - > aligned_paired_reads.bam`

-----------------------------

`-M` : Mark shorter split hits as secondary (for Picard compatibility).<br>
`-t` : Number of threads<br>
`-R` : Complete read group header line. '\t' can be used in STR and will be converted to a TAB in the output SAM.
       &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The read group ID will be attached to every read in the output. An example is ’@RG\tID:foo\tSM:bar’.<br>
`samtools view -1 -bS`: to sort and compress your sam file to the bam format

</font>

<font style="font-family: Verdana; font-size:1.2em;">

#### Validate BAM file

`java -jar picard.jar ValidateSamFile \
    I=aligned_paired_reads.bam \
    MODE=SUMMARY`

</font>

<font style="font-family: Verdana; font-size:1.2em;">

### 6. MergeBamAlignments
`java -jar picard.jar MergeBamAlignment \
       ALIGNED=aligned_paired_reads.bam \ 
       UNMAPPED=unmapped.bam \ 
       O=merge_alignments.bam \
       R=hg19.fasta`
       
------------------------

Merge alignment data from a `SAM` or `BAM` with data in an `unmapped BAM` file. This tool produces a new `SAM` or `BAM` file that includes all aligned and unaligned reads and also carries forward additional read attributes from the `unmapped BAM` (attributes that are otherwise lost in the process of alignment). The purpose of this tool is to use information from the `unmapped BAM` to fix up aligner output. The resulting file will be valid for use by other Picard tools. For simple BAM file merges, use `MergeSamFiles`. Note that `MergeBamAlignment` expects to find a sequence dictionary in the same directory as `REFERENCE_SEQUENCE` and expects it to have the same base name as the reference FASTA except with the extension ".dict". If the output sort order is not coordinate, then reads that are clipped due to adapters or overlapping will not contain the `NM`, `MD`, or `UQ` tags
<br>[Link](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_sam_MergeBamAlignment.php)
    
</font>

<font style="font-family: Verdana; font-size:1.2em;">

### 7. MarkDuplicates

`java -jar picard.jar MarkDuplicates \
      I=merge_alignments.bam \
      O=marked_duplicates.bam \
      M=marked_dup_metrics.txt`
      
---------------------------
[Link](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_sam_markduplicates_MarkDuplicates.php)

</font>

<font style="font-family: Verdana; font-size:1.2em;">

### 8. SortSam
`java -jar picard.jar SortSam \
      I=marked_duplicates.bam \
      O=sorted.bam \
      SORT_ORDER=coordinate`

---------------------

[Link](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_sam_SortSam.php)
</font>

<font style="font-family: Verdana; font-size:1.2em;">

### 9. Get stats
#### Alignment Summary Matrix

`java -jar picard.jar CollectAlignmentSummaryMetrics \
    REFERENCE=hg19.fa \
    INPUT=sorted.bam \
    OUTPUT=alignment_summary_matirx.txt`

[Link](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_analysis_CollectAlignmentSummaryMetrics.php)

---------------------------
#### Insert Size Matrix

`java -jar picard.jar CollectInsertSizeMetrics \
      I=sorted.bam \
      O=insert_size_metrics.txt \
      H=insert_size_histogram.pdf \
      M=0.5`

[Link](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_analysis_CollectInsertSizeMetrics.php)

---------------------------
#### Depth

`samtools depth -a sorted_reads.bam > depth_out.txt`

[Link](http://www.htslib.org/doc/samtools-depth.1.html)

---------------------------    
    
</font>


```python

```

<font style="font-family: Verdana; font-size:1.2em;">

### 10. Base (Quality Score) Recalibration
    
`java -jar gatk.jar BaseRecalibrator \
   -I sorted.bam \
   -R hg19.fa \
   --known-sites sites_of_variation.vcf \
   --known-sites another/optional/setOfSitesToMask.vcf \
   -O recal_data.table`
 <br><br>parameter
`--known-sites` : vcf of `ExAc`, `gnomAD`, or `dbSNP` etc
    
-----------------------------
    
[Link](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_bqsr_BaseRecalibrator.php)

#### Apply Recalibration

`gatk ApplyBQSR \
   -R hg19.fa \
   -I sorted.bam \
   --bqsr-recal-file recal_data.table \
   -O recalibrated.bam`

------------------------------

[Link](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_bqsr_ApplyBQSR.php)
    

</font>

<font style="font-family: Verdana; font-size:1.2em;">

### 11. AnalyzeCovariates

[Link](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_bqsr_AnalyzeCovariates.php)

#### Generate the first pass recalibration table file

`java -jar GenomeAnalysisTK.jar \
      -T BaseRecalibrator \
      -R hg19.fa \
      -I sorted.bam \
      -knownSites my-trusted-snps.vcf \  
      -knownSites my-trusted-indels.vcf \ 
      -o firstpass.table`

-----------------------------------
#### Generate the second pass recalibration table file

`java -jar GenomeAnalysisTK.jar \
      -T BaseRecalibrator \
      -R hg19.fa \
      -I sorted.bam \
      -knownSites my-trusted-snps.vcf \
      -knownSites my-trusted-indels.vcf \
      -BQSR firstpass.table \
      -o secondpass.table`

-----------------------------------
#### Finally generate the plots and also keep a copy of the csv (optional)

`java -jar GenomeAnalysisTK.jar \
      -T AnalyzeCovariates \
      -R hg19.fa \
      -before firstpass.table \
      -after secondpass.table \
      -csv BQSR.csv \
      -plots BQSR.pdf`

-----------------------------------
[Link](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_bqsr_AnalyzeCovariates.php)

</font>

<font style="font-family: Verdana; font-size:1.2em;">

### 12. Build BAM Index

`java -jar picard.jar \
    BuildBamIndex INPUT=recalibrated.bam`

---------------------------------------
    
[Link](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_sam_BuildBamIndex.php)

</font>

<font color=red>=======================================================================</font>

<font style="font-family: Verdana; font-size:1.2em;color:red">

#### RealignerTargetCreator

`java -jar GenomeAnalysisTK.jar \
   -T RealignerTargetCreator \
   -R hg19.fa \
   -I recalibrated.bam \
   --known indels.vcf \
   -o forIndelRealigner.intervals`
    
----------------------------    
[Link](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_indels_RealignerTargetCreator.php)
<br>**This is removed in latest GATK** <br> [Link](https://www.biostars.org/p/309094/)
</font>

<font style="font-family: Verdana; font-size:1.2em;color:red">

#### IndelRealigner
`java -jar GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R hg19.fa -I recalibrated.bam \
    -targetIntervals realignment_targets.list \
    -o realigned_reads.bam`

----------------------------  
[Link](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_indels_IndelRealigner.php)
    
**This is removed in latest GATK** <br>[Link](https://www.biostars.org/p/309094/)
</font>

<font color=red>=======================================================================</font>

<font style="font-family: Verdana; font-size:1.2em;">
    
### 13. Call Variants (HaplotypeCaller)

#### Single-sample GVCF calling (outputs intermediate GVCF)
    
`java -jar GenomeAnalysisTK.jar \
   -T HaplotypeCaller  \
   -R hg19.fa \
   -I recalibrated.bam \
   -O raw_variants.g.vcf.gz \
   -ERC GVCF`
    
------------------------------
#### Single-sample GVCF calling with allele-specific annotations

`java -jar GenomeAnalysisTK.jar \
   -T HaplotypeCaller  \
   -R hg19.fa \
   -I recalibrated.bam \
   -O raw_variants.g.vcf.gz \
   -ERC GVCF \
   -G Standard \
   -G AS_Standard`

------------------------------
#### Variant calling with bamout to show realigned reads
 
`java -jar GenomeAnalysisTK.jar \
   -T HaplotypeCaller \
   -R hg19.fa \
   -I recalibrated.bam \
   -O raw_variants.vcf.gz \
   -bamout bamout.bam`
    
<br>[Link](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.php)
<br>[Link2](https://software.broadinstitute.org/gatk/documentation/article?id=9622)
<br>[link3](https://software.broadinstitute.org/gatk/documentation/article?id=5484)
</font>

<font style="font-family: Verdana; font-size:1.2em;">

### 14. Extract SNPs & Indels

#### SNP    
`java -jar GenomeAnalysisTK.jar \
    -T SelectVariants \
    -R hg19.fa \
    -V raw_variants.vcf \
    -selectType SNP \
    -O raw_variants_snp.vcf`

------------------------------------------
#### Indels

`java -jar GenomeAnalysisTK.jar \
    -T SelectVariants \
    -R hg19.fa \
    -V raw_variants.vcf \
    -selectType INDEL \
    -O raw_variants_indel.vcf`

-----------------------------------------

[Link](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_variantutils_SelectVariants.php)    
</font>

<font style="font-family: Verdana; font-size:1.2em;">

### 15. Filter SNPs & Indels
#### SNPs

`java -jar GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R hg19.fa \
    -V raw_variants_snp.vcf \
    --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' \
    --filterName "basic_snp_filter" -o filtered_snps.vcf`
                                                           
------------------------------------------------------
[Link](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_filters_VariantFiltration.php)                                                           

`QD`,`FS`,`MQ` [Link](https://software.broadinstitute.org/gatk/documentation/article.php?id=1255)                               

------------------------------------------------------
#### Indels

`java -jar GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R hg19.fa \
    -V raw_variants_indel.vcf \
    --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' \
    --filterName "basic_indel_filter" -o filtered_indels.vcf`
    
----------------------------------------------------    



</font>


```python

```

<font style="font-family: Verdana; font-size:1.2em;">

### 16. Annotate using snpEff


`java -jar snpEff.jar \
    -v snpeff_db filtered_snps.vcf > filtered_snps_final.ann.vcf`

</font>

### 17. Coverage statistics
`bedtools genomecov -bga -ibam recal_reads.bam > genomecov.bedgraph`


```python

```
