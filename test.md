```
exec &> $outputDir/$sampleID.StdOp.Err.txt; #Redirect the stdout and stderror to a text file

readGroup="@RG\tID:${sampleID}\tSM:${sampleID}\tLB:${sampleID}\tPL:Illumina"

echo "BWA mem runtime start" `date +%d/%m/%Y\ %H:%M:%S` 2>> $outputDir/$sampleID.StdOp.Err.txt
bwa mem -M -t 16 -R "$readGroup" $reference -o $outputDir/$sampleID/${sampleID}_raw.sam $trimmedR1 $trimmedR2 2>> $outputDir/$sampleID.StdOp.Err.txt

echo "BWA mem runtime end" `date +%d/%m/%Y\ %H:%M:%S` 2>> $outputDir/$sampleID.StdOp.Err.txt

echo "sambamba view and sort runtime start" `date +%d/%m/%Y\ %H:%M:%S` 2>> $outputDir/$sampleID.StdOp.Err.txt
#sambamba to convert sam to sorted bam file
sambamba view --sam-input --nthreads=16 --format=bam --compression-level=0 $outputDir/$sampleID/${sampleID}_raw.sam | sambamba sort --nthreads=16 --out=$outputDir/$sampleID/${sampleID}_sorted.bam /dev/stdin 2>> $outputDir/$sampleID.StdOp.Err.txt
echo "sambamba view and sort runtime end" `date +%d/%m/%Y\ %H:%M:%S` 2>> $outputDir/$sampleID.StdOp.Err.txt

echo "GATK mark Dups runtime start" `date +%d/%m/%Y\ %H:%M:%S` 2>> $outputDir/$sampleID.StdOp.Err.txt
# GATK to mark PCR duplicates
gatk MarkDuplicates -I $outputDir/$sampleID/${sampleID}_sorted.bam --METRICS_FILE $outputDir/$sampleID/${sampleID}_markdups_gatk_metrics.txt -O $outputDir/$sampleID/${sampleID}_sorted_gatkMarkDups.bam --VALIDATION_STRINGENCY LENIENT 2>> $outputDir/$sampleID.StdOp.Err.txt
echo "GATK mark Dups runtime end" `date +%d/%m/%Y\ %H:%M:%S` 2>> $outputDir/$sampleID.StdOp.Err.txt


echo "GATK Base recalibrator runtime start" `date +%d/%m/%Y\ %H:%M:%S` 2>> $outputDir/$sampleID.StdOp.Err.txt
#Recalibrating bases based on known Indels, snp data
gatk BaseRecalibrator --input $outputDir/$sampleID/${sampleID}_sorted_gatkMarkDups.bam --known-sites $MillsIndels --known-sites $Indels1kGenome --known-sites $DBSNP --output $outputDir/$sampleID/${sampleID}_recal_data.table --reference $reference 2>> $outputDir/$sampleID.StdOp.Err.txt
echo "GATK Base recalibrator runtime end" `date +%d/%m/%Y\ %H:%M:%S` 2>> $outputDir/$sampleID.StdOp.Err.txt

echo "GATK apply BQSR runtime start" `date +%d/%m/%Y\ %H:%M:%S` 2>> $outputDir/$sampleID.StdOp.Err.txt
# Applying the recalibrated base info on Bam
gatk ApplyBQSR --input $outputDir/$sampleID/${sampleID}_sorted_gatkMarkDups.bam --reference $reference --bqsr-recal-file $outputDir/$sampleID/${sampleID}_recal_data.table --output $outputDir/$sampleID/${sampleID}_sorted.deduped.bqsr.bam 2>> $outputDir/$sampleID.StdOp.Err.txt
echo "GATK apply BQSR runtime end" `date +%d/%m/%Y\ %H:%M:%S` 2>> $outputDir/$sampleID.StdOp.Err.txt
```
