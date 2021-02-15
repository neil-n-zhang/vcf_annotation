# Human VCF file annotation

Input file: vcf file using human_g1k_v37.fasta as the reference genome, the vcf file should contain the following columns as the specific order: CHROM POS ID REF ALT QUAL FILTER INFO

Output file: annotated variants in csv contains the following information:
chr: Chromosome 
pos: Chromosome position of the variant 
ref: Reference sequence
alt: Variant sequence	
dp: Total read depth at the locus	
qual: variant quality score	
var_reads_fq: The percentage of the variant reads	
exac_fq: The variant frequency from the Exac database
type: The type of the variant	
eff: The effect of the variant predicted by snpEff

## 1. Download snpEff for variant effects annotation

wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip

unzip snpEff_latest_core.zip

## 2. Configuration of vcfAnnotation.sh

Open vcfAnnotation.sh with the text editor and set up the snpEff.jar location

## 3. Make vcfAnnotation.sh executable and run vcfAnnotation.sh

chmod +x vcfAnnotation.sh

bash vcfAnnotation.sh Challenge_data.vcf
