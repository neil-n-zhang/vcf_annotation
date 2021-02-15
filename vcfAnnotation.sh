#!/usr/bin/env bash 

inputfile=(${1//./ })

java -jar ~/bioinfo_tools/snpEff/snpEff.jar GRCh37.p13.RefSeq ${inputfile[0]}.vcf > ${inputfile[0]}_snpEff.vcf

chmod +x vcfAnnotation.R

./vcfAnnotation.R ${inputfile[0]}_snpEff.vcf
