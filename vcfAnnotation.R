#!/usr/bin/env Rscript

library(httr)
library(jsonlite)
library(dplyr)

#Read Variant Call Format (vcf) file.  
#Input: file path. Output: a list contains the raw file ($vcf), and the information column ($info) of the file
readVcf <-function(vcfpath)
{
  vcf <- read.table(vcfpath, sep = "\t", stringsAsFactors = F)
  info=strsplit(vcf$V8,';')
  #For loss-of-function variants, SNPeff add one or two extra columns, therefore we only keep the firt 42 columns
  info=lapply(info,head,n=42)
  info=as.data.frame(info)
  info=as.data.frame(t(info))
  rownames(info)=NULL
  #DP is the total depth, RO is the number of reads for the reference, AO is the number of reads for variant
  colnames(info)=c("AB","ABP","AC","AF","AN","AO","CIGAR","DP","DPB","DPRA","EPP","EPPR","GTI","LEN","MEANALT","MQM",
                   "MQMR","NS","NUMALT","ODDS","PAIRED","PAIREDR","PAO","PQA","PQR","PRO","QA","QR","RO","RPL","RPP",
                   "RPPR","RPR","RUN","SAF","SAP","SAR","SRF","SRP","SRR","TYPE","EFF")
  #remove the head part (XX=) for each values
  info=data.frame(lapply(info, gsub, pattern=".*=",replace=""), stringsAsFactors = F)
  return(list("vcf" = vcf, "info" = info))
}

#Process the vcf file, extract and calculate useful information  
#Input: vcf list from readVcf. Output: a data.frame contains the variants calling and all relevant features
processVcf<-function(vcf,info)
{
  vcf_anno=data.frame(chr=vcf$V1,pos=vcf$V2,ref=vcf$V4,alt=vcf$V5,qual=vcf$V6,type=info$TYPE,dp=info$DP,ro=info$RO,ao=info$AO, stringsAsFactors = F)
  #SnpEff already put the more deleterious effects first, select the first effect annotation
  eff=strsplit(info$EFF,"\\|")
  vcf_anno$eff=sapply(eff,`[`,2)
  
  #Add effect annotations for sites contain to two different variants
  #output example: CATATATATATA:5_prime_UTR_variant|CATATATATATATATA:5_prime_UTR_variant
  alt_list=strsplit(vcf_anno$alt,',')
  alt_length=sapply(alt_list,length)
  multi_alt_index=alt_length>1
  eff_list=strsplit(info$EFF[multi_alt_index],',')
  multi_alt_n=length(eff_list)
  eff_new=c()
  for (i in 1:multi_alt_n){
    eff=strsplit(eff_list[[i]],"\\|")
    eff=as.data.frame(lapply(eff,head,n=2))
    eff=as.data.frame(t(eff), stringsAsFactors = F)
    rownames(eff)=NULL
    colnames(eff)=c("alt","eff")
    
    alt_list=unique(eff$alt)
    row_index=c()
    
    for (alt in alt_list){
      row_index=c(row_index,grep(alt,eff$alt)[1])
    }
    eff=eff[row_index,]
    eff$summary=paste(eff$alt,eff$eff,sep = ':')
    eff_summary=paste(eff$summary,collapse ='|')
    eff_new=c(eff_new,eff_summary)
  }
  vcf_anno$eff[multi_alt_index]=eff_new
  
  #Calculate the percentage of reads supporting the variant var_reads_fq=n_variantreads/(n_variantreads+n_refreads)
  vcf_anno$var_reads_fq[!multi_alt_index]=as.numeric(vcf_anno$ao[!multi_alt_index])/
    (as.numeric(vcf_anno$ro[!multi_alt_index])+as.numeric(vcf_anno$ao[!multi_alt_index]))
  vcf_anno$var_reads_fq=round(vcf_anno$var_reads_fq, digits=3)
  vcf_anno$var_reads_fq=as.character(vcf_anno$var_reads_fq)
  
  #Calculate the percentage of reads supporting the variant for sites contain two different variants
  var_reads_fq_multi_alt=strsplit(vcf_anno$ao[multi_alt_index],',')
  var_reads_fq_multi_alt=lapply(var_reads_fq_multi_alt, as.numeric)
  dp_multi_alt=sapply(var_reads_fq_multi_alt,sum)+as.numeric(vcf_anno$ro[multi_alt_index])
  var_reads_fq_multi_alt = mapply(FUN = `/`, var_reads_fq_multi_alt, dp_multi_alt, SIMPLIFY = FALSE)
  var_reads_fq_multi_alt=lapply(var_reads_fq_multi_alt,round, digits=3)
  var_reads_fq_multi_alt=as.character(var_reads_fq_multi_alt)
  var_reads_fq_multi_alt=sub("c\\(",'',var_reads_fq_multi_alt)
  var_reads_fq_multi_alt=sub(')','',var_reads_fq_multi_alt)
  vcf_anno$var_reads_fq[multi_alt_index]=var_reads_fq_multi_alt
  return(vcf_anno)
}

#Query allele frequency of variants from ExAC API
#Input: a data.frame from processVcf. Output: a data.frame contains ExAC allele frequency
addExacfq <-function(vcf_anno)
{
  #bulk query using an JSON array
  Exaccall=paste0(vcf_anno$chr,"-",vcf_anno$pos,"-",vcf_anno$ref,"-",vcf_anno$alt)
  Exacquery=POST(url="http://exac.hms.harvard.edu/rest/bulk/variant/variant", body=toJSON(Exaccall), encode = "json")
  jsoncontent <- content(Exacquery)
  Exacfreq=sapply(jsoncontent,`[`,'allele_freq')
  Exacfreq_num=c()
  for (fq in Exacfreq){
    if (!is.null(fq)){
      Exacfreq_num=c(Exacfreq_num,round(fq,3))
    }
    else{
      Exacfreq_num=c(Exacfreq_num,'NA')
    }
  }
  
  #The output doesn't follow the same order as the input, reorder the result
  Exacfreq_num=as_tibble(Exacfreq_num)
  record=strsplit(names(Exacfreq),'-')
  chr=sapply(record,`[`,1)
  pos=sapply(record,`[`,2)
  Exacfreq_num$chr=chr
  Exacfreq_num$pos=as.numeric(pos)
  Exacfreq_num=arrange(Exacfreq_num,chr,pos)
  vcf_anno$exac_fq=Exacfreq_num$value
  
  #Query allele frequency for sites contain two different variants
  alt_list=strsplit(vcf_anno$alt,',')
  alt_length=sapply(alt_list,length)
  multi_alt_index=which(alt_length>1)
  for (i in multi_alt_index){
    alt_list=strsplit(vcf_anno$alt[i],',')[[1]]
    Exaccall=paste0(vcf_anno$chr[i],"-",vcf_anno$pos[i],"-",vcf_anno$ref[i],"-",alt_list)
    Exacquery=POST(url="http://exac.hms.harvard.edu/rest/bulk/variant/variant", body=toJSON(Exaccall), encode = "json")
    jsoncontent <- content(Exacquery)
    Exacfreq=sapply(jsoncontent,`[`,'allele_freq')
    Exacfreq_num=c()
    for (fq in Exacfreq){
      if (!is.null(fq)){
        Exacfreq_num=c(Exacfreq_num,round(fq,3))
      }
      else{
        Exacfreq_num=c(Exacfreq_num,'NA')
      }
    }
    vcf_anno$exac_fq[i]=paste(Exacfreq_num,collapse=',')
  }
  return(vcf_anno)
}

#save the annotated variants into a csv file
saveVcfanno <-function(vcf_anno)
{
  vcf_anno_filter=vcf_anno[,c(1,2,3,4,7,5,11,12,6,10)]
  write.csv(vcf_anno_filter,'Challenge_data_annotated.csv')
}


###Main script
args = commandArgs(trailingOnly = TRUE)
vcf_info=readVcf(args[1L])
vcf_processed=processVcf(vcf_info$vcf,vcf_info$info)
vcf_anno=addExacfq(vcf_processed)
saveVcfanno(vcf_anno)