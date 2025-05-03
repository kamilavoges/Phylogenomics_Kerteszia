#VCF pipeline
mkdir /draft5/kamila/vcf_cruzii_nano
cd /draft5/kamila/vcf_cruzii_nano
mkdir 0_index   1_data   2_fastqc   3_mapping   4_processing  5_freebayes  6_filtering


step=reducing the coverage of the fastq files to 30x
{

# Subsample 10000 read pairs from two large paired FASTQ files (remember to use the same random seed to keep pairing):

  # seqtk sample -s100 read1.fq 10000 > sub1.fq
  # seqtk sample -s100 read2.fq 10000 > sub2.fq


mkdir /draft5/kamila/vcf_cruzii_nano/1_data/raw

cd /draft5/kamila/vcf_cruzii_nano/1_data/raw

#soft-linking all Illumina  reads

ln -s /data2/HN00148433_abr2021/boc104_[12].fastq.gz ./
ln -s /data2/HN00148436_abr2021/boc_140_[12].fastq.gz ./
ln -s /data2/HN00115205/boc56_[12].fastq.gz ./
ln -s /data2/HN00115205/Boc146_[12].fastq.gz ./
ln -s /data2/HN00148433_abr2021/boc160_[12].fastq.gz ./
ln -s /data2/HN00148438_abr2021/boc66_[12].fastq.gz ./
ln -s /data2/HN00148438_abr2021/boc85_[12].fastq.gz ./
ln -s /data2/HN00148438_abr2021/boc86_[12].fastq.gz ./
ln -s /data2/HN00203827_set2023/BOC_3F_[12].fastq.gz ./
ln -s /data2/HN00148438_abr2021/UCAD7_5_[12].fastq.gz ./
ln -s /data2/HN00115205/UCAD13-7m_[12].fastq.gz ./
ln -s /data2/HN00148438_abr2021/UCAD10_11_[12].fastq.gz ./
ln -s /data2/HN00148438_abr2021/UCAD6_2_[12].fastq.gz ./
ln -s /data2/HN00125086_mar2020/Guapi_F_[12].fastq.gz ./
ln -s /data2/HN00176340_jul2022/Guapi2F_[12].fastq.gz ./
ln -s /data2/HN00125086_mar2020/M_lago_8_1_F_[12].fastq.gz ./
ln -s /data2/HN00125086_mar2020/CdoPesq1_1_M_[12].fastq.gz ./
ln -s /data2/HN00125086_mar2020/M_lago_11_1_M_[12].fastq.gz ./
ln -s /data2/HN00148438_abr2021/M_lago2_1_[12].fastq.gz ./
ln -s /data2/HN00148438_abr2021/M_lago12_1_[12].fastq.gz ./
ln -s /data2/HN00148438_abr2021/Pousada_1_[12].fastq.gz ./
ln -s /data2/HN00125086_mar2020/ResSLBf3_1_M_[12].fastq.gz ./
ln -s /data2/HN00203827_set2023/CJ4F_[12].fastq.gz ./
ln -s /data2/HN00203827_set2023/CJ6M_[12].fastq.gz ./
ln -s /data2/HN00203827_set2023/CJ7M_[12].fastq.gz ./
ln -s /data2/HN00203827_set2023/CJ9M_[12].fastq.gz ./
ln -s /data5/HN00223709_jul2024/Pcaba1_[12].fastq.gz ./
ln -s /data5/HN00223709_jul2024/Pcaba2_[12].fastq.gz ./ 
ln -s /data5/HN00223709_jul2024/Pcaba3_[12].fastq.gz ./
ln -s /data5/HN00223709_jul2024/Pcaba4_[12].fastq.gz ./
ln -s /data5/HN00223709_jul2024/Pcaba5_[12].fastq.gz ./

# 
cd /draft5/kamila/vcf_cruzii_nano/1_data/

#read are 150bp paired end , so each pair has 300bp. For 30x raw coverage we need 170000000 * 30 / 300 =  17,000,000 pairs
cd /draft5/kamila/vcf_cruzii_nano/1_data/raw
ls *.gz | sed 's/\.fastq\.gz//' > raw_file.list


wc -l raw_file.list # Ok


#  nohup bash subsample_30x.sh &
#  subsample_30x.sh
cd /draft5/kamila/vcf_cruzii_nano/1_data
while read -r file
do
    echo "starting process "  $file
    seqtk sample -s100 raw/$file.fastq.gz 17000000 | pigz -c > ${file}_30x.fastq.gz
done < /draft5/kamila/vcf_cruzii_nano/1_data/raw/raw_file.list


nohup bash subsample_30x.sh > nohup_subsample_30x.txt & 



step= 0_index (index reference genome )
{
cd /draft5/kamila/vcf_cruzii_nano/0_index
ln -s  /draft5/luisa/vcf/GCA_943734635.1_idAnoCruzAS_RS32_06_genomic.fna
bwa index -p ref_cruzii   GCA_943734635.1_idAnoCruzAS_RS32_06_genomic.fna
samtools faidx GCA_943734635.1_idAnoCruzAS_RS32_06_genomic.fna

}



step=3_mapping
{
cd /draft5/kamila/vcf_cruzii_nano/3_mapping
sed 's/_[12]$//' /draft5/kamila/vcf_cruzii_nano/1_data/raw/raw_file.list | uniq >  raw_file.list2

wc -l raw_file.list2 # ok



step=3_mapping

cd /draft5/kamila/vcf_cruzii_nano/3_mapping
while read -r file_prefix
do
    date
    fq1=/draft5/kamila/vcf_cruzii_nano/1_data/$file_prefix"_1_30x.fastq.gz" ; fq2=/draft5/kamila/vcf_cruzii_nano/1_data/$file_prefix"_2_30x.fastq.gz"
    ls -alth $fq2
    echo "start processing fastq files " $file_prefix
    bwa mem -M -t 100  /draft5/kamila/vcf_cruzii_nano/1_data/ref_cruzii   $fq1  $fq2  2>/dev/null  | samtools view -buS - >   ${file_prefix}_30x.bam 
    echo "finished processing fastq files " $file_prefix
done < /draft5/kamila/vcf_cruzii_nano/3_mapping/raw_file.list2


nohup bash bwa_mem_cruzii_30x.sh  > nohup_bwa_mem_so_cruzii_nano_30x.txt &
#in /draft5/kamila/vcf_cruzii_nano/4_processing  



cp -p /draft6/vcf_pipeline/4_processing/bwa_mem_cruzii_30x.sh  /draft5/kamila/vcf_cruzii_nano/4_processing

step= 4_processing  (picard) 
{


# creating names.txt or Anopheles data:
cd /draft5/kamila/vcf_cruzii_nano/3_mapping
ls *.bam | awk '{file_id = gensub(/\.bam/,"",1,$0); print file_id , file_id }' > /draft5/kamila/vcf_cruzii_nano/4_processing/names.txt
# boc104_30x boc104_30x
# boc_140_30x boc_140_30x
# Boc146_30x Boc146_30x

cd /draft5/kamila/vcf_cruzii_nano/4_processing


nohup bash picard2.sh > nohup_picard2.txt & 



cd /draft5/kamila/vcf_cruzii_nano/5_freebayes
cp -p /draft5/kamila/vcf_cruzii_nano/0_index/GCA* /draft5/kamila/vcf_cruzii_nano/5_freebayes


python2.7 /home6/tools/Solexa/scripts/fasta_generate_regions.py GCA_943734635.1_idAnoCruzAS_RS32_06_genomic.fna.fai 100000 > regions_list.txt
wc -l regions_list.txt 


ln -s /draft5/kamila/vcf_cruzii_nano/4_processing/*sorted-md-rg.bam ./
ln -s /draft5/kamila/vcf_cruzii_nano/4_processing/*sorted-md-rg.bam.bai ./

conda activate freebayes_condaenv

nohup bash run_freebayes.sh 30 > nohup_run_freebayes.txt & 

###rodar como:
#nohup bash indexar_cruzii.sh > nohup_indexar_cruzii_nano.txt &


nohup /home6/tools/bcftools/bcftools merge BOC_3F_raw.sorted.vcf.gz Boc146_raw.sorted.vcf.gz CJ4F_raw.sorted.vcf.gz CJ6M_raw.sorted.vcf.gz CJ7M_raw.sorted.vcf.gz CJ9M_raw.sorted.vcf.gz CdoPesq1_raw.sorted.vcf.gz Guapi2F_raw.sorted.vcf.gz Guapi_F_raw.sorted.vcf.gz M_lago12_1_raw.sorted.vcf.gz M_lago2_1_raw.sorted.vcf.gz M_lago_11_1_M_raw.sorted.vcf.gz M_lago_8_1_F_raw.sorted.vcf.gz Pcaba1_raw.sorted.vcf.gz Pcaba2_raw.sorted.vcf.gz Pcaba3_raw.sorted.vcf.gz Pcaba4_raw.sorted.vcf.gz Pcaba5_raw.sorted.vcf.gz Pousada_1_raw.sorted.vcf.gz ResSLBf3_1_M_raw.sorted.vcf.gz UCAD10_11_raw.sorted.vcf.gz UCAD13-7_raw.sorted.vcf.gz UCAD6_2_raw.sorted.vcf.gz UCAD7_5_raw.sorted.vcf.gz boc104_raw.sorted.vcf.gz boc140_raw.sorted.vcf.gz boc160_raw.sorted.vcf.gz boc56_raw.sorted.vcf.gz boc66_raw.sorted.vcf.gz boc85_raw.sorted.vcf.gz boc86_raw.sorted.vcf.gz -o raw_cruzii_spp_sorted_merged_files.vcf.gz &
 

gunzip raw_cruzii_spp_sorted_merged_files.vcf.gz

#####extracting chrx
nohup vcftools --vcf raw_cruzii_spp_sorted_merged_files.vcf --chr OX030881.1 --out chrx_cruzii_merged --recode & 

#####extracting chr2
nohup vcftools --vcf raw_cruzii_spp_sorted_merged_files.vcf  --chr OX030879.1 --out chr2_cruzii_merged --recode & 

#####extracting chr3
nohup vcftools --vcf raw_cruzii_spp_sorted_merged_files.vcf  --chr OX030880.1 --out chr3_cruzii_merged --recode & 



# *I made separate vcf for each sample ok
# *I informed frebays with cnv.bed who was M and F ok
# *I joined everything with merge ok
# *I will separate by chromosome ok
# *do the fix jar by chromosome



find ./ -name "*.bam" > bams.list

#make a fixjar for each chromosome

#fixjar chrx
/home6/tools/bcftools/bcftools  annotate -x 'FORMAT/GL' chrx_cruzii_merged.recode.vcf   > chrx_cruzii_merged.recode.noGL.vcf 

nohup java -jar /home/scripts/fixvcfmissinggenotypes.jar --depth 5 --output chrx_cruzii_merged_recode_fixed_vcf.gz -B bams.list < chrx_cruzii_merged.recode.noGL.vcf > nohup_chrx_cruzii_merged_recode_fixed.out & 


##fixjar chr 2
/home6/tools/bcftools/bcftools  annotate -x 'FORMAT/GL' chr2_cruzii_merged.recode.vcf > chr2_cruzii_merged_recode.vcf.noGL.vcf 

nohup java -jar /home/scripts/fixvcfmissinggenotypes.jar --depth 5 --output chr2_cruzii_merged_recode_fixed.vcf  -B bams.list < chr2_cruzii_merged_recode.vcf.noGL.vcf   > nohup_chr2_cruzii_merged_recode_fixed.out & 



##fixjar chr 3
/home6/tools/bcftools/bcftools  annotate -x 'FORMAT/GL' chr3_cruzii_merged.recode.vcf > chr3_cruzii_merged_recode_noGL.vcf 

nohup java -jar /home/scripts/fixvcfmissinggenotypes.jar --depth 5 --output chr3_cruzii_merged_recode_fixed.vcf.gz   -B bams.list < chr3_cruzii_merged_recode_noGL.vcf   > nohup_chr3_cruzii_merged_recode_fixed.out & 



######################################pipeline estimate Fst#############################################

nohup sh separar_genes_chr2_cruzii08jan.sh & #(I put the .sh separately on github)

nohup sh separar_genes_chr3_cruzii08jan.sh &

nohup sh separar_genes_chrx_cruzii08jan.sh & 

mkdir /draft5/kamila/vcf_cruzii_nano/5_freebayes/pipe_fst_cruzii08jan/chr2

mkdir /draft5/kamila/vcf_cruzii_nano/5_freebayes/pipe_fst_cruzii08jan/chr3

mkdir /draft5/kamila/vcf_cruzii_nano/5_freebayes/pipe_fst_cruzii08jan/chrx

gunzip chr2_cruzii_merged_recode_fixed.vcf.gz #unzip all vcf


# Inside the folder there are several .vcfs ending in “.recode.vcf” (e.g.: 128279028.recode.vcf).
# In this same folder, make a list with the names of all the genes you want to analyze.
# To make this list, use the command below (do this inside the folder):

ls -1 *.vcf | cut -c 1-9 > gene_list.txt

#This command will take the first 9 characters of each .vcf and place them in this list, as each gene ID has only 9 characters.

#fst
#below is an example of how to calculate the Fst for chromosome 2, in the comparison cruzii_Bocaina1 vs cruzii_Guapimirim
##############I did this for all comparisons

#cruzii_Bocaina1 vs cruzii_Guapimirim

mkdir /draft5/kamila/vcf_cruzii_nano/5_freebayes/pipe_fst_cruzii08jan/chr2/output_fst_boc1_Gua


./chr2_boc1Gua_script_fst_all_together.sh ./chr2/ ./chr2/gene_list.txt ./chr2/output_fst_boc1_Gua

# ###########################chr2_boc1Gua_script_fst_all_together.sh

# #!/bin/bash

# # Verifique se todos os argumentos necessários foram fornecidos
# if [ "$#" -ne 3 ]; then
  # echo "Uso: $0 <gene_ids_directory> <gene_ids_file> <output_directory>"
  # exit 1
# fi

# # Atribuir argumentos a variáveis
# input_directory="/draft5/kamila/vcf_cruzii_nano/5_freebayes/pipe_fst_cruzii08jan/chr2/"
# gene_ids_directory="/draft5/kamila/vcf_cruzii_nano/5_freebayes/pipe_fst_cruzii08jan/chr2/"
# output_directory="/draft5/kamila/vcf_cruzii_nano/5_freebayes/pipe_fst_cruzii08jan/chr2/output_fst_boc1_Gua"

# # Obter o caminho completo do arquivo de gene_ids
# gene_ids=$gene_ids_directory/gene_list.txt

# # Executar o comando 'cat' e colocar toda a saída em uma nova variável, chamada gene_ids
# gene_ids=$(cat $gene_ids)

# # Criar o diretório de saída 
# mkdir -p $output_directory

# # Este loop percorre as 'palavras' na variável $gene_ids
# for gene_id in $gene_ids
# do
  # vcf_file=$input_directory/$gene_id.recode.vcf
  # Rscript chr2_boc1Gua_script_fst_all_together.R $vcf_file $output_directory/$gene_id
# done



# ##################################chr2_boc1Gua_script_fst_all_together.R:


# library("gdsfmt")
# library("SNPRelate")

# # take the commandline arguments
# args <- commandArgs(trailingOnly = TRUE)

# # die if not enough args
# if (length(args) < 2) {
  # stop("not enough command line arguments")
# }

# # first argument = INPUT FILENAME (full path from current directory)
# vcf.fn <- args[1]

# # second argument = OUTPUT FILE PREFIX 
# # all output filenames will start with this prefix
# output.prefix <- args[2]

# # Loading the VCF file
# # create a temporary filename for the gds intermediate file
# # this uses the Linux process ID of the running script itself so that
# # we can run it in parallel without each process fighting over the same file
# gds.fn <- paste("temp", Sys.getpid(), "gds", sep=".")

# # Reformating to GDS : two method options "biallelic-only" or "copy.num.of.ref"
# snpgdsVCF2GDS(vcf.fn, gds.fn, method="biallelic.only")
# snpgdsSummary(gds.fn)

# # Data Analysis
# # load gds file into memory
# genofile <- snpgdsOpen(gds.fn)
# print("read genofile")

# # clean up temporary file
# file.remove(gds.fn)
# ##########################

# # Fst calculation:

# sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

# pop_code <- scan("filename_cruzii", what=character())

# flag <- pop_code %in% c("cruzii_Bocaina1","cruzii_Guapimirim")

# samp.sel <- sample.id[flag]

# pop.sel <- pop_code[flag]

# fst <- snpgdsFst(genofile, autosome.only=F, sample.id=samp.sel, population=as.factor(pop.sel),method="W&C84")

# #########################
# # OUTPUT FILE:
# write.table(fst$Fst, paste(output.prefix, "fst_boc1_gua_chr2", sep="."), sep="\t",col.names = FALSE,row.names = FALSE)
# #########################