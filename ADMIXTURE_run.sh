ADMIXTURE run (example for chrom 2)

#  plink2 filtering and conversion to bed format :
  plink2 --vcf chr2_cruzii_merged_recode_fixed2.vcf --make-bed --out
./chr2_cruzii_merged_recode_fixed2_geno10  --vcf-half-call m
--max-alleles 2 --maf 0.01 --geno 0.10 --var-min-qual 30 --threads 100

for K in 1 2 3 4 5 6 7 8 9
do
   admixture -j200   --cv chr2_cruzii_merged_recode_fixed2_geno10.bed
$K > log_K${K}_chr2_fixed2_geno10.out
   echo "finished K=" $K
done


#########

for K in 1 2 3 4 5 6 7 8 9
do
  awk '{print $2}' ../chrx_cruzii_merged_recode_fixed2_geno10.fam > simple.fam
  paste simple.fam   ../chrx_cruzii_merged_recode_fixed2_geno10.$K.Q
> temp1.txt
  awk 'BEGIN{printf "%s ","sample"}; (NR==1){for(i=1; i<NF;
i++){printf "V%s ",i}; print ""}; {print $0}' temp1.txt > temp2.txt
  casador_flexivel.awk  add_file_output_fields=no_index
add_missing_char="missing"
add_file=/draft6/bernardo/An_cruzii/admixture_oct2025/admixture_cruzii_samples.txt
 temp2.txt >  chrx_cruzii_fixed2_geno10.$K.out
  plot_ADMIXTURE_v0.py --input chrx_cruzii_fixed2_geno10.$K.out
--group species,population1 --group_sep --output
chrx_cruzii_fixed2_geno10.$K  --graph_style 2 --bar_label ML_ID_simple
--height 2.5 --Xlabel_size 12
done