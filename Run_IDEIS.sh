# Runs all IDEIS analysis
# Must have run all preprocessing scripts beforehand!

# GSE187515
# run this together with subsampling
python GSE187515_subset.py

# GSE212998
python IDEIS/IDEIS_main.py Homo_Sapiens Benchmarks/subset_bams/GSE212998_GFP-pos-w1_possorted_genome_bam.bam GSE212998/IDEIS/GFP-pos-w1 -g Homo-Sapiens --whitelist GSE212998/Whitelists/Whitelist_GFP-pos-w1.csv
python IDEIS/IDEIS_main.py Homo_Sapiens Benchmarks/subset_bams/GSE212998_GFP-pos-w2_possorted_genome_bam.bam GSE212998/IDEIS/GFP-pos-w2 -g Homo-Sapiens --whitelist GSE212998/Whitelists/Whitelist_GFP-pos-w2.csv
python IDEIS/IDEIS_main.py Homo_Sapiens Benchmarks/subset_bams/GSE212998_GFP-pos-w1-expand_possorted_genome_bam.bam GSE212998/IDEIS/GFP-pos-w1-expand -g Homo-Sapiens --whitelist GSE212998/Whitelists/Whitelist_GFP_pos-w1-expand.csv
python IDEIS/IDEIS_main.py Homo_Sapiens Benchmarks/subset_bams/GSE212998_GFP-pos-w2-expand_possorted_genome_bam.bam GSE212998/IDEIS/GFP-pos-w2-expand -g Homo-Sapiens --whitelist GSE212998/Whitelists/Whitelist_GFP_pos-w2-expand.csv
python IDEIS/IDEIS_main.py Homo_Sapiens Benchmarks/subset_bams/GSE212998_GFP-neg-expand_possorted_genome_bam.bam GSE212998/IDEIS/GFP-neg-expand -g Homo-Sapiens --whitelist GSE212998/Whitelists/Whitelist_GFP_neg_inf.csv

mkdir GSE212998/Ptprc
cp GSE212998/IDEIS/GFP-pos-w1/Results/Iso_Counts/raw_matrix.rds GSE212998/Ptprc/Ptprc_GFP-pos-w1.rds
cp GSE212998/IDEIS/GFP-pos-w2/Results/Iso_Counts/raw_matrix.rds GSE212998/Ptprc/Ptprc_GFP-pos-w2.rds
cp GSE212998/IDEIS/GFP-pos-w1-expand/Results/Iso_Counts/raw_matrix.rds GSE212998/Ptprc/Ptprc_GFP-pos-w1-expand.rds
cp GSE212998/IDEIS/GFP-pos-w2-expand/Results/Iso_Counts/raw_matrix.rds GSE212998/Ptprc/Ptprc_GFP-pos-w2-expand.rds
cp GSE212998/IDEIS/GFP-neg-expand/Results/Iso_Counts/raw_matrix.rds GSE212998/Ptprc/Ptprc_GFP-neg-expand.rds

# GSE155006
python IDEIS/IDEIS_main.py Mus_Musculus GSE155006/Bam/Liver_Aged.bam GSE155006/IDEIS/Liver_Aged --whitelist  GSE155006/Whitelists/Whitelist_Liver_Aged.csv
python IDEIS/IDEIS_main.py Mus_Musculus GSE155006/Bam/Liver_Young.bam GSE155006/IDEIS/Liver_Young --whitelist  GSE155006/Whitelists/Whitelist_Liver_Young.csv
python IDEIS/IDEIS_main.py Mus_Musculus GSE155006/Bam/Lung_Aged.bam GSE155006/IDEIS/Lungs_Aged --whitelist  GSE155006/Whitelists/Whitelist_Lungs_Aged.csv
python IDEIS/IDEIS_main.py Mus_Musculus GSE155006/Bam/Lung_Young.bam GSE155006/IDEIS/Lungs_Young --whitelist  GSE155006/Whitelists/Whitelist_Lungs_Young.csv
python IDEIS/IDEIS_main.py Mus_Musculus GSE155006/Bam/Peritoneal_cells_Aged.bam GSE155006/IDEIS/Peritoneal_Cells_Aged --whitelist  GSE155006/Whitelists/Whitelist_Peritoneal_Cells_Aged.csv
python IDEIS/IDEIS_main.py Mus_Musculus GSE155006/Bam/Peritoneal_cells_Young.bam GSE155006/IDEIS/Peritoneal_Cells_Young --whitelist  GSE155006/Whitelists/Whitelist_Peritoneal_Cells_Young.csv
python IDEIS/IDEIS_main.py Mus_Musculus GSE155006/Bam/Spleen_Aged.bam GSE155006/IDEIS/Spleen_Aged --whitelist GSE155006/Whitelists/Whitelist_Spleen_Aged.csv
python IDEIS/IDEIS_main.py Mus_Musculus GSE155006/Bam/Spleen_Young.bam GSE155006/IDEIS/Spleen_Young --whitelist GSE155006/Whitelists/Whitelist_Spleen_Young.csv

mkdir GSE155006/Ptprc
cp GSE155006/IDEIS/Liver_Aged/Results/Iso_Counts/raw_matrix.rds GSE155006/Ptprc/Ptprc_Liver_Aged.rds
cp GSE155006/IDEIS/Liver_Young/Results/Iso_Counts/raw_matrix.rds GSE155006/Ptprc/Ptprc_Liver_Young.rds
cp GSE155006/IDEIS/Lungs_Aged/Results/Iso_Counts/raw_matrix.rds GSE155006/Ptprc/Ptprc_Lungs_Aged.rds
cp GSE155006/IDEIS/Lungs_Young/Results/Iso_Counts/raw_matrix.rds GSE155006/Ptprc/Ptprc_Lungs_Young.rds
cp GSE155006/IDEIS/Peritoneal_Cells_Aged/Results/Iso_Counts/raw_matrix.rds GSE155006/Ptprc/Ptprc_Peritoneal_Cells_Aged.rds
cp GSE155006/IDEIS/Peritoneal_Cells_Young/Results/Iso_Counts/raw_matrix.rds GSE155006/Ptprc/Ptprc_Peritoneal_Cells_Young.rds
cp GSE155006/IDEIS/Spleen_Aged/Results/Iso_Counts/raw_matrix.rds GSE155006/Ptprc/Ptprc_Spleen_Aged.rds
cp GSE155006/IDEIS/Spleen_Young/Results/Iso_Counts/raw_matrix.rds GSE155006/Ptprc/Ptprc_Spleen_Young.rds

# PRJEB40376
# here we need to divide BAM by length (under and over 100 bp for R2) because of division

mkdir PRJEB40376/Bams_split
samtools view -bS -e 'length(seq)<=100' -@ 32 PRJEB40376/Pool_1/outs/possorted_genome_bam.bam > PRJEB40376/Bams_split/Pool_1_possorted_genome_bam_100minus.bam
samtools view -bS -e 'length(seq)>100' -@ 32 PRJEB40376/Pool_1/outs/possorted_genome_bam.bam > PRJEB40376/Bams_split/Pool_1_possorted_genome_bam_100plus.bam
samtools view -bS -e 'length(seq)<=100' -@ 32 PRJEB40376/Pool_2/outs/possorted_genome_bam.bam > PRJEB40376/Bams_split/Pool_2_possorted_genome_bam_100minus.bam
samtools view -bS -e 'length(seq)>100' -@ 32 PRJEB40376/Pool_2/outs/possorted_genome_bam.bam > PRJEB40376/Bams_split/Pool_2_possorted_genome_bam_100plus.bam
samtools view -bS -e 'length(seq)<=100' -@ 32 PRJEB40376/Pool_3/outs/possorted_genome_bam.bam > PRJEB40376/Bams_split/Pool_3_possorted_genome_bam_100minus.bam
samtools view -bS -e 'length(seq)>100' -@ 32 PRJEB40376/Pool_3/outs/possorted_genome_bam.bam > PRJEB40376/Bams_split/Pool_3_possorted_genome_bam_100plus.bam
samtools view -bS -e 'length(seq)<=100' -@ 32 PRJEB40376/Pool_4/outs/possorted_genome_bam.bam > PRJEB40376/Bams_split/Pool_4_possorted_genome_bam_100minus.bam
samtools view -bS -e 'length(seq)>100' -@ 32 PRJEB40376/Pool_4/outs/possorted_genome_bam.bam > PRJEB40376/Bams_split/Pool_4_possorted_genome_bam_100plus.bam
samtools view -bS -e 'length(seq)<=100' -@ 32 PRJEB40376/Pool_5/outs/possorted_genome_bam.bam > PRJEB40376/Bams_split/Pool_5_possorted_genome_bam_100minus.bam
samtools view -bS -e 'length(seq)>100' -@ 32 PRJEB40376/Pool_5/outs/possorted_genome_bam.bam > PRJEB40376/Bams_split/Pool_5_possorted_genome_bam_100plus.bam
samtools view -bS -e 'length(seq)<=100' -@ 32 PRJEB40376/Pool_6/outs/possorted_genome_bam.bam > PRJEB40376/Bams_split/Pool_6_possorted_genome_bam_100minus.bam
samtools view -bS -e 'length(seq)>100' -@ 32 PRJEB40376/Pool_6/outs/possorted_genome_bam.bam > PRJEB40376/Bams_split/Pool_6_possorted_genome_bam_100plus.bam
samtools view -bS -e 'length(seq)<=100' -@ 32 PRJEB40376/Pool_7/outs/possorted_genome_bam.bam > PRJEB40376/Bams_split/Pool_7_possorted_genome_bam_100minus.bam
samtools view -bS -e 'length(seq)>100' -@ 32 PRJEB40376/Pool_7/outs/possorted_genome_bam.bam > PRJEB40376/Bams_split/Pool_7_possorted_genome_bam_100plus.bam
samtools view -bS -e 'length(seq)<=100' -@ 32 PRJEB40376/Pool_8/outs/possorted_genome_bam.bam > PRJEB40376/Bams_split/Pool_8_possorted_genome_bam_100minus.bam
samtools view -bS -e 'length(seq)>100' -@ 32 PRJEB40376/Pool_8/outs/possorted_genome_bam.bam > PRJEB40376/Bams_split/Pool_8_possorted_genome_bam_100plus.bam
samtools view -bS -e 'length(seq)<=100' -@ 32 PRJEB40376/Pool_9/outs/possorted_genome_bam.bam > PRJEB40376/Bams_split/Pool_9_possorted_genome_bam_100minus.bam
samtools view -bS -e 'length(seq)>100' -@ 32 PRJEB40376/Pool_9/outs/possorted_genome_bam.bam > PRJEB40376/Bams_split/Pool_9_possorted_genome_bam_100plus.bam
samtools view -bS -e 'length(seq)<=100' -@ 32 PRJEB40376/Pool_10/outs/possorted_genome_bam.bam > PRJEB40376/Bams_split/Pool_10_possorted_genome_bam_100minus.bam
samtools view -bS -e 'length(seq)>100' -@ 32 PRJEB40376/Pool_10/outs/possorted_genome_bam.bam > PRJEB40376/Bams_split/Pool_10_possorted_genome_bam_100plus.bam

samtools index PRJEB40376/Bams_split/Pool_1_possorted_genome_bam_100minus.bam -@ 32
samtools index PRJEB40376/Bams_split/Pool_1_possorted_genome_bam_100plus.bam -@ 32
samtools index PRJEB40376/Bams_split/Pool_2_possorted_genome_bam_100minus.bam -@ 32
samtools index PRJEB40376/Bams_split/Pool_2_possorted_genome_bam_100plus.bam -@ 32
samtools index PRJEB40376/Bams_split/Pool_3_possorted_genome_bam_100minus.bam -@ 32
samtools index PRJEB40376/Bams_split/Pool_3_possorted_genome_bam_100plus.bam -@ 32
samtools index PRJEB40376/Bams_split/Pool_4_possorted_genome_bam_100minus.bam -@ 32
samtools index PRJEB40376/Bams_split/Pool_4_possorted_genome_bam_100plus.bam -@ 32
samtools index PRJEB40376/Bams_split/Pool_5_possorted_genome_bam_100minus.bam -@ 32
samtools index PRJEB40376/Bams_split/Pool_5_possorted_genome_bam_100plus.bam -@ 32
samtools index PRJEB40376/Bams_split/Pool_6_possorted_genome_bam_100minus.bam -@ 32
samtools index PRJEB40376/Bams_split/Pool_6_possorted_genome_bam_100plus.bam -@ 32
samtools index PRJEB40376/Bams_split/Pool_7_possorted_genome_bam_100minus.bam -@ 32
samtools index PRJEB40376/Bams_split/Pool_7_possorted_genome_bam_100plus.bam -@ 32
samtools index PRJEB40376/Bams_split/Pool_8_possorted_genome_bam_100minus.bam -@ 32
samtools index PRJEB40376/Bams_split/Pool_8_possorted_genome_bam_100plus.bam -@ 32
samtools index PRJEB40376/Bams_split/Pool_9_possorted_genome_bam_100minus.bam -@ 32
samtools index PRJEB40376/Bams_split/Pool_9_possorted_genome_bam_100plus.bam -@ 32
samtools index PRJEB40376/Bams_split/Pool_10_possorted_genome_bam_100minus.bam -@ 32
samtools index PRJEB40376/Bams_split/Pool_10_possorted_genome_bam_100plus.bam -@ 32

python IDEIS/IDEIS_main.py Homo_Sapiens PRJEB40376/Bams_split/Pool_1_possorted_genome_bam_100minus.bam PRJEB40376/IDEIS/Pool_1_minus -g Homo-Sapiens --sequencing-type 3-prime --whitelist PRJEB40376/Whitelists/Whitelist_Pool_1.csv
python IDEIS/IDEIS_main.py Homo_Sapiens PRJEB40376/Bams_split/Pool_2_possorted_genome_bam_100minus.bam PRJEB40376/IDEIS/Pool_2_minus -g Homo-Sapiens --sequencing-type 3-prime --whitelist PRJEB40376/Whitelists/Whitelist_Pool_2.csv
python IDEIS/IDEIS_main.py Homo_Sapiens PRJEB40376/Bams_split/Pool_3_possorted_genome_bam_100minus.bam PRJEB40376/IDEIS/Pool_3_minus -g Homo-Sapiens --sequencing-type 3-prime --whitelist PRJEB40376/Whitelists/Whitelist_Pool_3.csv
python IDEIS/IDEIS_main.py Homo_Sapiens PRJEB40376/Bams_split/Pool_4_possorted_genome_bam_100minus.bam PRJEB40376/IDEIS/Pool_4_minus -g Homo-Sapiens --sequencing-type 3-prime --whitelist PRJEB40376/Whitelists/Whitelist_Pool_4.csv
python IDEIS/IDEIS_main.py Homo_Sapiens PRJEB40376/Bams_split/Pool_5_possorted_genome_bam_100minus.bam PRJEB40376/IDEIS/Pool_5_minus -g Homo-Sapiens --sequencing-type 3-prime --whitelist PRJEB40376/Whitelists/Whitelist_Pool_5.csv
python IDEIS/IDEIS_main.py Homo_Sapiens PRJEB40376/Bams_split/Pool_6_possorted_genome_bam_100minus.bam PRJEB40376/IDEIS/Pool_6_minus -g Homo-Sapiens --sequencing-type 3-prime --whitelist PRJEB40376/Whitelists/Whitelist_Pool_6.csv
python IDEIS/IDEIS_main.py Homo_Sapiens PRJEB40376/Bams_split/Pool_7_possorted_genome_bam_100minus.bam PRJEB40376/IDEIS/Pool_7_minus -g Homo-Sapiens --sequencing-type 3-prime --whitelist PRJEB40376/Whitelists/Whitelist_Pool_7.csv
python IDEIS/IDEIS_main.py Homo_Sapiens PRJEB40376/Bams_split/Pool_8_possorted_genome_bam_100minus.bam PRJEB40376/IDEIS/Pool_8_minus -g Homo-Sapiens --sequencing-type 3-prime --whitelist PRJEB40376/Whitelists/Whitelist_Pool_8.csv
python IDEIS/IDEIS_main.py Homo_Sapiens PRJEB40376/Bams_split/Pool_9_possorted_genome_bam_100minus.bam PRJEB40376/IDEIS/Pool_9_minus -g Homo-Sapiens --sequencing-type 3-prime --whitelist PRJEB40376/Whitelists/Whitelist_Pool_9.csv
python IDEIS/IDEIS_main.py Homo_Sapiens PRJEB40376/Bams_split/Pool_10_possorted_genome_bam_100minus.bam PRJEB40376/IDEIS/Pool_10_minus -g Homo-Sapiens --sequencing-type 3-prime --whitelist PRJEB40376/Whitelists/Whitelist_Pool_10.csv

python IDEIS/IDEIS_main.py Homo_Sapiens PRJEB40376/Bams_split/Pool_1_possorted_genome_bam_100plus.bam PRJEB40376/IDEIS/Pool_1_plus -g Homo-Sapiens --sequencing-type 3-prime --whitelist PRJEB40376/Whitelists/Whitelist_Pool_1.csv
python IDEIS/IDEIS_main.py Homo_Sapiens PRJEB40376/Bams_split/Pool_2_possorted_genome_bam_100plus.bam PRJEB40376/IDEIS/Pool_2_plus -g Homo-Sapiens --sequencing-type 3-prime --whitelist PRJEB40376/Whitelists/Whitelist_Pool_2.csv
python IDEIS/IDEIS_main.py Homo_Sapiens PRJEB40376/Bams_split/Pool_3_possorted_genome_bam_100plus.bam PRJEB40376/IDEIS/Pool_3_plus -g Homo-Sapiens --sequencing-type 3-prime --whitelist PRJEB40376/Whitelists/Whitelist_Pool_3.csv
python IDEIS/IDEIS_main.py Homo_Sapiens PRJEB40376/Bams_split/Pool_4_possorted_genome_bam_100plus.bam PRJEB40376/IDEIS/Pool_4_plus -g Homo-Sapiens --sequencing-type 3-prime --whitelist PRJEB40376/Whitelists/Whitelist_Pool_4.csv
python IDEIS/IDEIS_main.py Homo_Sapiens PRJEB40376/Bams_split/Pool_5_possorted_genome_bam_100plus.bam PRJEB40376/IDEIS/Pool_5_plus -g Homo-Sapiens --sequencing-type 3-prime --whitelist PRJEB40376/Whitelists/Whitelist_Pool_5.csv
python IDEIS/IDEIS_main.py Homo_Sapiens PRJEB40376/Bams_split/Pool_6_possorted_genome_bam_100plus.bam PRJEB40376/IDEIS/Pool_6_plus -g Homo-Sapiens --sequencing-type 3-prime --whitelist PRJEB40376/Whitelists/Whitelist_Pool_6.csv
python IDEIS/IDEIS_main.py Homo_Sapiens PRJEB40376/Bams_split/Pool_7_possorted_genome_bam_100plus.bam PRJEB40376/IDEIS/Pool_7_plus -g Homo-Sapiens --sequencing-type 3-prime --whitelist PRJEB40376/Whitelists/Whitelist_Pool_7.csv
python IDEIS/IDEIS_main.py Homo_Sapiens PRJEB40376/Bams_split/Pool_8_possorted_genome_bam_100plus.bam PRJEB40376/IDEIS/Pool_8_plus -g Homo-Sapiens --sequencing-type 3-prime --whitelist PRJEB40376/Whitelists/Whitelist_Pool_8.csv
python IDEIS/IDEIS_main.py Homo_Sapiens PRJEB40376/Bams_split/Pool_9_possorted_genome_bam_100plus.bam PRJEB40376/IDEIS/Pool_9_plus -g Homo-Sapiens --sequencing-type 3-prime --whitelist PRJEB40376/Whitelists/Whitelist_Pool_9.csv
python IDEIS/IDEIS_main.py Homo_Sapiens PRJEB40376/Bams_split/Pool_10_possorted_genome_bam_100plus.bam PRJEB40376/IDEIS/Pool_10_plus -g Homo-Sapiens --sequencing-type 3-prime --whitelist PRJEB40376/Whitelists/Whitelist_Pool_10.csv

mkdir PRJEB40376/Ptprc
cp PRJEB40376/IDEIS/Pool_1_minus/Results/Iso_Counts/raw_matrix.rds PRJEB40376/Ptprc/Ptprc_Pool_1_100minus_raw_matrix.rds
cp PRJEB40376/IDEIS/Pool_2_minus/Results/Iso_Counts/raw_matrix.rds PRJEB40376/Ptprc/Ptprc_Pool_2_100minus_raw_matrix.rds
cp PRJEB40376/IDEIS/Pool_3_minus/Results/Iso_Counts/raw_matrix.rds PRJEB40376/Ptprc/Ptprc_Pool_3_100minus_raw_matrix.rds
cp PRJEB40376/IDEIS/Pool_4_minus/Results/Iso_Counts/raw_matrix.rds PRJEB40376/Ptprc/Ptprc_Pool_4_100minus_raw_matrix.rds
cp PRJEB40376/IDEIS/Pool_5_minus/Results/Iso_Counts/raw_matrix.rds PRJEB40376/Ptprc/Ptprc_Pool_5_100minus_raw_matrix.rds
cp PRJEB40376/IDEIS/Pool_6_minus/Results/Iso_Counts/raw_matrix.rds PRJEB40376/Ptprc/Ptprc_Pool_6_100minus_raw_matrix.rds
cp PRJEB40376/IDEIS/Pool_7_minus/Results/Iso_Counts/raw_matrix.rds PRJEB40376/Ptprc/Ptprc_Pool_7_100minus_raw_matrix.rds
cp PRJEB40376/IDEIS/Pool_8_minus/Results/Iso_Counts/raw_matrix.rds PRJEB40376/Ptprc/Ptprc_Pool_8_100minus_raw_matrix.rds
cp PRJEB40376/IDEIS/Pool_9_minus/Results/Iso_Counts/raw_matrix.rds PRJEB40376/Ptprc/Ptprc_Pool_9_100minus_raw_matrix.rds
cp PRJEB40376/IDEIS/Pool_10_minus/Results/Iso_Counts/raw_matrix.rds PRJEB40376/Ptprc/Ptprc_Pool_10_100minus_raw_matrix.rds

cp PRJEB40376/IDEIS/Pool_1_plus/Results/Iso_Counts/raw_matrix.rds PRJEB40376/Ptprc/Ptprc_Pool_1_100plus_raw_matrix.rds
cp PRJEB40376/IDEIS/Pool_2_plus/Results/Iso_Counts/raw_matrix.rds PRJEB40376/Ptprc/Ptprc_Pool_2_100plus_raw_matrix.rds
cp PRJEB40376/IDEIS/Pool_3_plus/Results/Iso_Counts/raw_matrix.rds PRJEB40376/Ptprc/Ptprc_Pool_3_100plus_raw_matrix.rds
cp PRJEB40376/IDEIS/Pool_4_plus/Results/Iso_Counts/raw_matrix.rds PRJEB40376/Ptprc/Ptprc_Pool_4_100plus_raw_matrix.rds
cp PRJEB40376/IDEIS/Pool_5_plus/Results/Iso_Counts/raw_matrix.rds PRJEB40376/Ptprc/Ptprc_Pool_5_100plus_raw_matrix.rds
cp PRJEB40376/IDEIS/Pool_6_plus/Results/Iso_Counts/raw_matrix.rds PRJEB40376/Ptprc/Ptprc_Pool_6_100plus_raw_matrix.rds
cp PRJEB40376/IDEIS/Pool_7_plus/Results/Iso_Counts/raw_matrix.rds PRJEB40376/Ptprc/Ptprc_Pool_7_100plus_raw_matrix.rds
cp PRJEB40376/IDEIS/Pool_8_plus/Results/Iso_Counts/raw_matrix.rds PRJEB40376/Ptprc/Ptprc_Pool_8_100plus_raw_matrix.rds
cp PRJEB40376/IDEIS/Pool_9_plus/Results/Iso_Counts/raw_matrix.rds PRJEB40376/Ptprc/Ptprc_Pool_9_100plus_raw_matrix.rds
cp PRJEB40376/IDEIS/Pool_10_plus/Results/Iso_Counts/raw_matrix.rds PRJEB40376/Ptprc/Ptprc_Pool_10_100plus_raw_matrix.rds

