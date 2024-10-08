This directory contains all required file to reproduce analyses presented in paper of Michalik et al.  

Link to paper:

## Before starting

Chceck that you have RStudio with R 4.2.1 installed along with required packages (mainly Seurat 4.1.1) and python 3.8.10. The software samtools of version at least 1.14 is also required. IDEIS itself requires Salmon Alevin of version .  

## Usage

1.) Clone this directory where you want to then go into it 

```
git clone  
cd  
```

2.) Clone IDEIS into this directory, then download bamtofastq into IDEIS directory and return to root of the project: 


```
git clone https://github.com/Lab-of-Adaptive-Immunity/IDEIS  
cd IDEIS  
wget https://github.com/10XGenomics/bamtofastq/releases/download/v1.4.1/bamtofastq_linux  
mv bamtofastq_linux bamtofastq-1.4.1
cd ..
```

3.) Clone CD45er project there as well:  

```
git clone https://github.com/getzlab/10x-cd45-isoform-quantification  
```

4.) Get subset-bam from 10X and download it to Benchmarks directory:  

```
cd Benchmarks  
wget https://github.com/10XGenomics/subset-bam/releases/download/v1.1.0/subset-bam_linux  
cd ..
```

5.) Run Data_downloader.sh, which will download all required data:  

```
bash Data_downloader.sh
```

6.) Run Path_injector.py. This will inject absolute path of root to library templates for Cell Ranger to run. 

```
python Path_injector.py
```

7.) Run Run_Cellranger.py, which will run Cell Ranger analyses.

```
python Run_Cellranger.py
```

8.) Perform preparations of data with GSE187515_preparation.Rmd, GSE212998_preparation.Rmd, GSE155006_preparation.Rmd and PRJEB40376_preparation.Rmd. This is optimally done in RStudio, using 'Run All' option from 'Run' menu of RStudio. 

9.) Subset created BAM files to barcodes in whitelists (which were created in previous step), then Run benchmarks.

```
cd Benchmarks
python subset_bams_by_whitelists.py
python benchmark_CD45er.py
python benchmark_IDEIS.py
cd ..
```

10.) Subset Pool 1 file from PRJEB40376 to cells in whitelist only for the subsampling.

```
cd PRJEB40376
python subset_bams_by_whitelists_PRJEB40376.py
cd ..
```

11.) Run IDEIS using Run_IDEIS.sh. Warning: In case of PRJEB40376 the BAM files are split by length beforehand, which takes a long time.

```
bash Run_IDEIS.sh
```

12.) Create subsampled data of GSE187515 and PRJE40376 with only certain percentages used, run IDEIS on them, then use next step to count reads that are mapping to used genome for each subsampling percentage.  Warning: This step can take a lot of time.

```
python GSE187515_subset.py
python PRJEB40376_subset.py
python count_mapped_reads.py
```

13.) Perform analyses of data with GSE187515_analysis.Rmd, GSE212998_analysis.Rmd, GSE155006_analysis.Rmd and PRJEB40376_analysis.Rmd. The best way is to run these files in RStudio with 'Run All' option.

14.) Run Benchmarks.Rmd. Again, the best way is to use RStudio and 'Run All' option.

15.) Run supplementary analyses done during Revision. These are in Revisions directory.

```
cd Revisions
```

Once there, run Figures_for_revision_tile_vln_plots.Rmd, Pseudobulk_R_regressions.Rmd and UMAP_isoforms_as_features.Rmd (in this order and ideally in R Studio with 'Run All' option). Return back in main project directory once finished.

```
cd ..
```

