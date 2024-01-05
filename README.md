This directory contains all required file to reproduce analyses presented in paper of Michalik et al.  

Link to paper:

## Before starting

Chceck that you have RStudio with R 4.2.1 installed along with required packages (mainly Seurat 4.1.1) and python 3.8.10. The software samtools is also required. IDEIS itself requires salmon.  

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

8.) Perform analyses with RStudio of files

