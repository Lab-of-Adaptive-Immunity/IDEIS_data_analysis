# Performs subsampling of bam file from GSE187515 and

import os, sys , glob

ncores = 32
  
bamtofastq_path = 'IDEIS/bamtofastq-1.4.1'
subsampled_path = 'GSE187515/Subset_bam'
source_exps = ['2004', '2023Rep1', '2023Rep2']
fastq_path = 'GSE187515/Subset_fastqs'
analysis_path = 'GSE187515/Subset_cellranger'

cellranger_path = '/mnt/scratch/beluga/michalik/cellranger3.0.2/cellranger-3.0.2/cellranger-cs/3.0.2/bin/cellranger'
path_to_reference_3 = 'Reference3'
transcriptome = os.path.abspath(os.path.join(path_to_reference_3, 'refdata-cellranger-GRCh38-3.0.0'))

if __name__ == "__main__":
  # converts subsampled bam to fastq
  if not os.path.isdir(fastq_path):
    os.makedirs(fastq_path)
  
  for source_exp in source_exps:
    bam_path_complete = os.path.join(subsampled_path, source_exp, "%s_possorted_genome_bam_10.bam"%(source_exp)) # taking 10% file
    fastq_path_complete = os.path.abspath(os.path.join(fastq_path, source_exp))
    print('%s %s %s'%(bamtofastq_path, bam_path_complete, fastq_path_complete))
    os.system('%s %s %s'%(bamtofastq_path, bam_path_complete, fastq_path_complete))

  if not os.path.isdir(analysis_path):
    os.makedirs(analysis_path)

  current_path = os.getcwd()

  # run cellranger
  for source_exp in source_exps:
    
    fastq_path_complete = os.path.abspath(os.path.join(fastq_path, source_exp))
    print(fastq_path_complete)
    
    os.chdir(analysis_path)
    print('%s count --id=%s_subset --transcriptome=%s --fastqs=%s --sample=bamtofastq --localcores=%d'%(cellranger_path, source_exp, transcriptome, fastq_path_complete, ncores))
    os.system('%s count --id=%s_subset --transcriptome=%s --fastqs=%s --sample=bamtofastq --localcores=%d'%(cellranger_path, source_exp, transcriptome, fastq_path_complete, ncores))

    os.chdir(current_path)
