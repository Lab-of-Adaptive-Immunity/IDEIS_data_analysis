import os, sys

if __name__ == '__main__':

  ###########################
  # Parameters              # 
  ###########################
  
  # Fill in configuration parameters
  
  ncores = 32
  path_to_cellranger_3 = '/mnt/scratch/beluga/michalik/cellranger3.0.2/cellranger-3.0.2/cellranger-cs/3.0.2/bin/cellranger' # location of cellranger3.0.2
  path_to_cellranger_5 = '/mnt/scratch/beluga/michalik/cellranger5.0.1/cellranger-5.0.1/cellranger' # location of cellranger5.0.1
  path_to_reference_3 = 'Reference3'  # where reference for Cell Ranger 3.0.2 will be downloaded
  path_to_ref_data_5 = 'Reference5'   # where files for references for Cell Ranger 5.0.1 will be downloaded; here alsoCell Ranger reference will be created
  
  # get absolute paths
  path_to_cellranger_3 = os.path.abspath(path_to_cellranger_3)
  path_to_cellranger_5 = os.path.abspath(path_to_cellranger_5)
  path_to_reference_3 = os.path.abspath(path_to_reference_3)
  path_to_ref_data_5 = os.path.abspath(path_to_ref_data_5)

  if ('' in [path_to_cellranger_3, path_to_cellranger_5, path_to_reference_3, path_to_ref_data_5]):
    sys.exit('One of required paths is undefined. Please fix this and run script again.')
  
  ##############################################
  # Download refererence for Cell Ranger 3.0.2 #
  ##############################################
  
  current_path = os.getcwd()
  if not os.path.exists(path_to_reference_3):
    os.makedirs(path_to_reference_3)
  os.chdir(path_to_reference_3)
  os.system('wget -nc "https://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-GRCh38-3.0.0.tar.gz"')
  if not os.path.exists('refdata-cellranger-GRCh38-3.0.0'):
    os.system('tar -xvzf refdata-cellranger-GRCh38-3.0.0.tar.gz')
  
  os.chdir(current_path)
  
  ##########################################
  # Create reference for Cell Ranger 5.0.1 #
  ##########################################
  
  # download required files from ensembl to path
  if not os.path.exists(path_to_ref_data_5):
    os.makedirs(path_to_ref_data_5)
    
  os.chdir(path_to_ref_data_5)
  os.system('wget -nc https://ftp.ensembl.org/pub/release-102/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz')
  os.system('wget -nc https://ftp.ensembl.org/pub/release-102/gtf/homo_sapiens/Homo_sapiens.GRCh38.102.gtf.gz')
  if(not os.path.exists('Homo_sapiens.GRCh38.dna.primary_assembly.fa')):
    os.system('gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz')
  
  if(not os.path.exists('Homo_sapiens.GRCh38.102.gtf')):  
    os.system('gunzip Homo_sapiens.GRCh38.102.gtf.gz')
  
  # filter gtf
  if not os.path.exists('Homo_sapiens.GRCh38.102.filtered.gtf'):
    os.system('%s mkgtf Homo_sapiens.GRCh38.102.gtf Homo_sapiens.GRCh38.102.filtered.gtf --attribute=gene_biotype:protein_coding --attribute=gene_biotype:lincRNA --attribute=gene_biotype:antisense --attribute=gene_biotype:IG_LV_gene --attribute=gene_biotype:IG_V_gene --attribute=gene_biotype:IG_V_pseudogene --attribute=gene_biotype:IG_D_gene --attribute=gene_biotype:IG_J_gene --attribute=gene_biotype:IG_J_pseudogene --attribute=gene_biotype:IG_C_gene --attribute=gene_biotype:IG_C_pseudogene --attribute=gene_biotype:TR_V_gene --attribute=gene_biotype:TR_V_pseudogene --attribute=gene_biotype:TR_D_gene --attribute=gene_biotype:TR_J_gene --attribute=gene_biotype:TR_J_pseudogene --attribute=gene_biotype:TR_C_gene'%path_to_cellranger_5)
  
  # make reference
  if not os.path.exists('GRCh38_102'): # skip if directory exists; remove if corrupted and retry
    os.system('%s mkref --genome=GRCh38_102 --fasta=Homo_sapiens.GRCh38.dna.primary_assembly.fa --genes=Homo_sapiens.GRCh38.102.filtered.gtf'%path_to_cellranger_5)
  os.chdir(current_path)
  
  ################################
  # Run all Cell Ranger analyses #
  ################################
  
  subdirs = {'GSE187515':['2004', '2023Rep1', '2023Rep2'], 
             'GSE212998':['GFP-pos-w1', 'GFP-pos-w2', 'GFP-pos-w1-expand', 'GFP-pos-w2-expand', 'GFP-neg-expand'], 
             'PRJEB40376':['Pool_%d'%i for i in range(1,11)]}
  
  for subdir in subdirs:
    cellranger_path = path_to_cellranger_5
    transcriptome = os.path.join(path_to_ref_data_5, 'GRCh38_102')
    if (subdir == 'GSE187515'):
      cellranger_path = os.path.join(path_to_cellranger_3)
      transcriptome = os.path.join(path_to_reference_3, 'refdata-cellranger-GRCh38-3.0.0')
    
    os.chdir(subdir)
    for id_run in subdirs[subdir]: 
      libpath = 'Library_%s.csv'%id_run
      print('%s count --id=%s --transcriptome=%s --feature-ref=FeatureReference.csv --libraries=%s --localcores=%d'%(cellranger_path, id_run, transcriptome, libpath, ncores))
      os.system('%s count --id=%s --transcriptome=%s --feature-ref=FeatureReference.csv --libraries=%s --localcores=%d'%(cellranger_path, id_run, transcriptome, libpath, ncores))
      
    os.chdir(current_path)
    
    
    
