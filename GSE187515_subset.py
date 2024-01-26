# Performs subsampling of bam file from GSE187515 and

import os, sys , glob

ncores = 20
bam_path = 'Benchmarks/subset_bams'
subsampled_path = 'GSE187515/Subset_bam'
whitelists = 'GSE187515/Whitelists'
source_exps = ['2004', '2023Rep1', '2023Rep2']
IDEIS_main = 'IDEIS/IDEIS_main.py'

target_path = 'GSE187515/IDEIS_subset'

if __name__ == "__main__":
  for source_exp in source_exps:
    bam_path_complete = os.path.join(bam_path, 'GSE187515_' + source_exp + '_possorted_genome_bam.bam')
    subsampled_path_complete = os.path.join(subsampled_path, source_exp)
    target_path_complete = os.path.join(target_path, '%s_subset'%source_exp)
    
    if not os.path.isdir(subsampled_path_complete):
      os.makedirs(subsampled_path_complete)
    
    if not os.path.isdir(target_path_complete):
      os.makedirs(target_path_complete)
      
    for perc in range(5, 105, 5):
      # first do subsampling and index new bamfile
      subsampled_file =  os.path.join(subsampled_path_complete, "%s_possorted_genome_bam_%d.bam"%(source_exp, perc))
      
      print('samtools view -s %f -b %s -@ %d > %s'%(0.01*perc, bam_path_complete, ncores, subsampled_file))
      os.system('samtools view --subsample %f --subsample-seed 42 -b %s -@ %d > %s'%(0.01*perc, bam_path_complete, ncores, subsampled_file))
      os.system('samtools index %s'%(subsampled_file))
      
      # now run IDEIS
      result_dir = os.path.join(target_path_complete, "%s_subset_%d"%(source_exp, perc))
      whitelist_path = os.path.join(whitelists, 'Whitelist_%s.csv'%source_exp)
      print('python IDEIS_main.py Homo-Sapiens %s %s -g Homo-Sapiens --whitelist %s'%(subsampled_file, result_dir, whitelist_path))
      os.system('python %s Homo-Sapiens %s %s -g Homo-Sapiens --whitelist %s'%(IDEIS_main, subsampled_file, result_dir, whitelist_path))
      
      auto_mat_path = os.path.join(result_dir, 'Results' , 'Iso_Counts', 'raw_matrix.rds')
      new_mat_path = os.path.join('GSE187515', 'Ptprc', 'Ptprc_%s_subset_%d.rds'%(source_exp, perc))
      new_mat_dir = os.path.dirname(new_mat_path)
      print(auto_mat_path, new_mat_path)
      if not os.path.isdir(new_mat_dir):
        os.makedirs(new_mat_dir)
      
      os.system('cp %s %s'%(auto_mat_path, new_mat_path))
