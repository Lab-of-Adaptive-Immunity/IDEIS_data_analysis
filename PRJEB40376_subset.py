# Performs subsampling of bam file from GSE187515 and IDEIS analysis

import os, sys , glob

ncores = 20
bam_path = 'PRJEB40376/subset_bams'
subsampled_path = 'PRJEB40376/Subsampled_bams'
whitelists = 'PRJEB40376/Whitelists'
source_exps = ['Pool_1']
IDEIS_main = 'IDEIS/IDEIS_main.py'

target_path = 'PRJEB40376/IDEIS_subset'

if __name__ == "__main__":
  for source_exp in source_exps:
    bam_path_complete = os.path.join(bam_path, 'PRJEB40376_' + source_exp + '_possorted_genome_bam.bam')
    subsampled_path_complete = os.path.join(subsampled_path, source_exp)
    target_path_complete = os.path.join(target_path, '%s_subset'%source_exp)
    
    if not os.path.isdir(subsampled_path_complete):
      os.makedirs(subsampled_path_complete)
    
    if not os.path.isdir(target_path_complete):
      os.makedirs(target_path_complete)
      
    for perc in range(25, 105, 25):
      # first do subsampling and index new bamfile
      subsampled_file =  os.path.join(subsampled_path_complete, "%s_possorted_genome_bam_%d.bam"%(source_exp, perc))
      
      print('samtools view -s %f -b %s -@ %d > %s'%(0.01*perc, bam_path_complete, ncores, subsampled_file))
      os.system('samtools view --subsample %f --subsample-seed 42 -b %s -@ %d > %s'%(0.01*perc, bam_path_complete, ncores, subsampled_file))
      os.system('samtools index %s'%(subsampled_file))
      
      # split by length
      subsampled_file_minus_100 = os.path.join(subsampled_path_complete, "%s_possorted_genome_bam_%d_minus100.bam"%(source_exp, perc))
      subsampled_file_plus_100 = os.path.join(subsampled_path_complete, "%s_possorted_genome_bam_%d_plus100.bam"%(source_exp, perc))
      print("samtools view -bS -e 'length(seq)<=100' -@ %d %s > %s"%(ncores, subsampled_file, subsampled_file_minus_100))
      print("samtools view -bS -e 'length(seq)>100' -@ %d %s > %s"%(ncores, subsampled_file, subsampled_file_plus_100))
      os.system("samtools view -bS -e 'length(seq)<=100' -@ %d %s > %s"%(ncores, subsampled_file, subsampled_file_minus_100))
      os.system("samtools view -bS -e 'length(seq)>100' -@ %d %s > %s"%(ncores, subsampled_file, subsampled_file_plus_100))
      os.system('samtools index %s'%(subsampled_file_minus_100))
      os.system('samtools index %s'%(subsampled_file_plus_100))
      
      # now run IDEIS
      result_dir_minus100 = os.path.join(target_path_complete, "%s_subset_%d_minus100"%(source_exp, perc))
      result_dir_plus100 = os.path.join(target_path_complete, "%s_subset_%d_plus100"%(source_exp, perc))
      whitelist_path = os.path.join(whitelists, 'Whitelist_%s.csv'%source_exp)
      print('python %s Homo-Sapiens %s %s -g Homo-Sapiens --sequencing-type 3-prime --whitelist %s'%(IDEIS_main, subsampled_file_minus_100, result_dir_minus100, whitelist_path))
      os.system('python %s Homo-Sapiens %s %s -g Homo-Sapiens --sequencing-type 3-prime --whitelist %s'%(IDEIS_main, subsampled_file_minus_100, result_dir_minus100, whitelist_path))
      print('python %s Homo-Sapiens %s %s -g Homo-Sapiens --sequencing-type 3-prime --whitelist %s'%(IDEIS_main, subsampled_file_plus_100, result_dir_plus100, whitelist_path))
      os.system('python %s Homo-Sapiens %s %s -g Homo-Sapiens --sequencing-type 3-prime --whitelist %s'%(IDEIS_main, subsampled_file_plus_100, result_dir_plus100, whitelist_path))
      
      # move results to final directory
      auto_mat_path_minus100 = os.path.join(result_dir_minus100, 'Results' , 'Iso_Counts', 'raw_matrix.rds')
      auto_mat_path_plus100 = os.path.join(result_dir_plus100, 'Results' , 'Iso_Counts', 'raw_matrix.rds')
      new_mat_path_minus100 = os.path.join('PRJEB40376', 'Ptprc', 'Ptprc_%s_subset_%d_minus100.rds'%(source_exp, perc))
      new_mat_path_plus100 = os.path.join('PRJEB40376', 'Ptprc', 'Ptprc_%s_subset_%d_plus100.rds'%(source_exp, perc))
      new_mat_dir_minus100 = os.path.dirname(new_mat_path_minus100)
      new_mat_dir_plus100 = os.path.dirname(new_mat_path_plus100)
      print(auto_mat_path_minus100, new_mat_path_minus100, auto_mat_path_plus100, new_mat_path_plus100)
      if not os.path.isdir(new_mat_dir_minus100):
        os.makedirs(new_mat_dir_minus100)

      if not os.path.isdir(new_mat_dir_plus100):
        os.makedirs(new_mat_dir_plus100)
      
      os.system('cp %s %s'%(auto_mat_path_minus100, new_mat_path_minus100))
      os.system('cp %s %s'%(auto_mat_path_plus100, new_mat_path_plus100))
