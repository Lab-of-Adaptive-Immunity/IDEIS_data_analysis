# Benchmarks CD45er
# Download CD45er from here first! https://github.com/getzlab/10x-cd45-isoform-quantification

import os, sys, glob
import subprocess
import multiprocessing

def run_CD45er(software_path, bam_path, reference_path, n_run, queue):
  proc = subprocess.Popen(['/usr/bin/time', '-f', '%e,%U,%S', 'python', software_path, bam_path, reference_path], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
  output, error = proc.communicate()

  time_results = str(error, 'utf-8') # real, user and system time; total cpu time is user + system time
  print('On time:')
  print(time_results.split('\n')[-2])
  time_results = time_results.split('\n')[-2].strip().split(',')
  res_line = '%s\t%d\t%s\t%s\t%s\t%s\t%s\n'%(paths[0].replace('/', '_'), n_run+1, 'CD45er', time_results[0], time_results[1], time_results[2], str(float(time_results[1]) + float(time_results[2])))
  queue.put(res_line)
  
if __name__ == "__main__":
  
  n_runs = 40 # number of benchmarked runs
  
  software_path = '../../../../10x-cd45-isoform-quantification/CD45er/run.py'
  reference_path = '../../../../10x-cd45-isoform-quantification/reference/CD45_exons_nochr.hg38.txt'
  
  paths_list = ['GSE187515/2004', 'GSE187515/2023Rep1', 'GSE187515/2023Rep2', 
                'GSE212998/GFP-neg-expand', 'GSE212998/GFP-pos-w1', 'GSE212998/GFP-pos-w2', 'GSE212998/GFP-pos-w1-expand', 'GSE212998/GFP-pos-w2-expand']
  bams_list = ['../../../subset_bams/GSE187515_2004_possorted_genome_bam.bam',
    '../../../subset_bams/GSE187515_2023Rep1_possorted_genome_bam.bam',
    '../../../subset_bams/GSE187515_2023Rep2_possorted_genome_bam.bam',
    '../../../subset_bams/GSE212998_GFP-neg-expand_possorted_genome_bam.bam',
    '../../../subset_bams/GSE212998_GFP-pos-w1_possorted_genome_bam.bam',
    '../../../subset_bams/GSE212998_GFP-pos-w2_possorted_genome_bam.bam',
    '../../../subset_bams/GSE212998_GFP-pos-w1-expand_possorted_genome_bam.bam',
    '../../../subset_bams/GSE212998_GFP-pos-w2-expand_possorted_genome_bam.bam',]

  current_path = os.getcwd()
  
  with open('CD45er_benchmark_report.csv', 'w') as report:
    report.write('data\trun\tanalysis\treal.time\tuser.time\tsystem.time\tCPU.time\n')
    for count, paths in enumerate(zip(paths_list, bams_list)):
      process_list = []
      report_lines = []
      queue = multiprocessing.Queue()
      
      for n_run in range(n_runs):
        print(n_run, paths[0])
        os.makedirs(os.path.join(paths[0], 'CD45er_Benchmark_%d'%n_run), exist_ok = True)
        abs_res_path = os.path.abspath(os.path.join(paths[0], 'CD45er_Benchmark_%d'%n_run))
        bam_path = paths[1]
        os.chdir(abs_res_path)
        print(os.getcwd())
        process = multiprocessing.Process(target = run_CD45er, args=(software_path, bam_path, reference_path, n_run, queue,))
        process_list.append(process)
        
        process.start()
        os.chdir(current_path)
        
      for process in process_list:
        report_line = queue.get()              
        report_lines.append(report_line)
              
      for process in process_list:
        process.join()

      report_lines = sorted(report_lines)
      for line in report_lines:
        report.write(line)   

