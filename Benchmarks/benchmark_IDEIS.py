# runs benchmark for IDEIS
# first, get IDEIS here!: https://github.com/Lab-of-Adaptive-Immunity/IDEIS

import os, sys, glob
import subprocess
import multiprocessing

def run_IDEIS(software_path, reference_path, bam_path, abs_res_path, cell_select_ops, n_run, option, queue):
  print(['/usr/bin/time', '-f', '%e,%U,%S', 'python', software_path, reference_path, bam_path, abs_res_path, '-g', 'Homo-Sapiens'] + cell_select_ops)
  proc = subprocess.Popen(['/usr/bin/time', '-f', '%e,%U,%S', 'python', software_path, reference_path, bam_path, abs_res_path, '-g', 'Homo-Sapiens'] + cell_select_ops, 
    stdout = subprocess.PIPE, stderr = subprocess.PIPE)
 
  output, error = proc.communicate()
  time_results = str(error, 'utf-8') # real, user and system time; total cpu time is user + system time
  #print(time_results.split('\n')[-2])
  time_results = time_results.split('\n')[-2].strip().split(',')
  res_line = '%s\t%d\t%s\t%s\t%s\t%s\t%s\n'%(paths[0].replace('/', '_'), n_run+1, 'IDEIS_' + option, time_results[0], time_results[1], time_results[2], str(float(time_results[1]) + float(time_results[2])))
  queue.put(res_line)

if __name__ == "__main__":
  
  n_runs = 40 # number of benchmarked runs
  
  software_path = '../IDEIS/IDEIS_main.py'
  
  paths_list = ['GSE187515/2004', 'GSE187515/2023Rep1', 'GSE187515/2023Rep2', 
                'GSE212998/GFP-neg-expand', 'GSE212998/GFP-pos-w1', 'GSE212998/GFP-pos-w2', 'GSE212998/GFP-pos-w1-expand', 'GSE212998/GFP-pos-w2-expand']

  # where directories with results will be stored                
  # bams subset to only reads present within qualifying cells   
  bams_list = ['subset_bams/GSE187515_2004_possorted_genome_bam.bam',
    'subset_bams/GSE187515_2023Rep1_possorted_genome_bam.bam',
    'subset_bams/GSE187515_2023Rep2_possorted_genome_bam.bam',
    'subset_bams/GSE212998_GFP-neg-expand_possorted_genome_bam.bam',
    'subset_bams/GSE212998_GFP-pos-w1_possorted_genome_bam.bam',
    'subset_bams/GSE212998_GFP-pos-w2_possorted_genome_bam.bam',
    'subset_bams/GSE212998_GFP-pos-w1-expand_possorted_genome_bam.bam',
    'subset_bams/GSE212998_GFP-pos-w2-expand_possorted_genome_bam.bam',]

  # bams subset to only reads present within qualifying cells       
  whitelist_paths = ['../GSE187515/Whitelists/Whitelist_2004.csv',
    '../GSE187515/Whitelists/Whitelist_2023Rep1.csv',
		'../GSE187515/Whitelists/Whitelist_2023Rep2.csv',
		'../GSE212998/Whitelists/Whitelist_GFP_neg_inf.csv',
		'../GSE212998/Whitelists/Whitelist_GFP-pos-w1.csv',
		'../GSE212998/Whitelists/Whitelist_GFP-pos-w2.csv',
		'../GSE212998/Whitelists/Whitelist_GFP_pos-w1-expand.csv',
		'../GSE212998/Whitelists/Whitelist_GFP_pos-w2-expand.csv']

  options = ['force_cells', 'whitelist']
  current_path = os.getcwd()

  with open('IDEIS_benchmark_report.csv', 'w') as report:
    report.write('data\trun\tanalysis\treal.time\tuser.time\tsystem.time\tCPU.time\n')
    for count, paths in enumerate(zip(paths_list, bams_list, whitelist_paths)):
      for option in options:
        process_list = []
        report_lines = []
        queue = multiprocessing.Queue()
        
        reference_path = 'Homo_Sapiens_Ref'
        transcriptome_path = 'Transcriptome'

        os.makedirs(reference_path, exist_ok = True)
        os.makedirs(transcriptome_path, exist_ok = True)
        
        bam_path = paths[1]

        cell_select_ops = []
        if option == 'force_cells':    
          cell_select_ops = ['--force-cells', '3000']
        elif option == 'whitelist':
          cell_select_ops = ['--whitelist', paths[2]]
      
        # make indexing run
        os.makedirs(os.path.join(paths[0], 'IDEIS_Benchmark_indexing_%s'%(option)), exist_ok = True)
        abs_res_path_indexing = os.path.join(paths[0], 'IDEIS_Benchmark_indexing_%s'%(option))
        
        proc = subprocess.Popen(['/usr/bin/time', '-f', '%e,%U,%S', 'python', software_path, reference_path, bam_path, abs_res_path_indexing, '-g', 'Homo-Sapiens'] + cell_select_ops, 
            stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        print('Indexing run (not counted)')
        output, error = proc.communicate()
          
        for n_run in range(n_runs):
          print('Benchmarking run %d ...'%n_run)
          os.makedirs(os.path.join(paths[0], 'IDEIS_Benchmark_%d_%s'%(n_run, option)), exist_ok = True)
          abs_res_path = os.path.join(paths[0], 'IDEIS_Benchmark_%d_%s'%(n_run, option))
          bam_path = paths[1]

          proc = multiprocessing.Process(target = run_IDEIS, args = (software_path, reference_path, bam_path, abs_res_path, cell_select_ops, n_run, option, queue,))
          process_list.append(proc)
          proc.start()

        for process in process_list:
          report_line = queue.get()              
          report_lines.append(report_line)
                
        for process in process_list:
          process.join()
  
        report_lines = sorted(report_lines)
        for line in report_lines:
          report.write(line)   
