import os, sys, glob
import subprocess

if __name__ == "__main__":
  with open('GSE187515/counted_bam.txt', 'w') as report:
    report.write('data\tpercent\tread.count\n')
    dircs = sorted(glob.glob('GSE187515/Subset_bam/20*'))
    
    print(dircs)
    for dirc in dircs:
      for perc in range(5, 105, 5):
        bam_path = os.path.join(dirc, '%s_possorted_genome_bam_%d.bam'%(os.path.basename(dirc.strip()), perc))
        print(bam_path)
        
        proc = subprocess.Popen(['samtools', 'view', '-c', '-F', '260', bam_path], stdout = subprocess.PIPE)
        (output, error) = proc.communicate()
        read_count = str(output, 'utf-8').strip()
        report.write('%s\t%d\t%s\n'%(dirc.strip(), perc, read_count))
        print(read_count)
