# subsets file to only reads from cells in analyzed objects
# you need subset-bam. Get it here (for Linux): https://github.com/10XGenomics/subset-bam/releases/download/v1.1.0/subset-bam_linux

import os

if __name__ == '__main__':	

	ncores = 32	

	bams_list = ['../GSE187515/2004/outs/possorted_genome_bam.bam',
		'../GSE187515/2023Rep1/outs/possorted_genome_bam.bam',
		'../GSE187515/2023Rep2/outs/possorted_genome_bam.bam',
		'../GSE212998/GFP-pos-w1/outs/possorted_genome_bam.bam',
		'../GSE212998/GFP-pos-w2/outs/possorted_genome_bam.bam',
		'../GSE212998/GFP-pos-w1-expand/outs/possorted_genome_bam.bam',
		'../GSE212998/GFP-pos-w2-expand/outs/possorted_genome_bam.bam',
		'../GSE212998/GFP-neg-expand/outs/possorted_genome_bam.bam']	

	whitelist = ['../GSE187515/Whitelists/Whitelist_2004.csv',
		'../GSE187515/Whitelists/Whitelist_2023Rep1.csv',
		'../GSE187515/Whitelists/Whitelist_2023Rep2.csv',
		'../GSE212998/Whitelists/Whitelist_GFP-pos-w1.csv',
		'../GSE212998/Whitelists/Whitelist_GFP-pos-w2.csv',
		'../GSE212998/Whitelists/Whitelist_GFP_pos-w1-expand.csv',
		'../GSE212998/Whitelists/Whitelist_GFP_pos-w2-expand.csv',
		'../GSE212998/Whitelists/Whitelist_GFP_neg_inf.csv']
	
	if not os.path.exists('subset_bams'):
		os.makedirs('subset_bams')
	
	for i in range(0, len(whitelist)):
		with open('whitelist.tmp', 'w') as tmplist:
			with open(whitelist[i], 'r') as permlist:
				for line in permlist:
					tmplist.write(line.strip() + '-1\n')
	
	  
		out_name = os.path.join('subset_bams', '_'.join(bams_list[i].split('/')[1:3]) + '_possorted_genome_bam.bam')
		
		print(out_name)
		print('./subset-bam_linux --bam %s --cell-barcodes whitelist.tmp --out-bam %s --cores %d'%(bams_list[i], out_name, ncores))
		os.system('./subset-bam_linux --bam %s --cell-barcodes whitelist.tmp --out-bam %s --cores %d'%(bams_list[i], out_name, ncores))
		os.system('samtools index %s'%out_name)
	
