import os, sys,glob

if __name__ == '__main__':
	script_path = os.path.dirname(os.path.realpath(__file__)) # get path of this script
	file_list = sorted(glob.glob(script_path + '/GSE*/Library*template*') + glob.glob(script_path + '/PRJEB*/Library*template*'))
	
	for path in file_list:
		with open(path, 'r') as template:
			contents = template.readlines()
			out_path = path.replace('_template','')
			with open(out_path, 'w') as injected:
				for line in contents:
					injected.write(line.replace('/path/to/IDEIS_data_analysis', script_path))
				
	
