# Bash script of all links to download files for paper
#
# It is easier and faster to get files from ENA, so that server will be used in priority
#
#
# NOTE: If file crashes remove last incomplete file and run it again. Otherwise incomplete file will be skipped and not re-downloaded.
#

#########################################
# GSE187515 - Data from Collora's paper #
#########################################

cd GSE187515

mkdir Fastqs
cd Fastqs

# Sample GSM5668725 - 2004 - GEX
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/024/SRR16763924/SRR16763924_1.fastq.gz -O 2004GEX_S17_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/024/SRR16763924/SRR16763924_2.fastq.gz -O 2004GEX_S17_L001_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/025/SRR16763925/SRR16763925_1.fastq.gz -O 2004GEX_S17_L001_R1_002.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/025/SRR16763925/SRR16763925_2.fastq.gz -O 2004GEX_S17_L001_R2_002.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/026/SRR16763926/SRR16763926_1.fastq.gz -O 2004GEX_S17_L001_R1_003.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/026/SRR16763926/SRR16763926_2.fastq.gz -O 2004GEX_S17_L001_R2_003.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/027/SRR16763927/SRR16763927_1.fastq.gz -O 2004GEX_S17_L001_R1_004.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/027/SRR16763927/SRR16763927_2.fastq.gz -O 2004GEX_S17_L001_R2_004.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/028/SRR16763928/SRR16763928_1.fastq.gz -O 2004GEX_S17_L001_R1_005.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/028/SRR16763928/SRR16763928_2.fastq.gz -O 2004GEX_S17_L001_R2_005.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/029/SRR16763929/SRR16763929_1.fastq.gz -O 2004GEX_S17_L001_R1_006.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/029/SRR16763929/SRR16763929_2.fastq.gz -O 2004GEX_S17_L001_R2_006.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/030/SRR16763930/SRR16763930_1.fastq.gz -O 2004GEX_S17_L001_R1_007.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/030/SRR16763930/SRR16763930_2.fastq.gz -O 2004GEX_S17_L001_R2_007.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/031/SRR16763931/SRR16763931_1.fastq.gz -O 2004GEX_S17_L001_R1_008.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/031/SRR16763931/SRR16763931_2.fastq.gz -O 2004GEX_S17_L001_R2_008.fastq.gz

# Sample GSM5668726 - 2023Rep1 - GEX
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/032/SRR16763932/SRR16763932_1.fastq.gz -O 2023Rep1GEX_S18_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/032/SRR16763932/SRR16763932_2.fastq.gz -O 2023Rep1GEX_S18_L001_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/033/SRR16763933/SRR16763933_1.fastq.gz -O 2023Rep1GEX_S18_L001_R1_002.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/033/SRR16763933/SRR16763933_2.fastq.gz -O 2023Rep1GEX_S18_L001_R2_002.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/034/SRR16763934/SRR16763934_1.fastq.gz -O 2023Rep1GEX_S18_L001_R1_003.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/034/SRR16763934/SRR16763934_2.fastq.gz -O 2023Rep1GEX_S18_L001_R2_003.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/035/SRR16763935/SRR16763935_1.fastq.gz -O 2023Rep1GEX_S18_L001_R1_004.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/035/SRR16763935/SRR16763935_2.fastq.gz -O 2023Rep1GEX_S18_L001_R2_004.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/036/SRR16763936/SRR16763936_1.fastq.gz -O 2023Rep1GEX_S18_L001_R1_005.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/036/SRR16763936/SRR16763936_2.fastq.gz -O 2023Rep1GEX_S18_L001_R2_005.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/037/SRR16763937/SRR16763937_1.fastq.gz -O 2023Rep1GEX_S18_L001_R1_006.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/037/SRR16763937/SRR16763937_2.fastq.gz -O 2023Rep1GEX_S18_L001_R2_006.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/038/SRR16763938/SRR16763938_1.fastq.gz -O 2023Rep1GEX_S18_L001_R1_007.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/038/SRR16763938/SRR16763938_2.fastq.gz -O 2023Rep1GEX_S18_L001_R2_007.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/039/SRR16763939/SRR16763939_1.fastq.gz -O 2023Rep1GEX_S18_L001_R1_008.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/039/SRR16763939/SRR16763939_2.fastq.gz -O 2023Rep1GEX_S18_L001_R2_008.fastq.gz

# Sample GSM5668727 - 2023Rep2 - GEX
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/040/SRR16763940/SRR16763940_1.fastq.gz -O 2023Rep2GEX_S19_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/040/SRR16763940/SRR16763940_2.fastq.gz -O 2023Rep2GEX_S19_L001_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/041/SRR16763941/SRR16763941_1.fastq.gz -O 2023Rep2GEX_S19_L001_R1_002.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/041/SRR16763941/SRR16763941_2.fastq.gz -O 2023Rep2GEX_S19_L001_R2_002.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/042/SRR16763942/SRR16763942_1.fastq.gz -O 2023Rep2GEX_S19_L001_R1_003.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/042/SRR16763942/SRR16763942_2.fastq.gz -O 2023Rep2GEX_S19_L001_R2_003.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/043/SRR16763943/SRR16763943_1.fastq.gz -O 2023Rep2GEX_S19_L001_R1_004.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/043/SRR16763943/SRR16763943_2.fastq.gz -O 2023Rep2GEX_S19_L001_R2_004.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/044/SRR16763944/SRR16763944_1.fastq.gz -O 2023Rep2GEX_S19_L001_R1_005.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/044/SRR16763944/SRR16763944_2.fastq.gz -O 2023Rep2GEX_S19_L001_R2_005.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/045/SRR16763945/SRR16763945_1.fastq.gz -O 2023Rep2GEX_S19_L001_R1_006.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/045/SRR16763945/SRR16763945_2.fastq.gz -O 2023Rep2GEX_S19_L001_R2_006.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/046/SRR16763946/SRR16763946_1.fastq.gz -O 2023Rep2GEX_S19_L001_R1_007.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/046/SRR16763946/SRR16763946_2.fastq.gz -O 2023Rep2GEX_S19_L001_R2_007.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/047/SRR16763947/SRR16763947_1.fastq.gz -O 2023Rep2GEX_S19_L001_R1_008.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/047/SRR16763947/SRR16763947_2.fastq.gz -O 2023Rep2GEX_S19_L001_R2_008.fastq.gz

# Sample GSM5668744 - 2004 - CITE-seq
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/092/SRR16763992/SRR16763992_1.fastq.gz -O 2004FB_S17_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/092/SRR16763992/SRR16763992_2.fastq.gz -O 2004FB_S17_L001_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/093/SRR16763993/SRR16763993_1.fastq.gz -O 2004FB_S17_L001_R1_002.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/093/SRR16763993/SRR16763993_2.fastq.gz -O 2004FB_S17_L001_R2_002.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/094/SRR16763994/SRR16763994_1.fastq.gz -O 2004FB_S17_L001_R1_003.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/094/SRR16763994/SRR16763994_2.fastq.gz -O 2004FB_S17_L001_R2_003.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/095/SRR16763995/SRR16763995_1.fastq.gz -O 2004FB_S17_L001_R1_004.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/095/SRR16763995/SRR16763995_2.fastq.gz -O 2004FB_S17_L001_R2_004.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/096/SRR16763996/SRR16763996_1.fastq.gz -O 2004FB_S17_L001_R1_005.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/096/SRR16763996/SRR16763996_2.fastq.gz -O 2004FB_S17_L001_R2_005.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/097/SRR16763997/SRR16763997_1.fastq.gz -O 2004FB_S17_L001_R1_006.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/097/SRR16763997/SRR16763997_2.fastq.gz -O 2004FB_S17_L001_R2_006.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/098/SRR16763998/SRR16763998_1.fastq.gz -O 2004FB_S17_L001_R1_007.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/098/SRR16763998/SRR16763998_2.fastq.gz -O 2004FB_S17_L001_R2_007.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/099/SRR16763999/SRR16763999_1.fastq.gz -O 2004FB_S17_L001_R1_008.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/099/SRR16763999/SRR16763999_2.fastq.gz -O 2004FB_S17_L001_R2_008.fastq.gz

# Sample GSM5668745 - 2023Rep1 - CITE-seq
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/000/SRR16764000/SRR16764000_1.fastq.gz -O 2023Rep1FB_S18_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/000/SRR16764000/SRR16764000_2.fastq.gz -O 2023Rep1FB_S18_L001_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/001/SRR16764001/SRR16764001_1.fastq.gz -O 2023Rep1FB_S18_L001_R1_002.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/001/SRR16764001/SRR16764001_2.fastq.gz -O 2023Rep1FB_S18_L001_R2_002.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/002/SRR16764002/SRR16764002_1.fastq.gz -O 2023Rep1FB_S18_L001_R1_003.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/002/SRR16764002/SRR16764002_2.fastq.gz -O 2023Rep1FB_S18_L001_R2_003.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/003/SRR16764003/SRR16764003_1.fastq.gz -O 2023Rep1FB_S18_L001_R1_004.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/003/SRR16764003/SRR16764003_2.fastq.gz -O 2023Rep1FB_S18_L001_R2_004.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/004/SRR16764004/SRR16764004_1.fastq.gz -O 2023Rep1FB_S18_L001_R1_005.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/004/SRR16764004/SRR16764004_2.fastq.gz -O 2023Rep1FB_S18_L001_R2_005.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/005/SRR16764005/SRR16764005_1.fastq.gz -O 2023Rep1FB_S18_L001_R1_006.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/005/SRR16764005/SRR16764005_2.fastq.gz -O 2023Rep1FB_S18_L001_R2_006.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/006/SRR16764006/SRR16764006_1.fastq.gz -O 2023Rep1FB_S18_L001_R1_007.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/006/SRR16764006/SRR16764006_2.fastq.gz -O 2023Rep1FB_S18_L001_R2_007.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/007/SRR16764007/SRR16764007_1.fastq.gz -O 2023Rep1FB_S18_L001_R1_008.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/007/SRR16764007/SRR16764007_2.fastq.gz -O 2023Rep1FB_S18_L001_R2_008.fastq.gz

# Sample GSM5668746 - 2023Rep2 - CITE-seq
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/008/SRR16764008/SRR16764008_1.fastq.gz -O 2023Rep2FB_S19_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/008/SRR16764008/SRR16764008_2.fastq.gz -O 2023Rep2FB_S19_L001_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/009/SRR16764009/SRR16764009_1.fastq.gz -O 2023Rep2FB_S19_L001_R1_002.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/009/SRR16764009/SRR16764009_2.fastq.gz -O 2023Rep2FB_S19_L001_R2_002.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/010/SRR16764010/SRR16764010_1.fastq.gz -O 2023Rep2FB_S19_L001_R1_003.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/010/SRR16764010/SRR16764010_2.fastq.gz -O 2023Rep2FB_S19_L001_R2_003.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/011/SRR16764011/SRR16764011_1.fastq.gz -O 2023Rep2FB_S19_L001_R1_004.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/011/SRR16764011/SRR16764011_2.fastq.gz -O 2023Rep2FB_S19_L001_R2_004.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/012/SRR16764012/SRR16764012_1.fastq.gz -O 2023Rep2FB_S19_L001_R1_005.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/012/SRR16764012/SRR16764012_2.fastq.gz -O 2023Rep2FB_S19_L001_R2_005.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/013/SRR16764013/SRR16764013_1.fastq.gz -O 2023Rep2FB_S19_L001_R1_006.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/013/SRR16764013/SRR16764013_2.fastq.gz -O 2023Rep2FB_S19_L001_R2_006.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/014/SRR16764014/SRR16764014_1.fastq.gz -O 2023Rep2FB_S19_L001_R1_007.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/014/SRR16764014/SRR16764014_2.fastq.gz -O 2023Rep2FB_S19_L001_R2_007.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/015/SRR16764015/SRR16764015_1.fastq.gz -O 2023Rep2FB_S19_L001_R1_008.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR167/015/SRR16764015/SRR16764015_2.fastq.gz -O 2023Rep2FB_S19_L001_R2_008.fastq.gz

cd ../..

####################################
# GSE212998 - Data from Yu's paper #
####################################

cd GSE212998

mkdir Fastqs
cd Fastqs

# GEX for all samples

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR214/005/SRR21497505/SRR21497505_1.fastq.gz -O GFP-pos-w1-GEX_S11_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR214/005/SRR21497505/SRR21497505_2.fastq.gz -O GFP-pos-w1-GEX_S11_L001_R2_001.fastq.gz

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR214/004/SRR21497504/SRR21497504_1.fastq.gz -O GFP-pos-w2-GEX_S12_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR214/004/SRR21497504/SRR21497504_2.fastq.gz -O GFP-pos-w2-GEX_S12_L001_R2_001.fastq.gz

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR214/093/SRR21497493/SRR21497493_1.fastq.gz -O GFP-pos-w1-expand-GEX_S21_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR214/093/SRR21497493/SRR21497493_2.fastq.gz -O GFP-pos-w1-expand-GEX_S21_L001_R2_001.fastq.gz 

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR214/092/SRR21497492/SRR21497492_1.fastq.gz -O GFP-pos-w2-expand-GEX_S22_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR214/092/SRR21497492/SRR21497492_2.fastq.gz -O GFP-pos-w2-expand-GEX_S22_L001_R2_001.fastq.gz

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR214/091/SRR21497491/SRR21497491_1.fastq.gz -O GFP-neg-w1-expand-GEX_S23_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR214/091/SRR21497491/SRR21497491_2.fastq.gz -O GFP-neg-w1-expand-GEX_S23_L001_R2_001.fastq.gz

# CITE-seq for all samples

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR214/099/SRR21497499/SRR21497499_1.fastq.gz -O GFP-pos-w1-CITE_S11_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR214/099/SRR21497499/SRR21497499_2.fastq.gz -O GFP-pos-w1-CITE_S11_L001_R2_001.fastq.gz

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR214/098/SRR21497498/SRR21497498_1.fastq.gz -O GFP-pos-w2-CITE_S12_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR214/098/SRR21497498/SRR21497498_2.fastq.gz -O GFP-pos-w2-CITE_S12_L001_R2_001.fastq.gz

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR214/087/SRR21497487/SRR21497487_1.fastq.gz -O GFP-pos-w1-expand-CITE_S21_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR214/087/SRR21497487/SRR21497487_2.fastq.gz -O GFP-pos-w1-expand-CITE_S21_L001_R2_001.fastq.gz

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR214/086/SRR21497486/SRR21497486_1.fastq.gz -O GFP-pos-w2-expand-CITE_S22_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR214/086/SRR21497486/SRR21497486_2.fastq.gz -O GFP-pos-w2-expand-CITE_S22_L001_R2_001.fastq.gz

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR214/085/SRR21497485/SRR21497485_1.fastq.gz -O GFP-neg-w1-expand-CITE_S23_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR214/085/SRR21497485/SRR21497485_2.fastq.gz -O GFP-neg-w1-expand-CITE_S23_L001_R2_001.fastq.gz

cd ../..

###########################################
# GSE155006 - Data from Mogilenko's paper #
###########################################

cd GSE155006

mkdir Bam
cd Bam

# download the BAM files
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/SRR111/SRR11114659/Liver_Aged.bam.1 -O Liver_Aged.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/SRR111/SRR11114659/Liver_Aged.bam.1.bai -O Liver_Aged.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/SRR111/SRR11114660/Liver_Young.bam.1 -O Liver_Young.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/SRR111/SRR11114660/Liver_Young.bam.1.bai -O Liver_Young.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/SRR111/SRR11114661/Lung_Aged.bam.1 -O Lung_Aged.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/SRR111/SRR11114661/Lung_Aged.bam.1.bai -O Lung_Aged.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/SRR111/SRR11114662/Lung_Young.bam.1 -O Lung_Young.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/SRR111/SRR11114662/Lung_Young.bam.1.bai -O Lung_Young.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/SRR111/SRR11114663/Peritoneal_cells_Aged.bam.1 -O Peritoneal_cells_Aged.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/SRR111/SRR11114663/Peritoneal_cells_Aged.bam.1.bai -O Peritoneal_cells_Aged.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/SRR111/SRR11114664/Peritoneal_cells_Young.bam.1 -O Peritoneal_cells_Young.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/SRR111/SRR11114664/Peritoneal_cells_Young.bam.1.bai -O Peritoneal_cells_Young.bam.bai  
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/SRR111/SRR11114665/Spleen_Aged.bam.1 -O Spleen_Aged.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/SRR111/SRR11114665/Spleen_Aged.bam.1.bai -O Spleen_Aged.bam.bai
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/SRR111/SRR11114666/Spleen_Young.bam.1 -O Spleen_Young.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/SRR111/SRR11114666/Spleen_Young.bam.1.bai -O Spleen_Young.bam.bai

cd ..
mkdir Count_matrices
cd Count_matrices

mkdir GSM4321524_Liver_Aged GSM4321525_Liver_Young GSM4321526_Lungs_Aged GSM4321527_Lungs_Young 
mkdir GSM4321528_Peritoneal_cells_Aged GSM4321529_Peritoneal_cells_Young GSM4321530_Spleen_Aged GSM4321531_Spleen_Young

wget -nc 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE155006&format=file' -O count_matrices.tar
tar -xvf count_matrices.tar

mv GSM4321524_*.gz GSM4321524_Liver_Aged
mv GSM4321525_*.gz GSM4321525_Liver_Young
mv GSM4321526_*.gz GSM4321526_Lungs_Aged
mv GSM4321527_*.gz GSM4321527_Lungs_Young
mv GSM4321528_*.gz GSM4321528_Peritoneal_cells_Aged
mv GSM4321529_*.gz GSM4321529_Peritoneal_cells_Young
mv GSM4321530_*.gz GSM4321530_Spleen_Aged
mv GSM4321531_*.gz GSM4321531_Spleen_Young

for f in GSM4321524_Liver_Aged/GSM4321524_Liver_Aged_*; do mv "$f" "${f/GSM4321524_Liver_Aged\/GSM4321524_Liver_Aged_/GSM4321524_Liver_Aged\/}";done
for f in GSM4321525_Liver_Young/GSM4321525_Liver_Young_*; do mv "$f" "${f/GSM4321525_Liver_Young\/GSM4321525_Liver_Young_/GSM4321525_Liver_Young\/}";done
for f in GSM4321526_Lungs_Aged/GSM4321526_Lung_Aged_*; do mv "$f" "${f/GSM4321526_Lungs_Aged\/GSM4321526_Lung_Aged_/GSM4321526_Lungs_Aged\/}";done
for f in GSM4321527_Lungs_Young/GSM4321527_Lung_Young_*; do mv "$f" "${f/GSM4321527_Lungs_Young\/GSM4321527_Lung_Young_/GSM4321527_Lungs_Young\/}";done
for f in GSM4321528_Peritoneal_cells_Aged/GSM4321528_Peritoneal_cells_Aged_*; do mv "$f" "${f/GSM4321528_Peritoneal_cells_Aged\/GSM4321528_Peritoneal_cells_Aged_/GSM4321528_Peritoneal_cells_Aged\/}";done
for f in GSM4321529_Peritoneal_cells_Young/GSM4321529_Peritoneal_cells_Young_*; do mv "$f" "${f/GSM4321529_Peritoneal_cells_Young\/GSM4321529_Peritoneal_cells_Young_/GSM4321529_Peritoneal_cells_Young\/}";done
for f in GSM4321530_Spleen_Aged/GSM4321530_Spleen_Aged_*; do mv "$f" "${f/GSM4321530_Spleen_Aged\/GSM4321530_Spleen_Aged_/GSM4321530_Spleen_Aged\/}";done
for f in GSM4321531_Spleen_Young/GSM4321531_Spleen_Young_*; do mv "$f" "${f/GSM4321531_Spleen_Young\/GSM4321531_Spleen_Young_/GSM4321531_Spleen_Young\/}";done

cd ../..

#########################################
# PRJEB40376 - Data from Lawlor's paper #
#########################################

cd PRJEB40376

mkdir Fastqs
cd Fastqs

# RNA
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448067/1-RNA_pooled_and_merged_S1_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448067/1-RNA_pooled_and_merged_S1_L001_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448092/2-RNA_pooled_and_merged_S1_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448092/2-RNA_pooled_and_merged_S1_L001_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448083/3-RNA_pooled_and_merged_S1_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448083/3-RNA_pooled_and_merged_S1_L001_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448080/4-RNA_pooled_and_merged_S1_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448080/4-RNA_pooled_and_merged_S1_L001_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448082/5-RNA_pooled_and_merged_S1_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448082/5-RNA_pooled_and_merged_S1_L001_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448074/6-RNA_pooled_and_merged_S1_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448074/6-RNA_pooled_and_merged_S1_L001_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448081/7-RNA_pooled_and_merged_S1_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448081/7-RNA_pooled_and_merged_S1_L001_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448086/8-RNA_pooled_and_merged_S1_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448086/8-RNA_pooled_and_merged_S1_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448089/9-RNA_pooled_and_merged_S1_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448089/9-RNA_pooled_and_merged_S1_L001_R1_001.fastq.gz 
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448066/10-RNA_pooled_and_merged_S1_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448066/10-RNA_pooled_and_merged_S1_L001_R2_001.fastq.gz

# ADT
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448068/1-ADT_pooled_and_merged_S1_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448068/1-ADT_pooled_and_merged_S1_L001_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448070/2-ADT_pooled_and_merged_S1_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448070/2-ADT_pooled_and_merged_S1_L001_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448064/3-ADT_pooled_and_merged_S1_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448064/3-ADT_pooled_and_merged_S1_L001_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448065/4-ADT_pooled_and_merged_S1_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448065/4-ADT_pooled_and_merged_S1_L001_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448072/5-ADT_pooled_and_merged_S1_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448072/5-ADT_pooled_and_merged_S1_L001_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448077/6-ADT_pooled_and_merged_S1_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448077/6-ADT_pooled_and_merged_S1_L001_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448084/7-ADT_pooled_and_merged_S1_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448084/7-ADT_pooled_and_merged_S1_L001_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448063/8-ADT_pooled_and_merged_S1_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448063/8-ADT_pooled_and_merged_S1_L001_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448090/9-ADT_pooled_and_merged_S1_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448090/9-ADT_pooled_and_merged_S1_L001_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448087/10-ADT_pooled_and_merged_S1_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448087/10-ADT_pooled_and_merged_S1_L001_R2_001.fastq.gz

# HTO
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448088/1-HTO_pooled_and_merged_S1_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448088/1-HTO_pooled_and_merged_S1_L001_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448078/2-HTO_pooled_and_merged_S1_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448078/2-HTO_pooled_and_merged_S1_L001_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448091/3-HTO_pooled_and_merged_S1_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448091/3-HTO_pooled_and_merged_S1_L001_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448085/4-HTO_pooled_and_merged_S1_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448085/4-HTO_pooled_and_merged_S1_L001_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448076/5-HTO_pooled_and_merged_S1_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448076/5-HTO_pooled_and_merged_S1_L001_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448075/6-HTO_pooled_and_merged_S1_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448075/6-HTO_pooled_and_merged_S1_L001_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448079/7-HTO_pooled_and_merged_S1_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448079/7-HTO_pooled_and_merged_S1_L001_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448073/8-HTO_pooled_and_merged_S1_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448073/8-HTO_pooled_and_merged_S1_L001_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448071/9-HTO_pooled_and_merged_S1_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448071/9-HTO_pooled_and_merged_S1_L001_R2_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448069/10-HTO_pooled_and_merged_S1_L001_R1_001.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR744/ERR7448069/10-HTO_pooled_and_merged_S1_L001_R2_001.fastq.gz

cd ../..
