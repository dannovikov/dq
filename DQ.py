#!/usr/bin/python3.8

# VERSION FROM JANUARY 25, 2017
# UPGRADED TO BE PYTHON3 COMPATIBLE ON JULY 25, 2022

# THIS SCRIPT WILL USE BEALIGN AND TN93DIST TO DETERMINE IF A SET OF QUERY SEQUENCES
# CONTAINS HIV PROTEASE, RT, AND/OR INTEGRASE REFERENCE SEQUENCES, TO BE WRITTEN LOCALLY AS SCRIPT RUNS 
# AND IT WILL DETERMINE IF ANY OF THESE SEQUENCES ARE IN REVERSE ORIENTATION
# IT WILL ALSO CHECK FOR SEQUENCES WITH STRINGS OF NNNs AND NON-IUPAC CHARACTERS
# IT WILL ALSO LOOK TO SEE IF RT AND PR ARE PRESENT AND IN CORRECT ORIENTATION, BUT IN THE WRONG ORDER
# IT WILL ALSO LOOK FOR SEQUENCES THAT LINK TO HXB2 IN THE PRRT REGION (1.5% FULLY RESOLVED TN93 DISTANCE)

# AND THIS SCRIPT TAKES A LONG TIME TO RUN WHEN THERE ARE LOTS OF SEQUENCES, SORRY

import os
import sys, getopt
import re
import time

# READ IN INPUT FASTA FILE NAME AND CREATE PATH FOR OUTPUT FILES
# INPUT FILE MUST END IN '.fa' OR '.fasta'
# EACH SEQUENCE MUST HAVE A UNIQUE ID
 
def usage():
	print('\nERROR: use -i to provide FULL path to input (.fa or .fasta) file\n')
	print('To skip bealign and TN93dist, use -s 1')
try:
	options, args = getopt.getopt(sys.argv[1:], "i:s:")
except getopt.GetoptError as err:
	print('\n'+str(err))
	usage()
	sys.exit(2)

skip = 0

for opt,arg in options:
	if opt == '-i': in_file = arg # input file path and name
	if opt == '-s': skip = int(arg)
if len(options) < 1: usage()

#in_file = '/Users/Joel/MHS/Superinfection/Test test/Test/Input.fa'

# DETERMINE PROPER SUFFIX FOR FASTA FILE
if '.fasta' in in_file: f = '.fasta'
else: f = '.fa'

in_path = ''
for AB in in_file.split('/')[:-1]: in_path = in_path + AB + '/'


# READ IN INITIAL SEQUENCE ALIGNMENT
print('\nREADING IN SEQUENCE ALIGNMENT\n')

seq_dict = {}
for line in open(str(in_file)):
	if ">" in line:
		seq_id = line.strip('>').rstrip()
		if seq_id in seq_dict: sys.exit('\nERROR: all sequences must have a unique ID\n\n'+str(seq_id)+' occurs at least twice in input file\n\nEXITING PROGRAM\n')
		seq_dict.update({seq_id:{'seq':''}})
	else: seq_dict[seq_id]['seq'] = seq_dict[seq_id]['seq'] + line.rstrip()

# IDENTIFY SEQUENCES WITH NON-IUPAC CHARACTERS
# IDENTIFY SEQUENCES WITH LONG STRINGS OF NNNs
print('CHECKING SEQUENCES FOR NON-IUPAC CHARACTERS, STRINGS OF NNNs, AND HIGH FREQUENCIES OF NT AMBIGUITIES (>5%)\n')
iupac_code = set(['A','C','G','T','U','R','Y','S','W','K','M','B','D','H','V','N','-','*'])
ambig_code = set(['R','Y','S','W','K','M','B','D','H','V','N'])
for seq_id in seq_dict:
	if set(seq_dict[seq_id]['seq']).issubset(iupac_code): seq_dict[seq_id]['UnknownCharacters'] = 0
	else:
		if '(' in seq_dict[seq_id]['seq']: seq_dict[seq_id]['seq'] = re.sub("\(.*?\)","",seq_dict[seq_id]['seq'])
		# REMOVE ALL NON-IUPAC, NON-GAP, AND NON-ASTERISK CHARACTERS FROM SEQUENCE
		# MIGHT NEED TO FIND A WAY TO HANDLE DIFFICULT CHARACTERS LIKE ? OR EVEN * 
			# TRY FOLLOWING: b = re.sub("\\"+str(char),"",b) PATTERN
		for char in set(seq_dict[seq_id]['seq']).difference(iupac_code.union('*')): seq_dict[seq_id]['seq'] = re.sub(char,"",seq_dict[seq_id]['seq'])
		seq_dict[seq_id]['UnknownCharacters'] = 1
	if 'NNNNN' in seq_dict[seq_id]['seq']: seq_dict[seq_id]['Nstring'] = 1
	else: seq_dict[seq_id]['Nstring'] = 0
	seq_dict[seq_id]['PR_Uni'] = 0
	seq_dict[seq_id]['PR_Bi'] = 0
	seq_dict[seq_id]['RT_Uni'] = 0
	seq_dict[seq_id]['RT_Bi'] = 0
	seq_dict[seq_id]['INT_Uni'] = 0
	seq_dict[seq_id]['INT_Bi'] = 0
	seq_dict[seq_id]['RTPR'] = 0
	seq_dict[seq_id]['PRRT_HXB2_Link'] = 0
	ambig_count = 0
	for nt in seq_dict[seq_id]['seq']:
		if nt in ambig_code: ambig_count += 1
		if float(ambig_count)/len(seq_dict[seq_id]['seq']) >= 0.05: seq_dict[seq_id]['Ambig5%'] = 1
		else: seq_dict[seq_id]['Ambig5%'] = 0

print('WRITING OUT NEW SEQUENCE FILES FOR FURTHER ANALYSIS (EXCLUDING SEQUENCES CONTAINING NON-IUPAC CHARACTERS)\n')

seqfile = open(str(in_file.split(f)[0])+'.NoUnkChar.fa','w')
for seq_id in seq_dict:
	if seq_dict[seq_id]['UnknownCharacters'] == 0 and len(seq_dict[seq_id]['seq']) > 100:
		seqfile.write('>'+str(seq_id)+'\n'+str(seq_dict[seq_id]['seq'])+'\n')
	del seq_dict[seq_id]['seq']

seqfile.close()

# WRITE HXB2 REFERENCE FASTA FILES TO LOCAL DIRECTORY
# PERFORM ALIGNMENT TO REFERENCES USING UNIDIRECTIONAL AND BIDIRECTIONAL PARAMETERS
# IGNORE SEQUENCES THAT ARE <10% IDENTICAL TO HXB2 REFERENCE STRAIN
# THIS IS A TIME SAVING STEP WHICH ALSO PREVENTS SEQUENCES WITH NO MEANINGFUL OVERLAP FROM BREAKING TN93
# PERFORM TN93 DISTANCE COMPARISON TO HXB2 REFERENCE FILES

# IF SEQS OVERLAP REFERENCE BY 250 NT, HAVE >10% IDENTIY, AND ARE WITHIN 20% TN93 DISTANCE AFTER ALIGNMENT
# THEN THESE SEQUENCES ARE CONSIDERED HOMOLOGOUS TO THE REFERENCE SEQUENCE

HXB2_ref = {'PR':['>HXB2_K03455_PR','CCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTT'],
'RT':['>HXB2_K03455_RT','CCCATTAGCCCTATTGAGACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAAAGTTAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTACAGAGATGGAAAAGGAAGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAAGACTTCTGGGAAGTTCAATTAGGAATACCACATCCCGCAGGGTTAAAAAAGAAAAAATCAGTAACAGTACTGGATGTGGGTGATGCATATTTTTCAGTTCCCTTAGATGAAGACTTCAGGAAGTATACTGCATTTACCATACCTAGTATAAACAATGAGACACCAGGGATTAGATATCAGTACAATGTGCTTCCACAGGGATGGAAAGGATCACCAGCAATATTCCAAAGTAGCATGACAAAAATCTTAGAGCCTTTTAGAAAACAAAATCCAGACATAGTTATCTATCAATACATGGATGATTTGTATGTAGGATCTGACTTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAGCTGAGACAACATCTGTTGAGGTGGGGACTTACCACACCAGACAAAAAACATCAGAAAGAACCTCCATTCCTTTGGATGGGTTATGAACTCCATCCTGATAAATGGACAGTACAGCCTATAGTGCTGCCAGAAAAAGACAGCTGGACTGTCAATGACATACAGAAGTTAGTGGGGAAATTGAATTGGGCAAGTCAGATTTACCCAGGGATTAAAGTAAGGCAATTATGTAAACTCCTTAGAGGAACCAAAGCACTAACAGAAGTAATACCACTAACAGAAGAAGCAGAGCTAGAACTGGCAGAAAACAGAGAGATTCTAAAAGAACCAGTACATGGAGTGTATTATGACCCATCAAAAGACTTAATAGCAGAAATACAGAAGCAGGGGCAAGGCCAATGGACATATCAAATTTATCAAGAGCCATTTAAAAATCTGAAAACAGGAAAATATGCAAGAATGAGGGGTGCCCACACTAATGATGTAAAACAATTAACAGAGGCAGTGCAAAAAATAACCACAGAAAGCATAGTAATATGGGGAAAGACTCCTAAATTTAAACTGCCCATACAAAAGGAAACATGGGAAACA'],
'INT':['>HXB2_K03455_INT','TTTTTAGATGGAATAGATAAGGCCCAAGATGAACATGAGAAATATCACAGTAATTGGAGAGCAATGGCTAGTGATTTTAACCTGCCACCTGTAGTAGCAAAAGAAATAGTAGCCAGCTGTGATAAATGTCAGCTAAAAGGAGAAGCCATGCATGGACAAGTAGACTGTAGTCCAGGAATATGGCAACTAGATTGTACACATTTAGAAGGAAAAGTTATCCTGGTAGCAGTTCATGTAGCCAGTGGATATATAGAAGCAGAAGTTATTCCAGCAGAAACAGGGCAGGAAACAGCATATTTTCTTTTAAAATTAGCAGGAAGATGGCCAGTAAAAACAATACATACTGACAATGGCAGCAATTTCACCGGTGCTACGGTTAGGGCCGCCTGTTGGTGGGCGGGAATCAAGCAGGAATTTGGAATTCCCTACAATCCCCAAAGTCAAGGAGTAGTAGAATCTATGAATAAAGAATTAAAGAAAATTATAGGACAGGTAAGAGATCAGGCTGAACATCTTAAGACAGCAGTACAAATGGCAGTATTCATCCACAATTTTAAAAGAAAAGGGGGGATTGGGGGGTACAGTGCAGGGGAAAGAATAGTAGACATAATAGCAACAGACATACAAACTAAAGAATTACAAAAACAAATTACAAAAATTCAAAATTTTCGGGTTTATTACAGGGACAGCAGAAATCCACTTTGGAAAGGACCAGCAAAGCTCCTCTGGAAAGGTGAAGGGGCAGTAGTAATACAAGATAATAGTGACATAAAAGTAGTGCCAAGAAGAAAAGCAAAGATCATTAGGGATTATGGAAAACAGATGGCAGGTGATGATTGTGTGGCAAGTAGACAGGATGAGGAT']}

for ref in HXB2_ref:
	if skip != 1:
		print('PERFORMING ALIGNMENT AND TN93 DISTANCE CALCULATIONS TO HXB2 REFERENCE FOR '+str(ref)+' REGION\n')
		outfile = open(str(in_file.split(f)[0])+'.'+str(ref)+'_HXB2.fa','w')
		outfile.write(str(HXB2_ref[ref][0])+'\n'+str(HXB2_ref[ref][1]))
		outfile.close()
		os.system('bealign "'+str(in_file.split(f)[0])+'.NoUnkChar.fa" "'+str(in_file.split(f)[0])+'.NoUnkChar.'+str(ref)+'_HXB2.Unidirectional.bam" -r "'+str(in_file.split(f)[0])+'.'+str(ref)+'_HXB2.fa" -m HIV_BETWEEN_F')
		time.sleep(1)
		os.system('bealign "'+str(in_file.split(f)[0])+'.NoUnkChar.fa" "'+str(in_file.split(f)[0])+'.NoUnkChar.'+str(ref)+'_HXB2.Bidirectional.bam" -r "'+str(in_file.split(f)[0])+'.'+str(ref)+'_HXB2.fa" -m HIV_BETWEEN_F -R')
		time.sleep(1)
		os.system('bam2msa "'+str(in_file.split(f)[0])+'.NoUnkChar.'+str(ref)+'_HXB2.Unidirectional.bam" "'+str(in_file.split(f)[0])+'.NoUnkChar.'+str(ref)+'_HXB2.Unidirectional.fa"')
		time.sleep(1)
		os.system('bam2msa "'+str(in_file.split(f)[0])+'.NoUnkChar.'+str(ref)+'_HXB2.Bidirectional.bam" "'+str(in_file.split(f)[0])+'.NoUnkChar.'+str(ref)+'_HXB2.Bidirectional.fa"')
		time.sleep(1)
		os.system('rm *bam*')
		# ENSURE EACH SEQUENCE HAS 25% NON-GAPPED CHARACTERS
		# DEALS WITH BUG/FEATURE of tn93 THAT DIES IF NOT SUFFICIENT DATA FOR COMPARISON
		print('\nREPLACING '+str(ref)+' SEQUENCE FILES WITH NEW SEQUENCES, EXCLUDING SEQUENCES WITH MORE THAN 75% GAPS\n')
		print('\nREMOVING '+str(ref)+' SEQUENCE FILES WITH ALIGNMENT ERRORS\n')
		for d in ['Uni','Bi']:
			for line in open(str(in_file.split(f)[0])+'.NoUnkChar.'+str(ref)+'_HXB2.'+str(d)+'directional.fa'):
				if '>' in line:
					seq_id = line.strip('>').rstrip()
					seq_dict[seq_id]['seq'] = ''
				else: seq_dict[seq_id]['seq'] = seq_dict[seq_id]['seq'] + line.rstrip()
			outfile = open(str(in_file.split(f)[0])+'.NoUnkChar.'+str(ref)+'_HXB2.'+str(d)+'directional.fa','w')
			for seq_id in seq_dict:
				if seq_dict[seq_id]['UnknownCharacters'] == 0 and 'seq' in seq_dict[seq_id]:
					if float(seq_dict[seq_id]['seq'].count('-'))/len(seq_dict[seq_id]['seq']) <= 0.75 and len(seq_dict[seq_id]['seq']) == len(HXB2_ref[ref][1]):
						outfile.write('>'+str(seq_id)+'\n'+str(seq_dict[seq_id]['seq'])+'\n')
			outfile.close()
		os.system('tn93 -t 1.0 -f csv -l 250 -o "'+str(in_file.split(f)[0])+'.NoUnkChar.'+str(ref)+'_HXB2.Unidirectional.csv" -s "'+str(in_file.split(f)[0])+'.'+str(ref)+'_HXB2.fa" "'+str(in_file.split(f)[0])+'.NoUnkChar.'+str(ref)+'_HXB2.Unidirectional.fa"')
		time.sleep(1)
		os.system('tn93 -t 1.0 -f csv -l 250 -o "'+str(in_file.split(f)[0])+'.NoUnkChar.'+str(ref)+'_HXB2.Bidirectional.csv" -s "'+str(in_file.split(f)[0])+'.'+str(ref)+'_HXB2.fa" "'+str(in_file.split(f)[0])+'.NoUnkChar.'+str(ref)+'_HXB2.Bidirectional.fa"')
		time.sleep(1)
		print('\n')
	else: print('SKIPPING ALIGNMENT AND TN93 DISTANCE CALCULATIONS FOR PR, RT, AND INT REGIONS\n')
	for d in ['Uni','Bi']:
		for line in open(str(in_file.split(f)[0])+'.NoUnkChar.'+str(ref)+'_HXB2.'+str(d)+'directional.csv','r'):
			if 'Distance' not in line:
				if len(line.split(',')) == 3:
					if float(line.split(',')[2].rstrip()) < 0.2: seq_dict[line.split(',')[0]][str(ref)+'_'+(d)] = 1
				# ADD REPORT THAT SPECIFIC SEQUENCE THAT THE ANALYSIS IS DYING ON 
				else: print('\n\t\t***ERROR***\n\n\tIMPROPERLY FORMATED TN93 DISTANCE FILE\n\nCHECK THE LINE CONTAINING: '+str(line)+'\n\nSKIPPING THE '+str(d)+'directional '+str(ref)+' ANALYSIS\
				\n\nTHIS TN93 DISTANCE ANALYSIS SHOULD BE RERUN MANUALLY USING CODE: tn93 -t 1.0 -f csv -l 250 -o "'+str(in_file.split(f)[0])+'.NoUnkChar.'+str(ref)+'_HXB2.'+str(d)+'directional.csv" -s "'+str(in_file.split(f)[0])+'.'+str(ref)+'_HXB2.fa" "'+str(in_file.split(f)[0])+'.NoUnkChar.'+str(ref)+'_HXB2.'+str(d)+'directional.fa"\n\n\
				DescribeSeq.py CAN BE THEN BE RE-RUN WITH THE -s 1 FLAG\n\n')
		
# CHECK FOR LARGE GAP IN SEQUENCES WITH BOTH PR AND RT WHEN ALIGNED TO PR/RT
# A LARGE GAP (>250 -) APPEARS IN THE FIRST 300 NT (PR) WHEN RT IS BEFORE PR IN SEQUENCE
# ALSO GENERATE BIDIRECTIONAL PRRT FOR HXB2 SCREENING (ADDED BY JOW ON 8/2/2018)
if skip != 1:
	print('CHECKING FOR INVERTED SEQUENCES IN PR/RT REGION\n')
	outfile = open(str(in_file.split(f)[0])+'.PRRT_HXB2.fa','w')
	outfile.write('>HXB2_K03455_PRRT\nCCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTTCCCATTAGCCCTATTGAGACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAAAGTTAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTACAGAGATGGAAAAGGAAGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAAGACTTCTGGGAAGTTCAATTAGGAATACCACATCCCGCAGGGTTAAAAAAGAAAAAATCAGTAACAGTACTGGATGTGGGTGATGCATATTTTTCAGTTCCCTTAGATGAAGACTTCAGGAAGTATACTGCATTTACCATACCTAGTATAAACAATGAGACACCAGGGATTAGATATCAGTACAATGTGCTTCCACAGGGATGGAAAGGATCACCAGCAATATTCCAAAGTAGCATGACAAAAATCTTAGAGCCTTTTAGAAAACAAAATCCAGACATAGTTATCTATCAATACATGGATGATTTGTATGTAGGATCTGACTTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAGCTGAGACAACATCTGTTGAGGTGGGGACTTACCACACCAGACAAAAAACATCAGAAAGAACCTCCATTCCTTTGGATGGGTTATGAACTCCATCCTGATAAATGGACAGTACAGCCTATAGTGCTGCCAGAAAAAGACAGCTGGACTGTCAATGACATACAGAAGTTAGTGGGGAAATTGAATTGGGCAAGTCAGATTTACCCAGGGATTAAAGTAAGGCAATTATGTAAACTCCTTAGAGGAACCAAAGCACTAACAGAAGTAATACCACTAACAGAAGAAGCAGAGCTAGAACTGGCAGAAAACAGAGAGATTCTAAAAGAACCAGTACATGGAGTGTATTATGACCCATCAAAAGACTTAATAGCAGAAATACAGAAGCAGGGGCAAGGCCAATGGACATATCAAATTTATCAAGAGCCATTTAAAAATCTGAAAACAGGAAAATATGCAAGAATGAGGGGTGCCCACACTAATGATGTAAAACAATTAACAGAGGCAGTGCAAAAAATAACCACAGAAAGCATAGTAATATGGGGAAAGACTCCTAAATTTAAACTGCCCATACAAAAGGAAACATGGGAAACA\n')
	outfile.close()
	os.system('bealign "'+str(in_file.split(f)[0])+'.NoUnkChar.fa" "'+str(in_file.split(f)[0])+'.NoUnkChar.PRRT_HXB2.Unidirectional.bam" -r "'+str(in_file.split(f)[0])+'.PRRT_HXB2.fa" -m HIV_BETWEEN_F')
	time.sleep(1)
	os.system('bam2msa "'+str(in_file.split(f)[0])+'.NoUnkChar.PRRT_HXB2.Unidirectional.bam" "'+str(in_file.split(f)[0])+'.NoUnkChar.PRRT_HXB2.Unidirectional.fa"')
	time.sleep(1)
	os.system('bealign "'+str(in_file.split(f)[0])+'.NoUnkChar.fa" "'+str(in_file.split(f)[0])+'.NoUnkChar.PRRT_HXB2.Bidirectional.bam" -r "'+str(in_file.split(f)[0])+'.PRRT_HXB2.fa" -m HIV_BETWEEN_F -R')
	time.sleep(1)
	os.system('bam2msa "'+str(in_file.split(f)[0])+'.NoUnkChar.PRRT_HXB2.Bidirectional.bam" "'+str(in_file.split(f)[0])+'.NoUnkChar.PRRT_HXB2.Bidirectional.fa"')
	time.sleep(1)
	os.system('rm *bam*')

prrt_dict = {}
for line in open(str(in_file.split(f)[0])+'.NoUnkChar.PRRT_HXB2.Unidirectional.fa','r'):
	if '>' in line:
		seq_id = line.strip('>').rstrip()
			[seq_id] = ''
	else: prrt_dict[seq_id] = prrt_dict[seq_id] + line.rstrip()

for seq_id in prrt_dict:
	if len(prrt_dict[seq_id]) >= 300:
		if '-'*250 in prrt_dict[seq_id][0:300]: seq_dict[seq_id]['RTPR'] = 1
	else: print('\n\t\t***ERROR***\n\n\tIMPROPERLY FORMATED PRRT ALIGNMENT FILE\n\n')

# CHECK FOR HXB2 LINKS IN PRRT REGION
# CHECK BIDIRECTIONAL PRRT FOR HXB2 SCREENING (ADDED BY JOW ON 8/2/2018)
print('CHECKING FOR HXB2 LINKED SEQUENCES IN PR/RT REGION\n')
os.system('tn93 -t 0.015 -a resolve -g 1.0 -f csv -l 250 -o "'+str(in_file.split(f)[0])+'.NoUnkChar.PRRT_HXB2.Bidirectional.csv" -s "'+str(in_file.split(f)[0])+'.PRRT_HXB2.fa" "'+str(in_file.split(f)[0])+'.NoUnkChar.PRRT_HXB2.Bidirectional.fa"')
for line in open(str(in_file.split(f)[0])+'.NoUnkChar.PRRT_HXB2.Bidirectional.csv','r'):
	if 'Distance' not in line:
		seq_dict[line.split(',')[0]]['PRRT_HXB2_Link'] = 1

# WRITE OUTPUT SUMMARY TO FILE

print('WRITING OUTPUT SUMMARY FILE TO:\t' + str(in_file.split(f)[0])+'.SequenceSummary.tab\n')

outfile = open(str(in_file.split(f)[0])+'.SequenceSummary.tab','w')
outfile.write('SeqID\tUnknownCharacters\tNString\tAmbig5%\tPRRT_HXB2_Link\tPR_Forward\tPR_Reverse\tRT_Forward\tRT_Reverse\tINT_Forward\tINT_Reverse\tRTPR_Swap\n')

for seq_id in seq_dict:
	outfile.write(str(seq_id)+'\t')
	if seq_dict[seq_id]['UnknownCharacters'] == 1: outfile.write('1\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\n')
	else:
		if seq_dict[seq_id]['Nstring'] == 1: outfile.write('0\t1\t')
		else: outfile.write('0\t0\t')
		if seq_dict[seq_id]['Ambig5%'] == 1: outfile.write('1\t')
		else: outfile.write('0\t')
		if seq_dict[seq_id]['PRRT_HXB2_Link'] == 1: outfile.write('1\t')
		else: outfile.write('0\t')
		if seq_dict[seq_id]['PR_Uni'] == 1 and seq_dict[seq_id]['PR_Bi'] == 1: outfile.write('1\t0\t')
		elif seq_dict[seq_id]['PR_Bi'] == 1: outfile.write('0\t1\t')
		else: outfile.write('0\t0\t')
		if seq_dict[seq_id]['RT_Uni'] == 1 and seq_dict[seq_id]['RT_Bi'] == 1: outfile.write('1\t0\t')
		elif seq_dict[seq_id]['RT_Bi'] == 1: outfile.write('0\t1\t')
		else: outfile.write('0\t0\t')
		if seq_dict[seq_id]['INT_Uni'] == 1 and seq_dict[seq_id]['INT_Bi'] == 1: outfile.write('1\t0\t')
		elif seq_dict[seq_id]['INT_Bi'] == 1: outfile.write('0\t1\t')
		else: outfile.write('0\t0\t')
		if seq_dict[seq_id]['RTPR'] == 1 and seq_dict[seq_id]['PR_Uni'] == 1 and seq_dict[seq_id]['RT_Uni'] == 1: outfile.write('1\n')
		else: outfile.write('0\n')
			
outfile.close()