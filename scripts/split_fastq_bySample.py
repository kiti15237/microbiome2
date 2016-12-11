import re
import getopt
import sys

#f = open("C:\Users\ctataru\Documents/Lab/data/raw_fastqs/test.fastq")
f = open(sys.argv[1])
i = 0
headers = []
seq = []
plus = []
qual = []
sids = []
lines_by_sample = {}
sid = ""
for line in f.readlines() :
	if i % 4 == 0:
		ind_dash = line.index('_')
		space = line.index(' ')
		sid = line[1:ind_dash]
		newLine = '@' + line[space + 1:]
		headers.append(newLine)
		if sid in lines_by_sample.keys() :
			lines_by_sample[sid] += (newLine)
		else:
			lines_by_sample[sid] = newLine
	else:
		lines_by_sample[sid] += line
	i += 1

for key in lines_by_sample.keys():
	if re.search('R1', sys.argv[1]):
		out = open(key + "_R1.fastq", 'w')
	elif re.search("R2", sys.argv[1]):
		out = open(key + "_R2.fastq", 'w')
	out.write(lines_by_sample[key])

