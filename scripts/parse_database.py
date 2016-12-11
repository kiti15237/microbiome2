import csv
import re

rows = []
f= open("C:/Users/Christine/Documents/Lab/16S/data_participantresponse.txt") 
s = f.read()
split = str.split(s, "),")

split_mod = []
for line in split:
	split_mod.append(re.sub("\(", "", line))

file = open("C:/Users/Christine/Documents/Lab/16S/data_participantresponse_split.txt", 'w')
for line in split_mod:
	file.write(line + '\n')
