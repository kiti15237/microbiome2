
import sys, os, re, math, csv
from sklearn.metrics.cluster import normalized_mutual_info_score


# MUTUAL INFORMATION CALCULATION
def calcMI(ab, aNotb, abNot, aNotbNot):
    MIFirstTerm = [1 for x in range(ab)] + [0 for x in range(aNotb)] + [1 for x in range(abNot)] + [0 for x in range(aNotbNot)]
    MISecondTerm = [1 for x in range(ab + aNotb)] + [0 for x in range(abNot + aNotbNot)]
    result = normalized_mutual_info_score(MIFirstTerm, MISecondTerm)
    return result


##### START: LOADING AND ACCESSING DATA FUNCTIONS #######
otuTable = []
with open('otu_table_normCSS.txt', 'r') as fileToGetFeatures:
    mastersToTrain = []
    allData = csv.reader(fileToGetFeatures, delimiter='\t')
    for lineOfData in allData:
        otuTable.append(lineOfData)

# Got this straight from the internet.
def xlsx(fname):
    import zipfile
    from xml.etree.ElementTree import iterparse
    z = zipfile.ZipFile(fname)
    strings = [el.text for e, el in iterparse(z.open('xl/sharedStrings.xml')) if el.tag.endswith('}t')]
    rows = []
    row = []
    value = ''
    for e, el in iterparse(z.open('xl/worksheets/sheet1.xml')):
        if el.tag.endswith('}v'): # <v>84</v>
            value = el.text
        if el.tag.endswith('}c'): # <c r="A3" t="s"><v>84</v></c>
            if el.attrib.get('t') == 's':
                value = strings[int(value)]
            letter = el.attrib['r'] # AZ22
            while letter[-1].isdigit():
                letter = letter[:-1]
            row.append(value)
            value = ''
        if el.tag.endswith('}row'):
            rows.append(row)
            row = []
    return rows

def isSampleAutistic(allRows, sampleToRetrieve):
	indexOfTreatment = allRows[0].index('Treatment ')
	indexOfSampleID = allRows[0].index('#SampleID')
	dataOfThisSampleID = [row for row in allRows if row[indexOfSampleID] == str(sampleToRetrieve)][0]
	return dataOfThisSampleID[indexOfTreatment]

def getValueOfOTU(otuTable, otuIDToRetrieve, sampleToRetrieve):
	indexOfSampleID = otuTable[0].index(sampleToRetrieve)
	indexOfOTUID = otuTable[0].index('OTUID')
	OTUDataOfThisSampleID = [row for row in otuTable if row[indexOfOTUID] == str(otuIDToRetrieve)][0]
	#print OTUDataOfThisSampleID[indexOfSampleID]
	return OTUDataOfThisSampleID[indexOfSampleID]

def getAllOTUs(otuTable):
	return [otudata[0] for otudata in otuTable[1:]]

def getAllSampleIDs(allRows):
	indexOfSampleID = allRows[0].index('#SampleID')
	return [sampleData[indexOfSampleID] for sampleData in allRows[1:]]

allRows = xlsx('mapping_merged_oct.xlsx')
#print otuTable

# sampleToRetrieve = '51' # ---> Aut
# sampleToRetrieve = '173.saliva' # --> Control
# print isSampleAutistic(allRows, sampleToRetrieve = '173.saliva')
# print getValueOfOTU(otuTable, otuIDToRetrieve='891031', sampleToRetrieve='175')

allOTUs = getAllOTUs(otuTable)
allSampleIDs = getAllSampleIDs(allRows)

##### END: LOADING AND ACCESSING DATA FUNCTIONS #######

dontDoTheseSamples = ['384', '172.saliva', '168.1', '367', '173', '173.saliva']


MIMapping = {}
for specificOTU in allOTUs:
	print "Cacluating for OTU %s" % specificOTU
	countAutisticAndOTU = 0
	countAutisticAndNotOTU = 0
	countNotAutisticAndOTU = 0
	countNotAutisticAndNotOTU = 0
	for specificSample in allSampleIDs:
		if specificSample in dontDoTheseSamples:
			# Data missing 
			continue
		autistic = isSampleAutistic(allRows, specificSample)
		OTUValue = getValueOfOTU(otuTable, specificOTU, specificSample)

		# These we can play with, i'm just thresholding otu value.
		if autistic == 'Aut' and float(OTUValue) > 0.0:
			countAutisticAndOTU += 1
		elif autistic == 'Control' and float(OTUValue) > 0.0:
			countNotAutisticAndOTU += 1
		elif autistic == 'Aut' and float(OTUValue) == 0.0:
			countAutisticAndNotOTU += 1
		elif autistic == 'Control' and float(OTUValue) == 0.0:
			countNotAutisticAndNotOTU += 1
		else:
			print autistic
			print OTUValue
			print "ERROR"
			sys.exit()

	MI = calcMI(countAutisticAndOTU, countNotAutisticAndOTU, countAutisticAndNotOTU, countNotAutisticAndNotOTU)
	metricsForCalculation = "(isAutistic & hasOTU %d, isAutistic' & hasOTU %d, isAutistic & hasOTU' %d, isAutistic' & hasOTU' %d)" % (countAutisticAndOTU, countNotAutisticAndOTU, countAutisticAndNotOTU, countNotAutisticAndNotOTU)
	MIMapping[specificOTU] = {'MutualInformation': MI, 'Metrics': metricsForCalculation}
	print MIMapping[specificOTU]

print sorted([v['MutualInformation'] for v in MIMapping.values()])
print "As you can see, very little MI individually for any of these."

MIMappingList = [[otu, data['MutualInformation']] for otu, data in MIMapping.items()]
sortedMIsForOTUs = sorted(MIMappingList, key=lambda x: x[1], reverse=True)
print sortedMIsForOTUs



