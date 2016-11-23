
import sys, os, re, math, csv, pickle
from sklearn.metrics.cluster import normalized_mutual_info_score

percentKeep = 0.25

# MUTUAL INFORMATION CALCULATION
def calcMI(tp, fp, fn, tn):
    firstTerm = [1 for x in range(tp)]+[0 for x in range(fp)]+[1 for x in range(fn)]+[0 for x in range(tn)]
    secondTerm = [1 for x in range(tp + fp)]+[0 for x in range(fn + tn)]
    result = normalized_mutual_info_score(firstTerm, secondTerm)
    return result

##### START: LOADING AND ACCESSING DATA FUNCTIONS #######
otuTable = []

# Which data do I want to load:
if True:
    with open('otu_table_normCSS.txt', 'r') as fileToGetFeatures:
        mastersToTrain = []
        allData = csv.reader(fileToGetFeatures, delimiter='\t')
        for lineOfData in allData:
            otuTable.append(lineOfData)
else:
    with open('unchangedOTUData.txt', 'r') as fileToGetFeatures:
        mastersToTrain = []
        allData = csv.reader(fileToGetFeatures, delimiter=' ')
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
    return float(OTUDataOfThisSampleID[indexOfSampleID])

def getAllOTUs(otuTable):
    return [otudata[0] for otudata in otuTable[1:]]

def getAllSampleIDs(allRows):
    indexOfSampleID = allRows[0].index('#SampleID')
    return [sampleData[indexOfSampleID] for sampleData in allRows[1:]]

def getSiblingSampleID(allRows, sampleID):
    indexOfSampleID = allRows[0].index('#SampleID')
    indexOfPairID = allRows[0].index('Pair')
    siblingPair = [sampleData[indexOfPairID] for sampleData in allRows[1:] if sampleData[indexOfSampleID] == sampleID][0]
    sampleIDsWithPair = [sampleData[indexOfSampleID] for sampleData in allRows[1:] if sampleData[indexOfPairID] == siblingPair]
    return sampleIDsWithPair

def getAllSampleIDsFromBatch(allRows, batchNumber):
    indexOfBatchNum = allRows[0].index('batch')
    indexOfSampleID = allRows[0].index('#SampleID')
    return [sampleData[indexOfSampleID] for sampleData in allRows[1:] if sampleData[indexOfBatchNum] == str(batchNumber)]    


def returnBatch(allRows):
    numberOfBatches = 5
    eachBatch = []
    for i in range(0,numberOfBatches):
        eachBatch.append(getAllSampleIDsFromBatch(allRows, i))
    return eachBatch

def returnSiblingPairs(allRows):
    allSamplesLookedAt = []
    result = []
    for specificSample in allSampleIDs:
        siblings = getSiblingSampleID(allRows, specificSample)
        appendThisPair = True
        for sibling in siblings:
            if sibling in allSamplesLookedAt:
                appendThisPair = False
                continue
        if appendThisPair:
            allSamplesLookedAt.extend(siblings)
            result.append(siblings)
    return result

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
# dropOutOnlyOne = True

def findMI(excludeSamples=None):
    MIMapping = {}
    allSamplesReferencedBuildingMI = []

    for specificOTU in allOTUs:
        # print "Cacluating for OTU %s" % specificOTU
        countAutisticAndOTU = 0
        countAutisticAndNotOTU = 0
        countNotAutisticAndOTU = 0
        countNotAutisticAndNotOTU = 0
        AutisticAndOTU = []
        AutisticAndNotOTU = []
        NotAutisticAndOTU = []
        NotAutisticAndNotOTU = []    
        for specificSample in allSampleIDs:
            if specificSample in dontDoTheseSamples or specificSample in excludeSamples:
                # Data missing or purposely being held out.
                continue
            autistic = isSampleAutistic(allRows, specificSample)
            OTUValue = getValueOfOTU(otuTable, specificOTU, specificSample)
            allSamplesReferencedBuildingMI.append(specificSample)

            # These we can play with, i'm just thresholding otu value.
            if autistic == 'Aut' and float(OTUValue) > 0.0:
                countAutisticAndOTU += 1
                AutisticAndOTU.append(specificSample)
            elif autistic == 'Control' and float(OTUValue) > 0.0:
                countNotAutisticAndOTU += 1
                NotAutisticAndOTU.append(specificSample)
            elif autistic == 'Aut':# and float(OTUValue) == 0.0:
                countAutisticAndNotOTU += 1
                AutisticAndNotOTU.append(specificSample)
            elif autistic == 'Control':# and float(OTUValue) == 0.0:
                countNotAutisticAndNotOTU += 1
                NotAutisticAndNotOTU.append(specificSample)
            else:
                print autistic
                print OTUValue
                print "ERROR"
                sys.exit()

        MI = calcMI(countAutisticAndOTU, countNotAutisticAndOTU, countAutisticAndNotOTU, countNotAutisticAndNotOTU)
        metricsForCalculation = "(isAutistic & hasOTU %d, isAutistic' & hasOTU %d, isAutistic & hasOTU' %d, isAutistic' & hasOTU' %d)" % (countAutisticAndOTU, countNotAutisticAndOTU, countAutisticAndNotOTU, countNotAutisticAndNotOTU)
        MIMapping[specificOTU] = {'MutualInformation': MI, 'Metrics': metricsForCalculation, 'ExtremelyVerbose':{'Aut_OTU':AutisticAndOTU, 'AutN_OTU':NotAutisticAndOTU, 'Aut_OTUN':AutisticAndNotOTU, 'AutN_OTUN':NotAutisticAndNotOTU}}
        # print MIMapping[specificOTU]

    print sorted([v['MutualInformation'] for v in MIMapping.values()])
    # print "As you can see, very little MI individually for any of these."

    if True:
        MIMappingList = [[otu, data['MutualInformation'], data['Metrics']] for otu, data in MIMapping.items()]
        sortedMIsForOTUs = sorted(MIMappingList, key=lambda x: x[1], reverse=True)
        excludeOTU = []
        for idx, sortedMI in enumerate(sortedMIsForOTUs):
            #print "OTU: %s, MI: %s, Calculated from: %s" % (sortedMI[0], sortedMI[1], sortedMI[2])
            if idx > len(sortedMIsForOTUs) * percentKeep:
                excludeOTU.append(sortedMI[0])
            pass

    #pickle.dump(sortedMIsForOTUs, open( "MI.p", "wb" ))

    return MIMapping, list(set(allSamplesReferencedBuildingMI)), excludeOTU


def runFullTestWithHoldout():
    rightPair = []
    wrongPair = []
    correctCount = 0
    incorrectCount = 0
    resultsList = []

    correctSamples = []
    wrontSamples = []

    runBatch = True
    if not runBatch:
        siblingPairs = returnSiblingPairs(allRows)
    else:
        siblingPairs = returnBatch(allRows)
    print "Batch Pairs:"
    print siblingPairs

    for batchIdx, siblingPair in enumerate(siblingPairs):

        if len(set(siblingPair) & set(dontDoTheseSamples)) != 0:
            print "Removing samples that didn't pass QA" 
            siblingPair = [i for i in siblingPair if i not in dontDoTheseSamples]
        if len(siblingPair) > 2:
            print "Batch size greater than two."
            #continue
        if len(siblingPair) < 2:
            print "This sample has no sibling, batch size one. %s" % siblingPair
            #continue

        MIMapping, allSamplesReferencedBuildingMI, excludeOTU = findMI(siblingPairs)

        #### TODO ### check if siblingPair in allSamplesReferencedBuildingMI. 
        for sibling in siblingPairs:
            if sibling in allSamplesReferencedBuildingMI:
                print "ERRRRRORORRRRR"
                sys.exit()

        autPair = []
        for sibling in siblingPair:
            autistic = isSampleAutistic(allRows, sibling)
            autPair.append(autistic)


        totalsAutistic = []
        totalsControl = []
        trueValues = []

        for idx, sampleToTest in enumerate(siblingPair):
            print "Predicting for sample: %s" % sampleToTest

            totalForAutistic = 0
            totalForControl = 0
            verboseInfoControl = []
            verboseInfoAut = []
            allTrainSamplesUsedUponComparisonTime = []

            countAut = 0
            countCont = 0
            for trainSample in allSampleIDs:
                if trainSample in dontDoTheseSamples or trainSample in siblingPair:
                    continue    
                autistic = isSampleAutistic(allRows, trainSample)
                if autistic == 'Aut':
                    countAut += 1
                else:
                    countCont += 1

            for specificOTU in allOTUs:
                OTUValue = getValueOfOTU(otuTable, specificOTU, sampleToTest)
                if OTUValue != 0 and specificOTU not in excludeOTU:
                    autCount = {'count': 0, 'otuCount': 0}
                    controlCount = {'count': 0, 'otuCount': 0}
                    for trainSample in allSampleIDs:
                        if trainSample in dontDoTheseSamples or trainSample in siblingPair:
                            continue    
                        else:
                            if trainSample not in allTrainSamplesUsedUponComparisonTime:
                                allTrainSamplesUsedUponComparisonTime.append(trainSample)

                        autistic = isSampleAutistic(allRows, trainSample)
                        OTUTrainValue = getValueOfOTU(otuTable, specificOTU, trainSample)

                        # So as to calculate the fractions.
                        if autistic == 'Aut':
                            autCount['count'] += 1
                        if autistic == 'Control':
                            controlCount['count'] += 1
                        if autistic == 'Aut' and float(OTUTrainValue) > 0.0:
                            autCount['otuCount'] += 1
                        if autistic == 'Control' and float(OTUTrainValue) > 0.0:
                            controlCount['otuCount'] += 1
                
                    FiAut = math.log(float(autCount['otuCount'] + 1)/float(autCount['count'] + len(allOTUs)))
                    FiControl = math.log(float(controlCount['otuCount'] + 1)/float(controlCount['count']+ len(allOTUs))) 
                    totalForAutistic += FiAut
                    totalForControl += FiControl


            totalForAutistic += math.log(float(countAut)/float(countAut + countCont))
            totalForControl += math.log(float(countCont)/float(countAut + countCont))
            if totalForAutistic > totalForControl:
                prediction = 'Aut'
            elif totalForAutistic < totalForControl:
                prediction = 'Control'

            if not runBatch:
                resultsList.append({'sampleToTest': sampleToTest, 'Aut': totalForAutistic, 'Control': totalForControl, 'trueValue': autPair[idx]})
            else:
                resultsList.append({'sampleToTest': sampleToTest, 'Aut': totalForAutistic, 'Control': totalForControl, 'trueValue': autPair[idx], 'batch':batchIdx})

            # Holdout verification, REALLY double check the OTU was indeed held out:
            print ""
            print "allTrainSamplesUsedUponComparisonTime:"
            print allTrainSamplesUsedUponComparisonTime
            if sampleToTest in allTrainSamplesUsedUponComparisonTime:
                print "ERROR: was used in tally!"
                sys.exit()
            print "It has been verified, the sampleToTest was not used in run-time comparison."

            ### PRINTING RESULTS:
            print ""
            if prediction == 'Unsure':
                print "UNSURE"
            elif autPair[idx] == prediction:
                correctCount += 1
                print "CCCCCOOOORRRRECT"
                correctSamples.append(sampleToTest)
            else:
                incorrectCount += 1
                print "INCORRECT"
                wrontSamples.append(sampleToTest)

            print resultsList
            print "Correct:"
            print correctCount
            print "Incorrect:"
            print incorrectCount
            print "SAMPLES"
            print correctSamples
            print wrontSamples

            totalsAutistic.append(totalForAutistic)
            totalsControl.append(totalForControl)
            trueValues.append(autPair[idx])

        if len(siblingPair) == 2:

            print "Comparing the siblings:"
            print totalsAutistic
            print totalsControl
            print trueValues
            if totalsAutistic[0] > totalsAutistic[1] and trueValues[0] == 'Aut':
                print 'GOT right on aut'
            if totalsControl[0] > totalsControl[1] and trueValues[0] == 'Control':
                print 'GOT right on cont'   

            if totalsAutistic[0] + totalsControl[1] > totalsAutistic[1] + totalsControl[0] and autPair[0] == 'Aut':
                print "RIGHTPAIR"
                rightPair.append({'sampleToTest': siblingPair, 'Aut-Cont': totalsAutistic[0] + totalsControl[1], 'Cont-Aut': totalsAutistic[1] + totalsControl[0], 'trueValue': autPair})
            else:
                print "WRONGPAIR"
                wrongPair.append({'sampleToTest': siblingPair, 'Aut-Cont': totalsAutistic[0] + totalsControl[1], 'Cont-Aut': totalsAutistic[1] + totalsControl[0], 'trueValue': autPair})

            print "rightPair:"
            print rightPair
            print "wrongPair:"
            print wrongPair

    print "resultsList"
    print resultsList
    print "Correct Count:"
    print correctCount
    print "Incorrect Count:"
    print incorrectCount

    print "SAMPLES"
    print correctSamples
    print wrontSamples
    print "percentKeep:"
    print percentKeep


# findMI([])
# sys.exit()

runFullTestWithHoldout()
# print getAllSampleIDsFromBatch(allRows, 0)
# runTestWithBatchHoldout(0)


# Try only doing using a single batch.
# dontDoTheseSamples += getAllSampleIDsFromBatch(allRows, 1) + getAllSampleIDsFromBatch(allRows, 2) + getAllSampleIDsFromBatch(allRows, 3) + getAllSampleIDsFromBatch(allRows, 4)
# runFullTestWithHoldout()










