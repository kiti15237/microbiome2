
import sys, os, re, math, csv, pickle
from sklearn.metrics.cluster import normalized_mutual_info_score

# This is an important metric! 3.25 works best for us.
p = 3.25
p = 2
p = 2.5

NormalizeTotals = True


# MUTUAL INFORMATION CALCULATION
def calcMI(tp, fp, fn, tn):
    firstTerm = [1 for x in range(tp)]+[0 for x in range(fp)]+[1 for x in range(fn)]+[0 for x in range(tn)]
    secondTerm = [1 for x in range(tp + fp)]+[0 for x in range(fn + tn)]
    result = normalized_mutual_info_score(firstTerm, secondTerm)
    return result


##### START: LOADING AND ACCESSING DATA FUNCTIONS #######
otuTable = []

# Which data do I want to load:
if False:
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

allRows = xlsx('mapping_merged_oct.xlsx')
#print otuTable

# sampleToRetrieve = '51' # ---> Aut
# sampleToRetrieve = '173.saliva' # --> Control
# print isSampleAutistic(allRows, sampleToRetrieve = '173.saliva')
# print getValueOfOTU(otuTable, otuIDToRetrieve='891031', sampleToRetrieve='175')

allOTUs = getAllOTUs(otuTable)
allSampleIDs = getAllSampleIDs(allRows)

##### END: LOADING AND ACCESSING DATA FUNCTIONS #######


# Newtons method, needed for normalization. To determine root.
def dx(f, x, weightsList, autTotalW):
    return abs(0-f(x, weightsList, autTotalW))
def newtons_method(f, df, x0, e, weightsList, autTotalW):
    delta = dx(f, x0, weightsList, autTotalW)
    while delta > e:
        x0 = x0 - f(x0, weightsList, autTotalW)/df(x0, weightsList)
        print x0
        delta = dx(f, x0, weightsList, autTotalW)
    # print 'Root is at: %s', x0
    # print 'f(x) at root is: %s', f(x0, weightsList, autTotalW)
    return x0
def f(x, weightsList, autTotalW):
    # print "fx:"
    # print "firstPart: %s" % sum([math.pow(w,x) for w in weightsList])
    # print sum([math.pow(w,x) for w in weightsList] + [-autTotalW])
    return sum([math.pow(w,x) for w in weightsList] + [-autTotalW])
def df(x, weightsList):
    # print "dfx"
    # print sum([math.pow(w,x) * math.log(w) for w in weightsList])
    return sum([math.pow(w,x) * math.log(w) for w in weightsList])
# newtons_method(f, df, p, 1e-7, weightsList, autTotalW)
#### Done Newtons Method.



dontDoTheseSamples = ['384', '172.saliva', '168.1', '367', '173', '173.saliva']

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

    if False:
        MIMappingList = [[otu, data['MutualInformation'], data['Metrics']] for otu, data in MIMapping.items()]
        sortedMIsForOTUs = sorted(MIMappingList, key=lambda x: x[1], reverse=True)
        for sortedMI in sortedMIsForOTUs:
            #print "OTU: %s, MI: %s, Calculated from: %s" % (sortedMI[0], sortedMI[1], sortedMI[2])
            pass


    # One possible method of normalizing. Different pControl. Uses newtons method to find.
    autOTUs = {}
    controlOTUs = {}
    for specificOTU in allOTUs:
        OTUValue = getValueOfOTU(otuTable, specificOTU, specificSample)
        autCount = {'count': 0, 'otuCount': 0}
        controlCount = {'count': 0, 'otuCount': 0}
        for trainSample in list(set(allSamplesReferencedBuildingMI)):
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
        fractionAut = float(autCount['otuCount'])/float(autCount['count'])
        fractionControl = float(controlCount['otuCount'])/float(controlCount['count'])    
        if fractionAut > 0 and fractionControl == 0:
            autOTUs[specificOTU] = math.pow(MIMapping[specificOTU]['MutualInformation'], p)
        if fractionControl > 0 and fractionAut == 0:
            controlOTUs[specificOTU] = MIMapping[specificOTU]['MutualInformation']#math.pow(MIMapping[specificOTU]['MutualInformation'], p)
    # print "Aut Total W:"
    # print sum(autOTUs.values())
    # print "Control Total W:"
    # print sum(controlOTUs.values())
    pControl = newtons_method(f, df, p, 1e-12, controlOTUs.values(), sum(autOTUs.values()))
    # Characterize the results.
    onlyAut = list(set(autOTUs.keys()) - set(controlOTUs.keys()))
    onlyControl = list(set(controlOTUs.keys()) - set(autOTUs.keys()))
    onlyAutW = [w for otu, w in autOTUs.items() if otu in onlyAut]
    onlyControlW = [math.pow(w, p) for otu, w in controlOTUs.items() if otu in onlyControl]
    onlyControlWUsingPControl= [math.pow(w, pControl) for otu, w in controlOTUs.items() if otu in onlyControl]
    print "What normalizing would do."
    print "Aut - Control: %s, weighting: %s" % (len(onlyAut), sum(onlyAutW))
    print "Control - Aut: %s, weighting: %s-->%s" % (len(onlyControl), sum(onlyControlW), sum(onlyControlWUsingPControl))

    #pickle.dump(sortedMIsForOTUs, open( "MI.p", "wb" ))

    return MIMapping, list(set(allSamplesReferencedBuildingMI)), pControl


def runFullTestWithHoldout():
    correctCount = 0
    incorrectCount = 0
    resultsList = []
    samplesDone = []
    for oneSample in allSampleIDs:

        # trueValue = isSampleAutistic(allRows, sampleToHoldOut)
        # if trueValue == 'Aut':
        #     continue

        if oneSample in dontDoTheseSamples:
            continue

        selfAndSiblings = getSiblingSampleID(allRows, oneSample)
        # selfAndSiblings = ['51', '52']
        print "Siblings and itself:"
        print selfAndSiblings
        if len(selfAndSiblings) != 2:
            continue


        totalsAutistic = []
        totalsControl = []
        trueValues = []
        for idx, sampleToHoldOut in enumerate(selfAndSiblings):
            print "Predicting for sample: %s" % sampleToHoldOut

            if idx == 0:
                MIMapping, allSamplesReferencedBuildingMI, pControl = findMI(selfAndSiblings)

            totalForAutistic = 0
            totalForControl = 0
            verboseInfoControl = []
            verboseInfoAut = []
            allTrainSamplesUsedUponComparisonTime = []

            for specificOTU in allOTUs:
                OTUValue = getValueOfOTU(otuTable, specificOTU, sampleToHoldOut)
                if OTUValue != 0:
                    autCount = {'count': 0, 'otuCount': 0}
                    controlCount = {'count': 0, 'otuCount': 0}
                    for trainSample in allSampleIDs:
                        if trainSample in dontDoTheseSamples or trainSample in selfAndSiblings:
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
                
                    fractionAut = float(autCount['otuCount'])/float(autCount['count'])
                    fractionControl = float(controlCount['otuCount'])/float(controlCount['count'])    

                    if fractionAut > 0:
                        verboseInfoAut.append({'specificOTU': specificOTU, 'contribution':math.pow(MIMapping[specificOTU]['MutualInformation'], p), 'verboseInfo': MIMapping[specificOTU]['Metrics'], 'ExtremelyVerbose': MIMapping[specificOTU]['ExtremelyVerbose']})
                        totalForAutistic += math.pow(MIMapping[specificOTU]['MutualInformation'], p)
                        print "totalForAutistic incrementing: %s" % math.pow(MIMapping[specificOTU]['MutualInformation'], p)
                    if fractionControl > 0:
                        verboseInfoControl.append({'specificOTU': specificOTU, 'contribution':math.pow(MIMapping[specificOTU]['MutualInformation'], pControl), 'verboseInfo': MIMapping[specificOTU]['Metrics'], 'ExtremelyVerbose': MIMapping[specificOTU]['ExtremelyVerbose']})
                        totalForControl +=  math.pow(MIMapping[specificOTU]['MutualInformation'], pControl)
                        print "totalForControl incrementing: %s" % math.pow(MIMapping[specificOTU]['MutualInformation'], pControl)


                    # totalForAutistic += MIMapping[specificOTU]['MutualInformation']*float(fractionAut)
                    # totalForControl +=  MIMapping[specificOTU]['MutualInformation']*float(fractionControl)

                    # print "fractionAut"
                    # print fractionAut
                    # print "fractionControl"
                    # print fractionControl
                    # print "autCount"
                    # print autCount
                    # print "controlCount"
                    # print controlCount


            if totalForAutistic > totalForControl:
                prediction = 'Aut'
            else:
                prediction = 'Control'

            trueValue = isSampleAutistic(allRows, sampleToHoldOut)
            resultsList.append({'SampleToHoldOUt': sampleToHoldOut, 'Aut': totalForAutistic, 'Control': totalForControl, 'trueValue': trueValue})

            # Holdout verification, REALLY double check the OTU was indeed held out:
            print ""
            print "allTrainSamplesUsedUponComparisonTime:"
            print allTrainSamplesUsedUponComparisonTime
            print "allSamplesReferencedBuildingMI:"
            print allSamplesReferencedBuildingMI
            if sampleToHoldOut in allTrainSamplesUsedUponComparisonTime:
                print "ERROR: was used in tally!"
                sys.exit()
            if sampleToHoldOut in allSamplesReferencedBuildingMI:
                print "ERROR: was used in MI."
            # Just to be super duper sure.
            for vInfo in verboseInfoAut + verboseInfoControl:
                allSamples = vInfo['ExtremelyVerbose']['Aut_OTU'] + vInfo['ExtremelyVerbose']['AutN_OTU'] + vInfo['ExtremelyVerbose']['Aut_OTUN'] + vInfo['ExtremelyVerbose']['AutN_OTUN']
                if sampleToHoldOut in allSamples:
                    print "ERRRRROR: was used in MI's, and still missed earlier."
                    sys.exit()
            print "It has been verified, the sampleToHoldOut was not used in feature-selection or the run-time comparison."

            ### PRINTING RESULTS:
            print ""
            if trueValue == prediction:
                correctCount += 1
                print "CCCCCOOOORRRRECT"
            else:
                incorrectCount += 1
                print "INCORRECT"
            print resultsList
            print "Correct:"
            print correctCount
            print "Incorrect:"
            print incorrectCount
            print "powerOf:"
            print p

            print ""
            print "Justfication:"

            verboseInfoAut = sorted(verboseInfoAut, key=lambda x: x['contribution'], reverse=True)
            verboseInfoControl = sorted(verboseInfoControl, key=lambda x: x['contribution'], reverse=True)

            print "Aut justification,"
            for vInfo in verboseInfoAut:
                if vInfo['contribution'] < 0.00001:
                    continue
                print 'SpecificOTU: %s, Contribution: %s' % (vInfo['specificOTU'], vInfo['contribution'])
                print 'Counts for MI: %s' % vInfo['verboseInfo']
                #print "Verbose Counts for MI: \nAO:%s,A'O:%s,AO':%s,A'O':%s" % (vInfo['ExtremelyVerbose']['Aut_OTU'], vInfo['ExtremelyVerbose']['AutN_OTU'], vInfo['ExtremelyVerbose']['Aut_OTUN'], vInfo['ExtremelyVerbose']['AutN_OTUN']) 
            print "Control justification,"
            for vInfo in verboseInfoControl:
                if vInfo['contribution'] < 0.00001:
                    continue
                print 'SpecificOTU: %s, Contribution: %s' % (vInfo['specificOTU'], vInfo['contribution'])
                print 'Counts for MI: %s' % vInfo['verboseInfo']
                #print "Verbose Counts for MI: \nAO:%s,A'O:%s,AO':%s,A'O':%s" % (vInfo['ExtremelyVerbose']['Aut_OTU'], vInfo['ExtremelyVerbose']['AutN_OTU'], vInfo['ExtremelyVerbose']['Aut_OTUN'], vInfo['ExtremelyVerbose']['AutN_OTUN']) 
            print ""

            totalsAutistic.append(totalForAutistic)
            totalsControl.append(totalForControl)
            trueValues.append(trueValue)


        print "DFSDFSDFSDFDFDSFDSFSDFDSFSFSDFSDFSDFDSFS"
        print totalsAutistic
        print totalsControl
        print trueValues
        if totalsAutistic[0] > totalsAutistic[0] and trueValues[0] == 'Aut':
            print 'GOT right on aut'
        if totalsControl[0] > totalsControl[0] and trueValues[0] == 'Control':
            print 'GOT right on cont'       


    print "resultsList"
    print resultsList
    print "Correct Count:"
    print correctCount
    print "Incorrect Count:"
    print incorrectCount
    print "p:"
    print p
    print "PControl:"
    print pControl





def runTestWithBatchHoldout(batchNumber):
    correctCount = 0
    incorrectCount = 0
    resultsList = []
    
    samplesToHoldOut = getAllSampleIDsFromBatch(allRows, 0)
    MIMapping, allSamplesReferencedBuildingMI, pControl = findMI(samplesToHoldOut)

    for sampleToHoldOut in samplesToHoldOut:

        # trueValue = isSampleAutistic(allRows, sampleToHoldOut)
        # if trueValue == 'Aut':
        #     continue

        verboseInfoControl = []
        verboseInfoAut = []

        if sampleToHoldOut in dontDoTheseSamples:
            continue
        print "Predicting for sample: %s" % sampleToHoldOut

        totalForAutistic = 0
        totalForControl = 0

        allTrainSamplesUsedUponComparisonTime = []

        for specificOTU in allOTUs:
            OTUValue = getValueOfOTU(otuTable, specificOTU, sampleToHoldOut)
            if OTUValue != 0:
                autCount = {'count': 0, 'otuCount': 0}
                controlCount = {'count': 0, 'otuCount': 0}
                for trainSample in allSampleIDs:
                    if trainSample in dontDoTheseSamples or trainSample in samplesToHoldOut:
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
            
                fractionAut = float(autCount['otuCount'])/float(autCount['count'])
                fractionControl = float(controlCount['otuCount'])/float(controlCount['count'])    

                if fractionAut > 0:# and fractionControl == 0:
                    verboseInfoAut.append({'specificOTU': specificOTU, 'contribution':math.pow(MIMapping[specificOTU]['MutualInformation'], p), 'verboseInfo': MIMapping[specificOTU]['Metrics'], 'ExtremelyVerbose': MIMapping[specificOTU]['ExtremelyVerbose']})
                    totalForAutistic += math.pow(MIMapping[specificOTU]['MutualInformation'], p)
                if fractionControl > 0:# and fractionAut == 0:
                    verboseInfoControl.append({'specificOTU': specificOTU, 'contribution':math.pow(MIMapping[specificOTU]['MutualInformation'], p), 'verboseInfo': MIMapping[specificOTU]['Metrics'], 'ExtremelyVerbose': MIMapping[specificOTU]['ExtremelyVerbose']})
                    totalForControl +=  math.pow(MIMapping[specificOTU]['MutualInformation'], p)

                # totalForAutistic += MIMapping[specificOTU]['MutualInformation']*float(fractionAut)
                # totalForControl +=  MIMapping[specificOTU]['MutualInformation']*float(fractionControl)

                # print "fractionAut"
                # print fractionAut
                # print "fractionControl"
                # print fractionControl
                # print "autCount"
                # print autCount
                # print "controlCount"
                # print controlCount


        if totalForAutistic > totalForControl:
            prediction = 'Aut'
        else:
            prediction = 'Control'

        trueValue = isSampleAutistic(allRows, sampleToHoldOut)
        resultsList.append({'SampleToHoldOUt': sampleToHoldOut, 'Aut': totalForAutistic, 'Control': totalForControl, 'trueValue': trueValue})

        # Holdout verification, REALLY double check the OTU was indeed held out:
        print ""
        print "allTrainSamplesUsedUponComparisonTime:"
        print allTrainSamplesUsedUponComparisonTime
        print "allSamplesReferencedBuildingMI:"
        print allSamplesReferencedBuildingMI
        if sampleToHoldOut in allTrainSamplesUsedUponComparisonTime:
            print "ERROR: was used in tally!"
            sys.exit()
        if sampleToHoldOut in allSamplesReferencedBuildingMI:
            print "ERROR: was used in MI."
        # Just to be super duper sure.
        for vInfo in verboseInfoAut + verboseInfoControl:
            allSamples = vInfo['ExtremelyVerbose']['Aut_OTU'] + vInfo['ExtremelyVerbose']['AutN_OTU'] + vInfo['ExtremelyVerbose']['Aut_OTUN'] + vInfo['ExtremelyVerbose']['AutN_OTUN']
            if sampleToHoldOut in allSamples:
                print "ERRRRROR: was used in MI's, and still missed earlier."
                sys.exit()
        print "It has been verified, the sampleToHoldOut was not used in feature-selection or the run-time comparison."

        ### PRINTING RESULTS:
        print ""
        if trueValue == prediction:
            correctCount += 1
            print "CCCCCOOOORRRRECT"
        else:
            incorrectCount += 1
            print "INCORRECT"
        print resultsList
        print "Correct:"
        print correctCount
        print "Incorrect:"
        print incorrectCount
        print "powerOf:"
        print p

        print ""
        print "Justfication:"

        verboseInfoAut = sorted(verboseInfoAut, key=lambda x: x['contribution'], reverse=True)
        verboseInfoControl = sorted(verboseInfoControl, key=lambda x: x['contribution'], reverse=True)

        print "Aut justification,"
        for vInfo in verboseInfoAut:
            if vInfo['contribution'] < 0.00001:
                continue
            print 'SpecificOTU: %s, Contribution: %s' % (vInfo['specificOTU'], vInfo['contribution'])
            print 'Counts for MI: %s' % vInfo['verboseInfo']
            #print "Verbose Counts for MI: \nAO:%s,A'O:%s,AO':%s,A'O':%s" % (vInfo['ExtremelyVerbose']['Aut_OTU'], vInfo['ExtremelyVerbose']['AutN_OTU'], vInfo['ExtremelyVerbose']['Aut_OTUN'], vInfo['ExtremelyVerbose']['AutN_OTUN']) 
        print "Control justification,"
        for vInfo in verboseInfoControl:
            if vInfo['contribution'] < 0.00001:
                continue
            print 'SpecificOTU: %s, Contribution: %s' % (vInfo['specificOTU'], vInfo['contribution'])
            print 'Counts for MI: %s' % vInfo['verboseInfo']
            #print "Verbose Counts for MI: \nAO:%s,A'O:%s,AO':%s,A'O':%s" % (vInfo['ExtremelyVerbose']['Aut_OTU'], vInfo['ExtremelyVerbose']['AutN_OTU'], vInfo['ExtremelyVerbose']['Aut_OTUN'], vInfo['ExtremelyVerbose']['AutN_OTUN']) 
        print ""


    print "resultsList"
    print resultsList
    print "Correct Count:"
    print correctCount
    print "Incorrect Count:"
    print incorrectCount
    print "p:"
    print p





# findMI([])
# sys.exit()

runFullTestWithHoldout()
# print getAllSampleIDsFromBatch(allRows, 0)
# runTestWithBatchHoldout(0)


# Try only doing using a single batch.
# dontDoTheseSamples += getAllSampleIDsFromBatch(allRows, 1) + getAllSampleIDsFromBatch(allRows, 2) + getAllSampleIDsFromBatch(allRows, 3) + getAllSampleIDsFromBatch(allRows, 4)
# runFullTestWithHoldout()










