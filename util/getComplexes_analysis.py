# Useful miscellaneous functions for analysing complexes

import glob, itertools
import numpy as np
import os, re
import pandas as pd

from util import getData
from util.getComplexes import buildComplexes

def getCORUM_baselines(limitingLabels, corumLoading_specs, corumData_key='rawData_xls'):
    limitingLabels = set(list(limitingLabels.values())[0])
    if not all(isinstance(ele, str) for ele in limitingLabels):
        print('Limiting labels passed include non-string entries. For comparison with CORUM, ensure all limiting '
              'are string type.')
        return

    else:
        corumData_limited = dict()
        corumData = getData.LoadDict().source(corumLoading_specs)
        corumProteins_limited = corumData['proteins'].intersection(limitingLabels)
        corumData_limited['corumProteins_limited'] = corumProteins_limited

        corumComplexes_limited = dict()
        corumComplexes_lengthsLimited = dict()
        corumComplexes_limitedOverlap = dict()
        for key, val in corumData['complexes'].items():
            corumComplexes_limited[key] = val.intersection(limitingLabels)
            corumComplexes_lengthsLimited[key] = len(corumComplexes_limited[key])
            corumComplexes_limitedOverlap[key] = \
                len(corumComplexes_limited[key]) / len(val)
        corumComplexes_limited = set([frozenset(ele) for ele in corumComplexes_limited.values()])

        corumData_limited['corumComplexes_limited'] = corumComplexes_limited
        corumData_limited['corumComplexes_lengthsLimited'] = corumComplexes_lengthsLimited
        corumData_limited['corumComplexes_limitedOverlap'] = corumComplexes_limitedOverlap

        print('The max complex limited length is ' +
              str(np.amax(np.array(list(set(list(
                  corumComplexes_lengthsLimited.values())))))) + '.')
        print('The mean complex limited length is ' +
              str(np.mean(np.array(list(set(list(
                  corumComplexes_lengthsLimited.values())))))) + '.')
        print('The median complex limited length is ' +
              str(np.median(np.array(list(set(list(
                  corumComplexes_lengthsLimited.values())))))) + '.')
        print('The min complex limited length is ' +
              str(np.amin(np.array(list(set(list(
                  corumComplexes_lengthsLimited.values())))))) + '.')
        print('____________________________________________________________')

        print('The input dataset enables an average of ' + str(
            np.mean(np.array(list(corumComplexes_limitedOverlap.values())))) +
              ' maximum recapitulation among CORUM complexes.')
        overlaps = np.array(list(corumComplexes_limitedOverlap.values()))
        print('The input dataset enables 50% recapitulation or greater of ' +
              str(len(overlaps[overlaps >= 0.5])) + ' CORUM complexes out of ' +
              str(len(corumData[corumData_key])))

    return corumData_limited

def getComplexes_setPlus_proteins(complexesDF, colName=0):
    complexesCounter = 0

    proteinsSet = set()
    complexesSet = dict()
    complexesLengths = dict()
    for idx in complexesDF.index:
        placeholder = \
          frozenset([ele for ele in re.split(',\(|\(|\),|\)|\t|\s|;|,',
                                             complexesDF.loc[idx, colName]) if ele])
        if len(placeholder) != 0:
          complexesSet[idx] = placeholder
          complexesLengths[idx] = len(complexesSet[idx])
          proteinsSet = proteinsSet.union(complexesSet[idx])
          complexesCounter+=1

    print('There is a total of ' + str(len(complexesSet)) +
          ' predicted complexes and a total of ' + str(len(proteinsSet)) + ' proteins.')
    print('The minimum and maximum lengths of predicted complexes were ' +
          str(np.amin(list(complexesLengths.values()))) + ' and ' +
          str(np.amax(list(complexesLengths.values()))) + ' respectively.')

    return proteinsSet, complexesSet

def condensePreds_proteinSet(predictedComplexes_dict, verbose=False):
    proteinSet = set()
    for complexSet in list(predictedComplexes_dict.values()):
        proteinSet = proteinSet.union(complexSet)

    if verbose:
        print('The length of the base protein set is: ' +
              str(len(proteinSet)))

    return proteinSet

def quickCompare_contrast(refComplex, predictedComplexes_dict, proteinSet_adjust=True, jindex=0.6, returnResults=False):
    predsProteins_base = condensePreds_proteinSet(predictedComplexes_dict)
    print('Proteins overlap with refComplex: ' +
          str(len(predsProteins_base.intersection(refComplex))) +
          ' out of ' + str(len(refComplex)) + '.')
    if proteinSet_adjust:
        refComplex = refComplex.intersection(predsProteins_base)

    supersetCount = 0
    subsetCount = 0
    predLens = np.empty((0,))
    predLen_refConstituents = np.empty((0,))
    jaccardIndices = np.empty((0,))
    for complexSet in list(predictedComplexes_dict.values()):
        # count supersets
        supersetCount += complexSet.issuperset(refComplex)

        # count subsets
        subsetCount += complexSet.issubset(refComplex)

        # track predicted complexes' lengths
        predLens = np.concatenate((predLens, [len(complexSet)]))

        # calculate J indices
        jaccardIndices = \
            np.concatenate((jaccardIndices,
                            [len(complexSet.intersection(refComplex))/len(complexSet.union(refComplex))]))

        # quantify overlap
        predLen_refConstituents = \
            np.concatenate((predLen_refConstituents,
                            [len(complexSet.intersection(refComplex))]))

    predsInterest = predLens[np.argwhere(jaccardIndices>=jindex)]
    overlapInterest = predLen_refConstituents[np.argwhere(jaccardIndices>=jindex)]

    print('# of predictions overlapping with reference complex: ' +
          str(len(predLen_refConstituents[predLen_refConstituents>0])))
    print('# of predictions for which refComplex is a superset: ' +
          str(supersetCount))
    print('# of predictions for which refComplex is a subset: ' +
          str(subsetCount))
    print('# of predictions with >= 0.6 Jaccard index: ' +
          str(len(jaccardIndices[jaccardIndices>=jindex])))
    try:
        print('Max overlap of predictions with reference complex having Jaccard index of >=' + str(jindex) + ': '
          + str(np.amax(overlapInterest)))
        print('Mean overlap of predictions with reference complex having Jaccard index of >=' + str(jindex) + ': '
          + str(np.mean(overlapInterest)))
    except ValueError:
        pass
    print('____________________________________________________')

    if returnResults:
        results = dict()
        results['refComplex_membersPresent'] = \
            str(len(predsProteins_base.intersection(refComplex))) + '/' + \
            str(len(refComplex))
        results['supersetCount'] = supersetCount
        results['subsetCount'] = subsetCount
        results['overlapCount_Jindex'] = \
            len(jaccardIndices[jaccardIndices>=jindex])
        results['jaccardIndices'] = jaccardIndices
        results['pred2ref'] = predLen_refConstituents
        try:
            results['overlapCount_Jindex_max'] = np.amax(overlapInterest)
            results['overlapCount_Jindex_mean'] = np.mean(overlapInterest)
            results['reconstitutionFraction_max'] = np.amax(overlapInterest)/len(refComplex)
            results['reconstitutionFraction_mean'] = np.mean(overlapInterest)/len(refComplex)
            results['referencePurity_maxFraction'] = np.amax(overlapInterest)/predsInterest[np.argmax(overlapInterest)][0]
            results['referencePurity_meanFraction'] = np.mean(overlapInterest)/np.mean(predsInterest)
            results['predsInterest'] = predsInterest
            results['overlapInterest'] = overlapInterest
        except ValueError:
            pass

        return results

def compareContrast_all(refComplex_set, predictedComplexes_dict, alignRef_2Pred=True, similarityMetric=0.6):

  reconFrac_max = np.empty((len(refComplex_set), 1))
  reconFrac_mean = np.empty((len(refComplex_set), 1))
  purityFrac_max = np.empty((len(refComplex_set), 1))
  purityFrac_mean = np.empty((len(refComplex_set), 1))
  for refIdx, refComplex in zip(np.arange(len(refComplex_set)), list(refComplex_set)):
    predRef_analysisResults = \
      quickCompare_contrast(refComplex, predictedComplexes_dict,
                            proteinSet_adjust=alignRef_2Pred, jindex=similarityMetric, returnResults=True)

    if 'reconstitutionFraction_max' in predRef_analysisResults:
      reconFrac_max[refIdx] = predRef_analysisResults['reconstitutionFraction_max']
    else:
      reconFrac_max[refIdx] = 0
    if 'reconstitutionFraction_mean' in predRef_analysisResults:
      reconFrac_mean[refIdx] = predRef_analysisResults['reconstitutionFraction_mean']
    else:
      reconFrac_mean[refIdx] = 0
    if 'referencePurity_maxFraction' in predRef_analysisResults:
      purityFrac_max[refIdx] = predRef_analysisResults['referencePurity_maxFraction']
    else:
      purityFrac_max[refIdx] = 0
    if 'referencePurity_meanFraction' in predRef_analysisResults:
      purityFrac_mean[refIdx] = predRef_analysisResults['referencePurity_meanFraction']
    else:
      purityFrac_mean[refIdx] = 0

  return reconFrac_max, reconFrac_mean, purityFrac_max, purityFrac_mean

def generateInputs_complexBuilder_testCase(proteinSets, outputPrefix,
                                           fakeConfidence_value=0.9,
                                           complexBuilder_threshold=0.9,
                                           runBuilder=False,
                                           *outputDir):

    allPairs_theoreticalCombined = np.empty((0, 2), dtype=str)
    for proteinSet in proteinSets:
        allPairs_theoreticalCombined = \
            np.vstack(
                (allPairs_theoreticalCombined,
                 np.vstack(list(itertools.combinations(proteinSet, 2))))
            )
        print('The new shape of iteratively generated pairs is: ' +
              np.vstack(list(itertools.combinations(proteinSet, 2))).shape)

    allPairs_theoreticalCombined = \
        np.insert(
            allPairs_theoreticalCombined, 0,
            np.arange(len(allPairs_theoreticalCombined)), axis=1)

    allPairs_theoreticalCombined_confidenceVals = \
        np.empty((len(allPairs_theoreticalCombined), 1))
    allPairs_theoreticalCombined_confidenceVals[:] = str(fakeConfidence_value)

    os.mkdir('./' + outputPrefix)
    np.savetxt('./' + outputPrefix + '/' + outputPrefix + '_pairsInfo.tsv',
               allPairs_theoreticalCombined, delimiter='\t', fmt='%s')
    np.savetxt('./' + outputPrefix + '/' + outputPrefix + '_pairsPredictions.txt',
               allPairs_theoreticalCombined_confidenceVals, delimiter=',',
               fmt='%.4f')

    # build complexes
    if runBuilder:
        buildComplexes(outputPrefix, complexBuilder_threshold, outputDir)

def generateInputs_complexBuilder_testCase_realData(pairsDF, proteinSets, idNames,
                                                    confidenceCol_name, outputPrefix,
                                                    complexBuilder_threshold=0.9,
                                                    runBuilder=False,
                                                    *outputDir):
    allProteins = set()
    for proteinSet in proteinSets:
        allProteins = allProteins.union(proteinSet)

    idNames_plusConf = list(set(idNames).union(set([confidenceCol_name])))
    pairsAvailable = pairsDF.loc[
        (
                (pairsDF[idNames[0]].isin(allProteins)) &
                (pairsDF[idNames[1]].isin(allProteins))
        ), idNames_plusConf]
    pairsAvailable[idNames[0]] = pairsAvailable[idNames[0]].astype('str')
    pairsAvailable[idNames[1]] = pairsAvailable[idNames[1]].astype('str')

    pairIDs = pairsAvailable.loc[:, idNames].to_numpy()
    pairIDs = np.insert(pairIDs, 0, np.arange(len(pairIDs)), axis=1)
    pairConfidences = pairsAvailable.loc[:, confidenceCol_name]

    os.mkdir('./' + outputPrefix)
    np.savetxt('./' + outputPrefix + '/' + outputPrefix + '_pairsInfo.tsv',
               pairIDs, delimiter='\t', fmt='%s')
    np.savetxt('./' + outputPrefix + '/' + outputPrefix + '_pairsPredictions.txt',
               pairConfidences, delimiter=',', fmt='%.4f')

    # build complexes
    if runBuilder:
        buildComplexes(outputPrefix, complexBuilder_threshold, outputDir)

def generateInputs_complexBuilder(pairsDF, idNames, confidenceCol_name, outputPrefix,
                                  complexBuilder_threshold=0.9, runBuilder=False, *outputDir):

    pairsDF[idNames[0]] = pairsDF[idNames[0]].astype('str')
    pairsDF[idNames[1]] = pairsDF[idNames[1]].astype('str')
    pairIDs = pairsDF.loc[:, idNames].to_numpy()
    pairIDs = np.insert(pairIDs, 0, np.arange(len(pairIDs)), axis=1)
    pairConfidences = pairsDF.loc[:, confidenceCol_name]

    os.mkdir('./' + outputPrefix)
    np.savetxt('./' + outputPrefix + '/' + outputPrefix + '_pairsInfo.tsv',
               pairIDs, delimiter='\t', fmt='%s')
    np.savetxt('./' + outputPrefix + '/' + outputPrefix + '_pairsPredictions.txt',
               pairConfidences, delimiter=',', fmt='%.4f')

    # build complexes
    if runBuilder:
        buildComplexes(outputPrefix, complexBuilder_threshold, outputDir)

def acquirePredictions(dataDir_patt, modelResults_dirPatt, outputPrefix, suffix='',
                       complexBuilder_threshold=0.9, runBuilder=False, *outputDir):
    numFolds = np.arange(5)
    allLabels = {fold: np.loadtxt(predictionsDir, delimiter='\t', skiprows=1, dtype=np.int) for fold, predictionsDir in
                 zip(numFolds, glob.glob(dataDir_patt + 'testLabels.tsv'))}
    probsPos_combinedWeighted_byFolds = \
        {fold: np.loadtxt(predictionsDir, delimiter=',') for fold, predictionsDir in
         zip(numFolds, glob.glob(modelResults_dirPatt + 'probsPos_weighted' + suffix + '.csv'))}
    labels_acrossFolds = []
    combinedWeighted_probsPos_kCV = []
    for combinedWeighted_probsPos, labels_byFolds in zip(probsPos_combinedWeighted_byFolds.values(), allLabels.values()):
        labels_acrossFolds.append(labels_byFolds)
        combinedWeighted_probsPos_kCV.append(combinedWeighted_probsPos)

    labels_acrossFolds = np.concatenate(labels_acrossFolds)
    combinedWeighted_probsPos_kCV = np.concatenate(combinedWeighted_probsPos_kCV)
    print('The lengths of the labels\'s and predictions\'s datasets are ' +
          str(len(labels_acrossFolds)) + ' and ' +
          str(len(combinedWeighted_probsPos_kCV)) + ' respectively.')

    os.mkdir('./' + outputPrefix)
    np.savetxt('./' + outputPrefix + '/' + outputPrefix + '_pairsInfo.tsv',
               labels_acrossFolds, delimiter='\t', fmt='%d')
    np.savetxt('./' + outputPrefix + '/' + outputPrefix + '_pairsPredictions.txt',
               combinedWeighted_probsPos_kCV, delimiter=',', fmt='%.4f')

    # build complexes
    if runBuilder:
        buildComplexes(outputPrefix, complexBuilder_threshold, outputDir)

    return labels_acrossFolds, combinedWeighted_probsPos_kCV

def getComplex_predictions(outputDir, outputPrefix, realDir=True):

  if not realDir:
    predictionsPath = \
      '/content/' + outputPrefix + '/' + outputPrefix + '_pairsPredictions*.txt'
  else:
    predictionsPath = \
      outputDir + outputPrefix + '/' + outputPrefix + '_pairsPredictions*.txt'

  predictedCliques_sets = dict()
  for filename in glob.glob(predictionsPath):
    filenameExtension = filename.split('.txt')[0].split('_pairsPredictions')[1]
    if filenameExtension.isnumeric():
      if int(filenameExtension) != 2:
        clique = dict()
        data = pd.read_csv(filename, sep='\t',
                           header=None, names=['complex', 'score'])
        for idx in data.index:
          line = [ele for ele in re.split(',', data['complex'].loc[idx]) if ele]
          clique[idx] = set(list(map(str, line)))

        predictedCliques_sets[int(filenameExtension)] = \
          [frozenset(ele) for ele in list(clique.values())]

  return predictedCliques_sets

def contrastCORUM(refComplexes, testComplexes):

    subsetsPlus_supersetsCounts = np.zeros((len(refComplexes), 2))
    subsetsCharacteristics = []
    supersetsCharacteristics = []
    for idx, refComplex in zip(np.arange(len(refComplexes)), list(refComplexes)):
        for testComplex in list(testComplexes):
            if (len(refComplex) != 0) & (len(testComplex) != 0):
                if refComplex.issubset(testComplex):
                    subsetsPlus_supersetsCounts[idx, 0] += 1
                    subsetsCharacteristics.append(int(len(testComplex)))
                if refComplex.issuperset(testComplex):
                    subsetsPlus_supersetsCounts[idx, 1] += 1
                    supersetsCharacteristics.append(int(len(testComplex)))
    refComparision = pd.DataFrame(data=subsetsPlus_supersetsCounts, columns=['subsets of', 'supersets of'])

    subsetsCharacteristics = np.array(subsetsCharacteristics)
    _, counts = np.unique(subsetsCharacteristics, return_counts=True)
    subsetsCharacteristics_modeValue = np.argwhere(counts == np.max(counts))

    supersetsCharacteristics = np.array(supersetsCharacteristics)
    _, counts = np.unique(supersetsCharacteristics, return_counts=True)
    supersetsCharacteristics_modeValue = np.argwhere(counts == np.max(counts))

    numReference_inImpure = len(refComparision.loc[refComparision['subsets of'] > 0, 'subsets of'])
    numImpure_mixes = refComparision.loc[refComparision['subsets of'] > 0, 'subsets of'].sum()
    impureComplex_meanLength = numImpure_mixes/numReference_inImpure

    numReference_inSubsets = len(refComparision.loc[refComparision['supersets of'] > 0, 'supersets of'])
    numReference_subsets = refComparision.loc[refComparision['supersets of'] > 0, 'supersets of'].sum()
    pureSubset_meanLength = numReference_subsets/numReference_inSubsets

    print('There were ' + str(numImpure_mixes) + ' impure assemblies from ' + str(numReference_inImpure) +
          ' reference complexes with an average complex length of ' + str(impureComplex_meanLength) +
          ', a max complex length of ' + str(subsetsCharacteristics.max()) + ', a min complex length of ' +
          str(subsetsCharacteristics.min()) + ', and finally a mode of ' + str(subsetsCharacteristics_modeValue) +
          '.')

    print('There were ' + str(numReference_subsets) + ' pure subsets from ' + str(numReference_inSubsets)
          + ' reference complexes with an average complex length of ' + str(pureSubset_meanLength) +
          ', a max complex length of ' + str(supersetsCharacteristics.max()) + ', a min complex length of ' +
          str(supersetsCharacteristics.min()) + ', and finally a mode of ' + str(supersetsCharacteristics_modeValue) +
          '.')

    subsetsStats = {'numImpure_mixes': numImpure_mixes, 'numReference_inImpure': numReference_inImpure,
                    'meanImpure_complexLength': impureComplex_meanLength,
                    'maxImpure_complexLength': subsetsCharacteristics.max(),
                    'minImpure_complexLength': subsetsCharacteristics.min(),
                    'modeImpure_complexLength': subsetsCharacteristics_modeValue}

    supersetsStats = {'numReference_subsets': numReference_subsets, 'numReference_inSubsets': numReference_inSubsets,
                    'meanSubset_complexLength': pureSubset_meanLength,
                    'maxSubset_complexLength': supersetsCharacteristics.max(),
                    'minSubset_complexLength': supersetsCharacteristics.min(),
                    'modeSubset_complexLength': supersetsCharacteristics_modeValue}

    return subsetsStats, supersetsStats