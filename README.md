proteinPairs_complexMaps
====================

###Scripts for modeling positive protein interactions and predicting complexes

####data acquisition
```
python ./proteinPairs_complexMaps/util/getData.py dataSource *otherParameters
```
######Bioplex 1.0 (...will be adding v.2 and 3) 
input: bioplex, include protein set flag

output: dictionary with data (pandas), minimally processed; total protein set, if specified

*Example*
```
python ./proteinPairs_complexMaps/util/getData.py bioplex True
```


######CORUM (various releases 2012--2018, updated 2022.01.30)
input: corum, release date, dataset type, processing flags

output: dictionary with data (pandas) in multiple forms after incremental processing changes, total protein set

*Example*
```
python ./proteinPairs_complexMaps/util/getData.py corum yyyy.mm.dd datasetType ribosomeOnly_flag missingGeneids_flag
 explicitNone_flag onlyHumans_flag

The following command reconstituted the exact set of CORUM proteins used in Drew 2017 to generate their training and 
testing pairs

python ./proteinPairs_complexMaps/util/getData.py corum 2012.02.17 allComplexesCore False True False False
```

*Notes*
```
Release dates include:
2012.02.17
2016.12.15 ~website maps to 2017.05.03 in error
2017.05.03
2017.07.02
2018.07.01
2018.09.03

Dataset types (worksheets) include: 
allComplexesCore
allComplexes
spliceComplexes (available only for 2018 releases)

Processing flags include:
ribosomeOnly_flag: returns only data for small and large ribosomal subunits. Default: False.
missingGeneid_flag: filters the raw CORUM data for any lines possessing the empty data, e.g. NaN. Default: False.
explicitNone_flag: accommodates filter in older formats in which "None" maps to empty data instead of NaN
onlyHumans_flag: filters the raw CORUM data for only human complexes
*May need to investigate the joint presence potential of "None" and NaN   

Data were imported from .xls files using the source files (.txt and .csv) main delimiters (; and tab) to facilitate 
fewer changes in code for acquisitions across release types. All columns were imported as text.

Noted remaining format inconsistencies include: delimiter type, case, spacing in column names, and empty/NaN 
designators.  
```


######Fantom (...under development)
input: ...

output: ...

*Example*
```
python ./proteinPairs_complexMaps/util/getData.py ...
```


######GO (Cellular Component) (under development...)
input: go

output: dataframe (pandas) with hpa and GO cellular component (main-, and other locations)

*Example*
```
python ./proteinPairs_complexMaps/util/getData.py go
```


######GTEx (...under development)
input: ...

output: ...

*Example*
```
python ./proteinPairs_complexMaps/util/getData.py ...
```


######Hein 2015 
input: hein2015 rawForm_flag ignoreIsoforms_flag mergeGeneids_flag manMerge_refPlus_intact2Geneids_flag 
outputProteins_flag

output: <= 4 dataframes corresponding to 1 of 4 flags progressively, e.g. 3rd flag corresponds to 3rd dataframe, 
proteins set

*Example*
```
python ./proteinPairs_complexMaps/util/getData.py hein2015 False True True True True
```

*Notes*
```
Flags include:
rawForm_flag: returns dataset with minimal processing
ignoreIsoforms_flag: some entries appended with '-' indicating isoforms; exclusion may increase available # of geneid 
identifier matches

mergeGeneids_flag: returns dataset with geneid identifiers
manMerge_refPlus_intact2Geneids_flag: returns dataset containing additional geneid matches from manual search of a small
number of entries pertaining to refseq and intact databases

outputProteins_flag: returns set of unique proteins constituting study
```


######huMAP 1.0
input: humap1, return pairs flag, return complexes flag, return feature matrix flag, return results flag

output: dictionary(ies) with data (pandas), corresponding in order with each respective, subsequent flag 

*Example*
```
python ./proteinPairs_complexMaps/util/getData.py humap1 True True True False
```

*Note*
```
Development of results feature still under development pending utility analysis for our study's purposes
```


######huMAP 2.0
input: humap2, return format, return pairs flag, return complexes flag, return feature matrix flag, return results flag

output: dictionary(ies) with data (pandas), corresponding in order with each respective, subsequent flag 

*Example*
```
python ./proteinPairs_complexMaps/util/getData.py humap2 geneid True True True False
```

*Note*
```
Return formats include: acc, geneid, and genename (only for results)
Development of results feature still under development pending utility analysis for our study's purposes
```


######Humphreys 2021
input: ...

output: ...

*Example*
```
python ./proteinPairs_complexMaps/util/getData.py ...
```


######\'Kaggle\' (...under development) 
input: ...

output: ...

*Example*
```
python ./proteinPairs_complexMaps/util/getData.py ...
```


######NCI-60 (Gholami, et al., 2013)
input: nci60, maximize found GeneIDs flag

output: dictionary(ies) with data (pandas) corresponding to cell and tissue protein expression and RNA expression; 
protein sets corresponding to each type of expression, the union of the protein expression sets, and the intersection of 
the protein expression sets with the RNA expression set

*Example*
```
python ./proteinPairs_complexMaps/util/getData.py nci60 True
```


######PROPER-seq (...under development)
input: ...

output: ...

*Example*
```
python ./proteinPairs_complexMaps/util/getData.py ...
```


######SubCellBarCode
input: scbc, reviewSCBC_flag, reviewMappings_flag, mappingFlag, mappingPath, mergeData_flag, dataPath, labelsPath, 
dataHeader_flag, labelsHeader_flag

output: dictionary(ies) with data (pandas) corresponding to scbc cell-specific expression data; mappings pertaining 
to GeneCard-Uniprot, Uniprot-GeneID, and GeneCard-GeneID respectively; dataframe with all SCBC protein pairs' 
correlation values (saved to disk), and the input data, labels merged with SCBC cell-specific protein co-expression

*Example*
```
python ./proteinPairs_complexMaps/util/getData.py scbc True True True None True <path-to-input-data> 
<path-to-input-labels> True True
```

*Notes*
```
Flags, other arguments include:
reviewSCBC_flag: returns dataframe for original SCBC cell-specific protein expression data
reviewMappings_flag: returns dataframes representing each step of the mapping from Uniprot IDs through to GeneIDs
mappingFlag: generates or loads existing correlation values for all SCBC protein pairs
mappingPath: path to alternate location of cell-specific correlation values for SCBC protein pairs (perhaps unnecessary)
mergeData_flag: generates merger between SCBC cell-specific protein pair correlations and input data, labels
dataPath: path to data for merging with SCBC cell-specific protein pair correlation data
labelsPath: path to labels for merging with SCBC cell-specific protein pair correlation data 
(should have both IDs and label)

dataHeader_flag: directs processing to acknowledge header
labelsHeader_flag: directs processing to acknowledge header

*May need to investigate the joint presence potential of "None" and NaN   

Data were imported from .xls files using the source files (.txt and .csv) main delimiters (; and tab) to facilitate 
fewer changes in code for acquisitions across release types. All columns were imported as text.

Noted remaining format inconsistencies include: delimiter type, case, spacing in column names, and empty/NaN 
designators.  
```


######STRING
input: ...

output: ...

*Example*
```
python ./proteinPairs_complexMaps/util/getData.py ...
```


######UniProt
input: uniprot, indEntry_flag, removeEmpty_geneidFlag, hpaData_flag, dropEmpty_locsFlag, 
hpaSeparate_geneidSingle_entriesFlag, hpaSeparate_genenameSingle_entriesFlag, hpaMerge_geneidUniprot_flag, 
hpaMerge_geneidGenename_flag, hpaMerge_geneidUniprot_plusGenename_flag, hpaOutput_proteinsFlag

output: dictionary(ies) with data (pandas) corresponding to UniProt human data with varying degrees of processing to 
identify missed GeneIDs, swissProt data (used in Lugo-Martinez 2019) featuring subcellular location data in <= 6 stages
of iterative processing pertaining to respective flags; total protein set

*Example*
```
(most conservative)
python ./proteinPairs_complexMaps/util/getData.py uniprot False False True True False False False False False True
(most relaxed)
python ./proteinPairs_complexMaps/util/getData.py uniprot True False True True True True True True True True
```

*Notes*
```
indEntry_flag: \'explode\' UniProt dataframe\'s GeneCard, GeneID, and Uniprot entries into individual entries retaining
original referenced attributes

removeEmpty_geneidFlag: drop entries for which there isn\'t a GeneID identifier
hpaData_flag: return swissProt data used in Lugo-Martinez\'s 2019 study
dropEmpty_locsFlag: drop entries for which subcellular location data is empty
hpaSeparate_geneidSingle_entriesFlag: \'explode\' swissProt dataframe\'s GeneID entries into individual entries 
retaining original referenced attributes

hpaSeparate_genenameSingle_entriesFlag: \'explode\' swissProt dataframe\'s GeneName entries into individual entries 
retaining original referenced attributes

hpaMerge_geneidUniprot_flag: identify missed GeneID matches by merging GeneIDs from UniProt data to entries in swissProt
with matching Uniprot IDs and missing GeneIDs

hpaMerge_geneidGenename_flag: identify missed GeneID matches by merging GeneIDs from UniProt data to entries in swissProt
with matching GeneNames and missing GeneIDs

hpaMerge_geneidUniprot_plusGenename_flag: identify missed GeneID matches by merging GeneIDs from UniProt data to entries in swissProt
with matching UniProt IDs and GeneNames that are missing GeneIDs

hpaOutput_proteinsFlag: return total set of proteins
```


######Wan 2015
input: wan2015, include protein set flag

output: dictionary(ies) with data (pandas) corresponding to Wan 2015 data; total protein set

*Example*
```
python ./proteinPairs_complexMaps/util/getData.py wan2015 True
```

