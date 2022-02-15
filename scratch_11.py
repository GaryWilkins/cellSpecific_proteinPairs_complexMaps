import importlib, itertools
import numpy as np
import pandas as pd
import re

from util import getData

# test CORUM acquisition code
'''
dirs = ['2018.09.03', '2018.07.01', '2017.07.02', '2017.05.03', '2016.12.15', '2012.02.17']
complexSet_types = ['allComplexesCore', 'allComplexes', 'spliceComplexes']
flagsOptions = [True, False]
possibleCORUM_loadingSpecs = list(itertools.product(['corum'], dirs, complexSet_types, flagsOptions, flagsOptions, flagsOptions, flagsOptions))
print(len(possibleCORUM_loadingSpecs))
for count, loadingSpec in zip(np.arange(len(possibleCORUM_loadingSpecs)), possibleCORUM_loadingSpecs):
    print('count: ' + str(count))
    print(loadingSpec)
    getData.LoadDict().source(loadingSpec)
'''

# reference for Drew 2017 recapitulation and pairs origin investigation
'''
drew2017_corum2012_recapitulationSpecs = ['corum', '2012.02.17', 'allComplexesCore', False, True, False, False]
drew2017_corum2012_recapitulation = getData.LoadDict().source(drew2017_corum2012_recapitulationSpecs)

## compared to original, found in Drew github
drew2017_geneidSet_primarySource = pd.read_csv('./sourceData/corum/2012.02.17/allComplexesCore_geneid.txt',
                                               sep='\t', header=None)
                                               
drew2017_geneidSet_primarySource_set = set()
for idx in drew2017_geneidSet_primarySource.index:
    complexLine = [ele for ele in re.split(',\(|\(|\),|\)|\t|\s|;|,', drew2017_geneidSet_primarySource.loc[idx, 0]) if ele]
    complexDelimited = set(list(map(str, complexLine)))
    drew2017_geneidSet_primarySource_set = drew2017_geneidSet_primarySource_set.union(complexDelimited)
'''

# test Hein2015 acquisition code
'''
flagsOptions = [True, False]
possibleHein_loadingSpecs = list(itertools.product(['hein2015'], flagsOptions, flagsOptions, flagsOptions, flagsOptions, flagsOptions))
print(len(possibleHein_loadingSpecs))
for count, loadingSpec in zip(np.arange(len(possibleHein_loadingSpecs)), possibleHein_loadingSpecs):
    print('count: ' + str(count))
    print(loadingSpec)
    getData.LoadDict().source(loadingSpec)
'''

# (?=<\d+>)
# (?<=<\d+>)
# ',|\t|\s|\(|\)|;'
# ',\(|\),|\t|\s|;|(?<=<\d>),(?=<\d>)'
#importlib.reload(getData)  #for re-loading module after in-time edit