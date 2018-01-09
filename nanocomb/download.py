import urllib.request
import os
import pandas as pd
from pathlib import Path
import numpy as np

"""
Download "all" entries of emble
"""

def get_accessions():
    """
    load locally saved accessions from tsv files
    :return: array with accession numbers
    """
    path = os.getcwd()+'/virus/'
    suffix = 'tsv'
    file = Path(path)
    accessions = np.array([])
    for f in file.glob('*'+suffix):
        test = pd.DataFrame.from_csv(f, sep='\t')
        accessions = np.append(accessions, test['GenBank Accession'].as_matrix())

    print(accessions)
    return np.unique(accessions)

IDArray = get_accessions()
l1 = len(IDArray)
length = 500

for i in range(0,len(IDArray),length):
    ID = IDArray[i:i + length]
    IDString = ','.join(ID)
    url = 'http://www.ebi.ac.uk/ena/data/view/' + IDString + '&display=text&download=txt&offset=i'
    response = urllib.request.urlopen(url)
    data = response.read()      # a `bytes` object
    # text = data.decode('utf-8')
    outfile = open(os.getcwd()+'/samples/from' + str(i) + 'to' + str(i+length-1)+'.txt', 'wb')
    outfile.write(data)
    outfile.close()
    i += length

    print(f"\t===== {i+1} / {l1} -- {int((i+1) / l1 * 100)}% =====") #, end='\r')
