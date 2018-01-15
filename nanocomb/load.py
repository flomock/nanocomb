#!/home/go96bix/my/programs/Python-3.6.1/bin/python3.6

from xmljson import parker
from xml.etree.ElementTree import fromstring
import json
import os
from pathlib import Path
from pymongo import MongoClient
from uuid import uuid4
from Bio import SeqIO
import xmltodict


def xml_to_json(path):
    """
    CAUTION: parses correct but result hard to query
    load xml file from path and parse to json
    :param path: path to xml file
    :return: json version of input file
    """
    file = Path(path)
    # if path is single file, than parse file to json
    # save filename.json
    if file.is_file():
        outputs = []
        with open(file, "r") as file:
            data = file.read()
            output = xmltodict.parse(data)
            print(output)

        with open(path[:-4]+".json", 'w+') as outjson:
            outjson.write(json.dumps(output) + '\n')
            outputs.append(json.dumps(output))
        return outputs

    # if path is directory, than get all xml files
    # parse files to json and save nanocomb.json
    elif file.is_dir():
        outputs = []
        with open(path + "nanocomb.json", 'w+') as outjson:
            for f in file.glob('*.xml'):
                with open(f, "r") as f:
                    data = f.read()
                    output = parker.data(fromstring(data))

                outjson.write(json.dumps(output) + '\n')
                outputs.append(json.dumps(output))
        return outputs

    else:
        raise KeyError('no valid file path')

def files_to_json(path, format = "embl"):
    """
    CAUTION: parses correct but result hard to query
    load xml file from path and parse to json
    :param path: path to txt file
    :return: json version of input file
    """
    file = Path(path)
    # if path is single file, than parse file to json
    # save filename.json
    if file.is_file():
        outputs = []
        a = SeqIO.read(path,format=format)
        output = a.annotations
        output.pop('references')
        host = a.features[0].qualifiers['host']
        output.update({'host':host})
        output.update({'seq': str(a.seq)})

        with open(path[:-4]+".json", 'w+') as outjson:
            outjson.write(json.dumps(output) + '\n')
            outputs.append(json.dumps(output))
        return outputs

    # if path is directory, than get all files
    # parse files to json and save nanocomb.json
    elif file.is_dir():
        outputs = []

        if format == 'embl':
            suffix = '.txt'

        with open(path + "nanocomb.json", 'w+') as outjson:
            # i = 0
            l1 = 1198
            for i, f in enumerate(file.glob('*'+suffix)):
                with open(f, "r") as f:
                    records = SeqIO.parse(f, format=format)
                    for record in records:
                        output = record.annotations
                        output.pop('references')
                        output.update({'parent': record.annotations['organism']})
                        try:
                            host = record.features[0].qualifiers['host'][0]

                        except:
                            # print(record)
                            # print()
                            continue

                        output.update({'host': host})
                        output.update({'seq': str(record.seq)})

                        outjson.write(json.dumps(output) + '\n')
                        outputs.append(json.dumps(output))

                print(f"\t===== {i+1} / {l1} -- {int((i+1) / l1 * 100)}% =====")
                i += 1
        return outputs

    else:
        raise KeyError('no valid file path')

def init(file, db, cell, client='localhost:27017'):
    '''
    Load json to mongodb and assign UUID.
    :param file: input json file
    :param db: the data base name
    :param cell: the name of the cell in the db
    :param client: where to find the db
    :return:
    '''
    # open db insert each single valid json file
    c = MongoClient(client)[db][cell]
    inserted = 0
    for line in file:
        d = json.loads(line.strip())
        d['_id'] = str(uuid4())
        c.insert_one(d)
        inserted += 1
    print(inserted, 'entries inserted into cell', '"' + cell + '".')
    print('Primary key assigned to field "_id".')

cwd = os.getcwd()
file = files_to_json(cwd + "/samples/")
init(file,"FinalDB","allData")