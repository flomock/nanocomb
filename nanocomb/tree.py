# stackoverflow, 26146623
import csv
from collections import defaultdict
from pprint import pprint
import os
from pymongo import MongoClient
import pandas as pd
import numpy as np

def new_tree(): return defaultdict(new_tree)


def tree_add(t, path):
    for node in path:
        t = t[node]


def pprint_tree(tree_instance):
    def dicts(t): return {k: dicts(t[k]) for k in t}
    pprint(dicts(tree_instance))




def tree_to_newick(root):
    items = []
    for k in iter(root.keys()):
        s = ''
        if len(root[k].keys()) > 0:
            sub_tree = tree_to_newick(root[k])
            if sub_tree != '':
                s += '(' + sub_tree + ')'
        s += k
        items.append(s)
    return ','.join(items)


# def csv_to_weightless_newick(input):
#     t = csv_to_tree(input)
#     #pprint_tree(t)
#     return tree_to_newick(t)


def csv2dict(fp):
    '''
    Tested on ICTV virus taxonomy.

    from io import StringIO
    from Bio import Phylo


    with open('tax.nwk', 'w+') as out:
        out.write(tree_to_newick(t) + '\n')

    new_tree = Phylo.read(StringIO(csv2newick(fp)), 'newick')
    '''

    t = new_tree()
    with open(fp, 'r') as file:
        csv_reader = csv.reader(file, delimiter='\t')
        _ = next(csv_reader)  # discard header
        # for row in itertools.islice(csv_reader, 1230, 1240):
        for row in csv_reader:
            tax = row[1:6]
            try:
                tax.remove('')  # some entries don't have a subfamily entry
            except ValueError:
                pass
            tree_add(t, tax)
    return t

def csv2newick(filepath):
    return tree_to_newick(csv2dict(filepath))

def add_all_leafs(file, collection):
    """
    1. go through tree leaves as potential parent parameter
    2. search for parent db and add _id as leaf
    3. return finished tree
    :param file: basic dict, no _ids included
    :param collection: the db entries, use _ids as new leafs
    :return: dict with _ids as leafs
    """
    for k in iter(list(file)):

        if len(list(file[k])) > 0:
            add_all_leafs(file[k], collection)
        else:
            species = k
            ids = []
            for element in collection.find({"parent": species}, {"_id": 1}):
                ids.append(element['_id'])
            if len(ids) > 0:

                tree = file[species]
                for id in ids:
                    tree[id]
                # Uncomment if list wanted
                # file.update({species:ids})

    return



def get_db(mongoclient="localhost:27017", db_name="testDB", collection_name="testCell2"):
    client = MongoClient(mongoclient)
    db = client[db_name]
    collection = db.get_collection(collection_name)
    return collection

def get_samples(virus_file, host_file, collection, host_clade="all", virus_clade="all"):
    """
    get samples which fulfill the clade conditions
    :param virus_file: dict tree for viruses
    :param host_file: dict tree for hosts
    :param collection: the database
    :param host_clade: which host clade should be used
    :param virus_clade: which virus clade should be used
    :return: samples as pandas array
    """
    #zugriff auf daten durch angabe von Host und Clade
    def get_leaves(file):
        """
        go to leafs of tree and put together in array
        :param file: dict tree representation
        :return: Array[Ids]
        """
        leave = np.array([])
        for k in iter(list(file)):

            if len(list(file[k])) > 0:
                leave = np.append(leave,get_leaves(file[k]))
            else:
                leave = np.append(leave,k)

        return leave

    def get_ids(file, filter):
        """
        get ids for clade
        :param file: tree of hosts or viruses
        :param filter: clade of interest
        :return: all leaves of interest
        """
        for k in iter(file.keys()):
            if k == filter:
                virus_leaves = get_leaves(file[k])
                return virus_leaves
            else:
                if len(file[k].keys()) > 0:
                    found = get_ids(file[k],filter)
                    if len(found) > 0:
                        return found
        return []

    if virus_clade == "all":
        virus_leaves = get_leaves(virus_file)
    else:
        virus_leaves = get_ids(virus_file, virus_clade)

    if host_clade == "all": #TODO delete hard coded example
        host_leaves = np.array(["ae02c884-187c-492a-a1be-880df33c3f51", "71057421-e62a-4633-b3ba-f2348f77c4bc"])
        # host_leaves = get_leaves(host_file)
    else:
        # host_leaves = get_ids(host_file,host)
        host_leaves = np.array(["ae02c884-187c-492a-a1be-880df33c3f51","71057421-e62a-4633-b3ba-f2348f77c4bc"])
    assert len(virus_leaves) > 0, "no valid virus clade"
    assert len(host_leaves) > 0, "no valid host clade"

    def compare_ids(virus,host):
        ids = set(virus).intersection(host)
        return ids

    ids = compare_ids(virus_leaves,host_leaves)

    def get_table_samples(ids):
        """
        parse ids to usable sample table
        :param ids: ids of interest
        :return: pandas dataframe with id, host, seq and parent of samples
        """
        l=[]
        for id in ids:
            for element in collection.find({"_id": id}):
                l.append(element)
        gen = ((i['_id'], i['host'], i['seq'], i['parent']) for i in l)
        df = pd.DataFrame.from_records(gen, columns=['id', 'host', 'seq', 'parent'])
        print(df)
        # test = pd.DataFrame()
        # for h in HOSTS:
        #     df_host = df[df.host == h]
        #     test_sub = df_host.sample(frac=SIZE, random_state=SEED)
        #     test = test.append(test_sub)

    get_table_samples(ids)


cwd = os.getcwd()
collection = get_db()
virus_tree = csv2dict(cwd + "/examples/small.csv")
add_all_leafs(virus_tree, collection)
print(tree_to_newick(virus_tree))
samples = get_samples(virus_tree, None, collection, virus_clade="all")
