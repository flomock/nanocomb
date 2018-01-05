# stackoverflow, 26146623
import csv
from collections import defaultdict
# import itertools
from pprint import pprint
import os
from pymongo import MongoClient


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



def get_db(mongoclient="localhost:27017", db_name="testDB", collection_name="testCell"):
    client = MongoClient(mongoclient)
    db = client[db_name]
    collection = db.get_collection(collection_name)
    return collection

cwd = os.getcwd()
collection = get_db()
file = csv2dict(cwd + "/examples/small.csv")
add_all_leafs(file,collection)
print(tree_to_newick(file))
# add_all_leafs(tree,collection)