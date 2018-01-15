# stackoverflow, 26146623
import csv
from collections import defaultdict
from pprint import pprint
import os
from pymongo import MongoClient
import pymongo
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
import json
import re

def new_tree(): return defaultdict(new_tree)


def tree_add(t, path):
    for node in path:
        t = t[node]


def pprint_tree(tree_instance):
    def dicts(t): return {k: dicts(t[k]) for k in t}

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


def csv2dict(fp):
    '''
    parse ICTV csv to virus tree
    :param fp: filepath to csv
    :return: dict tree
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


def ncbi_tax2dict(data_dir):
    """
    get ncbi tax data and parses it to phylogenetic tree dict
    :param data_dir: path to names.dmp and nodes.dmp
    :return: phylo dict, and dict to translate id in name
    """

    col_delimiter = '\t|\t'
    row_delimiter = '\t|\n'
    print('Getting names...')
    scientific_names = {}
    common_names = {}
    with open(os.path.join(data_dir, 'names.dmp')) as names_file:
        for line in names_file:
            line = line.rstrip(row_delimiter)
            values = line.split(col_delimiter)
            tax_id, name_txt, _, name_type = values[:4]
            if name_type == 'scientific name':
                scientific_names[int(tax_id)] = name_txt
            elif name_type == 'common name':
                common_names[int(tax_id)] = name_txt

    print('Reading taxonomy...')
    nodes = new_tree()
    with open(os.path.join(data_dir, 'nodes.dmp')) as nodes_file:
        for line in nodes_file:
            line = line.rstrip(row_delimiter)
            values = line.split(col_delimiter)
            tax_id, parent_id = values[:2]
            tax_id = int(tax_id)
            parent_id = int(parent_id)
            # this_node = Newick.Clade(name=scientific_names[tax_id])
            this_node = nodes[parent_id]
            if len(this_node) == 0:
                nodes.update({parent_id: [tax_id]})
            else:
                names = this_node
                names.append(tax_id)

    def build_dict_tree(nodes_dict, file, parent='root'):

        file = file[parent]
        for child in nodes_dict[parent]:
            file[child]
            build_dict_tree(nodes_dict, file, parent=child)

    def prepare_dict(nodes_dict, parent='root'):
        dict_tree = new_tree()
        dict_tree['root']
        print("building dict tree...")
        dict_tree = dict_tree[parent]
        for child in nodes_dict[1]:
            # print(scientific_names[child])
            if scientific_names[child] != 'cellular organisms':
                continue
            build_dict_tree(nodes_dict, dict_tree, parent=child)
        return dict_tree

    dict_tree = prepare_dict(nodes)
    return dict_tree, scientific_names, common_names


def add_all_leaves_virus(file, collection,include_inner_nodes = False):
    """
    1. go through tree leaves as potential parent parameter
    2. search for parent db and add _id as leaf
    3. return finished tree
    :param file: basic dict, no _ids included
    :param collection: the db entries, use _ids as new leafs
    :return: dict with _ids as leafs
    """
    species_list = []
    for k in iter(list(file)):

        if len(list(file[k])) > 0 and not include_inner_nodes:
            add_all_leaves_virus(file[k], collection)
        else:
            if len(list(file[k])):
                add_all_leaves_virus(file[k], collection)

            species = k
            ids = []

            # for element in collection.find({"parent": species}, {"_id": 1}):
            # for element in collection.find({"parent": {'$in': [re.compile(species)]}}, {"_id": 1}):


            sub = species.split(" ")
            species_list.append(sub[0]+" "+sub[1])
            # for element in collection.find({"parent": re.compile(sub[0]+" "+sub[1])}, {"_id": 1}):
            #     ids.append(element['_id'])
            # if len(ids) > 0:
            #
            #     tree = file[species]
            #     for id in ids:
            #         tree[id]
                # Uncomment if list wanted
                # file.update({species:ids})

    for species in set(species_list):
        for element in collection.find({"parent": re.compile(species)}, {"_id": 1}):
            ids.append(element['_id'])
        if len(ids) > 0:

            tree = file[species]
            for id in ids:
                tree[id]

    return


def add_all_leafs_host(file, collection, scien_dict, common_dict, include_inner_nodes = True):
    """
    1. go through tree leaves as potential parent parameter
    2. search for parent db and add _id as leaf
    3. return finished tree
    :param file: basic dict, no _ids included
    :param collection: the db entries, use _ids as new leafs
    :return: dict with _ids as leafs
    """
    for k in iter(list(file)):

        if len(list(file[k])) > 0 and not include_inner_nodes:
            add_all_leafs_host(file[k], collection, scien_dict, common_dict)

        else:
            if len(list(file[k])) > 0:
                add_all_leafs_host(file[k], collection, scien_dict, common_dict)

            species_s = scien_dict[k]
            try:
                species_c = common_dict[k]
            except:
                species_c = -1
                pass
            ids = []

            for element in collection.find({"host": species_s}, {"_id": 1}):
                ids.append(element['_id'])

            if species_c != -1:
                for element in collection.find({"host": species_c}, {"_id": 1}):
                    ids.append(element['_id'])

            if len(ids) > 0:
                tree = file[k]
                for id in ids:
                    tree[id]
                # Uncomment if list wanted
                # file.update({species:ids})


    return


def get_db(mongoclient="localhost:27017", db_name="testDB", collection_name="testCell2"):
    """
    get mongoDB collection of interest
    :param mongoclient: which port to speak mongoDB
    :param db_name: name of Database
    :param collection_name: name of collection
    :return: collection
    """
    client = MongoClient(mongoclient)
    db = client[db_name]
    collection = db.get_collection(collection_name)
    collection.create_index([("host", pymongo.DESCENDING)])
    collection.create_index([("parent", pymongo.DESCENDING)])
    return collection


def get_samples(virus_file, host_file, collection, host_clade="all", virus_clade="all", min_samples=100):
    """
    get samples which fulfill the clade conditions
    :param virus_file: dict tree for viruses
    :param host_file: dict tree for hosts
    :param collection: the database
    :param host_clade: which host clade should be used
    :param virus_clade: which virus clade should be used
    :return: samples as pandas array
    """

    # zugriff auf daten durch angabe von Host und Clade
    def get_leaves(file):
        """
        go to leafs of tree and put together in array
        :param file: dict tree representation
        :return: Array[Ids]
        """
        leaf = np.array([])
        for k in iter(list(file)):

            if len(list(file[k])) > 0:
                leaf = np.append(leaf, get_leaves(file[k]))
            else:
                leaf = np.append(leaf, k)

        return leaf

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
                    found = get_ids(file[k], filter)
                    if len(found) > 0:
                        return found
        return []
    print("getting samples")
    if virus_clade == "all":
        virus_leaves = get_leaves(virus_file)
    else:
        virus_leaves = get_ids(virus_file, virus_clade)

    if host_clade == "all":
        # host_leaves = np.array(["ae02c884-187c-492a-a1be-880df33c3f51", "71057421-e62a-4633-b3ba-f2348f77c4bc"])
        host_leaves = get_leaves(host_file)
    else:
        host_leaves = get_ids(host_file,host_clade)
        # host_leaves = np.array(["ae02c884-187c-492a-a1be-880df33c3f51", "71057421-e62a-4633-b3ba-f2348f77c4bc"])
    assert len(virus_leaves) > 0, "no valid virus clade"
    assert len(host_leaves) > 0, "no valid host clade"

    def compare_ids(virus, host):
        ids = set(virus).intersection(host)
        return ids

    def get_table_samples(ids):
        """
        parse ids to usable sample table
        :param ids: ids of interest
        :return: pandas dataframe with id, host, seq and parent of samples
        """
        l = []
        for id in ids:
            for element in collection.find({"_id": id}):
                l.append(element)
        gen = ((i['_id'], i['host'], i['seq'], i['parent']) for i in l)
        df = pd.DataFrame.from_records(gen, columns=['id', 'host', 'seq', 'parent'])

        return df

    ids = compare_ids(virus_leaves, host_leaves)
    df = get_table_samples(ids)

    def get_training_sets(df, min_samples):
        output_samples = []
        output_df = pd.DataFrame()
        for host in df.host.unique():
            df_host = df[df.host == host]
            # print(df_host.host.count())
            if df_host.host.count() >= min_samples:
                samples = df_host.sample(n=min_samples)
                output_samples.append(samples)

        for i in output_samples:
            output_df = output_df.append(i)

        return output_df

    return get_training_sets(df, min_samples)


def save_set(samples, dir):
    """
    save samples as train and test set
    :param samples: dataframe with seq and host
    :param dir: where to save output
    :return: saves csv files for training and test
    """
    X_train, X_test, Y_train, Y_test = train_test_split(samples.seq, samples.host, test_size=0.2, random_state=0)
    X_train.to_csv(dir + '/X_train.csv', sep='\t', encoding='utf-8')
    X_test.to_csv(dir + '/X_test.csv', sep='\t', encoding='utf-8')
    Y_train.to_csv(dir + '/Y_train.csv', sep='\t', encoding='utf-8')
    Y_test.to_csv(dir + '/Y_test.csv', sep='\t', encoding='utf-8')

def get_virus_tree(collection, path_csv= os.getcwd() + "/examples/ICTV_Master_Species_List_2016v1.3.csv", path_output_json= os.getcwd() + "/examples/" + "virus_tree.json"):
    print("getting virus tree")
    virus_tree = csv2dict(path_csv)
    add_all_leaves_virus(virus_tree, collection)
    with open(path_output_json, 'w+') as outjson:
        output = json.dumps(virus_tree)
        outjson.write(output)
    return virus_tree

def get_host_tree(collection, path_ncbi_files= os.getcwd() + "/examples", path_output_json= os.getcwd() + "/examples/" + "host_tree.json"):
    print("getting host tree")
    host_tree, scientific_names, common_names = ncbi_tax2dict(path_ncbi_files)
    print("add leaves to host tree")
    add_all_leafs_host(host_tree, collection, scientific_names, common_names)
    with open(path_output_json, 'w+') as outjson:
        output = json.dumps(host_tree)
        outjson.write(output)

collection = get_db(db_name="FinalDB", collection_name="allData")
virus_tree = get_virus_tree(collection)
exit()
# virus_tree = json.load(open(os.getcwd() + "/examples/" + "virus_tree.json"))
host_tree = json.load(open(os.getcwd() + "/examples/" + "host_tree.json"))

samples = get_samples(virus_tree, host_tree, collection, virus_clade="all", host_clade="all", min_samples=100)
save_set(samples, os.getcwd())

# Todo Caution, use kfold for training.
