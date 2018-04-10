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
from math import gcd
import pickle
import editdistance

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
    phylo_names = {}

    with open(os.path.join(data_dir, 'nodes.dmp')) as nodes_file:
        for line in nodes_file:
            line = line.rstrip(row_delimiter)
            values = line.split(col_delimiter)
            tax_id, parent_id, phylo_rank = values[:3]
            tax_id = int(tax_id)
            parent_id = int(parent_id)
            this_node = nodes[parent_id]
            phylo_names.update({tax_id: str(phylo_rank)})

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
    return dict_tree, scientific_names, common_names, phylo_names


def add_all_leaves_virus(file, collection, include_inner_nodes=False):
    """
    1. go through tree leaves as potential parent parameter
    2. search for parent db and add _id as leaf
    3. return finished tree
    :param file: basic dict, no _ids included
    :param collection: the db entries, use _ids as new leafs
    :return: dict with _ids as leafs
    """
    # species_list = []
    for k in iter(list(file)):

        if len(list(file[k])) > 0 and not include_inner_nodes:
            add_all_leaves_virus(file[k], collection)
        else:
            if len(list(file[k])) > 0:
                add_all_leaves_virus(file[k], collection)

            species = k
            ids = []

            # print(species)

            # sub = species.split(" ")
            # try:
            #     species_list.append(sub[0]+" "+sub[1])
            # except:
            #     species_list.append(sub[0])
            for element in collection.find({"parent": species}, {"_id": 1}):
                # for element in collection.find({"parent": re.compile(sub[0]+" "+sub[1])}, {"_id": 1}):
                # for element in collection.find({"parent": {'$in': [re.compile(species)]}}, {"_id": 1}):
                ids.append(element['_id'])
            if len(ids) > 0:
                tree = file[species]
                for id in ids:
                    tree[id]
                # Uncomment if list wanted
                # file.update({species:ids})
    # return
    # for species in set(species_list):
    #     for element in collection.find({"parent": re.compile(species)}, {"_id": 1}):
    #         ids.append(element['_id'])
    #     if len(ids) > 0:
    #
    #         tree = file[species]
    #         for id in ids:
    #             tree[id]

    return


def pre_host_leafs(collection, scien_dict, common_dict):
    from multiprocessing.dummy import Pool as ThreadPool
    pool = ThreadPool(7)
    # pool.
    id_dict = {}
    id_tuple_arr = []

    def get_entries(item):
        ids = []
        key, value = item[0], item[1]
        try:
            for element in collection.find({"host": {'$in': [re.compile(value)]}}, {"_id": 1}):
                ids.append(element['_id'])
        except:
            print(value)
            ids = []
        try:
            value_common = common_dict[key]
        except:
            value_common = -1
        if value_common != -1:
            for element in collection.find({"host": {'$in': [re.compile(value_common)]}}, {"_id": 1}):
                ids.append(element['_id'])

        return (key, ids)

    # for item in scien_dict.items():
    #     id_tuple_arr.append(get_entries(item))
    # items = [scien_dict.items()]
    list_key_value = [[k, v] for k, v in scien_dict.items()]
    id_tuple_arr = pool.map(get_entries, list_key_value)
    for i in id_tuple_arr:
        id_dict.update({i[0]: i[1]})
    return id_dict


def add_all_leafs_host(file, id_dict, include_inner_nodes=True):
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
            add_all_leafs_host(file[k], id_dict)

        else:
            if len(list(file[k])) > 0:
                add_all_leafs_host(file[k], id_dict)

            # species_s = scien_dict[k]
            # try:
            #     species_c = common_dict[k]
            # except:
            #     species_c = -1
            #     pass
            # ids = []

            # for element in collection.find({"host": {'$in': [re.compile(species_s)]}}, {"_id": 1}):
            # # for element in collection.find({"host": species_s}, {"_id": 1}):
            # # for element in collection.find({"parent": {'$in': [re.compile(species)]}}, {"_id": 1}):
            #     ids.append(element['_id'])
            #
            # if species_c != -1:
            #     for element in collection.find({"host": {'$in': [re.compile(species_c)]}}, {"_id": 1}):
            #     # for element in collection.find({"host": species_c}, {"_id": 1}):
            #         ids.append(element['_id'])
            try:
                ids = id_dict[k]
            except:
                print("missing entry" + str(k))
                ids = []

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


def get_samples(virus_file, host_file, collection, path, host_clade="all", virus_clade="all", min_samples=100,
                del_Clade=[], repeats=False, tax_subgroup=True, use_old_data=False, filter_out_small=False):
    """
    get samples which fulfill the clade conditions
    :param virus_file: dict tree for viruses
    :param host_file: dict tree for hosts
    :param collection: the database
    :param host_clade: which host clade should be used
    :param virus_clade: which virus clade should be used
    :param min_samples: minimal number of samples per host
    :param del_Clade: exclude this clade from train data
    :param repeats: allow multiple use of each sample
    :param tax_subgroup: make different sets for each taxonomial category
    :return: samples as pandas array
    """

    # zugriff auf daten durch angabe von Host und Clade
    def get_leaves(file, parent=""):
        """
        go to leafs of tree and put together in array
        :param file: dict tree representation
        :return: Array[Ids]
        """
        leaf = np.array([])
        for k in iter(list(file)):

            if len(list(file[k])) > 0:
                leaf = np.append(leaf, get_leaves(file[k], k))
            else:
                leaf = np.append(leaf, [k, parent])

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

    def filter_out_non_species_samples(leaves):
        """
        some ids are multiple times in the tree, filters out only ids which belong to a species and no other taxon
        :param leaves: list of ids / leaves to filter
        :return: clean id list
        """
        host_tree, scientific_names, common_names, phylo_names = ncbi_tax2dict(os.getcwd() + "/examples")
        filtered = []
        for i in leaves:
            parent = i[1]
            taxon = phylo_names[int(parent)]
            if taxon == 'species':
                filtered.append(i)
        return np.array(filtered)

    def compare_ids(virus, host):
        ids = set(virus[:, 0]).intersection(host[:, 0])
        virus_sub = [i for i in virus if i[0] in ids]
        host_sub = [i for i in host if i[0] in ids]
        virus_sub_dict = {i[0]: i[1] for i in virus_sub}
        host_sub_dict = {i[0]: i[1] for i in host_sub}

        return ids, virus_sub_dict, host_sub_dict

    def get_table_samples(ids, virus_sub, host_sub):
        """
        parse ids to usable sample table
        :param ids: ids of interest
        :return: pandas dataframe with id, host, seq and parent of samples
        """
        l = []
        # for index in range(len(ids)):
        #     id = ids[index]

        for id in ids:
            for element in collection.find({"_id": id}):
                l.append(element)
        # gen = ((i['_id'], i['host'], i['seq'], i['parent']) for i in l)

        # speichere den parent im Baum als Host bzw Virus, nicht den eintrag in der DB
        gen = ((i['_id'], host_sub[i['_id']], i['seq'], virus_sub[i['_id']]) for i in l)
        df = pd.DataFrame.from_records(gen, columns=['id', 'host', 'seq', 'parent'])

        return df

    def get_training_sets(df, rank, min_samples):
        """
        get samples of host with at least min_samples examples
        :param df: dataframe with raw data
        :param min_samples: per host
        :return: clean training set
        """
        output_samples = []
        output_df = pd.DataFrame()
        for host in df[rank].unique():
            df_host = df[df[rank] == host]
            # print(df_host.host.count())
            if df_host.host.count() >= min_samples:
                samples = df_host.sample(n=min_samples)
                output_samples.append(samples)

        for i in output_samples:
            output_df = output_df.append(i)

        return output_df

    def remove_clade(tree, filter, min_samples, path):
        """
        remove clade from tree entries and return all clade entries for test-set use
        :param tree: tree containing clade
        :param filter: which clade to "delete"
        :return:
        """
        if type(filter) == str:
            clade_leaves = get_ids(tree, filter)
            assert len(clade_leaves) > 0, "no valid deleting clade\n" + did_u_mean(filter)
        else:
            clade_leaves = []
            for clade in filter:
                leaves = get_ids(tree,clade)
                assert len(leaves) > 0, "no valid deleting clade\n" + did_u_mean(clade)
                clade_leaves.append(leaves)
            clade_leaves = np.array(clade_leaves)

        clade_leaves = clade_leaves.reshape((clade_leaves.shape[0] // 2, 2))
        clade_ids, virus_sub, host_sub = compare_ids(clade_leaves, host_leaves)
        df = get_table_samples(clade_ids, virus_sub, host_sub)
        # samples = get_training_sets(df, 1)
        df = add_host_taxonomy(df)
        random_order_df = df.sample(frac=1)
        X_test = random_order_df.seq
        Y_test = random_order_df["species"]

        ids, virus_sub, host_sub = compare_ids(virus_leaves, host_leaves)
        # del test set entries from training set
        ids = set(ids).symmetric_difference(clade_ids)
        df = get_table_samples(ids, virus_sub, host_sub)
        df = add_host_taxonomy(df)
        samples = get_training_sets(df, "species", min_samples)
        random_order_df = samples.sample(frac=1)
        X_train = random_order_df.seq
        Y_train = random_order_df["species"]

        dir = os.getcwd()
        X_test.to_csv(dir +path+ '/X_test.csv', sep='\t', encoding='utf-8')
        X_train.to_csv(dir +path+ '/X_train.csv', sep='\t', encoding='utf-8')
        Y_test.to_csv(dir +path+ '/Y_test.csv', sep='\t', encoding='utf-8')
        Y_train.to_csv(dir +path+ '/Y_train.csv', sep='\t', encoding='utf-8')
        exit()

    def add_host_taxonomy(df):
        """
        1. get dataframe
        2. make new columns (genus,family etc)
        3. add for every id taxonomy
        4. return bigger df
        :param df:
        :return:
        """
        host_tree, scientific_names, common_names, phylo_names = ncbi_tax2dict(os.getcwd() + "/examples")
        host_tree = json.load(open(os.getcwd() + "/examples/" + "host_tree_big.json"))

        columns = set(phylo_names.values())
        for col in columns:
            df[col] = pd.Series(index=df.index)

        for host in df.host.unique():
            # name = str(scientific_names[int(host)])
            # print(name)
            taxonomy = get_tax_path(host, host_tree)

            for taxon in taxonomy:
                # df[genus]= homo for all entries with host
                df.loc[df.host == host, phylo_names[int(taxon)]] = scientific_names[int(taxon)]

        return df

    if not use_old_data:
        print("getting samples")
        if virus_clade == "all":
            virus_leaves = get_leaves(virus_file)
        else:
            virus_leaves = get_ids(virus_file, virus_clade)

        # make shape [[id,parent],[id,parent]...]
        assert len(virus_leaves) > 0, "no valid virus clade\n"+did_u_mean(virus_clade)
        virus_leaves = virus_leaves.reshape((virus_leaves.shape[0] // 2, 2))
        print(f"number of all virus samples {len(set(virus_leaves[:, 0]))}")
        # virus_leaves = filter_out_non_species_samples(virus_leaves)
        # print(len(set(virus_leaves[:,0])))

        # f = open(os.getcwd() + '/host_leaves.pkl', 'rb')
        # host_leaves = pickle.load(f)

        if host_clade == "all":
            host_leaves = get_leaves(host_file)
        else:
            host_leaves = get_ids(host_file, host_clade)
        host_leaves = host_leaves.reshape((host_leaves.shape[0] // 2, 2))

        # print(len(set(host_leaves[:, 0])))
        assert len(host_leaves) > 0, "no valid host clade\n"
        host_leaves = filter_out_non_species_samples(host_leaves)
        print(f"number of all host samples {len(set(host_leaves[:, 0]))}")

        if del_Clade != []:
            print(f"deleting {del_Clade}")
            remove_clade(virus_file, del_Clade, min_samples, path)

        # speichere lang zu berechnendes object

        # with open(str(os.getcwd())+'/host_leaves.pkl', 'wb') as f:
        #     pickle.dump(host_leaves, f)
        #
        # with open(str(os.getcwd())+'/virus_leaves.pkl', 'wb') as f:
        #     pickle.dump(virus_leaves, f)


        # f = open(os.getcwd() + '/virus_leaves.pkl', 'rb')
        # virus_leaves = pickle.load(f)
        #
        ids, virus_sub, host_sub = compare_ids(virus_leaves, host_leaves)
        print(f"number of useable samples: {len(ids)}")

        df = get_table_samples(ids, virus_sub, host_sub)

        df = add_host_taxonomy(df)
        df.to_csv(os.getcwd() +path+ "/sample_df.csv")
        df.to_pickle(os.getcwd() + path + "/sample_df")

    else:
        df = pd.read_pickle(os.getcwd() + "/sample_df")

    if filter_out_small:
        # df = df[len(df.seq)
        sizes = [len(i) for i in df.seq]
        sizes = np.array(sizes)
        cutoff = np.percentile(sizes,5)
        print(f"cutting out sequences smaller {cutoff}")
        df = df[df.seq.str.len() > cutoff]

    if repeats:
        if tax_subgroup:
            for index, rank in enumerate(["species", "genus", "family", "order", "class", "phylum"]):
                print(rank)
                get_training_set_with_repeats(df, rank, min_samples * (index + 1), path=path)
        else:
            get_training_set_with_repeats(df, "species", min_samples, path=path)
    else:
        if tax_subgroup:
            for index, rank in enumerate(["species", "genus", "family", "order", "class", "phylum"]):
                samples = get_training_sets(df, rank, min_samples * (index + 1))
                save_set(samples, rank, os.getcwd() + path)
        else:
            samples = get_training_sets(df, "species", min_samples)
            save_set(samples, "species", os.getcwd() + path)


def save_set(samples, rank, dir):
    """
    save samples as train and test set
    :param samples: dataframe with seq and host
    :param dir: where to save output
    :return: saves csv files for training and test
    """
    X_train, X_test, Y_train, Y_test = train_test_split(samples.seq, samples[rank], test_size=0.2, random_state=0,
                                                        stratify=samples[rank])
    # try:
    #     print(Y_test[1].unique())
    # except:
    #     print(Y_test.host.unique())
    X_train.to_csv(dir + '/X_train_' + str(rank) + '.csv', sep='\t', encoding='utf-8')
    X_test.to_csv(dir + '/X_test_' + str(rank) + '.csv', sep='\t', encoding='utf-8')
    Y_train.to_csv(dir + '/Y_train_' + str(rank) + '.csv', sep='\t', encoding='utf-8')
    Y_test.to_csv(dir + '/Y_test_' + str(rank) + '.csv', sep='\t', encoding='utf-8')


def get_virus_tree(collection, path_csv=os.getcwd() + "/examples/ICTV_Master_Species_List_2016v1.3.csv",
                   path_output_json=os.getcwd() + "/examples/" + "virus_tree.json"):
    print("getting virus tree")
    virus_tree = csv2dict(path_csv)
    add_all_leaves_virus(virus_tree, collection)
    with open(path_output_json, 'w+') as outjson:
        output = json.dumps(virus_tree)
        outjson.write(output)
    return virus_tree


def get_host_tree(collection, path_ncbi_files=os.getcwd() + "/examples",
                  path_output_json=os.getcwd() + "/examples/" + "host_tree.json"):
    print("getting host tree")
    host_tree, scientific_names, common_names, phylo_names = ncbi_tax2dict(path_ncbi_files)
    print("add leaves to host tree")
    id_dict = pre_host_leafs(collection, scientific_names, common_names)
    add_all_leafs_host(host_tree, id_dict)
    with open(path_output_json, 'w+') as outjson:
        output = json.dumps(host_tree)
        outjson.write(output)


def get_training_sets(df, min_samples):
    """

    :param df:
    :param min_samples: per host
    :return:
    """
    output_samples = []
    output_df = pd.DataFrame()
    for host in df.host.unique():
        df_host = df[df.host == host]
        # print(df_host.host.count())
        # if df_host.host.count() >= min_samples:
        samples = df_host.sample(n=min_samples)
        output_samples.append(samples)

    for i in output_samples:
        output_df = output_df.append(i)

    return output_df


def get_training_set_with_repeats(df, rank, min_samples, path, test_size=0.2, val_size=0.2):
    """
    get specific number of samples per host
    take from every virus which is infecting the host #Samples_needed/#virus-species_in_host
    if not possible, draw multiple times
    CAUTION: its not allowed to have an identical sample in training - val -test set

    1. selected host
    2. calc how many samples needed per virus
    3. draw samples from each virus
    4. append samples to sets
    5. go to new host, till all handled
    6. save sets
    :param df: input dataframe, with possible samples
    :param min_samples: per host
    :param test_size: size of the test set
    :param val_size: size of the validation set
    :return:
    """

    test = pd.DataFrame()
    val = pd.DataFrame()
    train = pd.DataFrame()
    foo = gcd(int(val_size * 100), int(test_size * 100))
    # logical min number samples per virus-host to divide fair
    logical_max_num_viruses = gcd(foo, int((1 - val_size - test_size) * 100))
    logical_min_samples_virus = min_samples // logical_max_num_viruses
    val_size = val_size * (1 / (1 - test_size))

    for host in df[rank].unique():
        df_host = df[df[rank] == host]

        # Note: True == 1, False == 0
        num_useable_viruses = sum(df_host.parent.value_counts() >= 3)

        # because min 3 samples per virus
        # need be caution not to take too much samples, if host has lots of viruses
        if num_useable_viruses > logical_max_num_viruses:
            num_useable_viruses = logical_max_num_viruses

        if num_useable_viruses == 0:
            continue
        min_samples_virus = min_samples // num_useable_viruses
        num_used_viruses = 0

        for virus in df_host.parent.unique():
            if num_used_viruses == num_useable_viruses:
                break

            df_virus = df_host[df_host.parent == virus]

            # if enough samples take #min_sample_virus
            # split into train, val, test set
            if df_virus.parent.count() >= min_samples_virus:
                samples = df_virus.sample(n=min_samples_virus)
                samples_test = samples.sample(frac=test_size)
                samples_val = samples.drop(samples_test.index).sample(frac=val_size)
                samples_train = samples.drop(samples_test.index)
                samples_train = samples_train.drop(samples_val.index)

                test = test.append(samples_test)
                val = val.append(samples_val)
                train = train.append(samples_train)

                num_used_viruses += 1

            # if not enough samples take all available
            # than take repeatably more samples till #min_samples_virus is reached
            elif df_virus.parent.count() < min_samples_virus and df_virus.parent.count() >= 3:
                n = int(len(df_virus) * test_size)
                samples_test = df_virus.sample(n=n)  # frac=test_size)
                if n == 0:
                    samples_test = df_virus.sample(n=1)
                n = int(len(df_virus.drop(samples_test.index)) * val_size)
                samples_val = df_virus.drop(samples_test.index).sample(n=n)  # frac=val_size)
                if n == 0:
                    samples_val = df_virus.drop(samples_test.index).sample(n=1)

                # have to remove the test and val samples to get training set
                samples_train = df_virus.drop(samples_test.index)
                samples_train = samples_train.drop(samples_val.index).sample(frac=1)

                num_samples_test = int(min_samples_virus * test_size)
                num_samples_val = int(min_samples_virus * (1 - test_size) * val_size)
                num_samples_train = min_samples_virus - num_samples_test - num_samples_val
                samples_test = samples_test.append(
                    samples_test.sample(n=num_samples_test - len(samples_test), replace=True))
                samples_val = samples_val.append(samples_val.sample(n=num_samples_val - len(samples_val), replace=True))
                samples_train = samples_train.append(
                    samples_train.sample(n=num_samples_train - len(samples_train), replace=True))

                test = test.append(samples_test)
                val = val.append(samples_val)
                train = train.append(samples_train)

                num_used_viruses += 1

            else:
                continue

    def save_3_sets(rank, path):
        test_shuffled = test.sample(frac=1)
        val_shuffled = val.sample(frac=1)
        train_shuffled = train.sample(frac=1)

        X_test = test_shuffled.seq
        Y_test = test_shuffled[rank]
        X_val = val_shuffled.seq
        Y_val = val_shuffled[rank]
        X_train = train_shuffled.seq
        Y_train = train_shuffled[rank]

        dir = os.getcwd() + path

        X_test.to_csv(dir + '/X_test_' + str(rank) + '.csv', sep='\t', encoding='utf-8')
        X_val.to_csv(dir + '/X_val_' + str(rank) + '.csv', sep='\t', encoding='utf-8')
        X_train.to_csv(dir + '/X_train_' + str(rank) + '.csv', sep='\t', encoding='utf-8')

        Y_test.to_csv(dir + '/Y_test_' + str(rank) + '.csv', sep='\t', encoding='utf-8')
        Y_val.to_csv(dir + '/Y_val_' + str(rank) + '.csv', sep='\t', encoding='utf-8')
        Y_train.to_csv(dir + '/Y_train_' + str(rank) + '.csv', sep='\t', encoding='utf-8')

    save_3_sets(rank, path=path)
    # exit()


def get_tax_path(species, file):
    for k in iter(file.keys()):
        if k == species:
            tax = []
            tax.insert(0, k)
            return tax
        else:
            if len(file[k].keys()) > 0:
                found = get_tax_path(species, file[k])
                if len(found) > 0:
                    found.insert(0, k)
                    return found
    return []


def virus_tree_pairs():
    """
    1. load all ids
    2. make df
    3. count & memories combinations
    4. change name to taxonomy
    5. save list
    """

    # 1 & 2
    df = pd.read_pickle(os.getcwd() + "/sample_df")
    foo = df.groupby(['host', 'parent']).count().reset_index()
    foo = foo[['host', 'parent', 'id']].rename(index=str, columns={"parent": "virus", "id": "number"})
    # print(foo)
    dict_tree, scientific_names, common_names, phylo_names = ncbi_tax2dict(os.getcwd() + "/examples")
    inv_map_scien = {v: k for k, v in scientific_names.items()}
    inv_map_com = {v: k for k, v in common_names.items()}

    # i = 0
    virus_tree = json.load(open(os.getcwd() + "/examples/" + "virus_tree.json"))

    print("parse viruses")

    for virus in foo.virus.unique():
        taxonomy = get_tax_path(virus, virus_tree)
        phyl = ""
        for taxon in taxonomy:
            phyl = phyl + "/" + str(taxon).replace(' ', '_')
        foo.virus[foo.virus == virus] = phyl

        # print(phyl)

    print(foo)

    print("parse hosts")

    host_tree = json.load(open(os.getcwd() + "/examples/" + "host_tree.json"))
    for host in foo.host.unique():
        try:
            name = str(inv_map_scien[host])
        except:
            name = str(inv_map_com[host])
        taxonomy = get_tax_path(name, host_tree)
        phyl = ""
        for taxon in taxonomy:
            phyl = phyl + "/" + scientific_names[int(taxon)].replace(' ', '_')

        foo.host[foo.host == host] = phyl

    print(foo)
    foo.to_csv(os.getcwd() + '/host_virus_pairs.csv', sep='\t', encoding='utf-8')


def pca_prep():
    df = pd.read_pickle(os.getcwd() + "/sample_df")
    X_test = df.seq
    Y_test = df.host

    dict_tree, scientific_names, common_names, phylo_names = ncbi_tax2dict(os.getcwd() + "/examples")
    inv_map_scien = {v: k for k, v in scientific_names.items()}
    inv_map_com = {v: k for k, v in common_names.items()}
    host_tree = json.load(open(os.getcwd() + "/examples/" + "host_tree.json"))

    for host in df.host.unique():
        # try:
        #     name = str(inv_map_scien[host])
        # except:
        #     name = str(inv_map_com[host])
        taxonomy = get_tax_path(host, host_tree)
        phyl = ""
        for taxon in taxonomy:
            phyl = phyl + "/" + scientific_names[int(taxon)].replace(' ', '_')

        Y_test[Y_test == host] = phyl

    X_test.to_csv(os.getcwd() + '/PCA_X_test.csv', sep='\t', encoding='utf-8')
    Y_test.to_csv(os.getcwd() + '/PCA_Y_test.csv', sep='\t', encoding='utf-8')

def did_u_mean(searchword):
    # virus_tree
    catalogue = []
    def all_viruses(file):
        for k in iter(file.keys()):
            catalogue.append(k)
            if len(file[k].keys()) > 0:
                all_viruses(file[k])

    # catalogue = scientific_names.values()
    all_viruses(virus_tree)
    min_score = 10000000
    candidate = []
    for word in catalogue:
        score = editdistance.eval(word,searchword)
        if score < min_score:
            min_score = score
            candidate = [word]
        elif score == min_score:
            candidate.append(word)
        else:
            continue

    string=""
    for word in candidate:
        string += f"Did you mean: {word}?\n"
    return string

collection = get_db(db_name="FinalDB", collection_name="allData")
# host_tree = get_host_tree(collection,path_output_json="/home/go96bix/projects/nanocomb/nanocomb/examples/host_tree_big.json")
# # virus_tree = get_virus_tree(collection)
# # exit()
# virus_tree = json.load(open(os.getcwd() + "/examples/" + "virus_tree.json"))
virus_tree = json.load(open(os.getcwd() + "/examples/" + "virus_tree_big.json"))
# host_tree = json.load(open(os.getcwd() + "/examples/" + "host_tree.json"))
host_tree = json.load(open(os.getcwd() + "/examples/" + "host_tree_big.json"))
# # exit()
# # ebola= 1570291
# # marburgvirus = 11269
# df = pd.read_pickle(os.getcwd()+"/sample_df")
# get_training_set_with_repeats(df,min_samples=100)
get_samples(virus_tree, host_tree, collection, path="/100Samples_noShortSeqs", virus_clade="all", host_clade="all",
            min_samples=100,repeats=True, tax_subgroup=True, use_old_data=True, filter_out_small=True)
#del_Clade=["Influenzavirus A","Influenzavirus B","Influenzavirus C","Influenzavirus D"])
# save_set(samples, os.getcwd())
# virus_tree_pairs()
# pca_prep()
