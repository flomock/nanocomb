# stackoverflow, 26146623
import csv
from collections import defaultdict
# import itertools
from pprint import pprint


def tree(): return defaultdict(tree)


def tree_add(t, path):
  for node in path:
    t = t[node]


def pprint_tree(tree_instance):
    def dicts(t): return {k: dicts(t[k]) for k in t}
    pprint(dicts(tree_instance))

# def csv_to_tree(input):
#     t = tree()
#     for row in csv.reader(input, quotechar='\''):
#         tree_add(t, row)
#     return t


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


def csv2newick(fp, out):
    '''
    Tested on ICTV virus taxonomy.

    from io import StringIO
    from Bio import Phylo


    with open('tax.nwk', 'w+') as out:
        out.write(tree_to_newick(t) + '\n')

    tree = Phylo.read(StringIO(csv2newick(fp)), 'newick')
    '''

    t = tree()
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
    return tree_to_newick(t)
