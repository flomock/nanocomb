# nanocomb
How to use:

1. get the data
-use download.py to get all data from embl, needs accession numbers as input

2. init database
-save output of step 1 in one folder eg. Samples
-use load.files_to_json(your directory or single file)
-this saves a json and returns a jsos like object
-than use this object to initialise the database, with load.init(object, mongo-DB, cell)

NOTE: you need to run a mongo server on your machine for the init method

3. get the host- and the viral- taxonomy tree
-from now use tree.py
-get the db
-parse the ICTV csv to virus_tree with csv2dict
-add sample Ids as leaves in tree (add_all_leaves_virus())
-parse the names.dmp and nodes.dmp to host_tree (ncbi_tax2dict())
-add sample Ids as leaves in host tree

4. select your samples of interest and get test-, train- set
-use get_samples with the clades and min samples per host defined
-save as train, test set with save_set
