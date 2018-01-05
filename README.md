# nanocomb
How to use:

save files of interest in one folder eg. examples

use load.files_to_json(your directory or single file)

this saves a json and returns a jsos like object

than use this object to initialise the database, with load.init(object, mongo-DB, cell)

NOTE: you need to run a mongo server on your machine for the init method