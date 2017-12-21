from xmljson import badgerfish as bf
from xml.etree.ElementTree import fromstring
from json import dumps
import os
from pathlib import Path

def xml_to_json(path):
    """
    load xml file from path and parse to json
    :param path:
    :return:
    """
    file = Path(path)
    if file.is_file():
        with open(file, "r") as file:
            data = file.read()
            output = bf.data(fromstring(data))
            return dumps(output)
    else:
        raise KeyError('no valid file path')


cwd = os.getcwd()
xml_to_json(cwd + "/examples/HQ186308.xml")