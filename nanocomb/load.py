from xmljson import badgerfish as bf
from xml.etree.ElementTree import fromstring
from json import dumps
import os

def xml_to_json():
    cwd = os.getcwd()
    with open(cwd + "/examples/HQ186308.xml","r") as file:
        data = file.read()
        output = bf.data(fromstring(data))
        print(dumps(output))


xml_to_json()