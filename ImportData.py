from xmljson import badgerfish as bf
from xml.etree.ElementTree import fromstring
from xml.etree.ElementTree import Element, tostring
from json import dumps

def xml_to_json():

    with open("/home/go96bix/projects/nanocomb/sampleData/HQ186308.xml","r") as file:
        data = file.read()
        print(data)
        output = bf.data(fromstring(data))
        dumps(bf.data(output))


xml_to_json()