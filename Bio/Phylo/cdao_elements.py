import xml.etree.ElementTree as ET
import os

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

cdao_elements = {}

with open(os.path.join(__location__, 'cdao.owl')) as owl_file:
    tree = ET.parse(owl_file)
    root = tree.getroot()
    for node_type in 'ObjectProperty', 'Class', 'DatatypeProperty':
        for element in root.findall('{http://www.w3.org/2002/07/owl#}%s' % node_type):
            obo = element.attrib['{http://www.w3.org/1999/02/22-rdf-syntax-ns#}about'].split('/')[-1]
            cdao = element.find('{http://www.w3.org/2000/01/rdf-schema#}label').text
            cdao_elements[cdao] = obo