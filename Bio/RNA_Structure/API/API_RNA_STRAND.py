# -*- coding: utf-8 -*-

# Copyright 2017 by Joanna Zbijewska, Agata Gruszczyńska, Michał Karlicki.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
# Note that:
# 1. RMSD scripts and license is available separately.
#    We added it in file: calculate_rmsd and README2.
#
# 2. RNAstructure license is available separately.
#    Please consult rna.urmc.rochester.edu .

"""
__author__ = "Joanna Zbijewska, Gruszka"
"""

import re
import itertools
import urllib
import requests
from lxml import html
from prettytable import PrettyTable
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO

def split_list(to_split, size):
    it = iter(to_split)
    item = list(itertools.islice(it, size))
    while item:
        yield item
        item = list(itertools.islice(it, size))

class RNA_STRAND():
    """
    This class is use to search secondary structure of RNA in online RNA Strand database.

    The input must be string consist of sequence without whitespaces, commas, enters or another.
    """
    _urlrna = "http://www.rnasoft.ca/strand/"

    def __init__(self, sequence, path=None):
        self.path = path
        if path is None:
            self.path = ""
        """Initializer of the type sequence argument for futher functions."""
        self.sequence = sequence

    def download_database(self):
         #Co tutaj?
        path = self.path
        """RNA Strand database downloader for users who think they may need it all."""
        urlrna = "http://www.rnasoft.ca/strand/"
        url = urlrna+"{path}RNA_STRAND_data.tar.gz".format(path=path)
        urllib.urlretrieve(url, "RNA_STRAND_data.tar.gz")
        return ("RNA STRAND database downloaded to given localizaton.")



    def search_by_sequence(self):
        """
        RNA Strand database search for desired sequence.
        Give a list consist of five-element lists.
        These elements are respectively name, type, source organism, length and ID of molecules include enter sequence.
        """
        urlrna = "http://www.rnasoft.ca/strand/"
        search = urlrna+"search_results.php?select%5B%5D=Any+type&org=&source%5B%5D=Any+source&source_id=&select2=Any+length&first1=&last1=&exp_proven=Any&select4=Any+number&first2=&last2=&select5=Any&select_duplicate=All+molecules&seq={}&abstractshapetype=ABSTRACT_SHAPE_5&abstractshape=&Submit=Perform+search&select6=Any+number&first3=&last3=&select7=Any+number&first4=&last4=&select8=Any+number&first5=&last5=&select12=Any+number&first7=&last7=&select13=Any+number&first8=&last8=&motif=&select19=Any+number&first10=&last10=&select20=Any+number&first11=&last11=&select25=Any+number&first13=&last13=&select26=Any+number&first14=&last14=&select27=absolute&select28=Any+number&first15=&last15=&select33=Any+number&first17=&last17=&select34=Any+number&first18=&last18=&select35=average&select36=absolute&select37=Any+number&first19=&last19=&select38=Any+number&first20=&last20=&select69=Any+number&first31=&last31=&select71=bands&select72=Any+number&first32=&last32=&select73=base+pairs&select74=Any+number&first33=&last33=&select75=Any&select76=Any+number&first34=&last34=&select55=any&select_ncbp_context_2=any&select_ncbp_context_4=any&select_ncbp_context_1=any&select_ncbp_context_5=any&select56=Any+number&first26=&last26=&start=0&limit=100&sort_by=length&order=ascending".format(self.sequence)
        #Domyślne ustawienia dla tego wyszukiwania - 100 cząsteczek i ułożone wedle długości sekwencji rosnąco
        found = requests.get(search)
        output = html.fromstring(found.content)
        some_data = output.xpath('//tr/td/a/text()')
        pattern = re.compile(r'(_0)')
        IDs = [y for y in some_data if pattern.search(y)!=None]
        output_list = output.xpath('//tr/td/text()')
        regex1 = re.compile(r'(Any)')
        regex2 = re.compile(r'(\n)')
        regex3 = re.compile(r'(\t)')
        regex4 = re.compile(r'\[')
        regex5 = re.compile(r'(In all)')
        regex6 = re.compile(r'\ ')
        regex7 = re.compile(r'(Remove outliers)')
        regex8 = re.compile(r'(Normalize by RNA type)')
        out = [regex6,regex7,regex8]
        result1=[]
        for i in output_list:
            if regex1.search(i)==None:
                result1.append(i)
        result2=[]
        for i in result1:
            if regex2.search(i)==None:
                result2.append(i)
        result3=[]
        for i in result2:
            if regex3.search(i)==None:
                result3.append(i)
        result4=[]
        for i in result3:
            if regex4.search(i)==None:
                result4.append(i)
        result5=[]
        for i in result4:
            if regex5.search(i)==None:
                result5.append(i)
        result = [x for x in result5 if not any (regex.match(x) for regex in out)]
        our_result = list(split_list(result,4))
        for num in range(len(our_result)):
            our_result[num].append(IDs[num])
        return(our_result)

    def print_results(self):
        """Printer of a table consist of ordinal number, name, type,
        source organism, length and ID of molecules include enter sequence.
        """
        our_result = self.search_by_sequence()
        t = PrettyTable(["No", "Name", "Type", "Organism", "Length", "ID"])
        for num in range(1,1+len(our_result)):
            x = our_result[num-1]
            t.add_row([num, x[0], x[1], x[2], x[3], x[4]])
        print(t)

    def choose_result(self):
        """Interaction with user.

        Returns RNA Strand database ID of molecule chosen by user.
        Choose is based on data about molecule include enter sequence (look print_result())

        """
        self.print_results()
        our_result = self.search_by_sequence()
        choose_molecule = input('Which one is your molecule? Choose number: ')
        choose_molecule = int(choose_molecule)
        mo = our_result[choose_molecule-1]
        id_mo = mo[4]
        return(id_mo)


    def metadata(self):
        """Metadata download for chosen sequence."""
        results = self.search_by_sequence()
        chosen_id = self.choose_result()
        for i in range(len(results)):
            l = results[i]
            if l[4] == chosen_id:
                meta_list = l[1:]
        urlrna = "http://www.rnasoft.ca/strand/"
        res = urlrna+"show_results.php?molecule_ID={}".format(chosen_id)
        content = requests.get(res)
        content = html.fromstring(content.content)
        table = content.xpath('//tr/td/text()')
        regex1 = re.compile(r'\t')
        regex2 = re.compile(r'\[')
        regex3 = re.compile(r'\ ')
        regex4 = re.compile(r'\]')
        table1 = [elem for elem in table if elem != '\n']
        table2 = [elem for elem in table1 if elem != ']:']
        table3 = [elem for elem in table2 if regex1.match(elem) == None]
        table4 = [elem for elem in table3 if regex2.match(elem) == None]
        table5 = [elem for elem in table4 if regex3.match(elem) == None]
        table6 = [elem for elem in table5 if regex4.search(elem) == None]
        data = table6[3:]
        data2 = []
        for string in data:
            string1 = string.lstrip('\n\t')
            string2 = string1.rstrip(' [')
            data2.append(string2)
        dat = [x for x in data2 if x != '']
        meta = [y for y in dat if re.search('Source',y) == None]
        return(meta)

    def metadata_to_file(self):
        meta = self.metadata()
        path = self.path
        ID = meta[1]
        metadata = list(split_list(meta,2))
        metadata = [" : ".join(metadata[i]) for i in range(len(metadata))]
        metadata = [dat + "\n" for dat in metadata]
        with open("{path}report_{}".format(ID,path=path), "w") as metafile:
            metafile.writelines(metadata)
        metafile.close()
        return("Metadata report is ready")

    def get_sequence(self):
        """Returns respectively chosen molecule's (look choose_result and print_results) ID and sequence in a list."""
        id_mo = self.choose_result()
        urlrna = "http://www.rnasoft.ca/strand/"
        res = urlrna+"show_file.php?format=FASTA&molecule_ID={}&check_out_the=View+the+RNA+sequence+and+secondary+structure+for+molecule+{}".format(id_mo, id_mo)
        content = requests.get(res)
        content = html.fromstring(content.content)
        raw_seq = content.xpath('//div/textarea/text()')
        letters = '([AUGC])\w+'
        raw_seq = raw_seq[0].split('\n')
        for y in raw_seq:
            match = re.match(letters,y)
            if match:
                seq = match.group(0)
        molecule = SeqRecord(Seq(seq,IUPAC.ambiguous_rna), id=id_mo, name = "RNA sequence", description="None")
        return(molecule)

    def download_fasta_sequence(self):
        """Saver chosen molecule's (look choose_result and print_results) sequence in fasta format.
        """
        to_save = self.get_sequence()
        id_mo = self.choose_result()
        path = self.path
        with open("{path}{}_sequence.fasta".format(id_mo,path=path),"w") as f:
                SeqIO.write(to_save, f, "fasta")
        return("Fasta file is ready.")

    def get_structure(self):
        """Returns chosen molecule's (look choose_result and print_results) ID and structure in a list.
        """
        id_mo = self.choose_result()
        urlrna = "http://www.rnasoft.ca/strand/"
        res = urlrna+"show_file.php?format=Bpseq&molecule_ID={}&check_out_the=View+the+RNA+sequence+and+secondary+structure+for+molecule+{}".format(id_mo,id_mo)
        content = requests.get(res)
        content = html.fromstring(content.content)
        raw_struct = content.xpath('//div/textarea/text()')
        raw_struct = raw_struct[0].split('\n')
        ID = raw_struct[1].strip('# File')
        ID = ID.strip('.ct')
        re_hash = re.compile(r'\#')
        raw_structure = [elem for elem in raw_struct if not re_hash.search(elem)]
        raw_structure = raw_structure[1:]
        raw_structure.insert(0,ID)
        return(raw_structure)

    def download_bpseq_structure(self):
        """Saver of chosen molecule's (look choose_result and print_results) structure in bpseq format.
        """
        struct_to_save = self.get_structure()
        with open('{path}{}_structure.bpseq'.format(struct_to_save[0],path=self.path), 'w') as bpseq:
            for ind in range(1,len(struct_to_save)):
                bpseq.write(struct_to_save[ind]+'\n')
        bpseq.close()
        return("Bpseq file is ready.")


#a = RNA_STRAND('UAAGCCCUA',path="/Users/michalkarlicki/Desktop/")
#a.metadata_to_file()
