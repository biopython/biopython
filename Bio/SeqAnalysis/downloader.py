import os
import re
import urllib
from multiprocessing.pool import ThreadPool
from os.path import isfile, join
from time import time as timer

from Bio import Entrez
from Bio import SeqIO


class SeqDownloader:
    def __init__(self, full_list, directory):
        if type(full_list) == str:
            clear_ids = [full_list]
        elif type(full_list) == list:
            clear_ids = full_list
        analysis_dir = directory

        print("Number of unique IDs: {0}".format(len(clear_ids)))
        print("Name of the current batch of files: {0}".format(analysis_dir))

        print("------------------------------------Directory check------------------------------------------")

        dirs = ['Downloads', 'Failed', 'Database']
        for single_dir in dirs:
            self.mk_dirs(single_dir)
        self.mk_subdirs(analysis_dir)

        obsolete = []
        http_error = []

        print("--------------------------------------Downloading--------------------------------------------")

        urls = self.paths_and_urls(clear_ids, analysis_dir)

        print("Downloading files. Please wait...")

        ThreadPool(20).imap_unordered(self.fetch_url, urls)
        start = timer()
        for entry in urls:
            self.fetch_url(entry, obsolete, http_error)

        Entrez.email = ''
        for idx in clear_ids:
            if re.match(r'([A-Z]){2,3}\w+', idx):
                self.download_ncbi(idx, analysis_dir, obsolete, http_error)

        obsolete = list(dict.fromkeys(obsolete))

        mypath = ('./Downloads/{0}'.format(analysis_dir))
        onlyfiles = [f for f in os.listdir(mypath) if isfile(join(mypath, f))]
        files_no_ext = [".".join(f.split(".")[:-1]) for f in os.listdir(mypath) if os.path.isfile(join(mypath, f))]

        total_number = len(onlyfiles)
        print("Download of " + str(total_number) + f" took {timer() - start} seconds")

        print("-------------------------------------Saving outputs------------------------------------------")
        print('Possibly obsolete entries were saved in Failed/{0}/Obsolete.txt'.format(analysis_dir))
        print('Number of obsolete entries: {0}'.format(str(len(obsolete))))

        with open('./Failed/{0}/Obsolete.txt'.format(analysis_dir), 'w+') as obs:
            for item in obsolete:
                obs.write(item + '\n')

        print('Entries with possible download error were saved in Failed/{0}/Error.txt'.format(analysis_dir))
        print('Number of download error entries: {0}'.format(str(len(http_error))))

        with open('./Failed/{0}/Error.txt'.format(analysis_dir), 'w+') as fail:
            for item in http_error:
                fail.write(item + '\n')

        print("Creating database. Please wait...")
        for file in onlyfiles:
            with open('./Downloads/{0}/{1}'.format(analysis_dir, file)) as big_fasta:
                one_fasta = big_fasta.read()
                with open('./Database/{0}.fasta'.format(analysis_dir), 'w+') as data:
                    data.write(one_fasta + '\n')

        print("Your sequences were saved in Database/{0}.fasta".format(analysis_dir))
        print("Done.")
        self.downloaded_files = files_no_ext

    def read_list_of_ids(self, file) -> list:
        try:
            data = ''
            with open(file, 'r') as f:
                for line in f.readlines():
                    line = re.sub('[!@/#$|:]', ' ', line).strip()
                    data = data + ' ' + line
                data = data.split()
                data = list(dict.fromkeys(data))
            return data
        except FileNotFoundError:
            raise Exception("!!!ERROR!!! \n"
                            "Please enter right name of file (with .txt).\n"
                            "Remember that your ID file should be in the same directory as this script.\n"
                            "If not please specify path to the ID file.\n"
                            "E.g. dir1/subdir1/list.txt")

    def mk_dirs(self, dirName):
        if not os.path.exists(dirName):
            os.mkdir(dirName)
            print("Directory {0} status: created".format(dirName))
        else:
            print("Directory {0} status: already exists".format(dirName))

    def mk_subdirs(self, dirName):
        if not os.path.exists('./Downloads/{0}/'.format(dirName)):
            os.mkdir('./Downloads/{0}/'.format(dirName))
            print("Directory Downloads/{0} status: created".format(dirName))
        else:
            print("Directory Downloads/{0} status: already exists".format(dirName))
        if not os.path.exists('./Failed/{0}/'.format(dirName)):
            os.mkdir('./Failed/{0}/'.format(dirName))
            print("Directory Failed/{0} status: created".format(dirName))
        else:
            print("Directory Failed/{0} status: already exists".format(dirName))

    def paths_and_urls(self, ids, cdir) -> list:
        pu = []
        for idx in ids:
            if re.match(r'([A-Z]){1}\d+', idx):
                path_url = tuple([('./Downloads/{0}/{1}.fasta'.format(cdir, idx)),
                                  ('https://www.uniprot.org/uniprot/' + idx + '.fasta')])
                pu.append(path_url)

        return pu

    def fetch_url(self, entry, obs, http) -> str:
        path, url = entry
        uniprot_id = re.sub('[/.]', ' ', path).strip().split()
        if not os.path.exists(path):
            try:
                fasta = urllib.request.urlopen(url)
                fasta = fasta.read().decode()

                if fasta == '':
                    obs.append(uniprot_id[2])
                else:
                    with open(path, 'w+') as f:
                        f.write(fasta)

            except urllib.error.HTTPError or TimeoutError or urllib.error.URLError:
                http.append(uniprot_id[2])
        return path

    def download_ncbi(self, ncbi_id, cdir, obs, http):
        try:
            handle = Entrez.efetch(db="protein", id=ncbi_id, rettype="fasta", retmode="text")
            record = handle.read()

            if record == '':
                obs.append(ncbi_id)
            else:
                with open('./Downloads/{0}/{1}.fasta'.format(cdir, ncbi_id), 'w+') as f:
                    f.write(record)

        except urllib.error.HTTPError or TimeoutError or urllib.error.URLError:
            http.append(ncbi_id)


class SequenceDatabase:
    def __init__(self, full_list, directory):
        download = SeqDownloader(full_list, directory)
        seqData = {}
        for index, record in zip(download.downloaded_files,
                                 SeqIO.parse("C:/Users/solak/Desktop/Docs/MGR/APB/Database/directory.fasta", "fasta")):
            if index not in seqData.keys():
                seqData[index] = record
        self.database = seqData

    def get(self, indexes=None) -> list:
        if indexes is None:
            records = list(self.database.values())
        else:
            if type(indexes) == str:
                keys = [indexes]
            elif type(indexes) == list:
                keys = indexes
            records = [self.database.get(key) for key in keys]
        return records
