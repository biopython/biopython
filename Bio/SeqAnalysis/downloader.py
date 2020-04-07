# Unique Sequence Downloader
import os
import re
import sys
import urllib
from multiprocessing.pool import ThreadPool
from os.path import isfile, join
from time import time as timer

from Bio import Entrez

try:
    filename = sys.argv[1]
    current_dir = sys.argv[2]


    def read_list_of_ids(filedir):
        try:
            data = ''
            with open(filedir, 'r') as f:
                for line in f.readlines():
                    line = re.sub('[!@/#$|:]', ' ', line).strip()
                    data = data + ' ' + line
            return data
        except FileNotFoundError:
            sys.exit("!!!ERROR!!! \n"
                     "Please enter right name of file (with .txt).\n"
                     "Remember that your ID file should be in the same directory as this script.\n"
                     "If not please specify path to the ID file.\n"
                     "E.g. dir1/subdir1/list.txt")


    ids = read_list_of_ids(filename)
    print("Name of the file with IDs: {0}".format(filename))
    clear_ids = ids.split()
    clear_ids = list(dict.fromkeys(clear_ids))
    print("Number of unique IDs in file: {0}".format(len(clear_ids)))
    print("Name of the current batch of files: {0}".format(current_dir))

    print("-----------------------------------------MENU-----------------------------------------------")

    print("Unique Sequence Downloader Menu:\n"
          "1)Download sequences and create database\n"
          "2)Only create database from files in directory Downloads/{0}\n"
          "0)Exit".format(current_dir))
    print("Please choose an action:", end="")
    choice = input()

    if choice == "1":

        print("------------------------------------Directory check------------------------------------------")


        def mk_dirs(dirName):
            if not os.path.exists(dirName):
                os.mkdir(dirName)
                print("Directory {0} status: created".format(dirName))
            else:
                print("Directory {0} status: already exists".format(dirName))


        def mk_subdirs(dirName):
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


        dirs = ['Downloads', 'Failed', 'Database']
        for single_dir in dirs:
            mk_dirs(single_dir)
        mk_subdirs(current_dir)

        obsolete = []
        http_error = []

        print("--------------------------------------Downloading--------------------------------------------")


        def paths_and_urls(ids, cdir):
            pu = []
            for idx in ids:
                if re.match(r'([A-Z]){1}\d+', idx):
                    path_url = tuple([('./Downloads/{0}/{1}.fasta'.format(cdir, idx)),
                                      ('https://www.uniprot.org/uniprot/' + idx + '.fasta')])
                    pu.append(path_url)

            return pu


        def fetch_url(entry):
            path, url = entry
            uniprot_id = re.sub('[/.]', ' ', path).strip().split()
            if not os.path.exists(path):
                try:
                    fasta = urllib.request.urlopen(url)
                    fasta = fasta.read().decode()

                    if fasta == '':
                        obsolete.append(uniprot_id[2])
                    else:
                        with open(path, 'w+') as f:
                            f.write(fasta)
                except urllib.error.HTTPError:
                    http_error.append(uniprot_id[2])
            return path


        def download_ncbi(ncbi_id, cdir):
            try:
                handle = Entrez.efetch(db="protein", id=ncbi_id, rettype="fasta", retmode="text")
                record = handle.read()

                if record == '':
                    obsolete.append(ncbi_id)
                else:
                    with open('./Downloads/{0}/{1}.fasta'.format(cdir, idx), 'w+') as f:
                        f.write(record)

            except urllib.error.HTTPError:
                http_error.append(ncbi_id)


        urls = paths_and_urls(clear_ids, current_dir)

        print("Downloading files. Please wait...")

        ThreadPool(20).imap_unordered(fetch_url, urls)
        start = timer()
        for entry in urls:
            fetch_url(entry)

        Entrez.email = ''
        for idx in clear_ids:
            if re.match(r'([A-Z]){2,3}\w+', idx):
                download_ncbi(idx, current_dir)

        obsolete = list(dict.fromkeys(obsolete))

        mypath = ('./Downloads/{0}'.format(current_dir))
        onlyfiles = [f for f in os.listdir(mypath) if isfile(join(mypath, f))]

        total_number = len(onlyfiles)
        print("Download of " + str(total_number) + f" took {timer() - start} seconds")

        print("-------------------------------------Saving outputs------------------------------------------")
        print('Possibly obsolete entries were saved in Failed/{0}/Obsolete.txt'.format(current_dir))
        print('Number of obsolete entries: {0}'.format(str(len(obsolete))))

        with open('./Failed/{0}/Obsolete.txt'.format(current_dir), 'w+') as obs:
            for item in obsolete:
                obs.write(item + '\n')

        print('Entries with possible download error were saved in Failed/{0}/Error.txt'.format(current_dir))
        print('Number of download error entries: {0}'.format(str(len(http_error))))

        with open('./Failed/{0}/Error.txt'.format(current_dir), 'w+') as fail:
            for item in http_error:
                fail.write(item + '\n')

        print("Creating database. Please wait...")
        for file in onlyfiles:
            with open('./Downloads/{0}/{1}'.format(current_dir, file)) as big_fasta:
                one_fasta = big_fasta.read()
                with open('./Database/{0}.fasta'.format(current_dir), 'a+') as data:
                    data.write(one_fasta + '\n')

        print("Your sequences were saved in Database/{0}.fasta".format(current_dir))

    elif choice == "2":
        try:
            mypath = ('./Downloads/{0}'.format(current_dir))
            onlyfiles = [f for f in os.listdir(mypath) if isfile(join(mypath, f))]
            print("------------------------------------Saving database------------------------------------------")

            print("Creating database. Please wait...")
            for file in onlyfiles:
                with open('./Downloads/{0}/{1}'.format(current_dir, file)) as big_fasta:
                    one_fasta = big_fasta.read()
                    with open('./Database/{0}.fasta'.format(current_dir), 'a+') as data:
                        data.write(one_fasta + '\n')

            print("Your sequences were saved in Database/{0}.fasta".format(current_dir))
        except FileNotFoundError:
            print("!!!ERROR!!!\nThere are no files in directory Downloads/{0} to create database ".format(current_dir))
    elif choice == "0":
        print("-----------------------------------------EXIT-----------------------------------------------")
        sys.exit("USD will stop working now.")
    else:
        print("!!!ERROR!!!\nYou can choose only listed options!")

except IndexError:
    sys.exit("!!!ERROR!!! \n"
             "Please pass arguments correctly.\n"
             "#1 argument: Name of the file with IDs.\n"
             "#2 argument: Name for the current batch of files.\n"
             "Remember that your ID file should be in the same directory as this script.\n"
             "If not, please specify path to the ID file.\n"
             "E.g. python this_script.py list.txt test")
