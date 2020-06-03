import os
import shutil
import sys

from Bio.SeqAnalysis.SeqDatabase import SeqDownloader, SeqDatabase

# Declaration of variables needed for running tests.
filename = "./test.txt"  # Additional file made only for reading function test.
dirs = ['Downloads', 'Failed', 'Database']
test_dir = "test"
test_tuple = ('./Downloads/test/A0A1A2.fasta', 'https://www.uniprot.org/uniprot/A0A1A2.fasta')
test_obs = []
test_http = []


def blockPrint():
    """
    Prevents some functions from printing their outputs.
    """
    sys.stdout = open(os.devnull, 'w')


def test_read_list_of_ids():
    """
    Check if the function reads the file correctly.
    """
    entries = SeqDownloader.read_list_of_ids(filename)
    assert entries == ['A0A1A2', 'B0B1B2', 'C0C1C2', 'D0D1D2', 'E0E1E2', 'F0F1F2', 'G0G1G2', 'H0H1H2',
                       'I0I1I2', 'J0J1J2'], \
        'Function execution failure. Check if "test.txt" is in the same directory.'


def test_mk_dirs():
    """
    Check if the function create directories correctly
    """
    for single_dir in dirs:
        SeqDownloader.mk_dirs(single_dir)
    for single_dir in dirs:
        assert os.path.exists(
            single_dir) is True, \
            "Function execution failure. Check if all test directories were deleted before running test."


def test_mk_subdirs():
    """
    Check if the function create directories correctly
    """
    SeqDownloader.mk_subdirs(test_dir)
    assert os.path.exists(
        "./Downloads/test") is True, \
        "Function execution failure. Check if all test directories were deleted before running test."
    assert os.path.exists(
        "./Failed/test") is True, \
        "Function execution failure. Check if all test directories were deleted before running test."


def test_paths_and_urls():
    """
    Check if the function create proper tuple of path and url.
    """
    urls = SeqDownloader.paths_and_urls(["A0A1A2"], test_dir)
    assert urls == [('./Downloads/test/A0A1A2.fasta',
                     'https://www.uniprot.org/uniprot/A0A1A2.fasta')], \
        "Function execution failure."


def test_fetch_url():
    """
    Check if the function returns file path, as well as ,downloads and saves file correctly.
    """
    fetched_url = SeqDownloader.fetch_url(test_tuple, test_obs, test_http)
    assert fetched_url == './Downloads/test/A0A1A2.fasta', \
        "Function execution failure."
    assert os.path.exists("./Downloads/test/A0A1A2.fasta") is True, \
        "Function execution failure."


def test_download_ncbi():
    """
    Check if the function downloads and saves NCBI file correctly using Entrez.
    """
    SeqDownloader.download_ncbi("YP_025292.1", test_dir, test_obs, test_http)
    assert os.path.exists("./Downloads/test/YP_025292.1.fasta") is True, \
        "Function execution failure."


def test_SeqDownloader_class():
    """
    Check if the SeqDownloader class downloads and saves files correctly.
    """
    SeqDownloader(["A0A1A2", "YP_025292.1"], test_dir)
    assert os.path.exists(
        "./Downloads/test/A0A1A2.fasta") is True, \
        "Class execution failure. Something wrong with one of the files."
    assert os.path.exists(
        "./Downloads/test/YP_025292.1.fasta") is True, \
        "Class execution failure. Something wrong with one of the files."


def test_get_and_SeqDatabase_class():
    """
    Check if the class SeqDatabase creates dict. Also checks if get() function works properly.
    """
    test_database = SeqDatabase(["A0A1A2", "YP_025292.1"], test_dir)
    assert type(test_database.database) == dict, \
        "Class execution failure."

    record = test_database.get()[0:1]
    same_record = test_database.get(["A0A1A2", "YP_025292.1"])
    assert record != same_record, \
        "Function execution failure. Records are not the same."


def rm_dirs():
    """
    Quick cleanup after certain tests and at the end of testing.
    """
    shutil.rmtree("./Database")
    shutil.rmtree("./Downloads")
    shutil.rmtree("./Failed")


if __name__ == '__main__':
    blockPrint()
    test_read_list_of_ids()
    test_mk_dirs()
    test_mk_subdirs()
    test_paths_and_urls()
    test_fetch_url()
    test_download_ncbi()
    rm_dirs()
    test_SeqDownloader_class()
    rm_dirs()
    test_get_and_SeqDatabase_class()
    rm_dirs()
