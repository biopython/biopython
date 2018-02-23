"""InterProScan."""
from .interpro_xml import InterproXmlParser


# if not used as a module, run the doctest
if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()
