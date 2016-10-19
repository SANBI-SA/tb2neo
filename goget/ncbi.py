"""
Interface to NCBI.
"""
import time
from urllib2 import HTTPError

from Bio import Entrez
from Bio import Medline


def fetch_publications(citation):
    """
    Fetch Publications.
    :param citation:
    :return:
    """
    print("=========================================")
    print("About to fetch Publication data from PubMed.")
    print("=========================================")
    time.sleep(2)
    Entrez.email = 'A.N.Other@example.com'
    try:
        h = Entrez.efetch(db='pubmed', id=citation, rettype='medline', retmode='text')
    except HTTPError:
        time.sleep(200)
        h = Entrez.efetch(db='pubmed', id=citation, rettype='medline', retmode='text')
    records = Medline.parse(h)
    return records
