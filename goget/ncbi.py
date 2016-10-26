"""
Interface to NCBI.
"""
import time
try:
  from urllib2 import HTTPError
except ImportError:
  from urllib.error import HTTPError

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

def fetch_publication_list(citations):
    """
    Fetch Publications.
    :param citations:
    :return:
    """
    print("=====================================================================")
    print("About to fetch Publication data for {} publications from PubMed.".format(len(citations)))
    print("=====================================================================")
    citation_string = ','.join(citations)
    Entrez.email = 'support@sanbi.ac.za'
    retries = 5
    failed = True
    for i in range(retries):
        try:
            h = Entrez.efetch(db='pubmed', id=citation_string, rettype='medline', retmode='text')
            failed = False
        except HTTPError:
            pass
        else:
            break
        finally:
            time.sleep(0.4) # we are not allowed to hit NCBI more than 3 times per second
    if failed:
        print("Retrieval from PubMed failed")
        records = []
    else:
        records = Medline.parse(h)
    return records
