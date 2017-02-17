from __future__ import print_function
import os.path
import sys
import json
import re
import requests
from tqdm import tqdm
from bs4 import BeautifulSoup
from combat_tb_model.model.core import Gene, Feature, FeatureSet
from .dbconn import graph, create_uniprot_nodes
from .uniprot import query_uniprot


def fetch_ortholog(locus_tag, ortholog_type='cdc1551'):
    ortholog_name = None
    if ortholog_type == 'cdc1551':
        tuberculist_url = 'http://tuberculist.epfl.ch/quicksearch.php?gene+name='
        response = requests.get(tuberculist_url + locus_tag)
        bs_obj = BeautifulSoup(response.text, 'html.parser')
        links = (bs_obj).find_all(href=re.compile('cmr.jcvi.org'))
        ortholog_name = links[0].text if links else None
    else:
        print('unknown ortholog type:', type, file=sys.stderr)
    return ortholog_name


def find_orthologs(ortholog_type='cdc1551', chunksize=400):
    # this yields lists of (locus_tag, gene_object) tuples
    # for known genes from the existing annotion, chunksize at a time
    cachefile = 'orthologs_cache.json'
    tmp_cachefile = cachefile + '.tmp'
    ortholog_cache = json.load(open(cachefile)) if os.path.exists(cachefile) else dict()
    ortholog_list = []
    # print(ortholog_cache)
    count = 1
    # TODO: find a better solution than this hack to identify
    # H37Rv genes
    for gene in Gene.select(graph).where('_.uniquename =~ "gene:Rv[0-9].*"'):
        # if count > 10:
        #     break
        # count += 1
        # print('.', end='', file=sys.stderr)
        locus_tag = list(gene.is_a)[0].uniquename.replace('gene:', '')
        if locus_tag in ortholog_cache:
            ortholog_name = ortholog_cache[locus_tag]
        else:
            ortholog_name = fetch_ortholog(locus_tag, ortholog_type)
            ortholog_cache[locus_tag] = ortholog_name

        json.dump(ortholog_cache, open(tmp_cachefile, 'w'))
        if os.path.exists(cachefile):
            os.unlink(cachefile)
        os.rename(tmp_cachefile, cachefile)
        if ortholog_name is not None:
            ortholog_list.append((ortholog_name, gene))
            count += 1
            if (count % chunksize) == 0:
                yield ortholog_list
                ortholog_list = []
    yield ortholog_list


def add_orthologs_to_db(ortholog_type='cdc1551'):
    ortholog_group = FeatureSet()
    ortholog_group.name = ortholog_type
    for chunk in find_orthologs(ortholog_type):
        ortholog_names = []
        for (name, gene) in tqdm(chunk):
            ortholog = Gene()
            ortholog.uniquename = ortholog.name = ("gene:"+name)
            ortholog_feature = Feature()
            ortholog_feature.uniquename = ortholog_feature.name = ("gene:"+name)
            ortholog.is_a.add(ortholog_feature)
            graph.create(ortholog_feature)
            graph.create(ortholog)
            ortholog_group.contains.add(ortholog_feature)
            try:
                feature = next(iter(gene.is_a))
            except StopIteration:
                sys.exit("Database error: gene {} has no Feature".format(gene.uniquename))
            feature.orthologous_to.add(ortholog_feature)
            graph.push(feature)
            ortholog_names.append(name)
            # this will create Protein / Polypeptide nodes for all
            # the proteins associated with these orthologs
            # commented out for new because we don't have enough structure (Gene->CDS->Transcript) to associate Protein with
            # create_uniprot_nodes(query_uniprot(ortholog_names, taxonomy='83331', proteome='UP000001020'), add_protein_interactions=False)
    graph.create(ortholog_group)
        # exit(0) # for debug purposes, end here
