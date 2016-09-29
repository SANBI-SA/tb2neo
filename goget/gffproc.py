import pprint
from time import time

from BCBio import GFF
from BCBio.GFF import GFFExaminer
from dbconn import *
from tqdm import tqdm


def examine(gff_file):
    """
    Examine GFF file
    :param gff_file:
    :return:
    """
    examiner = GFFExaminer()
    in_file = open(gff_file)
    pprint.pprint(examiner.available_limits(in_file))
    in_file.close()


def parse_gff(gff_file):
    """
    Parse GFF file
    :return:
    """

    create_organism_nodes()
    limits = [["gene", "pseudogene", "exon", "tRNA_gene", "ncRNA_gene", "rRNA_gene"], ["transcript"], ["CDS"]]
    for limit in limits:
        print("Loading", limit, "...")
        load_gff_data(gff_file, limit)
    print("Done.")


def get_locus_tags(gff_file, chunk):
    """
    Return a list of locus tags from gff_file
    :param gff_file:
    :param chunk
    :return:
    """
    print("Getting locus_tags...")
    start = time()
    count = 0
    locus_tags = []
    for rec in GFF.parse(gff_file, limit_info=dict(gff_type=['gene'])):
        for gene in rec.features:
            locus_tag = gene.qualifiers["gene_id"][0]
            count += 1
            locus_tags.append(locus_tag.lower())
            if count == chunk:
                yield locus_tags
                locus_tags = []
                count = 0
    end = time()
    print('Done getting locus_tags in', end - start, 'secs.')
    yield locus_tags


def load_gff_data(gff_file, limit):
    """
    Extract and load features to Neo4j
    :param gff_file:
    :param limit:
    :return:
    """
    in_file = open(gff_file)
    limit_info = dict(gff_type=limit)
    for rec in GFF.parse(gff_file, limit_info=limit_info):
        for feature in tqdm(rec.features):
            rna = ["tRNA_gene", "ncRNA_gene", "rRNA_gene"]
            create_feature_nodes(feature)
            if feature.type == 'gene':
                create_gene_nodes(feature)
            elif feature.type == 'pseudogene':
                create_pseudogene_nodes(feature)
            elif feature.type == 'exon':
                create_exon_nodes(feature)
            elif feature.type in rna:
                create_rna_nodes(feature)
            elif feature.type == 'transcript':
                create_transcript_nodes(feature)
            elif feature.type == 'CDS':
                create_cds_nodes(feature)
            create_featureloc_nodes(feature)
    in_file.close()
