"""
Interface for GFF processing.
"""
from __future__ import print_function
import pprint

from BCBio import GFF
from BCBio.GFF import GFFExaminer
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


def parse_gff(db, gff_file, featureset_name=None, featureset_description=None):
    """
    Parse GFF file
    :return:
    """

    db.create_organism_nodes()
    # we are not interested in exons as this is a bacterial genome
    limits = [["gene", "pseudogene", "tRNA_gene", "ncRNA_gene", "rRNA_gene"],
              ["transcript"], ["CDS"]]
    for limit in limits:
        print("Loading", limit, "...")
        load_gff_data(db, gff_file, limit)
    print("Done.")


def get_locus_tags(gff_file, chunk):
    """
    Return a list of locus tags from gff_file
    :param gff_file:
    :param chunk
    :return:
    """
    print("Getting locus_tags...")
    count = 0
    locus_tags = []
    for rec in GFF.parse(gff_file, limit_info=dict(gff_type=['gene'])):
        for gene in rec.features:
            locus_tag = gene.qualifiers["gene_id"][0]
            count += 1
            locus_tags.append(locus_tag)
            if count == chunk:
                yield locus_tags
                locus_tags = []
                count = 0
    yield locus_tags


def load_gff_data(db, gff_file, limit):
    """
    Extract and load features to Neo4j
    :param gff_file:
    :param limit:
    :return:
    """
    in_file = open(gff_file)
    limit_info = dict(gff_type=limit)
    transaction = db.graph.begin()
    for rec in GFF.parse(gff_file, limit_info=limit_info):
        for feature in tqdm(rec.features):
            rna = ["tRNA_gene", "ncRNA_gene", "rRNA_gene"]
            db.create_feature_nodes(feature, transaction)
            db.create_featureloc_nodes(feature, transaction)
            if feature.type == 'gene':
                db.create_gene_nodes(feature, transaction)
            elif feature.type == 'pseudogene':
                db.create_pseudogene_nodes(feature, transaction)
            elif feature.type == 'exon':
                db.create_exon_nodes(feature, transaction)
            elif feature.type in rna:
                db.create_rna_nodes(feature, transaction)
            elif feature.type == 'transcript':
                db.create_transcript_nodes(feature, transaction)
            elif feature.type == 'CDS':
                db.create_cds_nodes(feature, transaction)

    transaction.commit()
    db.create_feature_indexes()
    in_file.close()
