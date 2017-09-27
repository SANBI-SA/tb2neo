"""
Interface for GFF processing.
"""
from __future__ import print_function

import pprint
import sys

import datetime
import pytz

from BCBio import GFF
from BCBio.GFF import GFFExaminer
from tqdm import tqdm

from combat_tb_model.model.core import *

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


def strip_id_colon(id):
    if ':' in id:
        parts = id.split(':')
        assert len(parts) == 2, \
            "Expected ID with 2 parts, like thing:ID, got " + id
        return parts[1]
    else:
        return id


def get_parent(feature, thing, lookup_dict):
    """Get parent object of a thing

    :param feature: Bio.SeqFeature.SeqFeature
    :param thing: combat_tb_neomodel.model.core.Feature
    :param lookup_dict: dict
    :return: combat_tb_neomodel.model.core.Feature
    """
    assert 'Parent' in feature.qualifiers,\
        "No Parent attribute found in feature for {} {}".format(
            thing.__class__.__name__,
            thing.uniquename
        )
    parent_id = feature.qualifiers['Parent'][0]
    assert parent_id in lookup_dict,\
        "Parent with ID {} not found for {} {} with biotype {}".format(
            parent_id,
            thing.__class__.__name__,
            thing.uniquename,
            thing.biotype
        )
    parent_thing = lookup_dict[parent_id]
    return parent_thing

# def has_qualifier(feature, key, value):
#     if key in feature.qualifiers and feature.qualifiers[key][0] == value:
#         return True
#     else:
#         return False

def set_if_has(thing, feature, key, transform=None):
    """ set_if_has - key is in feature's qualifiers, set that attribute on thing, optionally transforming with transform

    :param thing: neomodel.StructuredNode
    :param feature: Bio.SeqFeature
    :param key: str
    :param transform: (str) -> str
    :return:
    """
    if key in feature.qualifiers:
        value = feature.qualifiers[key][0]
        if transform is not None:
            value = transform(value)
        setattr(thing, key.lower(), value)

def import_gff(graph, gff_file):
    """import_gff(graph, gff_file) - import features from GFF3 format gff_file to graph db

    :param graph: GraphDb
    :param gff_file: str
    :return:
    """
    in_handle = open(gff_file)
    gff_limits = dict(gff_type=('gene', 'transcript', 'pseudogene',
                                'CDS', 'tRNA_gene', 'ncRNA_gene',
                                'rRNA_gene', 'rRNA', 'repeat_region'))
    with graph.graph.transaction:
        tb = Organism(genus='Mycobacterium',
                      species='tuberculosis',
                      strain='H37Rv',
                      common_name='M. tuberculosis').save()

        record = next(GFF.parse(in_handle, target_lines=100))
        assert 'sequence-region' in record.annotations,\
            'Could not find sequence-region in GFF3 annotations'
        (name, start, end) = record.annotations['sequence-region'][0]
        chrom_length = end
        chrom = Chromosome(seqlen=chrom_length, uniquename='1')
        chrom.save()
        chrom.belongs_to.connect(tb)
        in_handle.seek(0)  # rewind file

        count = 0
        old_count = 1000
        quit = False
        gene_dict = dict()
        pseudogene_dict = dict()
        transcript_dict = dict()
        trna_dict = dict()
        ncrna_dict = dict()
        rrna_dict = dict()
        cds_dict = dict()
        repeat_region_count = 0
        while not quit:
            if count == old_count:
                quit = True
                break
            old_count = count
            location_dict = dict()
            for record in GFF.parse(in_handle, limit_info=gff_limits,
                                    target_lines=1):
                count += 1
                feature = record.features[0]
                f_type = feature.type
                if 'ID' in feature.qualifiers:
                    id = feature.qualifiers['ID'][0]
                elif f_type == 'repeat_region':
                    repeat_region_count += 1
                    id = 'repeat:' + str(repeat_region_count)
                else:
                    exit("Feature found with no ID and not a repeat_region: " + repr(feature))
                location_key = '_'.join([str(feature.location.start),
                                         str(feature.location.end),
                                         str(feature.location.strand)])
                if location_key in location_dict:
                    location = location_dict[location_key]
                else:
                    location = Location(location_key=location_key,
                                        start=feature.location.start,
                                        end=feature.location.end,
                                        strand=str(feature.location.strand))
                    location.save()
                    location.feature.connect(chrom)
                    location_dict[location_key] = location
                feature_length = end - start
                assert feature_length >= 1,\
                    "Feature found with length < 1 after processing {} records".format(count)
                if f_type == 'gene':
                    thing = Gene(uniquename=id)
                    set_if_has(thing, feature, 'biotype')
                    set_if_has(thing, feature, 'description')
                    gene_dict[id] = thing
                elif f_type == 'pseudogene':
                    if 'transcript_id' in feature.qualifiers:
                        thing = Transcript(uniquename=id)
                        set_if_has(thing, feature, 'biotype')
                    else:
                        thing = PseudoGene(uniquename=id)
                        set_if_has(thing, feature, 'description')
                        pseudogene_dict[id] = thing
                elif (f_type == 'transcript'):
                    thing = Transcript(uniquename=id)
                    set_if_has(thing, feature, 'biotype')
                    transcript_dict[id] = thing
                elif f_type == 'CDS':
                    thing = CDS(uniquename=id)
                    cds_dict[id] = thing
                elif f_type == 'tRNA_gene':
                    thing = Trna(uniquename=id)
                    trna_dict[id] = thing
                elif f_type == 'ncRNA_gene':
                    thing = NCrna(uniquename=id)
                    ncrna_dict[id] = thing
                elif f_type == 'rRNA_gene':
                    thing = Rrna(uniquename=id)
                    rrna_dict[id] = thing
                elif f_type == 'rRNA':
                    thing = Transcript(uniquename=id)
                    set_if_has(thing, feature, 'biotype')
                elif f_type == 'repeat_region':
                    thing = RepeatRegion(uniquename=id)
                    set_if_has(thing, feature, 'description')
                else:
                    print("Unknown feature type:", f_type, file=sys.stderr)
                    continue

                set_if_has(thing, feature, 'Name', strip_id_colon)

                thing.seqlen = feature_length
                thing.timelastmodified = datetime.datetime.now(tz=pytz.utc)
                thing.save()
                thing.location.connect(location)
                if isinstance(thing, Transcript):
                    assert hasattr(thing, 'biotype'),\
                        "Transcript with ID {} has not biotype".format(thing.uniquename)
                    if thing.biotype == 'protein_coding':
                        lookup_dict = gene_dict
                    elif thing.biotype == 'tRNA':
                        lookup_dict = trna_dict
                    elif thing.biotype == 'ncRNA':
                        lookup_dict = ncrna_dict
                    elif thing.biotype == 'pseudogene':
                        lookup_dict = pseudogene_dict
                    elif thing.biotype == 'rRNA':
                        lookup_dict = rrna_dict
                    else:
                        exit("unknown biotype:" + thing.biotype)
                    parent_thing = get_parent(feature, thing, lookup_dict)
                    try:
                        thing.part_of.connect(parent_thing)
                    except ValueError:
                        exit("Failed to connect to a: " + parent_thing.__class__.__name__ + " " + thing.__class__.__name__)
                elif f_type == 'CDS':
                    parent_transcript = get_parent(feature, thing, transcript_dict)
                    thing.part_of.connect(parent_transcript)

                print('.', end='', file=sys.stderr)
                # print(thing.uniquename)
                # location = thing.location.connect(chrom)
                # location.start = start
                # location.end = end
                # location.strand = strand
                # location.save()

                thing.belongs_to.connect(tb)

def parse_gff(db, gff_file, featureset_name=None, featureset_description=None):
    """
    Parse GFF file
    :return:
    """
    sys.stdout.write("Parsing GFF...")
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
    sys.stdout.write("Extract and load features to Neo4j.")
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
