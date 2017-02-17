#!/usr/bin/env pytest

import pytest
from py2neo import Graph, getenv
from combat_tb_model.model import core

# NOTE: TODO
# 45 tRNAs
# 30 ncRNA_gene
# 3 rRNAs

@pytest.fixture
def graph():
    handle = Graph(host=getenv('DB', '127.0.0.1'), bolt=True, http_port=7474, password=getenv("NEO4J_PASSWORD", ""))
    # watch("neo4j.bolt")
    return handle

def test_rv0001(graph):
    gene = core.Gene.select(graph).where("_.uniquename = 'gene:Rv0001'").first()
    assert gene.name == 'dnaA', "Expected name 'dnaA'"
    feature = next(iter(gene.is_a))
    assert isinstance(feature, core.Feature), "feature is {}, not core.Feature".format(type(feature))

def test_rv0001_feature_loc(graph):
    feature = core.Feature.select(graph).where("_.parent = 'gene:Rv0001'").first()
    location = next(iter(feature.location))
    assert isinstance(location, core.FeatureLoc), "location is {}, not core.FeatureLoc".format(type(feature))
    assert location.fmin == 0 and location.fmax == 1524 and location.strand == 1, (
        "Expected strand 1, 0..1524, got strand {} {}..{}".format(location.strand, location.fmin, location.fmax))

def test_CCP46753_feature_loc(graph):
    feature = core.Feature.select(graph).where("_.parent = 'transcript:CCP46753'").first()
    location = next(iter(feature.location))
    assert isinstance(location, core.FeatureLoc), "location is {}, not core.FeatureLoc".format(type(feature))
    assert location.fmin == 4410785 and location.fmax == 4410929 and location.strand == -1, (
        "Expected strand -1, 4410785..4410929, got strand {} {}..{}".format(location.strand, location.fmin, location.fmax))

H37RV_GENE_COUNT = 4018
H37RV_PSEUDOGENE_COUNT = 13
CDC1551_ORTHOLOG_COUNT = 3785
TOTAL_GENE_COUNT = H37RV_GENE_COUNT + CDC1551_ORTHOLOG_COUNT
# total Polypeptide count is 3979 but 2 are not clearly associated with genes
TOTAL_CONNECTED_POLYPEPTIDE_COUNT = 3977
def test_gene_count(graph):
    count = graph.evaluate("MATCH (g:Gene) WHERE g.uniquename =~ 'gene:Rv[0-9].*' RETURN count(g)")
    assert count == H37RV_GENE_COUNT, "Expected {} genes in H37Rv".format(H37RV_GENE_COUNT)
    count = graph.evaluate("MATCH (pg:PseudoGene) WHERE pg.uniquename =~ 'gene:Rv[0-9].*' RETURN count(pg)")
    assert count == H37RV_PSEUDOGENE_COUNT, "Expected {} pseudogenes in H37Rv".format(H37RV_PSEUDOGENE_COUNT)

def test_ortholog_count(graph):
    count = graph.evaluate("MATCH (g:Gene) -[:IS_A]-> (f:Feature) -[:ORTHOLOGOUS_TO]-> (f2:Feature) <-[:IS_A]- (g2:Gene) WHERE g.uniquename =~ 'gene:Rv[0-9].*' RETURN count(distinct(g2))")
    assert count == CDC1551_ORTHOLOG_COUNT, "Expected {} CDC1551 orthologs".format(CDC1551_ORTHOLOG_COUNT)
    count = graph.evaluate("MATCH (g:Gene) RETURN count(g)")
    assert count == TOTAL_GENE_COUNT, "Expected {} total genes".format(TOTAL_GENE_COUNT)

def test_polypeptide_count(graph):
    count = graph.evaluate("MATCH (p:Polypeptide) -[:DERIVES_FROM]->(CDS) -[:PART_OF]->(Transcript)-[:PART_OF]->(Gene) RETURN count(distinct(p.uniquename))")
    assert count == TOTAL_CONNECTED_POLYPEPTIDE_COUNT, "Expected {} total polypeptides (connected to genes)".format(TOTAL_CONNECTED_POLYPEPTIDE_COUNT)

TRNA_COUNT = 45
NCRNA_COUNT = 30
RRNA_COUNT = 3
def test_noncoding_rna_count(graph):
    count = graph.evaluate("MATCH (t:TRna) RETURN count(t)")
    assert count == TRNA_COUNT, "Expected {} tRNAs".format(TRNA_COUNT)
    count = graph.evaluate("MATCH (nc:NCRna) RETURN count(nc)")
    assert count == NCRNA_COUNT, "Expected {} ncRNAs".format(NCRNA_COUNT)
    count = graph.evaluate("MATCH (r:RRna) RETURN count(r)")
    assert count == RRNA_COUNT, "Expected {} rRNAs".format(RRNA_COUNT)

def test_cds_connections(graph):
    count = graph.evaluate("MATCH (c:CDS) -[:PART_OF]-> (Transcript) -[:PART_OF]-> (Gene) RETURN count(c)")
    assert count == H37RV_GENE_COUNT, "Expected {} CDSs connected to genes".format(H37RV_GENE_COUNT)

PUB_COUNT = 1608
PP_PUBLISHED_COUNT = 3979
def test_publication_count(graph):
    count = graph.evaluate("MATCH (p:Publication) RETURN count(p)")
    assert count == PUB_COUNT, "Expected {} publications".format(PUB_COUNT)
    count = graph.evaluate("MATCH (p:Publication)<-[:PUBLISHED_IN]-(pp:Polypeptide) RETURN count(distinct(pp))")
    assert count == PP_PUBLISHED_COUNT, "Expected {} polypeptides linked to publications".format(PP_PUBLISHED_COUNT)

INTERPRO_DBXREF_COUNT = 7612
def test_dbxref_count(graph):
    count = graph.evaluate("MATCH (d:DbXref) RETURN count(d)")
    assert count == INTERPRO_DBXREF_COUNT, "Expected {} DbXrefs (Interpro terms)".format(INTERPRO_DBXREF_COUNT)
    connected_count = graph.evaluate("MATCH (Polypeptide) -[:XREF]-> (d:DbXref) RETURN count(distinct(d))")
    assert connected_count == count, "Expected connected Xref count to match Xref count: {}".format(count)

GO_CVTERM_COUNT = 2231
def test_cvterm_count(graph):
    count = graph.evaluate("MATCH (c:CvTerm) RETURN count(c)")
    assert count == GO_CVTERM_COUNT, "Expected {} CvTerms (GO terms)".format(GO_CVTERM_COUNT)
    connected_count = graph.evaluate("MATCH (Polypeptide) -[:ASSOC_WITH]-> (c:CvTerm) RETURN count(distinct(c))")
    assert connected_count == count, "Expected connected CvTerm count to match CvTerm count: {}".format(count)
