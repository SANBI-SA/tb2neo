"""
Interface to the Neo4j Database
"""
from urllib2 import HTTPError

from Bio import Entrez
from Bio import Medline
from combat_tb_model.model import *
from py2neo import watch, Graph, getenv
from quickgo import fetch_quick_go_data

graph = Graph(host=getenv("DB", "localhost"), bolt=True, password=getenv("NEO4J_PASSWORD", ""))
watch("neo4j.bolt")


def create_organism_nodes():
    """
    Create Organism Nodes
    :return:
    """
    abbrev = "H37Rv"
    genus = "Mycobacterium"
    species = "M. tuberculosis"
    common_name = "TB"

    organism = Organism(abbreviation=abbrev, genus=genus, species=species, common_name=common_name)
    graph.create(organism)


def create_gene_nodes(feature):
    """
    Create Gene Nodes
    :param feature:
    :return:
    """
    names = get_feature_name(feature)
    name = names.get("Name", names.get("UniqueName"))
    unique_name = names.get("UniqueName", name)

    gene = Gene()
    gene.ontology_id = gene.so_id
    gene.name = name
    gene.uniquename = unique_name
    graph.create(gene)


def create_transcript_nodes(feature):
    """
    Create Transcipt Nodes
    :param feature:
    :return:
    """
    names = get_feature_name(feature)
    name = names.get("Name", names.get("UniqueName"))
    unique_name = names.get("UniqueName", name)

    transcript = Transcript()
    transcript.ontology_id = transcript.so_id
    transcript.name = name
    transcript.uniquename = unique_name
    graph.create(transcript)


def create_pseudogene_nodes(feature):
    """
    Create Pseudogene Nodes
    :param feature:
    :return:
    """
    names = get_feature_name(feature)
    name = names.get("Name", names.get("UniqueName"))
    unique_name = names.get("UniqueName", name)

    pseudogene = PseudoGene()
    pseudogene.ontology_id = pseudogene.so_id
    pseudogene.name = name
    pseudogene.uniquename = unique_name
    graph.create(pseudogene)


def create_exon_nodes(feature):
    """
    Create Exon Nodes
    :param feature:
    :return:
    """
    names = get_feature_name(feature)
    name = names.get("Name", names.get("UniqueName"))
    unique_name = names.get("UniqueName", name)

    exon = Exon()
    exon.ontology_id = exon.so_id
    exon.name = name
    exon.uniquename = unique_name
    graph.create(exon)


def create_rna_nodes(feature):
    """
    Create RNA Nodes
    :param feature:
    :return:
    """
    names = get_feature_name(feature)
    name = names.get("Name", names.get("UniqueName"))
    unique_name = names.get("UniqueName", name)

    if feature.type == 'tRNA_gene':
        trna = TRna()
        trna.ontology_id = trna.so_id
        trna.name = name
        trna.uniquename = unique_name
        graph.create(trna)
    if feature.type == 'ncRNA_gene':
        ncrna = NCRna()
        ncrna.ontology_id = ncrna.so_id
        ncrna.name = name
        ncrna.uniquename = unique_name
        graph.create(ncrna)
    if feature.type == 'rRNA_gene':
        rrna = RRna()
        rrna.ontology_id = rrna.so_id
        rrna.name = name
        rrna.uniquename = unique_name
        graph.create(rrna)


def create_cds_nodes(feature):
    """
    Create CDS Nodes
    :param feature:
    :return:
    """
    names = get_feature_name(feature)
    name = names.get("Name", names.get("UniqueName"))
    unique_name = names.get("UniqueName", name)

    cds = CDS()
    cds.name = name
    cds.uniquename = unique_name
    cds.ontology_id = cds.so_id
    graph.create(cds)


def create_feature_nodes(feature):
    """
    Create Feature Nodes
    :param feature:
    :return:
    """
    names = get_feature_name(feature)
    name = names.get("Name", names.get("UniqueName"))
    unique_name = names.get("UniqueName", name)

    if feature.qualifiers.get('Parent'):
        parent = feature.qualifiers['Parent'][0]
    # [feature.qualifiers['Parent'][0].find(":") + 1:]
    else:
        parent = None

    _feature = Feature()
    _feature.name = name
    _feature.parent = parent
    _feature.uniquename = unique_name
    graph.create(_feature)


def create_featureloc_nodes(feature):
    """
    Create FeatureLoc Nodes
    :param feature:
    :return:
    """
    srcfeature_id = get_feature_name(feature).get("UniqueName")
    feature_loc = FeatureLoc(srcfeature_id=srcfeature_id, fmin=feature.location.start, fmax=feature.location.end,
                             strand=feature.location.strand)
    graph.create(feature_loc)


def get_feature_name(feature):
    """
    Get Feature Name and UniqueName
    :param feature:
    :return:
    """
    names = dict()
    if feature.qualifiers.get("Name"):
        names["Name"] = feature.qualifiers["Name"][0]
        names["UniqueName"] = feature.id
        # [feature.id.find(":") + 1:]
    else:
        names["Name"] = names["UniqueName"] = feature.id
        # [feature.id.find(":") + 1:]
    return names


def build_relationships():
    """
    Build relationships
    :return:
    """
    # TODO: Try optimize this
    print("Building Relationships...")
    features = Feature.select(graph)
    for feature in features:
        # Find organism via __primarykey__ and link with feature via BELONGS_TO
        org = Organism.select(graph, 'Mycobacterium').first()
        feature.belongs_to.add(org)

        # Find feature with a parent attr. matching this features uniquename and link them via RELATED_TO
        _feature = Feature.select(graph).where("_.parent = '{}'".format(feature.uniquename)).first()
        if _feature:
            feature.related_to.add(_feature)

        # Building is_a relationships
        gene = Gene.select(graph, feature.uniquename).first()
        if gene:
            gene.is_a.add(feature)
            graph.push(gene)
            # Find feature with this gene's uniquename as a parent
            _feature = Feature.select(graph).where("_.parent = '{}'".format(gene.uniquename)).first()
            if _feature:
                # Find transcript: A gene is a parent to it.
                transcript = Transcript.select(graph, _feature.uniquename).first()
                if transcript:
                    transcript.part_of.add(gene)
                    graph.push(transcript)

        p_gene = PseudoGene.select(graph, feature.uniquename).first()
        if p_gene:
            p_gene.is_a.add(feature)
            graph.push(p_gene)
        transcript = Transcript.select(graph, feature.uniquename).first()
        if transcript:
            transcript.is_a.add(feature)
            graph.push(transcript)
            # Find feature with this transcript's uniquename as a parent
            _feature = Feature.select(graph).where("_.parent = '{}'".format(transcript.uniquename)).first()
            if _feature:
                # Find exon: A transcript is a parent to it
                exon = Exon.select(graph, _feature.uniquename).first()
                if exon:
                    exon.part_of.add(transcript)
                    graph.push(exon)
                # Find cds: A transcript is a parent to it
                cds = CDS.select(graph, _feature.uniquename).first()
                if cds:
                    cds.part_of.add(transcript)
                    graph.push(cds)

        exon = Exon.select(graph, feature.uniquename).first()
        if exon:
            exon.is_a.add(feature)
            graph.push(exon)
        cds = CDS.select(graph, feature.uniquename).first()
        if cds:
            cds.is_a.add(feature)
            graph.push(cds)

        # Find feature location with a srcfeature_id attr. matching this features uniquename and link them via
        # LOCATED_AT
        _feature = FeatureLoc.select(graph, feature.uniquename).first()
        if _feature:
            feature.location.add(_feature)
        graph.push(feature)


def create_cv_term_nodes(polypeptide, bp, cc, mf):
    """
    Create CvTerm Nodes.
    :param polypeptide:
    :param bp:
    :param cc:
    :param mf:
    :return:
    """
    # go(biological process)
    go_bp_ids = [t[t.find('G'):-1] for t in bp.split('; ') if t is not '']
    go_bp_defs = [t[:t.find('[') - 1] for t in bp.split('; ') if t is not '']
    # go(cellular component)
    go_cc_ids = [t[t.find('G'):-1] for t in cc.split('; ') if t is not '']
    go_cc_defs = [t[:t.find('[') - 1] for t in cc.split('; ') if t is not '']
    # go(molecular function)
    go_mf_ids = [t[t.find('G'):-1] for t in mf.split('; ') if t is not '']
    go_mf_defs = [t[:t.find('[') - 1] for t in mf.split('; ') if t is not '']

    # TODO: Find a way to refactor this.
    for _id in go_bp_ids:
        cv = CvTerm()
        for _def in go_bp_defs:
            cv.name = _id
            cv.definition = _def
            cv.namespace = "biological process"
            graph.create(cv)
            polypeptide.cvterm.add(cv)
            graph.push(polypeptide)

    for _id in go_mf_ids:
        cv = CvTerm()
        for _def in go_mf_defs:
            cv.name = _id
            cv.definition = _def
            cv.namespace = "cellular component"
            graph.create(cv)
            polypeptide.cvterm.add(cv)
            graph.push(polypeptide)
    for _id in go_cc_ids:
        cv = CvTerm()
        for _def in go_cc_defs:
            cv.name = _id
            cv.definition = _def
            cv.namespace = "molecular function"
            graph.create(cv)
            polypeptide.cvterm.add(cv)
            graph.push(polypeptide)


def create_interpro_term_nodes(polypeptide, entry):
    """
    Create InterPro Term Nodes.
    :param polypeptide:
    :param entry:
    :return:
    """
    terms = [t for t in entry.split("; ") if t is not '']
    for interpro in terms:
        import time
        dbxref = DbXref(db="InterPro", accession=interpro, version=time.time())
        graph.create(dbxref)
        polypeptide.dbxref.add(dbxref)
        graph.push(polypeptide)


def fetch_publications(citation):
    """
    Fetch Publications.
    :param citation:
    :return:
    """
    Entrez.email = 'A.N.Other@example.com'
    try:
        h = Entrez.efetch(db='pubmed', id=citation, rettype='medline', retmode='text')
    except HTTPError:
        import time
        time.sleep(200)
        h = Entrez.efetch(db='pubmed', id=citation, rettype='medline', retmode='text')
    records = Medline.parse(h)
    return records


def create_author_nodes(publication, full_author):
    """
    Create Author Nodes.
    :param publication:
    :param full_author:
    :return:
    """
    # TODO: Get more info about Authors
    if full_author:
        for au in full_author:
            _author = Author()
            _author.givennames = au
            graph.create(_author)
            _author.wrote.add(publication)
            publication.author.add(_author)
            graph.push(_author)
            graph.push(publication)


def create_pub_nodes(polypeptide, pubs):
    """
    Create Publication Nodes
    :param polypeptide:
    :param pubs:
    :return:
    """
    citations = [c for c in pubs.split("; ") if c is not '']
    records = fetch_publications(citations)
    # https://www.nlm.nih.gov/bsd/mms/medlineelements.html
    for record in records:
        pm_id = record.get('PMID', None)
        title = record.get('TI', None)
        volume = record.get('VI', None)
        issue = record.get('IP', None)
        pages = record.get('PG', None)
        date_of_pub = record.get('DP', None)
        pub_place = record.get('PL', None)
        publisher = record.get('SO', None)
        author = record.get('AU', None)
        full_author = record.get('FAU', None)

        pub = Publication()
        pub.pmid = pm_id
        pub.title = title
        pub.volume = volume
        pub.issue = issue
        pub.pages = pages
        pub.year = date_of_pub
        pub.pubplace = pub_place
        pub.publisher = publisher
        graph.create(pub)

        polypeptide.published_in.add(pub)
        graph.push(polypeptide)
        create_author_nodes(pub, full_author)


def create_is_a_cv_term_rel():
    """
    Creating IS_A relationships between CVTerms
    :return:
    """
    cv = CvTerm.select(graph)
    for cv.name in cv:
        is_a_list = fetch_quick_go_data(cv.name)
        # cv = CvTerm.select(graph, _id).first()
        for go in is_a_list:
            goid = go[go.find('G'):go.find('!')].strip()
            cv_term = CvTerm.select(graph, goid).first()
            if cv_term:
                cv.is_a.add(cv_term)
                graph.push(cv)


def create_uniprot_nodes(uniprot_data):
    """
    Build DbXref nodes from UniProt results.
    :param uniprot_data:
    :return:
    """
    count = 0
    for entry in uniprot_data:
        count += 1

        dbxref = DbXref(db="UniProt", accession=entry[1], version=entry[0])
        graph.create(dbxref)

        polypeptide = Polypeptide()
        polypeptide.name = entry[9]
        polypeptide.uniquename = entry[0]
        polypeptide.ontology_id = polypeptide.so_id
        polypeptide.seqlen = entry[16]
        polypeptide.residues = entry[14]
        polypeptide.parent = entry[2]
        polypeptide.family = entry[17]
        polypeptide.function = entry[13]
        graph.create(polypeptide)

        gene = Gene.select(graph, "gene:" + entry[2]).first()
        if gene:
            _feature = Feature.select(graph).where("_.parent = '{}'".format(gene.uniquename)).first()
            if _feature:
                transcript = Transcript.select(graph, _feature.uniquename).first()
                if transcript:
                    cds = CDS.select(graph, "CDS" + transcript.uniquename[transcript.uniquename.find(":"):]).first()
                    if cds:
                        # Polypetide-derives_from->CDS
                        polypeptide.derives_from.add(cds)
                        cds.polypeptide.add(polypeptide)
                        graph.push(polypeptide)
                        graph.push(cds)

        polypeptide.dbxref.add(dbxref)
        graph.push(polypeptide)

        create_cv_term_nodes(polypeptide, entry[18], entry[19], entry[20])
        create_interpro_term_nodes(polypeptide, entry[5])
        create_pub_nodes(polypeptide, entry[11])
    create_is_a_cv_term_rel()
    print ("TOTAL:", count)
