"""
Interface to the Neo4j Database
"""

from __future__ import print_function
import re
import sys
from tqdm import tqdm
from combat_tb_model.model.core import *
from goget.ncbi import fetch_publication_list
from py2neo import Graph, getenv, watch
from .quickgo import fetch_quick_go_data
from .uniprot import *

graph = Graph(host=getenv("DB", "localhost"), bolt=True, password=getenv("NEO4J_PASSWORD", ""))

# watch("neo4j.bolt")


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


def create_gene_nodes(feature, transaction=None):
    """
    Create Gene Nodes
    :param feature:
    :return:
    """
    names = get_feature_name(feature)
    name = names.get("Name", names.get("UniqueName"))
    unique_name = names.get("UniqueName", name)
    description = feature.qualifiers["description"]
    biotype = feature.qualifiers['biotype'][0]

    gene = Gene()
    gene.ontology_id = gene.so_id
    gene.name = name
    gene.uniquename = unique_name
    gene.biotype = biotype
    gene.description = description
    if transaction is not None:
        transaction.create(gene)
    else:
        graph.create(gene)


def create_transcript_nodes(feature, transaction=None):
    """
    Create Transcipt Nodes
    :param feature:
    :return:
    """
    names = get_feature_name(feature)
    name = names.get("Name", names.get("UniqueName"))
    unique_name = names.get("UniqueName", name)
    biotype = feature.qualifiers['biotype'][0]

    transcript = Transcript()
    transcript.ontology_id = transcript.so_id
    transcript.name = name
    transcript.uniquename = unique_name
    transcript.biotype = biotype
    if transaction is not None:
        transaction.create(transcript)
    else:
        graph.create(transcript)


def create_pseudogene_nodes(feature, transaction=None):
    """
    Create Pseudogene Nodes
    :param feature:
    :return:
    """
    names = get_feature_name(feature)
    name = names.get("Name", names.get("UniqueName"))
    unique_name = names.get("UniqueName", name)
    description = feature.qualifiers["description"][0]
    biotype = feature.qualifiers['biotype'][0]

    pseudogene = PseudoGene()
    pseudogene.ontology_id = pseudogene.so_id
    pseudogene.name = name
    pseudogene.uniquename = unique_name
    pseudogene.description = description
    pseudogene.biotype = biotype
    if transaction is not None:
        transaction.create(pseudogene)
    else:
        graph.create(pseudogene)


def create_exon_nodes(feature, transaction=None):
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
    if transaction is not None:
        transaction.create(exon)
    else:
        graph.create(exon)


def create_rna_nodes(feature, transaction=None):
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
        if transaction is not None:
            transaction.create(trna)
        else:
            graph.create(trna)
    if feature.type == 'ncRNA_gene':
        ncrna = NCRna()
        ncrna.ontology_id = ncrna.so_id
        ncrna.name = name
        ncrna.uniquename = unique_name
        if transaction is not None:
            transaction.create(ncrna)
        else:
            graph.create(ncrna)
    if feature.type == 'rRNA_gene':
        rrna = RRna()
        rrna.ontology_id = rrna.so_id
        rrna.name = name
        rrna.uniquename = unique_name
        if transaction is not None:
            transaction.create(rrna)
        else:
            graph.create(rrna)


def create_cds_nodes(feature, transaction=None):
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
    if transaction is not None:
        transaction.create(cds)
    else:
        graph.create(cds)

def create_feature_nodes(feature, transaction=None):
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
    if transaction is not None:
        transaction.create(_feature)
    else:
        graph.create(_feature)

def create_featureloc_nodes(feature, graph):
    """
    Create FeatureLoc Nodes
    :param feature:
    :return:
    """
    srcfeature_id = get_feature_name(feature).get("UniqueName")
    feature_loc = FeatureLoc(srcfeature_id=srcfeature_id, fmin=feature.location.start, fmax=feature.location.end,
                             strand=feature.location.strand)
    if graph is not None:
        graph.create(feature_loc)
    else:
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
        names["UniqueName"] = feature.id if feature.id != '' else None
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
    feature_count = graph.run("MATCH (f:Feature) RETURN count(f) AS count").next()['count']
    features = Feature.select(graph).where("_.uniquename = 'transcript:CCP46448'")
    featureset = FeatureSet().select(graph).where("_.name = 'h37rv'").first()
    if featureset is None:
        featureset = FeatureSet()
        featureset.name = 'h37rv'
        featureset.description = 'M. tuberculosis H37Rv'
        print("Creating Featureset for {}".format(featureset.name))
        graph.create(featureset)
    org = Organism.select(graph, 'Mycobacterium').first()
    for feature in tqdm(features, total=feature_count):
        # Find organism via __primarykey__ and link with feature via BELONGS_TO
        if org is not None:
            feature.belongs_to.add(org)
        featureset.contains.add(feature)

        # Find feature with a parent attr. matching this features uniquename and link them via RELATED_TO
        print('ping', file=sys.stderr)
        _feature = Feature.select(graph).where("_.parent = '{}'".format(feature.uniquename)).first()
        print(feature, _feature, file=sys.stderr)

        if _feature:
            feature.related_to.add(_feature)

        # Building is_a relationships
        gene = Gene.select(graph, feature.uniquename).first()
        if gene:
            gene.is_a.add(feature)
            featureset.contains.add(gene)
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
            featureset.contains.add(p_gene)
            p_gene.is_a.add(feature)
            graph.push(p_gene)
        transcript = Transcript.select(graph, feature.uniquename).first()
        if transcript:
            featureset.contains.add(transcript)
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
            featureset.contains.add(exon)
            exon.is_a.add(feature)
            graph.push(exon)
        cds = CDS.select(graph, feature.uniquename).first()
        if cds:
            featureset.contains.add(cds)
            cds.is_a.add(feature)
            graph.push(cds)

        # Find feature location with a srcfeature_id attr. matching this features uniquename and link them via
        # LOCATED_AT
        _feature = FeatureLoc.select(graph, feature.uniquename).first()
        if _feature:
            feature.location.add(_feature)
        graph.push(feature)
    if featureset is not None:
        graph.push(featureset)

def create_cv_term_nodes(polypeptide, bp, cc, mf):
    """
    Create CvTerm Nodes and build Polypetide relationships.
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

    if "name" not in graph.schema.get_indexes("CvTerm"):
        graph.schema.create_index("CvTerm", "name")


def create_interpro_term_nodes(polypeptide, entry):
    """
    Create InterPro Term Nodes.
    :param polypeptide:
    :param entry:
    :return:
    """
    # http://generic-model-organism-system-database.450254.n5.nabble.com/Re-GMOD-devel-Storing-Interpro-domains-in-Chado-td459778.html
    terms = [t for t in entry.split("; ") if t is not '']
    for interpro in terms:
        import time
        dbxref = DbXref(db="InterPro", accession=interpro, version=time.time())
        graph.create(dbxref)
        polypeptide.dbxref.add(dbxref)
        graph.push(polypeptide)


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


# TODO: Fetch data from PubMed

def update_pub_nodes():
    publications = Publication.select(graph)
    print(len(list(publications)))
    pmids = [publication.pmid for publication in publications]
    publication_by_id = dict(zip(pmids, publications))
    num_ids = len(pmids)
    chunksize = 500
    records = []
    for start in range(0, num_ids, chunksize):
        subset = pmids[start:start + chunksize]
        records.extend(fetch_publication_list(subset))
    record_loaded_count = 0
    for record in tqdm(records):
        if len(record) < 2:
            pm_id = record['id:'][0][record['id:'][0].find('able: ') + 6:]
            print("PMID: {}".format(pm_id))
            record = fetch_publication_list(pm_id, rettype='xml')
            rec = next(record)
            article = rec['MedlineCitation']['Article']
            title = article['ArticleTitle']
            pages = article['Pagination']['MedlinePgn']
            volume = article['Journal']['JournalIssue']['Volume']
            issue = article['Journal']['JournalIssue']['Issue']
            date_of_pub = article['Journal']['JournalIssue']['PubDate']['Month'] + " " + \
                          article['Journal']['JournalIssue']['PubDate']['Year']
            pub_place = rec['MedlineCitation']['MedlineJournalInfo']['Country']
            publisher = None
            author = None
            # full_author = article['AuthorList']
            full_author = None
        else:
            # https://www.nlm.nih.gov/bsd/mms/medlineelements.html
            pm_id = record['PMID']
            # there is internal caching so using a dictionary here doesn't
            # actually seem to save any time - pvh
            title = record.get('TI', None)
            volume = record.get('VI', None)
            issue = record.get('IP', None)
            pages = record.get('PG', None)
            date_of_pub = record.get('DP', None)
            pub_place = record.get('PL', None)
            publisher = record.get('SO', None)
            author = record.get('AU', None)
            full_author = record.get('FAU', None)

        publication = publication_by_id[pm_id]  # Publication.select(graph, pm_id).first()
        publication.title = title
        publication.volume = volume
        publication.issue = issue
        publication.pages = pages
        publication.year = date_of_pub
        publication.pubplace = pub_place
        publication.publisher = publisher
        graph.push(publication)
        create_author_nodes(publication, full_author)
        record_loaded_count += 1


def create_pub_nodes(polypeptide, pubs):
    """
    Create Publication Nodes
    :param polypeptide:
    :param pubs:
    :return:
    """
    citations = [c for c in pubs.split("; ") if c is not '']

    for citation in citations:
        pub = Publication()
        pub.pmid = citation

        polypeptide.published_in.add(pub)
        graph.push(polypeptide)


def create_is_a_cv_term_rel():
    """
    Creating IS_A relationships between CVTerms
    :return:
    """
    cv_terms = CvTerm.select(graph)
    for cv in cv_terms:
        is_a_list = fetch_quick_go_data(cv.name)
        # cv = CvTerm.select(graph, _id).first()
        for go in is_a_list:
            goid = go[go.find('G'):go.find('!')].strip()
            cv_term = CvTerm.select(graph, goid).first()
            if cv_term:
                cv.is_a.add(cv_term)
                graph.push(cv)


def build_protein_interaction_rels(protein_interaction_dict):
    """
    Build protein-protein interactions
    :param protein_interaction_dict:
    :return:
    """
    for uni_id, interactors in protein_interaction_dict.items():
        if len(interactors) > 0:
            poly = Polypeptide.select(graph, uni_id).first()
            interactors = interactors.split('; ')
            for interactor in interactors:
                if interactor == 'Itself':
                    interactor = poly.uniquename
                _poly = Polypeptide.select(graph, interactor).first()
                if _poly is None:
                    print("No Polypeptide with uniquename: {}".format(interactor))
                    # time.sleep(2)
                else:
                    poly.interacts_with.add(_poly)
                    graph.push(poly)


def build_gene_protein_relationship(locus_tag, polypeptide):
    """
    Build a connection between CDS and Polypeptide
    :param locus_tag:
    :param polypeptide:
    :return:
    """
    try:
        cds_name = graph.run(
            """MATCH (c:CDS) -[:PART_OF]-> (Transcript) -[:PART_OF]-> (g:Gene)
            WHERE g.uniquename = 'gene:{}' RETURN c.uniquename AS name""".format(locus_tag)).next()['name']
        # print("CDS name:", cds_name, file=sys.stderr)
    except StopIteration:
        cds_name = None
    else:
        cds = CDS.select(graph, cds_name).first()
        if cds:
            locus_found = True
            print("found CDS", cds_name, file=sys.stderr)
            # Polypetide-derives_from->CDS
            polypeptide.derives_from.add(cds)
            cds.polypeptide.add(polypeptide)
            graph.push(polypeptide)
            graph.push(cds)


def create_uniprot_nodes(uniprot_data, add_protein_interactions=True, proteome_name='h37rv', locus_tag_re=re.compile('(Rv\d+\w?(?:\.\d)?\w?)')):
    """
    Build DbXref nodes from UniProt results.
    :param uniprot_data:
    :return:
    """
    print("=========================================")
    print("About to create Nodes from UniProt data.")
    print("=========================================")
    # time.sleep(2)
    count = 0
    protein_interaction_dict = dict()
    featureset = FeatureSet().select(graph).where("_.name = '{}'".format(proteome_name)).first()
    for entry in tqdm(uniprot_data):
        protein_interaction_dict[entry[0]] = entry[6]
        count += 1

        dbxref = DbXref(db="UniProt", accession=entry[1], version=entry[0])
        graph.create(dbxref)
        pdb_id = map_ue_to_pdb(entry[0])
        polypeptide = Polypeptide()
        if featureset is not None:
            featureset.contains.add(polypeptide)
        polypeptide.name = entry[9]
        polypeptide.uniquename = entry[0]
        polypeptide.ontology_id = polypeptide.so_id
        polypeptide.seqlen = entry[16]
        polypeptide.residues = entry[14]
        polypeptide.parent = entry[2]
        polypeptide.family = entry[17]
        polypeptide.function = entry[13]
        polypeptide.pdb_id = pdb_id
        polypeptide.mass = entry[15]
        polypeptide.three_d = entry[12]
        graph.create(polypeptide)
        genes_oln = entry[2]
        if genes_oln.strip() == '':
            # fallback to the genes column
            genes_oln = entry[3]
        # print("genes_oln", genes_oln, file=sys.stderr)
        locus_found = False
        for locus_tag in locus_tag_re.findall(genes_oln):
            print("locus_tag", locus_tag, file=sys.stderr)
            build_gene_protein_relationship(locus_tag, polypeptide)
        if not locus_found:
            for locus_tag in locus_tag_re.findall(entry[3]):
                print("second try, locus_tag", locus_tag, file=sys.stderr)
                build_gene_protein_relationship(locus_tag, polypeptide)

            # gene = Gene.select(graph, "gene:" + locus_tag).first()
            # if gene:
            #     try:
            #         _feature = next(iter(gene.is_a))
            #     except StopIteration:
            #         pass
            #     else:
            #         print("found feature", _feature.uniquename, file=sys.stderr)
            #         try:
            #             transcript_feature = next(iter(_feature.related_to))
            #         except StopIteration:
            #             pass
            #         else:
            #             transcript = Transcript.select(graph, transcript_feature.uniquename).first()
            #             if transcript:
            #                 print("found transcript", file=sys.stderr)
            #                 cds = CDS.select(graph, "CDS" + transcript.uniquename[transcript.uniquename.find(":"):]).first()
            #                 if cds:
            #                     print("found CDS", file=sys.stderr)
            #                     # Polypetide-derives_from->CDS
            #                     polypeptide.derives_from.add(cds)
            #                     cds.polypeptide.add(polypeptide)
            #                     graph.push(polypeptide)
            #                     graph.push(cds)
            # else:
            #     print("Gene gene:{} not found".format(locus_tag), file=sys.stderr)

        polypeptide.dbxref.add(dbxref)
        graph.push(polypeptide)

        create_cv_term_nodes(polypeptide, entry[18], entry[19], entry[20])
        create_interpro_term_nodes(polypeptide, entry[5])
        create_pub_nodes(polypeptide, entry[11])
    if featureset is not None:
        graph.push(featureset)
    if "uniquename" not in graph.schema.get_indexes("Polypeptide"):
        graph.schema.create_index("Polypeptide", "uniquename")
    if add_protein_interactions:
        build_protein_interaction_rels(protein_interaction_dict)
    print ("TOTAL:", count)

def create_chromosome(seqrecord, name, featureset_name):
    chromosome = Chromosome()
    chromosome.residues = str(seqrecord.seq)
    chromosome.seqlen = len(seqrecord.seq)
    chromosome.name = chromosome.uniquename = name
    graph.create(chromosome)
    featureset = FeatureSet().select(graph).where("_.name = '{}'".featureset_name).first()
    if featureset is not None:
        featureset.contains.add(chromosome)
        graph.push(featureset)

def create_feature_indexes():
    for label in ("CDS", "Transcript", "Gene", "Exon", "PseudoGene"):
        if "uniquename" not in graph.schema.get_indexes(label):
            graph.schema.create_index(label, "uniquename")
    if "srcfeature_id" not in graph.schema.get_indexes("FeatureLoc"):
        graph.schema.create_index("FeatureLoc", "srcfeature_id")
    if "parent" not in graph.schema.get_indexes("Feature"):
        graph.schema.create_index("Feature", "parent")
