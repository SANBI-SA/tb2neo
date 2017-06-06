"""
Interface to the Neo4j Database
"""

from __future__ import print_function
import re
import socket
import logging
import time
import os.path
import json
from tqdm import tqdm
from combat_tb_model.model.core import Organism, Gene, Exon, \
    Transcript, PseudoGene, TRna, NCRna, RRna, CDS, Polypeptide, Feature, \
    FeatureLoc, FeatureSet, CvTerm, DbXref, Publication, Chromosome, \
    Author
from goget.ncbi import fetch_publication_list
from py2neo import Graph, watch
from .quickgo import fetch_quick_go_data
from .uniprot import map_ue_to_pdb
from .orthologs import fetch_ortholog

# graph = Graph(host=getenv("DB", "localhost"), bolt=True,
#               password=getenv("NEO4J_PASSWORD", ""))


class GraphDb(object):
    def __init__(self, host, password=None,
                 bolt_port=7687, http_port=7474,
                 debug=False):
        if password is None:
            password = ''
        self.debug = debug
        self.graph = self.connect(host, password, bolt_port, http_port)

    def connect(self, host, password, bolt_port, http_port, timeout=30):
        """connect - make connection to Neo4j DB

        :type host: str - hostname or IP of Neo4j database server
        :type password: str - password for Neo4j database server
        :type bolt_port: int - port for Neo4j Bolt protocol
        :type http_port: int - port for Neo4j HTTP protocol
        :type timeout: int - timeout for waiting for the Neo4j connection"""

        connected = False
        while timeout > 0:
            try:
                socket.create_connection((host, http_port), 1)
            except socket.error:
                timeout -= 1
                time.sleep(1)
            else:
                connected = True
                break
        if not connected:
            raise socket.timeout('timed out trying to connect to {}'.format(
                     host, http_port))
        logging.debug(
            "connecting graph to http port: {} bolt_port: {} host: {}".format(
                http_port, bolt_port, host))
        uri = 'bolt://{}:{}'.format(host, bolt_port)
        graph = Graph(uri, bolt=True, password=password,
                      bolt_port=bolt_port, http_port=http_port)
        if self.debug:
            watch("neo4j.bolt")
        logging.debug("connected", graph)
        return graph

    def create_organism_nodes(self):
        """
        Create Organism Nodes
        :return:
        """
        abbrev = "H37Rv"
        genus = "Mycobacterium"
        species = "M. tuberculosis"
        common_name = "TB"

        organism = Organism(abbreviation=abbrev, genus=genus,
                            species=species, common_name=common_name)
        self.graph.create(organism)

    def create_gene_nodes(self, feature, transaction=None):
        """
        Create Gene Nodes
        :param feature:
        :return:
        """
        names = self.get_feature_name(feature)
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
            self.graph.create(gene)

    def create_transcript_nodes(self, feature, transaction=None):
        """
        Create Transcipt Nodes
        :param feature:
        :return:
        """
        names = self.get_feature_name(feature)
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
            self.graph.create(transcript)

    def create_pseudogene_nodes(self, feature, transaction=None):
        """
        Create Pseudogene Nodes
        :param feature:
        :return:
        """
        names = self.get_feature_name(feature)
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
            self.graph.create(pseudogene)

    def create_exon_nodes(self, feature, transaction=None):
        """
        Create Exon Nodes
        :param feature:
        :return:
        """
        names = self.get_feature_name(feature)
        name = names.get("Name", names.get("UniqueName"))
        unique_name = names.get("UniqueName", name)

        exon = Exon()
        exon.ontology_id = exon.so_id
        exon.name = name
        exon.uniquename = unique_name
        if transaction is not None:
            transaction.create(exon)
        else:
            self.graph.create(exon)

    def create_rna_nodes(self, feature, transaction=None):
        """
        Create RNA Nodes
        :param feature:
        :return:
        """
        names = self.get_feature_name(feature)
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
                self.graph.create(trna)
        if feature.type == 'ncRNA_gene':
            ncrna = NCRna()
            ncrna.ontology_id = ncrna.so_id
            ncrna.name = name
            ncrna.uniquename = unique_name
            if transaction is not None:
                transaction.create(ncrna)
            else:
                self.graph.create(ncrna)
        if feature.type == 'rRNA_gene':
            rrna = RRna()
            rrna.ontology_id = rrna.so_id
            rrna.name = name
            rrna.uniquename = unique_name
            if transaction is not None:
                transaction.create(rrna)
            else:
                self.graph.create(rrna)

    def create_cds_nodes(self, feature, transaction=None):
        """
        Create CDS Nodes
        :param feature:
        :return:
        """
        names = self.get_feature_name(feature)
        name = names.get("Name", names.get("UniqueName"))
        unique_name = names.get("UniqueName", name)

        cds = CDS()
        cds.name = name
        cds.uniquename = unique_name
        cds.ontology_id = cds.so_id
        if transaction is not None:
            transaction.create(cds)
        else:
            self.graph.create(cds)

    def create_feature_nodes(self, feature, transaction=None):
        """
        Create Feature Nodes
        :param feature:
        :return:
        """
        names = self.get_feature_name(feature)
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
            self.graph.create(_feature)

    def create_featureloc_nodes(self, feature, transaction=None):
        """
        Create FeatureLoc Nodes
        :param feature:
        :return:
        """
        srcfeature_id = self.get_feature_name(feature).get("UniqueName")
        feature_loc = FeatureLoc(srcfeature_id=srcfeature_id,
                                 fmin=feature.location.start,
                                 fmax=feature.location.end,
                                 strand=feature.location.strand)
        if transaction is not None:
            transaction.create(feature_loc)
        else:
            self.graph.create(feature_loc)

    def get_feature_name(self, feature):
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

    def build_relationships(self):
        """
        Build relationships
        :return:
        """
        # TODO: Try optimize this
        logging.info("Building Relationships...")
        feature_count = self.graph.run(
            "MATCH (f:Feature) RETURN count(f) AS count").next()['count']
        # THIS was the code before, not sure why the select
        # was there - pvh 6/6/2017
        # features = Feature.select(self.graph).where(
        #     "_.uniquename = 'transcript:CCP46448'")
        features = Feature.select(self.graph)
        featureset = FeatureSet().select(self.graph).where(
            "_.name = 'h37rv'").first()
        if featureset is None:
            featureset = FeatureSet()
            featureset.name = 'h37rv'
            featureset.description = 'M. tuberculosis H37Rv'
            logging.info("Creating Featureset for {}".format(featureset.name))
            self.graph.create(featureset)
        org = Organism.select(self.graph, 'Mycobacterium').first()
        for feature in tqdm(features, total=feature_count):
            # Find organism via __primarykey__ and link with feature via
            # BELONGS_TO
            if org is not None:
                feature.belongs_to.add(org)
            featureset.contains.add(feature)

            # Find feature with a parent attr. matching this feature's
            # uniquename and link them via RELATED_TO
            _feature = Feature.select(self.graph).where(
                "_.parent = '{}'".format(feature.uniquename)).first()
            logging.debug(feature, _feature)

            if _feature:
                feature.related_to.add(_feature)

            # Building is_a relationships
            gene = Gene.select(self.graph, feature.uniquename).first()
            if gene:
                gene.is_a.add(feature)
                featureset.contains.add(gene)
                self.graph.push(gene)
                # Find feature with this gene's uniquename as a parent
                _feature = Feature.select(self.graph).where(
                    "_.parent = '{}'".format(gene.uniquename)).first()
                if _feature:
                    # Find transcript: A gene is a parent to it.
                    transcript = Transcript.select(self.graph,
                                                   _feature.uniquename).first()
                    if transcript:
                        transcript.part_of.add(gene)
                        self.graph.push(transcript)

            p_gene = PseudoGene.select(self.graph, feature.uniquename).first()
            if p_gene:
                featureset.contains.add(p_gene)
                p_gene.is_a.add(feature)
                self.graph.push(p_gene)
            transcript = Transcript.select(self.graph,
                                           feature.uniquename).first()
            if transcript:
                featureset.contains.add(transcript)
                transcript.is_a.add(feature)
                self.graph.push(transcript)
                # Find feature with this transcript's uniquename as a parent
                _feature = Feature.select(self.graph).where(
                    "_.parent = '{}'".format(transcript.uniquename)).first()
                if _feature:
                    # Find exon: A transcript is a parent to it
                    exon = Exon.select(self.graph, _feature.uniquename).first()
                    if exon:
                        exon.part_of.add(transcript)
                        self.graph.push(exon)
                    # Find cds: A transcript is a parent to it
                    cds = CDS.select(self.graph, _feature.uniquename).first()
                    if cds:
                        cds.part_of.add(transcript)
                        self.graph.push(cds)

            exon = Exon.select(self.graph, feature.uniquename).first()
            if exon:
                featureset.contains.add(exon)
                exon.is_a.add(feature)
                self.graph.push(exon)
            cds = CDS.select(self.graph, feature.uniquename).first()
            if cds:
                featureset.contains.add(cds)
                cds.is_a.add(feature)
                self.graph.push(cds)

            # Find feature location with a srcfeature_id attr. matching this
            # features uniquename and link them via LOCATED_AT
            _feature = FeatureLoc.select(self.graph,
                                         feature.uniquename).first()
            if _feature:
                feature.location.add(_feature)
            self.graph.push(feature)
        if featureset is not None:
            self.graph.push(featureset)

    def create_cv_term_nodes(self, polypeptide, bp, cc, mf):
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
        go_bp_defs = [t[:t.find('[') - 1] for t in bp.split('; ')
                      if t is not '']
        # go(cellular component)
        go_cc_ids = [t[t.find('G'):-1] for t in cc.split('; ') if t is not '']
        go_cc_defs = [t[:t.find('[') - 1] for t in cc.split('; ')
                      if t is not '']
        # go(molecular function)
        go_mf_ids = [t[t.find('G'):-1] for t in mf.split('; ') if t is not '']
        go_mf_defs = [t[:t.find('[') - 1] for t in mf.split('; ')
                      if t is not '']

        # TODO: Find a way to refactor this.
        for _id in go_bp_ids:
            cv = CvTerm()
            for _def in go_bp_defs:
                cv.name = _id
                cv.definition = _def
                cv.namespace = "biological process"
                self.graph.create(cv)
                polypeptide.cvterm.add(cv)
                self.graph.push(polypeptide)

        for _id in go_mf_ids:
            cv = CvTerm()
            for _def in go_mf_defs:
                cv.name = _id
                cv.definition = _def
                cv.namespace = "cellular component"
                self.graph.create(cv)
                polypeptide.cvterm.add(cv)
                self.graph.push(polypeptide)
        for _id in go_cc_ids:
            cv = CvTerm()
            for _def in go_cc_defs:
                cv.name = _id
                cv.definition = _def
                cv.namespace = "molecular function"
                self.graph.create(cv)
                polypeptide.cvterm.add(cv)
                self.graph.push(polypeptide)

        if "name" not in self.graph.schema.get_indexes("CvTerm"):
            self.graph.schema.create_index("CvTerm", "name")

    def create_interpro_term_nodes(self, polypeptide, entry):
        """
        Create InterPro Term Nodes.
        :param polypeptide:
        :param entry:
        :return:
        """
        # http://generic-model-organism-system-database.450254.n5.nabble.com/Re-GMOD-devel-Storing-Interpro-domains-in-Chado-td459778.html
        terms = [t for t in entry.split("; ") if t is not '']
        for interpro in terms:
            dbxref = DbXref(db="InterPro", accession=interpro,
                            version=time.time())
            self.graph.create(dbxref)
            polypeptide.dbxref.add(dbxref)
            self.graph.push(polypeptide)

    def create_author_nodes(self, publication, full_author):
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
                self.graph.create(_author)
                _author.wrote.add(publication)
                publication.author.add(_author)
                self.graph.push(_author)
                self.graph.push(publication)

    # TODO: Fetch data from PubMed

    def update_pub_nodes(self):
        publications = Publication.select(self.graph)
        logging.info(len(list(publications)))
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
                logging.info("PMID: {}".format(pm_id))
                record = fetch_publication_list(pm_id, rettype='xml')
                rec = next(record)
                article = rec['MedlineCitation']['Article']
                title = article['ArticleTitle']
                pages = article['Pagination']['MedlinePgn']
                volume = article['Journal']['JournalIssue']['Volume']
                issue = article['Journal']['JournalIssue']['Issue']
                date_of_pub = \
                    article['Journal']['JournalIssue']['PubDate']['Month'] \
                    + " " \
                    + article['Journal']['JournalIssue']['PubDate']['Year']
                pub_place = \
                    rec['MedlineCitation']['MedlineJournalInfo']['Country']
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

            publication = publication_by_id[pm_id]
            publication.title = title
            publication.volume = volume
            publication.issue = issue
            publication.pages = pages
            publication.year = date_of_pub
            publication.pubplace = pub_place
            publication.publisher = publisher
            self.graph.push(publication)
            self.create_author_nodes(publication, full_author)
            record_loaded_count += 1

    def create_pub_nodes(self, polypeptide, pubs):
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
            self.graph.push(polypeptide)

    def create_is_a_cv_term_rel(self):
        """
        Creating IS_A relationships between CVTerms
        :return:
        """
        cv_terms = CvTerm.select(self.graph)
        for cv in cv_terms:
            is_a_list = fetch_quick_go_data(cv.name)
            # cv = CvTerm.select(graph, _id).first()
            for go in is_a_list:
                goid = go[go.find('G'):go.find('!')].strip()
                cv_term = CvTerm.select(self.graph, goid).first()
                if cv_term:
                    cv.is_a.add(cv_term)
                    self.graph.push(cv)

    def build_protein_interaction_rels(self, protein_interaction_dict):
        """
        Build protein-protein interactions
        :param protein_interaction_dict:
        :return:
        """
        for uni_id, interactors in protein_interaction_dict.items():
            if len(interactors) > 0:
                poly = Polypeptide.select(self.graph, uni_id).first()
                interactors = interactors.split('; ')
                for interactor in interactors:
                    if interactor == 'Itself':
                        interactor = poly.uniquename
                    _poly = Polypeptide.select(self.graph, interactor).first()
                    if _poly is None:
                        logging.info("No Polypeptide with uniquename: {}".
                                     format(interactor))
                        # time.sleep(2)
                    else:
                        poly.interacts_with.add(_poly)
                        self.graph.push(poly)

    def build_gene_protein_relationship(self, locus_tag, polypeptide):
        """
        Build a connection between CDS and Polypeptide
        :param locus_tag:
        :param polypeptide:
        :return:
        """
        try:
            cds_name = self.graph.run(
                """MATCH (c:CDS) -[:PART_OF]-> (Transcript)
                -[:PART_OF]-> (g:Gene) WHERE g.uniquename =
                'gene:{}' RETURN c.uniquename
                AS name""".format(locus_tag)).next()['name']
            # print("CDS name:", cds_name, file=sys.stderr)
        except StopIteration:
            cds_name = None
        else:
            cds = CDS.select(self.graph, cds_name).first()
            if cds:
                # locus_found = True
                logging.debug("found CDS", cds_name)
                # Polypetide-derives_from->CDS
                polypeptide.derives_from.add(cds)
                cds.polypeptide.add(polypeptide)
                self.graph.push(polypeptide)
                self.graph.push(cds)

    def create_uniprot_nodes(self, uniprot_data,
                             add_protein_interactions=True,
                             proteome_name='h37rv',
                             locus_tag_re=re.compile(
                                '(Rv\d+\w?(?:\.\d)?\w?)')):
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
        featureset = FeatureSet().select(self.graph).\
            where("_.name = '{}'".format(proteome_name)).first()
        for entry in tqdm(uniprot_data):
            protein_interaction_dict[entry[0]] = entry[6]
            count += 1

            dbxref = DbXref(db="UniProt", accession=entry[1], version=entry[0])
            self.graph.create(dbxref)
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
            self.graph.create(polypeptide)
            genes_oln = entry[2]
            if genes_oln.strip() == '':
                # fallback to the genes column
                genes_oln = entry[3]
            # print("genes_oln", genes_oln, file=sys.stderr)
            locus_found = False
            for locus_tag in locus_tag_re.findall(genes_oln):
                logging.debug("locus_tag", locus_tag)
                self.build_gene_protein_relationship(locus_tag, polypeptide)
            if not locus_found:
                for locus_tag in locus_tag_re.findall(entry[3]):
                    logging.debug("second try, locus_tag", locus_tag)
                    self.build_gene_protein_relationship(locus_tag,
                                                         polypeptide)

                # gene = Gene.select(graph, "gene:" + locus_tag).first()
                # if gene:
                #     try:
                #         _feature = next(iter(gene.is_a))
                #     except StopIteration:
                #         pass
                #     else:
                #         print("found feature", _feature.uniquename,
                #               file=sys.stderr)
                #         try:
                #             transcript_feature = \
                #                 next(iter(_feature.related_to))
                #         except StopIteration:
                #             pass
                #         else:
                #             transcript = Transcript.select(graph,
                #                        transcript_feature.uniquename).first()
                #             if transcript:
                #                 print("found transcript", file=sys.stderr)
                #                 cds = CDS.select(graph, "CDS" + \
                #                           transcript.uniquename[
                #                               transcript.uniquename.find(
                #                                    ":"):]).first()
                #                 if cds:
                #                     print("found CDS", file=sys.stderr)
                #                     # Polypetide-derives_from->CDS
                #                     polypeptide.derives_from.add(cds)
                #                     cds.polypeptide.add(polypeptide)
                #                     graph.push(polypeptide)
                #                     graph.push(cds)
                # else:
                #     print("Gene gene:{} not found".format(locus_tag),
                #                                           file=sys.stderr)

            polypeptide.dbxref.add(dbxref)
            self.graph.push(polypeptide)

            self.create_cv_term_nodes(polypeptide, entry[18],
                                      entry[19], entry[20])
            self.create_interpro_term_nodes(polypeptide, entry[5])
            self.create_pub_nodes(polypeptide, entry[11])
        if featureset is not None:
            self.graph.push(featureset)
        if "uniquename" not in self.graph.schema.get_indexes("Polypeptide"):
            self.graph.schema.create_index("Polypeptide", "uniquename")
        if add_protein_interactions:
            self.build_protein_interaction_rels(protein_interaction_dict)
        logging.info("TOTAL:", count)

    def create_chromosome(self, seqrecord, name, featureset_name):
        chromosome = Chromosome()
        chromosome.residues = str(seqrecord.seq)
        chromosome.seqlen = len(seqrecord.seq)
        chromosome.name = chromosome.uniquename = name
        self.graph.create(chromosome)
        featureset = FeatureSet().select(self.graph).\
            where("_.name = '{}'".featureset_name).first()
        if featureset is not None:
            featureset.contains.add(chromosome)
            self.graph.push(featureset)

    def create_feature_indexes(self):
        for label in ("CDS", "Transcript", "Gene", "Exon", "PseudoGene"):
            if "uniquename" not in self.graph.schema.get_indexes(label):
                self.graph.schema.create_index(label, "uniquename")
        if "srcfeature_id" not in self.graph.schema.get_indexes("FeatureLoc"):
            self.graph.schema.create_index("FeatureLoc", "srcfeature_id")
        if "parent" not in self.graph.schema.get_indexes("Feature"):
            self.graph.schema.create_index("Feature", "parent")

    def delete_data(self):
        """
        Delete existing data
        :return:
        """
        logging.info("Deleting all nodes and relationships in {}".
                     format(self.graph))
        self.graph.delete_all()
        for label in ["CDS", "Transcript", "Gene", "Exon", "PseudoGene",
                      "Feature", "FeatureLoc", "Polypeptide", "CvTerm"]:
            for property_name in self.graph.schema.get_indexes(label):
                self.graph.schema.drop_index(label, property_name)

    def find_orthologs(self, ortholog_type='cdc1551', chunksize=400):
        # this yields lists of (locus_tag, gene_object) tuples
        # for known genes from the existing annotion, chunksize at a time
        cachefile = 'orthologs_cache.json'
        tmp_cachefile = cachefile + '.tmp'
        ortholog_cache = json.load(open(cachefile)) \
            if os.path.exists(cachefile) else dict()
        ortholog_list = []
        # print(ortholog_cache)
        count = 1
        # TODO: find a better solution than this hack to identify
        # H37Rv genes
        for gene in Gene.select(self.graph).where(
                '_.uniquename =~ "gene:Rv[0-9].*"'):
            # if count > 10:
            #     break
            # count += 1
            # print('.', end='', file=sys.stderr)
            associated_features = list(gene.is_a)
            if len(associated_features) == 0:
                logging.debug("No IS_A related features found for {}".
                              format(gene.uniquename))
            else:
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

    def add_orthologs_to_db(self, ortholog_type='cdc1551'):
        ortholog_group = FeatureSet()
        ortholog_group.name = ortholog_type
        for chunk in self.find_orthologs(ortholog_type):
            ortholog_names = []
            for (name, gene) in tqdm(chunk):
                ortholog = Gene()
                ortholog.uniquename = ortholog.name = ("gene:"+name)
                ortholog_feature = Feature()
                ortholog_feature.uniquename =\
                    ortholog_feature.name = ("gene:"+name)
                ortholog.is_a.add(ortholog_feature)
                self.graph.create(ortholog_feature)
                self.graph.create(ortholog)
                ortholog_group.contains.add(ortholog_feature)
                try:
                    feature = next(iter(gene.is_a))
                except StopIteration:
                    raise ValueError("Database error: gene {} has no Feature".
                                     format(gene.uniquename))
                feature.orthologous_to.add(ortholog_feature)
                self.graph.push(feature)
                ortholog_names.append(name)
                # this will create Protein / Polypeptide nodes for all
                # the proteins associated with these orthologs
                # commented out for new because we don't have enough
                # structure (Gene->CDS->Transcript) to associate Protein with
                # create_uniprot_nodes(query_uniprot(ortholog_names,
                # taxonomy='83331', proteome='UP000001020'),
                # add_protein_interactions=False)
        self.graph.create(ortholog_group)
        # exit(0) # for debug purposes, end here
