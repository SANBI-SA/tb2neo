from combat_tb_model.model import Organism, Feature, FeatureLoc
from py2neo import watch, Graph, getenv

graph = Graph(host=getenv("DB", "localhost"), bolt=True, password=getenv("NEO4J_PASSWORD", ""))
watch("neo4j.bolt")


def create_organism_nodes():
    abbrev = "H37Rv"
    genus = "Mycobacterium"
    species = "M. tuberculosis"
    common_name = "TB"

    organism = Organism(abbreviation=abbrev, genus=genus, species=species, common_name=common_name)
    graph.create(organism)


def create_feature_nodes(feature):
    """
    Create Feature Nodes
    :param feature:
    :return:
    """
    exclude = ['gene', 'pseudogene', 'tRNA_gene', 'ncRNA_gene', 'rRNA_gene']
    if feature.type in exclude:
        parent = None
    else:
        parent = feature.qualifiers['Parent'][0][feature.qualifiers['Parent'][0].find(":") + 1:]

    names = get_feature_name(feature)
    name = names.get("Name", names.get("UniqueName"))
    unique_name = names.get("UniqueName", name)

    _feature = Feature()
    _feature.name = name
    _feature.type = feature.type
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
        names["UniqueName"] = feature.id[feature.id.find(":") + 1:]
    else:
        names["Name"] = names["UniqueName"] = feature.id[feature.id.find(":") + 1:]
    return names


def build_relationships():
    """
    Build relationships
    :return:
    """
    watch("neo4j.bolt")
    print("Building Relationships...")
    features = Feature.select(graph)
    for feature in features:
        # Find organism via __primarykey__ and link with feature via BELONGS_TO
        org = Organism.select(graph, 'Mycobacterium').first()
        feature.belongs_to.add(org)
        # Find feature with a parent attr. matching this features uniquename and link them via RELATED_TO
        _feature = Feature.select(graph).where("_.parent = '{}'".format(feature.uniquename)).first()
        if _feature and feature.type is not _feature.type:
            feature.related.add(_feature)
        # Find feature location with a srcfeature_id attr. matching this features uniquename and link them via
        # LOCATED_AT
        _feature = FeatureLoc.select(graph, feature.uniquename).first()
        if _feature:
            feature.location.add(_feature)
        graph.push(feature)
