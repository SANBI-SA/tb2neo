from py2neo.ogm import GraphObject, Property, RelatedTo, RelatedFrom

# class Phenotype(GraphObject):
#     __primarykey__ = 'type'
#     # type {XDR, DR, MDR, SUS}
#     _type = Property()
#     has_var = RelatedFrom("Variant", "HAS_VAR")
# Adapting GA4GH Variant Data Model
# https://ga4gh-schemas.readthedocs.io/en/latest/api/variants.html
# VariantSet = Phenotype
# TODO: Dataset and ReferenceSet?
class VariantSet(GraphObject):
    __primarykey__ = 'name'
    name = Property()

    has_var = RelatedTo("VariantSite", "HAS_VAR")
    has_call = RelatedTo("Call", "HAS_CALL")

    def __init__(self, name):
        self.name = name


# VariantSite = Variant
class VariantSite(GraphObject):
    ### NOTE: relies on FeatureLoc from core.py
    # __primarykey__ = 'pos'

    pos = Property()
    feature_id = Property()
    biotype = Property()
    chrom = Property()
    ref_allele = Property()
    alt_allele = Property()
    quality = Property()
    depth = Property()
    gene = Property()
    gene_id = Property()
    known = Property()
    novel = Property()

    occurs_in = RelatedTo("Gene", "OCCURS_IN")
    location = RelatedTo("FeatureLoc", "LOCATED_AT")

    has_call = RelatedTo("Call", "HAS_CALL")
    belongs_to_vset = RelatedTo("VariantSet", "BELONGS_TO_VSET")

    def __init__(self, chrom, pos, ref_allele, alt_allele, gene=None):
        self.chrom = chrom
        self.pos = pos
        self.ref_allele = ref_allele
        self.alt_allele = alt_allele
        self.gene = gene


# CallSet = VCF file
class CallSet(GraphObject):
    __primarykey__ = 'name'
    name = Property()

    has_call = RelatedTo("Call", "HAS_CALL")
    has_calls_in = RelatedTo("VariantSet", "HAS_CALLS_IN")

    def __init__(self, name):
        self.name = name


class Call(GraphObject):
    # __primarykey__ = 'pos'
    genotype = Property()
    ref_allele = Property()
    alt_allele = Property()
    gene = Property()
    pos = Property()

    associated_with = RelatedTo("VariantSite", "ASSOC_WITH_VARIANT")
    belongs_to_cset = RelatedTo("CallSet", "BELONGS_TO_SET")

    def __init__(self, pos, ref_allele, alt_allele, gene=None):
        self.pos = pos
        self.ref_allele = ref_allele
        self.alt_allele = alt_allele
        self.gene = gene
