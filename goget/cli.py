"""
Interface CLI commands.
"""

import click
from dbconn import build_relationships, update_pub_nodes, graph, create_is_a_cv_term_rel
from gffproc import examine, parse_gff, get_locus_tags
from uniprot import query_uniprot


def delete_data():
    """
    Delete existing data
    :return:
    """
    print("Deleting all nodes and relationships in {}".format(graph))
    graph.delete_all()


@click.group()
def cli():
    """
    This script parses a GFF file and builds a Neo4j Graph database.
    """
    pass


@cli.command()
@click.argument('gff_file', type=click.Path(exists=True, file_okay=True))
@click.option('--d', '--delete_all', default=False, is_flag=True, prompt='Delete existing database?',
              help='Delete existing data.')
@click.option('--r', '--relationships', default=False, is_flag=True, prompt='Build node relationships?',
              help='Build node relationships.')
@click.option('--u', '--uniprot', default=True, is_flag=True, prompt='Query UniProt?',
              help='Query UniProt using locus tags.')
@click.option('--p', '--publications', default=True, is_flag=True, prompt='Query PubMed?',
              help='Query PubMed using pmid.')
@click.option('--g', '--map_go', default=True, is_flag=True, prompt='Map GO Terms?',
              help='Using QuickGo to map GO is_a relationships.')
def init(gff_file, delete_all, relationships, uniprot, publications, map_go):
    """
    Load features from GFF file.
    :param map_go:
    :param publications:
    :param gff_file:
    :param delete_all:
    :param relationships:
    :param uniprot:
    :return:
    """
    if delete_all and relationships and not uniprot and not publications and not map_go:
        # Deleting existing data, load features and build relationships
        delete_data()
        parse_gff(gff_file)
        build_relationships()
    elif delete_all and not relationships and not uniprot and not publications and not map_go:
        # Deleting existing data, load features
        delete_data()
        parse_gff(gff_file)
    elif delete_all and relationships and uniprot and not publications and not map_go:
        # Deleting existing data, load features, build relationships, fetch data from UniProt and create nodes,
        # and build relationships then update Publication nodes with data from PubMed
        delete_data()
        parse_gff(gff_file)
        build_relationships()
        query_uniprot(get_locus_tags(gff_file, 400))
    elif delete_all and relationships and uniprot and publications and map_go:
        # Deleting existing data, load features, build relationships, fetch data from UniProt and create nodes,
        # build relationships then update Publication nodes with data from PubMed and map GO
        delete_data()
        parse_gff(gff_file)
        build_relationships()
        query_uniprot(get_locus_tags(gff_file, 400))
        update_pub_nodes()
        create_is_a_cv_term_rel()
    elif delete_all and relationships and uniprot and not publications and map_go:
        # Deleting existing data, load features, build relationships, fetch data from UniProt and create nodes,
        # build relationships then map GO
        delete_data()
        parse_gff(gff_file)
        build_relationships()
        query_uniprot(get_locus_tags(gff_file, 400))
        create_is_a_cv_term_rel()
    elif not delete_all and not uniprot and not publications and not map_go and relationships:
        # Build relationships from existing data
        build_relationships()
    elif not delete_all and not relationships and not publications and not map_go and uniprot:
        # Fetch data from UniProt and create nodes
        query_uniprot(get_locus_tags(gff_file, 400))
    elif not delete_all and not relationships and not uniprot and not map_go and publications:
        # Update Publication nodes with data from PubMed
        update_pub_nodes()
    elif not delete_all and not relationships and not uniprot and not publications and map_go:
        # Create is_a relationship between GO using QuickGo
        create_is_a_cv_term_rel()
    else:
        click.Abort()


@cli.command()
@click.argument('gff_file', type=click.Path(exists=True, file_okay=True))
def inspect(gff_file):
    """
    Examine GFF file.
    :param gff_file:
    :return:
    """
    examine(gff_file)


if __name__ == '__main__':
    cli()
