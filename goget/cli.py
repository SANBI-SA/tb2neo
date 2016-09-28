import click
from dbconn import *
from gffproc import examine, parse_gff, get_locus_tags
from uniprot import query_uniprot
from time import time


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
@click.option('--d', '--delete', default=False, is_flag=True, prompt='Delete existing database?',
              help='Delete existing data.')
@click.option('--r', '--relationships', default=False, is_flag=True, prompt='Build node relationships?',
              help='Build node relationships.')
@click.option('--u', '--uniprot', default=True, is_flag=True, prompt='Query UniProt?',
              help='Query UniProt using locus tags.')
def init(gff_file, delete, relationships, uniprot):
    """
    Load features from GFF file.
    :param gff_file:
    :param delete:
    :param relationships:
    :param uniprot:
    :return:
    """
    # Deleting existing data, load features and build relationships or
    # build relationships from existing data
    if delete and relationships:
        start = time()
        delete_data()
        parse_gff(gff_file)
        build_relationships()
        end = time()
        print("Done in ", end - start, "secs.")
    elif not delete and relationships:
        start = time()
        build_relationships()
        end = time()
        print("Done in ", end - start, "secs.")
    elif delete and not relationships:
        start = time()
        delete_data()
        parse_gff(gff_file)
        end = time()
        print("Done in ", end - start, "secs.")
    elif delete and relationships and uniprot:
        start = time()
        delete_data()
        parse_gff(gff_file)
        build_relationships()
        query_uniprot(get_locus_tags(gff_file, 400))
        end = time()
        print("Done in ", end - start, "secs.")
    elif not delete and not relationships and uniprot:
        start = time()
        query_uniprot(get_locus_tags(gff_file, 400))
        end = time()
        print("Done in ", end - start, "secs.")
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
