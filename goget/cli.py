import time

import click
from createdb import *
from parsegff import examine, parse_gff


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
    This script parses in features from a GFF file and builds a Neo4j Graph database.
    """
    pass


@cli.command()
@click.argument('gff_file', type=click.Path(exists=True, file_okay=True))
@click.option('-rel', '--relations', default=False, is_flag=True, help='Build node relationships.')
def gff(relations, gff_file):
    """
    Load features from GFF file.
    :param gff_file:
    :param relations:
    :return:
    """
    # Deleting existing data
    delete_data()
    parse_gff(gff_file)
    if relations:
        build_relationships()


@cli.command()
@click.argument('gff_file', type=click.Path(exists=True, file_okay=True))
def inspect(gff_file):
    """
    Examine GFF file for the type of data present.
    :param gff_file:
    :return:
    """
    examine(gff_file)


@cli.command()
def relationships():
    """
    Build relationships between loaded features.
    :return:
    """
    build_relationships()


if __name__ == '__main__':
    time.sleep(20)
    cli()
