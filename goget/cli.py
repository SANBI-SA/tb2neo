import click
from dbconn import *
from gffproc import examine, parse_gff


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
@click.option('--d', '--delete', default=True, is_flag=True, prompt='Delete existing database?',
              help='Delete existing data.')
@click.option('--r', '--relationships', default=False, is_flag=True, prompt='Build node relationships?',
              help='Build node relationships.')
def init(gff_file, delete, relationships):
    """
    Load features from GFF file.
    :param gff_file:
    :param delete:
    :param relationships:
    :return:
    """
    # Deleting existing data, load features and build relationships or
    # build relationships from existing data
    if delete and relationships:
        delete_data()
        parse_gff(gff_file)
        build_relationships()
    elif not delete and relationships:
        build_relationships()
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
