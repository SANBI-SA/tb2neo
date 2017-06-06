"""
Interface CLI commands.
"""

from __future__ import print_function
import click
import atexit
import os
import logging
import Bio.SeqIO
from .dbconn import GraphDb
from .gffproc import examine, parse_gff, get_locus_tags
from .uniprot import query_uniprot
from .docker import launch_neo4j_docker, find_docker_portmapping, kill_docker


class Context(object):
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


@click.group()
@click.option('--docker/--no-docker', default=False)
@click.option('--outputdir', type=click.Path(dir_okay=True))
@click.option('--dbhost', type=str, default='127.0.0.1')
@click.option('--dbpassword', type=str, default='')
@click.pass_context
def cli(ctx, docker, outputdir, dbhost, dbpassword):
    """
    This script parses a GFF file and builds a Neo4j Graph database.
    """

    logging.basicConfig(level=os.environ.get('GOGET_LOGLEVEL', 'INFO'))
    if docker:
        if outputdir is None:
            outputdir = os.getcwd()
        (proc, container_name) = launch_neo4j_docker(outputdir)
        atexit.register(kill_docker, proc, container_name)
        port_mapping = find_docker_portmapping(container_name)
        bolt_port = port_mapping[7687]
        http_port = port_mapping[7474]
        print("docker mode http port: {} bolt port: {}".
              format(http_port, bolt_port))
        db = GraphDb(host='localhost', password='',
                     bolt_port=bolt_port, http_port=http_port)
    else:
        # do "default" connection
        print("non-docker mode, connecting to host:", dbhost)
        db = GraphDb(host=dbhost, password=dbpassword)
    ctx.db = db


@cli.command()
@click.argument('gff_file', type=click.Path(exists=True, file_okay=True))
@click.option('-d/-D', '--delete_all/--no_delete_all', default=False,
              prompt='Delete existing database?',
              help='Delete existing data.')
@click.option('-r/-R', '--relationships/--no_relationships', default=False,
              prompt='Build node relationships?',
              help='Build node relationships.')
@click.option('-u/-U', '--uniprot/--no_uniprot', default=True,
              prompt='Query UniProt?',
              help='Query UniProt using locus tags.')
@click.option('-p/-P', '--publications/--no_publications', default=True,
              prompt='Query PubMed?',
              help='Query PubMed using pmid.')
@click.option('-g/-G', '--map_go/--no_map_go', default=True, is_flag=True,
              prompt='Map GO Terms?',
              help='Query QuickGo to map GO is_a relationships.')
@click.option('-o/-O', '--add_cdc1551_orthologs/--no_add_cdc1551_orthologs',
              default=True, is_flag=True, prompt='Add CDC1551 ortologs?',
              help='Add CDC1551 orthologs identified by Tuberculist')
@click.pass_context
def init(ctx, gff_file, delete_all, relationships, uniprot, publications,
         map_go, add_cdc1551_orthologs):
    """
    Load GFF features to Neo4j Graph database.
    :param map_go:
    :param publications:
    :param gff_file:
    :param delete_all:
    :param relationships:
    :param uniprot:
    :param add_cdc1551_orthologs:
    :return:
    """
    db = ctx.parent.db
    if (delete_all and relationships and not uniprot and
       not publications and not map_go):
        # Deleting existing data, load features and build relationships
        db.delete_data()
        parse_gff(db, gff_file)
        db.build_relationships()
    elif (delete_all and not relationships and
          not uniprot and not publications and not map_go):
        # Deleting existing data, load features
        db.delete_data()
        parse_gff(db, gff_file)
    elif (delete_all and relationships and uniprot and
          not publications and not map_go):
        # Deleting existing data, load features, build relationships,
        # fetch data from UniProt and create nodes, and build
        # relationships then update Publication nodes with data from PubMed
        db.delete_data()
        parse_gff(db, gff_file)
        db.build_relationships()
        db.create_uniprot_nodes(query_uniprot(get_locus_tags(gff_file, 400)))
    elif delete_all and relationships and uniprot and publications and map_go:
        # Deleting existing data, load features, build relationships,
        # fetch data from UniProt and create nodes, build relationships
        # then update Publication nodes with data from PubMed and map GO
        db.delete_data()
        parse_gff(db, gff_file)
        db.build_relationships()
        db.create_uniprot_nodes(query_uniprot(get_locus_tags(gff_file, 400)))
        db.update_pub_nodes()
        db.create_is_a_cv_term_rel()
    elif (delete_all and relationships and uniprot and
          publications and not map_go):
        # Deleting existing data, load features, build relationships,
        # fetch data from UniProt and create nodes, build relationships
        # then update Publication nodes with data  from PubMed
        db.delete_data()
        parse_gff(db, gff_file)
        db.build_relationships()
        db.create_uniprot_nodes(query_uniprot(get_locus_tags(gff_file, 400)))
        db.update_pub_nodes()
    elif (delete_all and relationships and uniprot and
          not publications and map_go):
        # Deleting existing data, load features, build relationships,
        # fetch data from UniProt and create nodes, build relationships
        # then map GO
        db.delete_data()
        parse_gff(db, gff_file)
        db.build_relationships()
        db.create_uniprot_nodes(query_uniprot(get_locus_tags(gff_file, 400)))
        db.create_is_a_cv_term_rel()
    elif (not delete_all and not relationships and
          not uniprot and publications and map_go):
        # Build relationships then map GO and update Publication
        # nodes with data from PubMed
        db.create_is_a_cv_term_rel()
        db.update_pub_nodes()
    elif (not delete_all and not uniprot and not publications and
          not map_go and relationships):
        # Build relationships from existing data
        db.build_relationships()
    elif (not delete_all and not relationships and not publications and
          not map_go and uniprot):
        # Fetch data from UniProt and create nodes
        db.create_uniprot_nodes(query_uniprot(get_locus_tags(gff_file, 400)))
    elif (not delete_all and not relationships and
          not uniprot and not map_go and publications):
        # Update Publication nodes with data from PubMed
        db.update_pub_nodes()
    elif (not delete_all and not relationships and not uniprot and
          not publications and map_go):
        # Create is_a relationship between GO using QuickGo
        db.create_is_a_cv_term_rel()
    elif (not delete_all and not relationships and uniprot and
          publications and map_go):
        # Build relationships, fetch data from UniProt and create nodes,
        # build relationships then map GO
        db.create_uniprot_nodes(query_uniprot(get_locus_tags(gff_file, 400)))
        db.update_pub_nodes()
        db.create_is_a_cv_term_rel()
    elif add_cdc1551_orthologs:
        db.add_orthologs_to_db()
    else:
        click.Abort()


@cli.command()
@click.pass_context
def delete_all(ctx):
    ctx.parent.db.delete_data()


@cli.command()
@click.argument('gff_file', type=click.Path(exists=True, file_okay=True))
@click.pass_context
def load_gff(ctx, gff_file):
    parse_gff(ctx.parent.db, gff_file)


@cli.command()
@click.pass_context
def build_rels(ctx):
    ctx.parent.db.build_relationships()


@cli.command()
@click.argument('gff_file', type=click.Path(exists=True, file_okay=True))
@click.pass_context
def load_uniprot(ctx, gff_file):
    ctx.parent.db.create_uniprot_nodes(query_uniprot(
                                       get_locus_tags(gff_file, 400)))


@cli.command()
@click.pass_context
def pub_nodes(ctx):
    ctx.parent.db.update_pub_nodes()


@cli.command()
@click.pass_context
def connect_go_terms(ctx):
    ctx.parent.db.create_is_a_cv_term_rel()


@cli.command()
@click.pass_context
def add_cdc1551_orthologs(ctx):
    ctx.parent.db.add_orthologs_to_db()


@cli.command()
@click.option('--name')
@click.option('--featureset')
@click.argument('fasta_file', type=click.File())
@click.pass_context
def add_chromosome(ctx, fasta_file, name=None, featureset='h37rv'):
    """
    Read in a FASTA file and add it as a Chromosome record in the database.
    :param fasta_file: click.File
    :param name: str
    :return:
    """
    seq = Bio.SeqIO.read(fasta_file, 'fasta')
    if name is None:
        name = seq.id
    ctx.parent.db.create_chromosome(seq, name, featureset)


@cli.command()
@click.argument('gff_file', type=click.Path(exists=True, file_okay=True))
def inspect(gff_file):
    """
    Inspect GFF file.
    :param gff_file:
    :return:
    """
    examine(gff_file)


if __name__ == '__main__':
    cli(Context())
