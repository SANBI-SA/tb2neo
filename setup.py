from setuptools import setup

setup(
    name='tb2neo',
    version='0.0.7',
    description='Parses M. tuberculosis annotation (in GFF file and online '
                'sources) and builds a Neo4j graph database'
                'storing the annotate features. It also maps these features '
                'to external services such as UniProt.',
    keywords='neo4j, bioservices, gff',
    packages=['tb2neo'],
    install_requires=[
        'click',
        'bioservices',
        'bcbio-gff',
        'biopython',
        'beautifulsoup4'
    ],
    entry_points={
        'console_scripts': ['tb2neo=tb2neo.cli:cli']
    },
)
