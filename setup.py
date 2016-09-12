from setuptools import setup

setup(
    name='goget',
    version='0.0.1',
    description='Parses GFF file and builds a graph database based on the features,'
                'It also maps these features to external services like UniProt using locus_tags.',
    keywords='neo4j, bioservices and gff',
    py_modules=['goget'],
    install_requires=[
        'click',
        'bioservices',
        'pandas',
        'bcbio-gff',
        'biopython',
        'beautifulsoup4'
    ],
    entry_points={
        'console_scripts': ['goget=goget:cli']
    },
)
