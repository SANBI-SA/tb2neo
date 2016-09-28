# **goget**

Parses GFF file and builds a Neo4j Graph database.

## Usage

**Clone this repository:**

```
$ git clone https://github.com/SANBI-SA/goget2neo.git
$ cd goget2neo
```
### With `docker-compose:`
**Run the following:**

```
$ docker-compose up -d
```

***Point your browser at [localhost:7474](http://localhost:7474]) .***

### Standalone

**Pull and run the [neo4j docker image](https://hub.docker.com/_/neo4j/):**

```
$ docker run -d -p 7474:7474 -p 7687:7687 --name neo -e NEO4J_AUTH=none -v=$HOME/neo4j/data:/data neo4j:3.0.4
```

**Create a virtual environment:**

```
$ virtualenv envname
$ source envname/bin/activate
$ pip install -r requirements.txt
$ pip install --editable .
$ goget --help
$ goget init --help
$ goget init --d --r --u data/MTB_H37rv.gff3
```
***Point your browser at [localhost:7474](http://localhost:7474]) .***

We used a [GFF file from EnsemblBacteria](ftp://ftp.ensemblgenomes.org/pub/bacteria/release-30/gff3/bacteria_0_collection/mycobacterium_tuberculosis_h37rv).

