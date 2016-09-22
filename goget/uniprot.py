from bioservices import UniProt

u = UniProt(verbose=False)

columns = "id, genes, go-id, interpro, interactor, genes(PREFERRED), feature(DOMAIN EXTENT), protein names, go, " \
          "citation, 3d, comment(FUNCTION), sequence, mass, length, families, go(biological process), " \
          "go(molecular function), go(cellular component)"


def query_uniprot(locus_tags, proteome='UP000001584'):
    result = None
    for tag_list in locus_tags:
        gene_query = '(' + '+OR+'.join(['gene:' + name for name in tag_list]) + ')'
        query = "proteome:{}+AND+{}".format(proteome, gene_query)

        result = u.search(query, frmt="tab", columns=columns, sort=None)
        print(result)

    return result
