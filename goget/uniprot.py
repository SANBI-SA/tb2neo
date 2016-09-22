from bioservices import UniProt

u = UniProt(verbose=False)


def search_uniprot(query, columns, proteome='UP000001584'):
    query = "proteome:{}+AND+{}".format(proteome, query)

    result = u.search(query=query, frmt="tab", columns=columns, sort=None)
    return result


def query_uniprot(locus_tags):
    """
    Get data from UniProt
    :param locus_tags:
    :return:
    """
    columns = "id, genes, go-id, interpro, interactor, genes(PREFERRED), feature(DOMAIN EXTENT), protein names, go, " \
              "citation, 3d, comment(FUNCTION), sequence, mass, length, families, go(biological process), " \
              "go(molecular function), go(cellular component)"
    result = None
    for tag_list in locus_tags:
        query = '(' + '+OR+'.join(['gene:' + name for name in tag_list]) + ')'
        result = search_uniprot(query, columns)
        print(result)

    return result
