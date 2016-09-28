import StringIO
import csv
from time import time

from bioservices import UniProt
from dbconn import create_uniprot_nodes
from tqdm import tqdm

u = UniProt(verbose=False)


def search_uniprot(query, columns, proteome='UP000001584'):
    """
    Search UniProt and return results as list
    :param query:
    :param columns:
    :param proteome:
    :return:
    """
    query = "proteome:{}+AND+{}".format(proteome, query)

    result = u.search(query=query, frmt="tab", columns=columns, sort=None)
    reader = csv.reader(StringIO.StringIO(result), delimiter='\t')
    try:
        reader.next()
    except StopIteration:
        return []
    else:
        return list(reader)


def query_uniprot(locus_tags):
    """
    Get data from UniProt
    :param locus_tags:
    :return:
    """
    start = time()
    columns = "id, genes, go-id, interpro, interactor, genes(PREFERRED), feature(DOMAIN EXTENT), protein names, go, " \
              "citation, 3d, comment(FUNCTION), sequence, mass, length, families, go(biological process), " \
              "go(molecular function), go(cellular component)"
    ambiguous_names = []
    uniprot_data = []
    for tag_list in tqdm(locus_tags):
        query = '(' + '+OR+'.join(['gene:' + name for name in tag_list]) + ')'
        result = search_uniprot(query, columns)
        locus_tag_set = set(tag_list)
        new_rows = []
        for row in result:
            names = row[1].replace('/', ' ').split()
            found_count = 0
            found_name = ''
            for name in names:
                if name.lower() in locus_tag_set:
                    found_name = name
                    found_count += 1
            if found_count == 0:
                for name in names:
                    if '.' in name:
                        name = name[:name.find('.')]
                        if name.lower() in locus_tag_set:
                            found_name = name
                            found_count += 1
            if found_name:
                assert found_count >= 1, "Problem with ({}); count ({});\n row ({})".format(found_name, found_count,
                                                                                            row)

            if found_count > 1:
                ambiguous_names.extend([name for name in names if name.startswith('Rv')])
            else:
                new_row = save_row(found_name, row)
                new_rows.append(new_row)
                uniprot_data.extend(new_rows)
        # print("ambiguous_names: ", len(ambiguous_names), ambiguous_names)
        for name in ambiguous_names:
            rows = search_uniprot(name, columns)
            assert len(rows) == 1, "Got multiple results for {}.".format(name)
            new_row = save_row(name, rows[0])
            uniprot_data.append(new_row)
    end = time()
    print("Done fetching data from UniProt in ", end - start, "secs.")
    create_uniprot_nodes(uniprot_data)
    return uniprot_data


def save_row(found_name, row):
    return [found_name, row[0], row[1], row[2].split('; '), row[3].split('; '), row[4], row[5], row[6].split('; '),
            row[7], row[8].split('; '), row[9].split('; '), row[10].split('; '), row[11], row[12], row[13], row[14],
            row[15], row[16], row[17], row[18]]
