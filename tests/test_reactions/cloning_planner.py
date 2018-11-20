from jdna import database
from jdna.core import Reaction

database.credentials['db_location'] = '/Volumes/120GB SSD/Google Drive/Klavinslab/DNAdb'

temp_database = {}

def capture(name_or_id):
    # try to convert to string
    name = None
    if isinstance(name_or_id, str):
        name = name_or_id
    else:
        name = database.get_sample(name_or_id)['name']

    # try temporary database first
    if name in temp_database:
        return temp_database[name]
    # else get from local database
    else:
        try:
            return database.get_jdna(name)
        except:
            pass

        try:
            return database.get_jdna_primer(name)
        except:
            pass

    # sample = database.get_sample(name_or_id)


def temp_fragment(template_name, p1_name, p2_name):
    template = capture(template_name, type='dna')
    p1 = capture(p1_name, type='primer')
    p2 = capture(p2_name, type='primer')
    print((Reaction.anneal_primer(template, p1)))
    print((Reaction.anneal_primer(template, p2)))
    print((template.is_cyclic()))
    products = Reaction.pcr(template, p1, p2)
    print(('{} products'.format(len(products))))
    if len(products) > 1:
        raise Exception('more than one product')
    #temp_seq(products[0])
    print((len(products[0])))
    return products[0]