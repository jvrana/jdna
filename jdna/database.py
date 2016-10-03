import mechanize
from bs4 import BeautifulSoup
import re
from aquariumapi import AquariumAPI
import shelve
from benchlingapi.convert import *
from benchlingapi import BenchlingAPI
import os
from jdna.core import Convert, Sequence, Reaction
import json
from copy import copy


credentials = None
with open('data/login.json', 'r') as f:
    credentials = json.load(f)

aqapi = AquariumAPI(credentials['aqapi_url'], credentials['aqapi_login'], credentials['aqapi_key'])
api = BenchlingAPI(credentials['benchling_api_key'])
DNAdb_location = credentials['db_location']

def login(url, login, password):
    br = mechanize.Browser()
    br.set_handle_robots(False)
    br.set_handle_refresh(False)
    response = br.open(url)
    br.form = list(br.forms())[0]
    control_login = br.form.find_control("session[login]")
    control_password = br.form.find_control("session[password]")
    control_login.value = login
    control_password.value = password
    response = br.submit()
    response = response.read()

    response = br.open(url).read()
    br.url = url
    return br, response


def goto(url, br):
    text = br.open(os.path.join(br.url, url)).read()
    return text


def get_sample_info(sample_id):
    br, response = login(credentials['aq_url'], credentials['aqapi_login'], credentials['aq_password'])
    txt = goto('samples/{}'.format(sample_id), br)
    soup = BeautifulSoup(txt)

    # Get field information
    items = soup.find('ul', {'class': 'list'}).findAll('li')
    field_info = {}

    for i in items:
        txt = i.text.replace('\n', '')
        field = i.b.text
        g = re.search('{}\s*(.+)'.format(field), txt)
        if g == None:
            field_info[field] = None
        else:
            field_info[field] = g.group(1).strip()
    return field_info


def get_sample(sample):
    s = None
    if isinstance(sample, basestring):
        s = aqapi.find('sample', {'name': sample})['rows'][0]
    elif isinstance(sample, int):
        s = aqapi.find('sample', {'id': sample})['rows'][0]
    elif isinstance(sample, dict):
        s = sample
    else:
        print 'Sample {} must be an int or string'.format(sample)
    s['fields'] = get_sample_info(s['id'])
    return s



########################


def db():
    return shelve.open(os.path.join(DNAdb_location, 'DNA.db'), writeback=False)


def db_get(key):
    dnadb = db()
    value = dnadb[key]
    dnadb.close()
    return value


def _save_to_fsa(item):
    item['id'] = str(item['id'])
    seqrec = benchling_to_seqrecord(item)
    filename = os.path.join(DNAdb_location, '{id}_{name}.fsa'.format(id=item['id'], name=item['name']))
    with open(filename, 'w') as handle:
        SeqIO.write(seqrec, handle, 'fasta')


def retrieve(key):
    dnadb = db()
    value = None
    try:
        value = dnadb[str(key)]
    except KeyError:
        pass
    finally:
        dnadb.close()
    return value


def register(item, savefsa=True):
    item = copy(item)
    dnadb = db()
    dnadb[str(item['id'])] = item
    dnadb[item['name']] = item['id']
    dnadb.close()
    if savefsa:
        _save_to_fsa(item)
    print 'Item {} {} registered.'.format(item['id'], item['name'])


def register_jdna(jseq, savefsa=True):
    bseq = Convert.to_benchling_json(jseq)
    bseq['id'] = jseq.id
    bseq['aliases'] = jseq.aliases
    register(bseq, savefsa=savefsa)


def get_seq_from_db(id):
    dnadb = db()
    return Convert.from_benchling(dnadb[str(id)])


###########

def construct_sequence(sample_name):
    print 'getting sequence for', sample_name
    s = get_sample(sample_name)
    sharelink = s['fields']['Sequence']
    bseq = api.getsequencefromsharelink(sharelink)
    bseq['name'] = s['name']
    bseq['id'] = s['id']
    bseq['aliases'] = [s['name'], s['id']]
    return bseq


def get_json_primer(primer_name):
    pid = retrieve(primer_name)
    p = retrieve(pid)
    if p is None:
        p = update_primer(primer_name)
    return p


def get_json_sequence(sample_name):
    sample_name = str(sample_name)
    bseq_id = retrieve(sample_name)
    bseq = retrieve(bseq_id)
    if bseq is None:
        bseq = update_sequence(sample_name)
    return bseq


def get_jdna(sample_name):
    sample_name = str(sample_name)
    bseq = get_json_sequence(sample_name)
    hseq = Convert.from_benchling(bseq)
    hseq.aliases = bseq['aliases']
    hseq.id = bseq['id']
    return hseq


def get_jdna_primer(primer_name):
    p = get_json_primer(primer_name)
    anneal = p['fields']['Anneal Sequence']
    overhang = p['fields']['Overhang Sequence']
    if overhang is None or 'one' in overhang:
        overhang = ''
    seq_str = str(overhang).lower().strip() + str(anneal).lower().strip()
    seq = Sequence(sequence=seq_str)
    seq.name = 'Primer {}_{}'.format(p['id'], p['name'])
    seq.description = p['description']
    # seq.create_feature(seq.name, 'primer', 0, len(seq)-1)
    return seq


def update_primer(primer_name):
    p = get_sample(primer_name)
    register(p, savefsa=False)
    return p


def update_sequence(sample_name):
    bseq = construct_sequence(sample_name)
    register(bseq)
    return bseq


def update_fragment_sequence(name_or_id):
    aqfragment = get_sample(name_or_id)
    fields = aqfragment['fields']
    template = get_jdna(fields['Template'])
    p1 = get_jdna_primer(fields['Forward Primer'])
    p2 = get_jdna_primer(fields['Reverse Primer'])
    products = Reaction.pcr(template, p1, p2)
    if len(products) > 1:
        raise Exception("More than one product found.")
    product = products[0]
    product.name = aqfragment['name']
    product.description = str(aqfragment['description'])
    product.id = aqfragment['id']
    product.aliases = [str(aqfragment['id']),
                       str(aqfragment['name']),
                       'Fragment {} {}'.format(aqfragment['id'], aqfragment['name'])]
    register_jdna(product)
    return product


def get_fragment_sequence(fragment_id):
    bfrag = retrieve(fragment_id)
    if bfrag is None:
        bfrag = update_fragment_sequence(fragment_id)
    else:
        bfrag = Convert.from_benchling(bfrag)
    return bfrag


def assembly_gibson(task_id):
    task = aqapi.find('task', {'id': task_id})['rows'][0]
    specs = task['specification']
    specs = json.loads(specs)
    fragments = specs['fragments Fragment']
    plasmid = aqapi.find('sample', {'id': specs['plasmid Plasmid']})['rows'][0]
    print 'Gibson Assembly: {} {}'.format(task_id, task['name'])
    print 'Plasmid: {} {}'.format(plasmid['id'], plasmid['name'])
    print 'Fragments: {}'.format(fragments)
    print 'OUT >',
    hfrags = [get_fragment_sequence(f) for f in fragments]
    products = Reaction.cyclic_assembly(hfrags, max_homology=61)
    if len(products) > 1:
        raise Exception
    result = products[0]
    result.name = str(plasmid['name'])
    result.description = str(plasmid['description'])
    result.id = plasmid['id']
    result.aliases = [str(result.id), str(result.name)]
    register_jdna(result)
    return result


def to_benchling(id, folder_name, overwrite=True):
    bseq = copy(retrieve(id))
    del bseq['id']
    bseq['folder'] = api.find_folder(folder_name, regex=True)['id']
    bseq['overwrite'] = overwrite
    seq = api.create_sequence(**bseq)
    print str('www.benchling.com' + seq['editURL'])


def gibson_workflow(task_id, benchling_folder_name='jdna', overwrite=True):
    result = assembly_gibson(task_id)
    to_benchling(result.id, benchling_folder_name, overwrite=overwrite)