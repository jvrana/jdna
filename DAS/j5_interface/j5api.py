'''
Project: jdna
File: j5_api
Author: Justin
Date: 2/9/17

Description: 

'''

import xmlrpclib
from zipfile import ZipFile
import json
import os
from glob import glob
import base64

directory = '/Users/Justin/Desktop/NewDesign/'

credentials = None
with open('j5_credentials.json') as handle:
    credentials = json.load(handle)


files = ['zipped_sequences', 'master_oligos', 'master_plasmids', 'master_direct_syntheses', 'parts_list',
         'j5_parameters', 'sequences_list', 'eugene_rules', 'target_part_order_list']

class J5Assembly(object):

    def __init__(self):
        self.zipped_sequences = None
        self.master_oligos = None
        self.master_plasmids = None
        self.master_direct_syntheses = None
        self.parts_list = None
        self.j5_parameters = None
        self.sequences_list = None
        self.eugene_rules = None
        self.target_part_order_list = None
        self.file_types = self.__dict__.keys()

    def encode(self, filepath):
        with open(filepath, 'rU') as handle:
            return base64.b64encode(handle.read())

    def assembly_params_from_dir(self, directory, session_id, assembly_method):

        def imply_file(x):
            g = glob(os.path.join(directory, x))
            assert len(g) == 1
            return g[0]

        params = {}
        for filetype in self.file_types:
            params['encoded_{}_file'.format(filetype)] = self.encode(imply_file('*{}*'.format(filetype)))
            params['reuse_{}_file'.format(filetype)] = False
        params['j5_session_id'] = session_id
        params['assembly_method'] = assembly_method
        self.params = params
        return params





class J5API(object):

    def __init__(self, username, password):
        print username, password
        self.proxy_home = 'https://j5.jbei.org/bin/j5_xml_rpc.pl'
        self.proxy = xmlrpclib.ServerProxy('https://j5.jbei.org/bin/j5_xml_rpc.pl')
        self.session = self.proxy.CreateNewSessionId({'username': username, 'password': password})
        self.session_id = self.session['j5_session_id']

    def open_design_file(self, design_path):
        self.assembly = J5Assembly()
        self.assembly.assembly_params_from_dir(design_path, self.session_id, 'SLIC/Gibson/CPEC')


    def run_assembly(self, design_path):
        self.open_design_file(design_path)
        return self.proxy.DesignAssembly(self.assembly.params)

#
#
# files = ['zipped_sequences', 'master_oligos', 'master_plasmids', 'master_direct_syntheses', 'parts_list',
#          'j5_parameters', 'sequences_list', 'eugene_rules', 'target_part_order_list']
#
# d = {}
# for filetype in files:
#     d['encoded_{}_file'.format(filetype)] = encode(f('*{}*'.format(filetype)))
#     d['reuse_{}_file'.format(filetype)] = False
# d['j5_session_id'] = session['j5_session_id']
# d['assembly_method'] = 'SLIC/Gibson/CPEC'
#
# out = proxy.DesignAssembly(d)
# with open('/Users/Justin/Desktop/error.zip', 'w') as handle:
#     handle.write(base64.decodestring(out['encoded_output_file']))
# print out
#
#
#
# with ZipFile(os.path.join(directory, 'zzzzipped_sequences.zip'), mode='w') as myzip:
#     for filename in glob(os.path.join('*')):
#         myzip.write(filename)