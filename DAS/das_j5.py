'''
Project: jdna
File: das_j5
Author: Justin
Date: 2/21/17

Description: 

'''
from das_assembly import *

import xmlrpc.client

class J5Part(Region):

    def __init__(self, name, source, five_internal, three_internal, region):
        super(J5Part, self).__init__(
            region.start,
            region.end,
            region.length,
            region.circular,
            direction=region.direction,
            start_index=region.bounds_start)
        self.name = name
        self.source = source
        self.rc = self.direction != Region.FORWARD
        self.five_internal = five_internal
        self.three_internal = three_internal



class J5Assembly(Assembly):
    PARTLABELS = ['Part Name', 'Part Source (Sequence Display ID)', 'Reverse Compliment?', 'Start (bp)', 'End (bp)',
                  'Five Prime Internal Preferred Overhangs?', 'Three Prime Internal Preferred Overhangs?']
    SEQLISTLABELS = ['Sequence File Name', 'Format']
    TARGETLABELS = ["(>Bin) or Part Name", "Direction", "Forced Assembly Strategy?",
                    "Forced Relative Overhang Position?", "Direct Synthesis Firewall?", "Extra 5' CPEC overlap bps",
                    "Extra 3' CPEC overlap bps"]
    MASTEROLIGOSLABELS = ['Oligo Name', 'Length', "Tm", "Tm (3' only)", "Sequence"]
    MASTERPLASMIDLABELS = ["Plasmid Name", "Alias", "Contents", "Length", "Sequence"]
    MASTERDIRECTLABELS = ["Direct Synthesis Name", "Alias", "Contents", "length", "Sequence"]

    # ZIP = None
    # EUGENE = None
    # PARAMS = None

    def __init__(self, assembly):
        super(J5Assembly, self).__init__(assembly.contigs, assembly.contig_container, assembly.primer_container)
        self.params = None
        self.proxy_home = 'https://j5.jbei.org/bin/j5_xml_rpc.pl'
        self.proxy = xmlrpc.client.ServerProxy('https://j5.jbei.org/bin/j5_xml_rpc.pl', verbose=True)
        self.session = None
        self.session_id = None

    def to_csv(self, labels, rows):
        csvrows = []
        for row in [labels] + rows:
            csvrows.append(','.join([str(x) for x in row]))
        return '\n'.join(csvrows)

    @staticmethod
    def encode64(filestring):
        return base64.encodestring(filestring)

    @staticmethod
    def decode64(string):
        return base64.decodestring(string)

    def all(self):
        self.process_assembly()
        self.get_target()
        self.get_parts()
        self.get_sequences()
        self.get_plasmids()
        self.get_primers()
        self.get_direct()
        self.get_eugene()
        self.get_parameters()

    def decode_all_to(self, destination):
        self.all()
        names = [
            'parts_csv',
            'target',
            'sequences',
            'plasmids',
            'direct',
            'eugene',
            'parameters',
        ]
        print self.__dict__.keys()
        for n in names:
            with open(os.path.join(destination, n + '.csv'), 'w') as handle:
                handle.write(J5Assembly.decode64(self.__dict__[n]))
        with open(os.path.join(destination, 'sequences.zip'), 'w') as zip:
            zip.write(J5Assembly.decode64(self.encoded_zip))

    def get_parameters(self):
        with open('j5_parameters.csv') as params:
            self.parameters = J5Assembly.encode64(params.read())

    def process_assembly(self):

        # Process Parts
        parts = []
        sequences = []

        pairs = zip(self.contigs[0:-1], self.contigs[1:])
        pairs.append((self.contigs[-1], self.contigs[0]))
        for l, r in pairs:
            _l = Region(l.query.start, l.query.end, l.query.length / 2.0, True,
                        direction=l.query.direction,
                        name=l.subject.name,
                        start_index=pairs[0][0].query.start)
            _r = Region(r.query.start, r.query.end, r.query.length / 2.0, True,
                        direction=l.query.direction,
                        name=r.subject.name,
                        start_index=pairs[0][0].query.start)
            part = J5Part(l.contig_id, l.filename,
                          '', '', _l)
            parts.append(part)
            d = _l.get_gap_degree(_r)
            if d > 0:
                gap_part = J5Part('{}_{}_gap'.format(_l.name, _r.name),
                                  self.meta.query_filename,
                                  '',
                                  '',
                                  _l.get_gap(_r))
                parts.append(gap_part)
        self.parts = parts

        # Process sequences


    def get_target(self):
        rows = []
        for p in self.parts:
            row = [
                p.name,
                'forward',
                '',
                '',
                '',
                '',
                ''
            ]
            rows.append(row)
        self.target = J5Assembly.encode64(self.to_csv(J5Assembly.TARGETLABELS, rows))

    def get_parts(self):
        rows = []
        for p in self.parts:
            row = [
                p.name,
                os.path.basename(p.source),
                str(p.direction != Region.FORWARD),
                p.start,
                p.end,
                '',
                ''
            ]
            rows.append(row)
        csv = self.to_csv(J5Assembly.PARTLABELS, rows)
        self.parts_csv = J5Assembly.encode64(csv)

    def get_direct(self):
        self.direct = J5Assembly.encode64(self.to_csv(J5Assembly.MASTERDIRECTLABELS, []))

    def get_plasmids(self):
        self.plasmids = J5Assembly.encode64(self.to_csv(J5Assembly.MASTERPLASMIDLABELS, []))

    def get_primers(self):
        # ['Oligo Name', 'Length', "Tm", "Tm (3' only)", "Sequence"]
        rows = []
        for p in self.primer_container.contigs:
            row = [
                p.subject.name,
                len(p.subject.seq),
                60,
                60,
                p.subject.seq
            ]
            rows.append(row)
        self.primers = J5Assembly.encode64(self.to_csv(J5Assembly.MASTEROLIGOSLABELS, rows))

    # TODO: add direct fragments to directmasterlist
    # TODO: if is direct, do something with j5 assembly parameters?
    def get_sequences(self):
        print 'SEQUENCE LIST'
        sequences = []
        for p in self.parts:
            print '\tSEQ: {}  @  {}'.format(os.path.basename(p.source), p.name)
            sequences.append((p.source, p.name,))
        sequences = list(set(sequences))

        # Save sequence_list

        csvrows = []
        for s in sequences:
            filename = os.path.basename(s[0])
            format = determine_format(s[0]).capitalize()
            csvrows.append((filename, format))
        csv = self.to_csv(J5Assembly.SEQLISTLABELS, csvrows)
        self.sequences = J5Assembly.encode64(csv)

        # TODO: move zip stuff to seqio
        # Save zipped sequence file
        temp = tempfile.TemporaryFile()  # open('/Users/Justin/Desktop/temporary.txt', 'w')#tempfile.TemporaryFile()
        with zipfile.ZipFile(temp, 'w') as zf:

            for s in sequences:
                filename = s[0]
                print '\tZipping {} as {}'.format(s[0], os.path.basename(filename))
                with open(filename, 'rU') as handle:
                    zf.writestr(os.path.join('ZippedTemplates', os.path.basename(filename)), handle.read())
        temp.seek(0)
        self.encoded_zip = J5Assembly.encode64(temp.read())
        temp.close()

    def get_eugene(self):
        self.eugene = J5Assembly.encode64('')

    def to_dict(self):
        d = {
            'zipped_sequences': self.encoded_zip,
            'master_oligos': self.primers,
            'master_plasmids': self.plasmids,
            'master_direct_syntheses': self.direct,
            'parts_list': self.parts_csv,
            'j5_parameters': self.parameters,
            'sequences_list': self.sequences,
            'eugene_rules': self.eugene,
            'target_part_order_list': self.target,
        }
        params = {}
        for k in d:
            params['encoded_{}_file'.format(k)] = d[k]
            params['reuse_{}_file'.format(k)] = d[k]
        return params

    def login(self, username, password):
        self.session = self.proxy.CreateNewSessionId({'username': username, 'password': password})
        self.session_id = self.session['j5_session_id']
        return self.session_id

    def submit(self, username, password):
        self.all()
        self.login = self.login(username, password)
        return self.run_assembly()

    def run_assembly(self, assembly_method='SLIC/Gibson/CPEC'):
        d = self.to_dict()
        print d
        d['j5_session_id'] = self.session_id
        d['assembly_method'] = assembly_method
        print d
        print self.to_dict().keys()
        return self.proxy.DesignAssembly(d)

    ## Not implemented, this remove the possibility of finding a cheaper or more successful solution
    # def remove_expensive_fragments(self):
    #     contigs_for_removal = []
    #     for c1 in self.contigs:
    #         for c2 in self.contigs:
    #             if c1 == c2:
    #                 continue
    #             if c1 in contigs_for_removal or c2 in contigs_for_removal:
    #                 continue
    #
    #             if c1.equivalent_location(c2):
    #                 c1_cost = self.get_pcr_cost(c1)
    #                 c2_cost = self.get_pcr_cost(c2)
    #                 if c1_cost > c2_cost



        # params['encoded_{}_file'.format(filetype)] = self.encode(imply_file('*{}*'.format(filetype)))
        #             params['reuse_{}_file'.format(filetype)] = False
