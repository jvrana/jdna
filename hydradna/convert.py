from hydradna import *
import re
class Convert:

    @staticmethod
    def from_benchling(data):

        def clean_data(dic):
            for key in dic:
                if isinstance(dic[key], basestring):
                    dic[key] = str(dic[key])

        clean_data(data)
        # Sequence, name
        seq = Sequence(sequence=data['bases'], name=data['name'])

        # Topology
        if data['circular']:
            seq.make_cyclic()

        def fix_end(end):
            end = end - 1
            if end < 0:
                end = len(seq) - 1
            return end

        for a in data['annotations']:
            # Feature name, type, strand, color

            name, start, end = Convert.parse_feature_name(a['name'])
            newf = Feature(name, a['type'], strand=a['strand'], color=a['color'])
            newf.start = start
            newf.end = end

            start = a['start']
            end = fix_end(a['end'])
            seq.add_feature(start, end, newf)
        return seq

    @staticmethod
    def to_benchling_json(seq):
        data = {}
        data['name'] = seq.name
        data['circular'] = seq.is_cyclic()
        data['annotations'] = []
        data['bases'] = str(seq)
        def add_annotation(f):
            data['annotations'].append(f)

        def fix_end(end):
            end += 1
            if end == len(seq):
                end = 0
            return end

        features = seq.get_feature_ranges()
        for f in features:
            ranges = features[f]
            if not len(ranges) == 1:
                raise Exception("Feature cannot be in more than one place.")
            start, end = ranges[0]

            add_annotation(dict(
                name=str(f),
                type=f.type,
                start=start,
                end=fix_end(end),
                strand=f.strand,
                color=f.color
            ))
        return data

    @staticmethod
    def parse_feature_name(fname):
        g = re.search('(.+)(\[(\d+)\.\.\.(\d*)\])|(.+)', fname)
        name, span, start, end, fullname = g.groups()
        if name is None:
            start = 0
            name = fullname
        return name, start, end
