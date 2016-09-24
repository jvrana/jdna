import re
from hydradna.feature import Feature
from hydradna.sequence import Sequence

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

            name, span_start, span_end, span_length = Convert.parse_feature_name(a['name'])
            newf = Feature(name, a['type'], strand=a['strand'], color=a['color'])
            start_index = 0
            if span_start is not None:
                start_index = int(span_start)
            seq.add_feature(a['start'], fix_end(a['end']), newf, start_index=start_index)
            if span_length is not None:
                newf.length = int(span_length)
        return seq

    @staticmethod
    def to_benchling_json(seq):
        data = {}
        data['name'] = seq.name
        data['circular'] = seq.is_cyclic()
        data['annotations'] = []
        data['bases'] = str(seq)
        data['description'] = seq.description
        def add_annotation(f):
            data['annotations'].append(f)

        def fix_end(end):
            end += 1
            if end == len(seq):
                end = 0
            return end

        features = seq.get_features()
        for f in features:
            range, span = features[f]
            start, end = range
            span = '[{},{},{}]'.format(span[0], span[1], f.length)
            add_annotation(dict(
                name=str(f.name) + ' ' + str(span),
                type=f.type,
                start=start,
                end=fix_end(end),
                strand=f.strand,
                color=f.color
            ))
        return data

    @staticmethod
    def parse_feature_name(fname):
        g = re.search('(.+)(\[(\d+),(\d*),(\d*)\])|(.+)', fname)
        name, span, span_start, span_end, span_length, fullname = g.groups()
        if name is None:
            name = fullname
        return name, span_start, span_end, span_length
