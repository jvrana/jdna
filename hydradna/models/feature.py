from hydradna.utilities import Utilities

class Feature(object):

    def __init__(self, name, type='misc feature', strand=1, color=None):
        self.name = name
        self.type = type
        self.start = 0
        self.end = None
        self.strand = 1
        if color is None:
            color = Utilities.random_color()
        self.color = color

    def __str__(self):
        return '{}{}'.format(self.name, self._feature_span())

    def _feature_span(self):
        e = self.end
        if e is None:
            e = ''
            if self.start == 0:
                return ''
        return '[{}...{}]'.format(self.start, e)

    def __repr__(self):
        return str(self)
