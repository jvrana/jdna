from hydradna.utilities import Utilities

class Feature(object):

    def __init__(self, name, type='misc feature', strand=1, color=None):
        self.name = name
        self.type = type
        self.strand = 1
        self._length = None
        if color is None:
            color = Utilities.random_color()
        self.color = color

    def __str__(self):
        return '{} {}'.format(self.name, self.type)

    def __repr__(self):
        return str(self)
