import random

class Feature(object):

    def __init__(self, name, type='misc feature', strand=1, color=None):
        self.name = name
        self.type = type
        self.strand = 1
        self._length = None
        if color is None:
            color = random_color()
        self.color = color

    def __str__(self):
        return '{} {}'.format(self.name, self.type)

    def __repr__(self):
        return str(self)

@staticmethod
def rgb_to_hex(r, g, b):
    def clamp(x):
        return max(0, min(x, 255))

    return "#{0:02x}{1:02x}{2:02x}".format(clamp(r), clamp(g), clamp(b))

@staticmethod
def random_color():
    rgb = [int(random.random()*255) for x in range(3)]
    return rgb_to_hex(*rgb)
