"""
misc utilities
"""

import random


def rgb_to_hex(r, g, b):
    def clamp(x):
        return max(0, min(x, 255))

    return "#{0:02x}{1:02x}{2:02x}".format(clamp(r), clamp(g), clamp(b))


def random_color():
    rgb = [int(random.random() * 255) for x in range(3)]
    return rgb_to_hex(*rgb)