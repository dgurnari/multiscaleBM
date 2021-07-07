# -*- coding: utf-8 -*-
"""
This is a module to be used as a reference for building other modules
"""

import matplotlib.pyplot as plt

from ._ballmapperinterfaces import BallMapperCppInterfacePython



class BallMapper():
    """BallMapper class
    """
    def __init__(self, points, values, epsilon):
        bm_output = BallMapperCppInterfacePython(points,
                                                 values,
                                                 epsilon)

        self.edges = bm_output[1]
