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
        (self.numer_of_covered_points , self.edges ,
        self.strength_of_edges , self.points_covered_by_landmarks ,
        self.landmarks , self.coloring ,
        self.coverage ) = BallMapperCppInterfacePython(points,
                                                       values,
                                                       epsilon)
