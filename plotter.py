# this plotter is used for plotting a series of analysis
# created by AvrinJoker 2016
# e-mail:avrinjoker@hotmail.com

import os

import matplotlib.pyplot as plt

import click

from methods_plotter import *

from overall_best import *

CONTEXT_SETTINGS = dict(help_option_names=['-h','--help'])

#methods_plotter('example_MAP4K4_RMSD.csv')

overall_best('example_methods_RMSD.csv','orange.png')
