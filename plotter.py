# this plotter is used for plotting a series of analysis
# created by AvrinJoker 2016
# e-mail:avrinjoker@hotmail.com

import os

import matplotlib.pyplot as plt

import click

import methods_plotter

CONTEXT_SETTINGS = dict(help_option_names=['-h','--help'])

methods_plotter.methods_plotter('example.csv')

