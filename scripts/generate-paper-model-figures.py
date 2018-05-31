#! /usr/bin/env python3

from __future__ import division
import numpy as np
import time, signal, sys, os, matplotlib, subprocess, re

if 'show' not in sys.argv:
    matplotlib.use('Agg')

import argparse
import datetime
import dynein.draw.cartoon as cartoon
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

poststroke = mpimg.imread('papers/paper/figures/schematic-poststroke.png')
prestroke = mpimg.imread('papers/paper/figures/schematic-prestroke.png')

plt.imshow(poststroke)
plt.axis('off')
plt.savefig("plots/poststroke-figure.pdf", bbox_inches='tight', format="pdf")

plt.imshow(prestroke)
plt.axis('off')
plt.savefig("plots/prestroke-figure.pdf", bbox_inches='tight', format="pdf")
