#!/usr/bin/env python3

import numpy as np
import sys

fileA, fileB = sys.argv[2:]
A, B = [np.genfromtxt(f, delimiter=1, dtype='S1') for f in (fileA, fileB)]
