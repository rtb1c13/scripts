#!/usr/bin/env python

# Script to perform fit of theoretical to experimental SAXS intensities while
# minimising Chi-squared. Constant offset (crysol -cst) currently not performed...

import numpy as np
from scipy.optimize import minimize

def chiSquare(
