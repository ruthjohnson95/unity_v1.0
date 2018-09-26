import logging
import numpy as np
import scipy.stats as st
import math
import scipy
from scipy.special import logsumexp
from optparse import OptionParser
import random
import sys
from scipy.special import logit

# global logging
logging.basicConfig(format='%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S', level=logging.INFO)


# global priors
beta_lam = 0.20 # prior for proportion p, p ~ Beta(beta_lam, beta_lam)

# global constants used to check for under/overflow
LOG_MIN = 1.7976931348623157e-308
LOG_MAX = 1.7976931348623157e+308
EXP_MAX = math.log(sys.float_info.max)
MAX = sys.float_info.max
MIN = sys.float_info.min

BURN = .25

# performs logsumexp calculation over a vector of values
def logsumexp_vector(a, axis=0):
    if axis is None:
        return logsumexp(a)
    a = np.asarray(a)
    shp = list(a.shape)
    shp[axis] = 1
    a_max = a.max(axis=axis)
    s = np.log(np.exp(a - a_max.reshape(shp)).sum(axis=axis))
    lse  = a_max + s
    return lse


# prints both to console and to outfile with file descriptor f
def print_func(line, f):
    print(line)
    sys.stdout.flush()
    f.write(line)
    f.write('\n')
    return
