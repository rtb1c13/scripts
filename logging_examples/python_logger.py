#!/usr/bin/env python

# Example use of logging module

import logging
import random

# Basic setup, logger name is current module name:
# Choose level of logging here too - maybe use an argument!
# 5 levels: CRITICAL, ERROR, WARNING, INFO, DEBUG
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG) # Will log DEBUG or higher...

# Set up outputs to file
handler = logging.FileHandler('python_logs.log')
# Will write out DEBUG or higher...
handler.setLevel(logging.DEBUG)
handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
logger.addHandler(handler)

# Do some work!

def multiply(a, b, c):
    logger.debug("Input arguments to multiply, step %s were: %s %s %s" % (i, a, b, c))
    return a*b*c

logger.info("Now doing some multiplication")
for i in range(1,11):
    x, y, z = random.sample(xrange(100), 3)
    result = multiply(x, y, z)
    logger.info("The answer for step %s was %s" % (i, result))
logger.info("Done!")

# Now do some work that will fail!
try:
    x, y, z = random.sample(xrange(100), 4)
    result = multiply(x, y, z)
except Exception, e:
    logger.error("An exception was raised:", exc_info=True)
