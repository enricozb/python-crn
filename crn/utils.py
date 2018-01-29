import os
import sys

from contextlib import contextmanager

@contextmanager
def quiet():
    sys.stdout = sys.stderr = open(os.devnull, "w")

    yield

    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__

def datadir(filename=""):
    """
    Returns the absolute path of python-crn's data directory. Creates the
    folder if it does not exist.
    """
    datadir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

    if not os.path.isdir(datadir):
        os.makedirs(datadir)

    if filename:
        return os.path.join(datadir, filename)
    return datadir

