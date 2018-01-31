import os
import sys
import matplotlib.pyplot as plt

from contextlib import contextmanager

@contextmanager
def no_output():
    sys.stdout = sys.stderr = open(os.devnull, "w")

    yield

    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__

@contextmanager
def plot_without_xserver():
    backend = plt.get_backend()
    plt.switch_backend("Svg")

    yield

    plt.switch_backend(backend)

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

