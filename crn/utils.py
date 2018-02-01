import matplotlib.pyplot as plt
import os
import pkgutil
import sys

from contextlib import contextmanager
from lib2to3.main import main as lib2to3_main
from types import ModuleType
from os.path import join, dirname

@contextmanager
def no_output():
    try:
        sys.stdout = open(os.devnull, "w")
        sys.stderr = open(os.devnull, "w")

        yield

        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__

    except:
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        raise


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


def stochpy_fix():
    with no_output():
        stochpy = pkgutil.get_loader("stochpy")

        if stochpy is None:
            import stochpy

        pysces_mini_model = join(dirname(stochpy.path),
            "modules", "PyscesMiniModel.py")

        lib2to3_main("lib2to3.fixes",
            f"-w -f has_key {pysces_mini_model}".split())

