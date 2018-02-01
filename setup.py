from setuptools import setup
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
long_description = open(path.join(here, 'README.md'), encoding='utf-8').read()
requirements = ['scipy', 'numpy', 'stochpy', 'sympy']

setup(
    name='crn',
    packages=['crn'],
    version='0.1.0a1',
    description='A Python3.6+ CRN Simulator',
    long_description=long_description,
    author='Enrico Borba',
    author_email='enricozb@gmail.com',
    url='https://github.com/enricozb/python-crn',
    download_url='https://github.com/enricozb/python-crn/archive/0.1.0a0.tar.gz',
    install_requires=requirements,
    python_requires='>=3.6',
    keywords=['crn', 'simulator', 'simulation'],
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'License :: OSI Approved :: MIT License',

        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ])

