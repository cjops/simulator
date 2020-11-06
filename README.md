# Exploring dynamic fitness landscapes

This repository is a bit chaotic, as the project has become a multi-year effort
and my programming/organizational skills have evolved along with it. I used
Python and Jupyter notebooks to do analysis and plotting. There is an
implementation of the simulator in C++ that must be built into a Python module
in order for the notebooks to run.

## Dependencies

The main external dependencies are [NumPy](https://numpy.org/install/),
[pandas](https://pandas.pydata.org/docs/getting_started/index.html#installation),
[Plotly](https://plotly.com/python/getting-started/#installation), and
[Jupyter](https://jupyter.org/install.html). My environment contains
`plotly==3.3.0` and Python 3.6.4. Newer versions of Plotly may break things.

## Getting started

1.  Make sure NumPy is installed in your environment.
2.  Inside `simulator-cpp`, run `python setup.py build`.
3.  Run `python setup.py install`, or move the resulting `.pyd` or `.so` from
    the `build/lib.system` directory into `notebooks-py`.
4.  Use Jupyter Notebook or JupyterLab to view the notebooks in `notebooks-py`.

## References

 -  Ogbunugafor, C., Eppstein, M. Competition along trajectories governs
    adaptation rates towards antimicrobial resistance. Nat Ecol Evol 1, 0007
    (2017). https://doi.org/10.1038/s41559-016-0007
 -  Mira, P. M. et al. Rational Design of Antibiotic Treatment Plans: A
    Treatment Strategy for Managing Evolution and Reversing Resistance. PLOS ONE
    10, e0122283 (2015). https://doi.org/10.1371/journal.pone.0122283
