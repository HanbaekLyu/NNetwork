# NNetwork

[![PyPI Version](https://img.shields.io/pypi/v/NNetwork.svg)](https://pypi.org/project/NNetwork/)
[![Supported Python Versions](https://img.shields.io/pypi/pyversions/NNetwork.svg)](https://pypi.org/project/NNetwork/)

Custom graph/network/multi-weighted network class based on storing list of neighbors for each nodes (as opposed to edge list) for scalable sampling and searching algorithms

---

## Installation

To install NNetwork, run this command in your terminal:

```bash
$ pip install -U NNetwork
```

This is the preferred method to install NNetwork, as it will always install the most recent stable release.

If you don't have [pip](https://pip.pypa.io) installed, these [installation instructions](http://docs.python-guide.org/en/latest/starting/installation/) can guide
you through the process.

## Quick Start
```python
>>> from NNetwork import Example
>>> a = Example()
>>> a.get_value()
10

```

## Citing
If you use our work in an academic setting, please cite our paper:



## Development
See [CONTRIBUTING.md](CONTRIBUTING.md) for information related to developing the code.

#### Suggested Git Branch Strategy
1. `master` is for the most up-to-date development, very rarely should you directly commit to this branch. Your day-to-day work should exist on branches separate from `master`. It is recommended to commit to development branches and make pull requests to master.4. It is recommended to use "Squash and Merge" commits when committing PR's. It makes each set of changes to `master`
atomic and as a side effect naturally encourages small well defined PR's.


#### Additional Optional Setup Steps:
* Create an initial release to test.PyPI and PyPI.
    * Follow [This PyPA tutorial](https://packaging.python.org/tutorials/packaging-projects/#generating-distribution-archives), starting from the "Generating distribution archives" section.

* Create a blank github repository (without a README or .gitignore) and push the code to it.

* Delete these setup instructions from `README.md` when you are finished with them.
