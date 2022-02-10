<p align="center">
<img width="700" src="https://github.com/HanbaekLyu/NNetwork/blob/master/nnetwork_logo.png?raw=true" alt="logo">
</p>

# NNetwork

[![PyPI Version](https://img.shields.io/pypi/v/NNetwork.svg)](https://pypi.org/project/NNetwork/)
[![Supported Python Versions](https://img.shields.io/pypi/pyversions/NNetwork.svg)](https://pypi.org/project/NNetwork/)

`NNetwork` is a Custom graph/network/multi-weighted network class optimized for scalable subgraph sampling and searching algorithms. NNetwork stores a dictionary that maps each node to a list of its neighbors to allow for O(1) access for finding neighbors. 

The efficiency of neighbor access is import for sampling algorithm such as random walks and Markov chain Monte Carlo motif sampling on graphs, which rely on accessing neighborhood information at every iteration of sampling. In comparison, many packages rely on calculations involving powers of adjacency matrices to calculate random walks of length k. 

The default class of `NNetwork` encodes a network with weighted edges, which can also have list-valued edge weights as its 'color'. 

Update for 0.1.0: 

Built-in functions contain sampling algorithms for mesoscale network patches using various MCMC motif sampling algorithms [1]. At stationary distribution, it computes a uniformly chosen k-walk in the graph, which can optionally enforced to be non-backtraking, and the induced adjacency pattern is returned as a k x k matrix. Algorithimically, a given k-walk is randomly updated using a suitable MCMC algorithm. The so-computed k x k `mesoscale patches` are basis of subgraph analysis and network dictionary learning in [2].  


### References:
[1] Hanbaek Lyu, Facundo Memoli, and David Sivakoff,
“Sampling random graph homomorphisms and applications to network data analysis.” (2019) https://arxiv.org/abs/1910.09483
[2] Hanbaek Lyu, Yacoub Kureh, Joshua Vendrow*, and Mason A. Porter,
“Learning low-rank mesoscale structures of networks” https://arxiv.org/abs/2102.06984 (2021)


By Josh Vendrow and Hanbaek Lyu

---

## Installation

To install NNetwork, run this command in your terminal:

```bash
$ pip install -U NNetwork
```

This is the preferred method to install NNetwork, as it will always install the most recent stable release.

If you don't have [pip](https://pip.pypa.io) installed, these [installation instructions](http://docs.python-guide.org/en/latest/starting/installation/) can guide
you through the process.

## Usage

**Undirected Graphs**

Create a simple undirected graph from an edgelist:
```python
>>> from NNetwork import NNetwork
>>> edgelist = [[1,2],[2,3],[3,4]]
>>> Network = NNetwork()
>>> Network.add_edges(edgelist)
>>> Network.has_edge(2,3)
True
```
Get the neighbors of a node:
```python
>>> Network.neighbors(3)
[2,4]
```

Find the intersection of edges with another network:
```python
>>> edgelist2 = [[2,3],[3,4],[5,7]]
>>> Network2 = NNetwork()
>>> Network2.add_edges(edgelist2)
>>> Network.intersection(Network2)
[[2,3],[3,4]]
```

**Weighted Graphs**

Create a weighted graph from an edgelist:
```python
>>> from NNetwork import Wtd_NNetwork
>>> edgelist = [[1,2,0.5],[2,3,0.8]]]
>>> Network = NNetwork()
>>> Network.add_wtd_edges(edgelist)
>>> Network.get_edge_weight([2,3])
0.8
```

Convert weighted graph to an unweighed graph by thresholding
```python
>>> Network_simple = Network.threshold2simple(0.7)
>>> Network_simple.edges()
[[2,3]]
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
