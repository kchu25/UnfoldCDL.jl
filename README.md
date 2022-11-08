# UnfoldCDL

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://kchu25.github.io/UnfoldCDL.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://kchu25.github.io/UnfoldCDL.jl/dev/)
[![Build Status](https://github.com/kchu25/UnfoldCDL.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/kchu25/UnfoldCDL.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/kchu25/UnfoldCDL.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/kchu25/UnfoldCDL.jl)

# What is UnfoldCDL?
UnfoldCDL (**Unfolded Convolutional Dictionary Learning**) is a method for [DNA sequence motif discovery](https://en.wikipedia.org/wiki/Sequence_motif). We first formulate a convolutional dictionary learning problem and then "unfold" its optimization algorithm into a neural network. The network is *fully interpretable, fast to train*, and outputs a *sparse representation* of the dataset. The sparse representation allows us to infer the motifs in the dataset efficiently.

## Why use UnfoldCDL for motif discovery?
Many methods can find statistically significant motifs, but they may be only "partially found" because motifs may have complex binding patterns. For example, motifs may have multiple modes: each mode may share similar patterns (multimeric binding, alternate structural conformations), have distinct "parts" (variable spacing), or have multiple motifs that look entirely dissimilar to each other (multiple DNA binding domains). Some traditional motif discovery method uses heuristics such as substring-masking to deal with the above scenarios, which leads to a sequential motif discovery method, and results in some secondary motif being masked and not revealed. The inference on the motifs using other black-box deep learning approaches is currently challenging.

The sparse representation we obtained from UnfoldCDL reveals where the enriched patterns are in the dataset, and we seek to use such representation to discover all the motifs simultaneously.


# How to Install
We are currently adding this package to the Julia registry. Once it's added, the user can simply install our package via the Julia's package manager:
```julia
pkg> add UnfoldCDL
```

## Software requirements
 This package requires [Weblogo](http://weblogo.threeplusone.com/manual.html#download). You need python3 and install Weblogo with following command:
 ```bash
 pip3 install weblogo
 ```

## Hardware requirements
We require the user to have an Nvidia GPU; we plan to implement a CPU version in the future.

# How to Use

In Julia, import the UnfoldCDL package first:
````julia
using UnfoldCDL
````
<br>


To do motif discovery on a single fasta file, execute
````julia
find_motif(<fasta-path>, <output-folder>)
````
- Perform motif discovery on the fasta file `<fasta-path>`, and 
- Output the result in a pre-specified folder `<output-folder>`. <br><br>



To do motif discovery on a batch of fasta files, execute
````julia
find_motif_fasta_folder(<fasta-folder-path>, <output-folder>)
````
- Perform motif discovery on all the fasta files in the `<fasta-folder-path>`, and 
- Output each of the results in a pre-specified folder `<output-folder>`.<br><br>

# Citation

The paper for UnfoldCDL is at [https://www.biorxiv.org/content/10.1101/2022.11.06.515322v2](https://www.biorxiv.org/content/10.1101/2022.11.06.515322v2). It can be cited using the following BibTex entry:
```
@article {Chu2022.11.06.515322,
	author = {Chu, Shane Kuei-Hsien and Stormo, Gary D},
	title = {Deep unfolded convolutional dictionary learning for motif discovery},
	elocation-id = {2022.11.06.515322},
	year = {2022},
	doi = {10.1101/2022.11.06.515322},
	publisher = {Cold Spring Harbor Laboratory},,
	URL = {https://www.biorxiv.org/content/early/2022/11/07/2022.11.06.515322},
	eprint = {https://www.biorxiv.org/content/early/2022/11/07/2022.11.06.515322.full.pdf},
	journal = {bioRxiv}
}
```

