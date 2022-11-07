# uCDL

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://kchu25.github.io/uCDL.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://kchu25.github.io/uCDL.jl/dev/)
[![Build Status](https://github.com/kchu25/uCDL.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/kchu25/uCDL.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/kchu25/uCDL.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/kchu25/uCDL.jl)

# What is uCDL?
uCDL (Unfolded Convolutional Dictionary Learning) is a method for motif discovery. We first formulate a convolutional dictionary learning problem, and then "unfold" its optimization method into a neural network. The training of the network is fast. The resulting a network is fully interpretable, and outputs a sparse representation of the dataset. The sparse representation allow us to carry out efficient inference to discover motifs.

## Why use uCDL for motif discovery?


# How to Install
We are currently adding this package to the Julia registry. Once it's added, the user can simply install our package via the Julia's package manager:
```
pkg> add uCDL
```

## Software requirements
 This package requires [Weblogo](http://weblogo.threeplusone.com/manual.html#download). You need python3 and install Weblogo with following command:
 ```bash
 pip3 install weblogo
 ```
You can check if you have weblogo installed by typing:
```bash
weblogo -h
```

## Hardware requirements
We require the user to have an Nvidia GPU; we plan to implement a CPU version in the future.

# How to Use

In Julia, import the uCDL package first:
````julia
using uCDL
````
<br>


To do motif discovery on a single fasta file, execute
````
find_motif(<fasta-path>, <output-folder>)
````
- Perform motif discovery on the fasta file `<fasta-path>`, and 
- Output the result in a pre-specified folder `<output-folder>`. <br><br>



To do motif discovery on a batch of fasta files, execute
````
find_motif_fasta_folder(<fasta-folder-path>, <output-folder>)
````
- Perform motif discovery on all the fasta files in the `<fasta-folder-path>`, and 
- Output each of the results in a pre-specified folder `<output-folder>`.<br><br>

# Citing uCDL

If you use `uCDL` in your work, please cite
```
Shane Kuei-Hsien Chu, Gary D Stormo, Deep unfolded convolutional dictionary learning for motif discovery, bioRxiv 2022.11.06.515322; doi: https://doi.org/10.1101/2022.11.06.515322
```
