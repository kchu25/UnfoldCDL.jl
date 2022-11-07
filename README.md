# uCDL

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://kchu25.github.io/uCDL.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://kchu25.github.io/uCDL.jl/dev/)
[![Build Status](https://github.com/kchu25/uCDL.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/kchu25/uCDL.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/kchu25/uCDL.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/kchu25/uCDL.jl)

# What is uCDL?

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


