# Generic-POA

[![Rust](https://github.com/anuradhawick/generic-poa/actions/workflows/rust.yml/badge.svg)](https://github.com/anuradhawick/generic-poa/actions/workflows/rust.yml)

Generic-POA is a command-line tool for performing Generalised Partial Order Alignment (POA). It offers flexibility in output formats and is designed for applications where alignment, visualisation and manipulation are necessary.

Core logic was adopted from [https://github.com/ljdursi/poapy](https://github.com/ljdursi/poapy). Following [Simpson Labs tutorial](http://simpsonlab.github.io/2015/05/01/understanding-poa/). Following are the original publications.

* [Christopher Lee, Catherine Grasso, Mark F. Sharlow, Multiple sequence alignment using partial order graphs , Bioinformatics, Volume 18, Issue 3, March 2002, Pages 452–464, https://doi.org/10.1093/bioinformatics/18.3.452](https://doi.org/10.1093/bioinformatics/18.3.452)
* [Christopher Lee, Generating consensus sequences from partial order multiple sequence alignment graphs, Bioinformatics, Volume 19, Issue 8, May 2003, Pages 999–1008, https://doi.org/10.1093/bioinformatics/btg109](https://doi.org/10.1093/bioinformatics/btg109)

This CLI program intends to generalise the POA from sequencing fragments to arbitrary objects in string form. Pontially aiming at applications where lenient matching of sequencial objects is needed (sentences, roads, etc)

## Installation

Clone this repository to your local machine using the following command:

```bash
git clone https://github.com/anuradhawick/generic-poa.git
cd generic-poa
```
Build the project using Cargo ([install](https://www.rust-lang.org/tools/install) rust if you already dont have it):

```bash
cargo build --release
```

## Usage

Use `gpoa --help` to get the complete command line help.

```bash
Generalised Partial Order Alignment

Usage: Generic-POA [OPTIONS] --input <INPUT> --output <OUTPUT>

Options:
  -i, --input <INPUT>    Input file path
  -o, --output <OUTPUT>  Input file path
      --html             Enable HTML output
      --graph            Enable graph output
  -h, --help             Print help
  -V, --version          Print version
```

## Examples

Consider the following file `examples/entries.tsv`.

```
seq_1	ACGT	ATTCC	ACGT	ACGT
seq_2	ACGT	ATTCC	ACGT
seq_3	ACGT	ATTCC	ACGT	TTGG	ACGT
```

You can run the following command.

```bash
mkdir test
./target/release/gpoa -i examples/entries.tsv  -o test/entries.tsv.aln --debug --graph
```

The output of this will be as follows. In the command line interface;

```
Graph    :  ACGT   ATTCC  ACGT   ACGT  
Match    :    |      |      |          
Alignment:  ACGT   ATTCC  ACGT     -   

Graph    :  ACGT   ATTCC  ACGT     -    ACGT  
Match    :    |      |      |             |   
Alignment:  ACGT   ATTCC  ACGT   TTGG   ACGT  
```

In `test/entries.tsv.aln` file;

```
seq_1  ACGT   ATTCC  ACGT     -    ACGT  
seq_2  ACGT   ATTCC  ACGT     -      -   
seq_3  ACGT   ATTCC  ACGT   TTGG   ACGT  
```

In `test/entries.tsv.aln.graph` file;

```
0:ACGT
	0 -> 1 ["seq_1", "seq_2", "seq_3"]
1:ATTCC
	1 -> 2 ["seq_1", "seq_2", "seq_3"]
2:ACGT
	2 -> 4 ["seq_3"]
	2 -> 3 ["seq_1"]
4:TTGG
	4 -> 3 ["seq_3"]
3:ACGT
```
