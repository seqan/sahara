<!--
    SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
    SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
    SPDX-License-Identifier: CC-BY-4.0
-->

# Sahara

Approximate searches using Optimum Search Schemes

## Usage
1. Create an index for a fasta file:
```bash
    $ sahara index somefastafile.fasta
```

2. Search inside the index and allow 2 errors accordingly to edit distance:
```bash
    $ sahara search --index somefastafile.fasta.idx --query queryfile.fasta --errors 2
```

## Compile from Source

To compile the source, download it through git and build it with cmake/make.
At the end we execute the help page.
```bash
git clone https://github.com/seqan/sahara.git
mkdir -p sahara/build && cd $_
cmake -DCMAKE_BUILD_TYPE=Release ..
make
./src/sahara/sahara --help
```
