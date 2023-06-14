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

