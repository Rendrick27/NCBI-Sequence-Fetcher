# NCBI Sequence Fetcher
This is a Python script that fetches sequences from a database using the Bio.Entrez library.

This work was proposed in the Curricular Unit of Biology Analysis and Sequences of the bioinformatics course https://stuntspt.gitlab.io/ASB/classes/class_05/index.html#/


## Requisites
* Python
* Bio.Entrez library

## Installation
### Biopython
```bash
pip install biopython
```

## Usage
```bash
python fetch_sequences.py [database name] [search term]
```
Where [database name] is the name of the database you want to search (e.g., "nucleotide" for the NCBI nucleotide database), and [search term] is the term you want to search for.

### Examples
To search the NCBI nucleotide database for sequences related to "E. coli", run:
```bash
python fetch_sequences.py nucleotide "Escherichia coli"
```

This will fetch the sequences in batches and print them to the console.

## Credits
This script uses the Bio.Entrez library from the Biopython project.

## License
MIT License
