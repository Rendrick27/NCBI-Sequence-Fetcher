#!/usr/bin/env python3

import sys
from Bio import Entrez

# Parse command line arguments
if len(sys.argv) != 3:
    print("Usage: python3 retrieve_sequences.py <database> <search_term>")
    sys.exit(1)

database = sys.argv[1]
search_term = sys.argv[2]

# Set email address (required by Entrez)
Entrez.email = "rendrickcarreira.social@gmail.com"

# Use Entrez API to search for sequences
handle = Entrez.esearch(db=database, term=search_term)
record = Entrez.read(handle)
handle.close()

# Retrieve sequences from search results
webenv = record["WebEnv"]
query_key = record["QueryKey"]
handle = Entrez.efetch(db=database, rettype="fasta", retmode="text", webenv=webenv, query_key=query_key, usehistory="y")
sequences = handle.read()
handle.close()

# Print sequences to stdout
print(sequences)
