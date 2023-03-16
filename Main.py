#Linux shebang line
#!/usr/bin/env python3

from Bio import Entrez # Import the Entrez module for accessing the NCBI databases
import sys # Import the sys module for handling command-line arguments

def fetch_sequences(database, term):

    # Email para Entrez API
    Entrez.email = "rendrickcarreira.social@gmail.com" # Replace with your email address
    Entrez.tool = "fetch_sequences.py"

    # Use esearch to search the database and save the results to a history
    search_results = Entrez.esearch(db=database, term=term, usehistory="y")
    search_results = Entrez.read(search_results)
    webenv = search_results["WebEnv"] # Save the WebEnv and QueryKey for use in later API calls
    query_key = search_results["QueryKey"]

    # Use efetch to fetch the results in batches and write to STDOUT
    batch_size = 10 # Number of records to fetch at a time
    for start in range(0, int(search_results["Count"]), batch_size):
        end = min(int(search_results["Count"]), start+batch_size) # Calculate the end of the current batch
        fetch_results = Entrez.efetch(db=database, rettype="fasta", retmode="text",
                                      retstart=start, retmax=batch_size, webenv=webenv, query_key=query_key)
        sys.stdout.write(fetch_results.read()) # Write the results to STDOUT
        fetch_results.close() # Close the handle to free up resources

if __name__ == '__main__':
    # Get the command-line arguments
    if len(sys.argv) != 3: # Check that there are exactly 2 arguments (the script name and the search term)
        print(__doc__) # If not, print the usage information and exit with an error code
        sys.exit(1)
    database = sys.argv[1] # Get the name of the database from the command line
    term = sys.argv[2] # Get the search term from the command line

    fetch_sequences(database, term) # Call the fetch_sequences function with the database and search term

##Link: https://stuntspt.gitlab.io/asb2023/classes/class_05/index.html#/
##Criado por
#David
#Helder
#Joana
#Rendrick