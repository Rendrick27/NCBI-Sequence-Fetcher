from Bio import Entrez
from Bio import SeqIO
import sys


def parse_command_line_arguments():
    """
    Parses commwand line arguments.
    Expects two arguments besides the script name: 
        the database name 
        the search term.
    Returns:
        db_name (str): The name of the database to query.
        search_term (str): The search term to use for querying the database.
    """
    if len(sys.argv) != 3:
        print("Usage: python script.py <database> <search_term>")
        sys.exit(1)
    db_name = sys.argv[1]
    search_term = sys.argv[2]
    return db_name, search_term


def search_database(db_name, search_term):
    """
    Uses the Entrez API to search for the given search term in the database.
    Utilizes the history feature to facilitate retrieval of search results.
    Arguments:
        db_name (str): The database to search.
        search_term (str): The term to search for.
    Returns:
        search_results: The results of the search, including WebEnv
                                                   and QueryKey.
    """
    Entrez.email = "rendrickcarreira@gmail.com"
    search_handle = Entrez.esearch(db=db_name, term=search_term,
                                   usehistory="y")
    search_results = Entrez.read(search_handle)
    search_handle.close()
    return search_results


def fetch_sequences(db_name, search_results):
    """
    Fetches the sequences from the database based on the search results.
    The results are fetched in FASTA format and printed to STDOUT.
    Arguments:
        db_name (str): The database from which to fetch the sequences.
        search_results: The search results containing WebEnv and QueryKey.
    """
    webenv = search_results["WebEnv"]
    query_key = search_results["QueryKey"]
    fetch_handle = Entrez.efetch(db=db_name, rettype="fasta", retmode="text",
                                 webenv=webenv, query_key=query_key)
    for seq_record in SeqIO.parse(fetch_handle, "fasta"):
        print(seq_record.format("fasta"))
    fetch_handle.close()


def main():
    """
    Main function to orchestrate the querying 
        and fetching of sequences from the Entrez database.
    """
    db_name, search_term = parse_command_line_arguments()
    search_results = search_database(db_name, search_term)
    fetch_sequences(db_name, search_results)


if __name__ == "__main__":
    main()
    # This is a line that runs the script.