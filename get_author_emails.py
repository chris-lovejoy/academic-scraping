from Bio import Entrez
import argparse
import config
import pandas as pd
from datetime import date


# (1) Define email for query based on config.py file
Entrez.email = config.email

# (2) Parse command line input prompts
parser = argparse.ArgumentParser(description='Get author emails from pubmed')
parser.add_argument('search_terms', type=str, nargs=1,
                    help='The main terms to query, provided as a single string.')
parser.add_argument('--affiliation', type=str, nargs='+',
                    help='The organisational affiliation to query. Will limit to papers where \
                    1+ authors have that affiliation. Can provide multiple - if so, will only \
                    return papers where ALL listed affiliations are present.')
# TODO: consider adding further search term options

args = parser.parse_args()


def create_search_term(args):
    query = args.search_terms[0] # Index into the search term itself
    if args.affiliation:
        query += " " + args.affiliation[0] + "[Affiliation]"    
    return query

def run_search(query):
    handle = Entrez.esearch(db='pubmed', 
                            sort='relevance', 
                            retmax='1000',
                            retmode='xml', 
                            term=query)

    search = Entrez.read(handle)
    pubmed_ids = search['IdList']

    ids_string = ','.join(pubmed_ids)
    handle = Entrez.efetch(db='pubmed',
                        retmode='xml',
                        id=ids_string)

    results = Entrez.read(handle)

    return results


def create_output(results):
    # TODO: convert search results appropriate file output.
    # (probably CSV, perhaps provide options)

    print("Search completed. {n} results found. File saved as: {file_name}")

    return None


if __name__ == "__main__":
    query = create_search_term(args)
    results = run_search(query)
    create_output(results)
