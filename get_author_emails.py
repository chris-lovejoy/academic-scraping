from Bio import Entrez
import argparse
import config
import pandas as pd
from datetime import date
import re

def create_search_term(args):
    query = args.search_terms[0] # Index into the search term itself
    if args.affiliation:
        for affil in args.affiliation:
            query += " " + affil + "[Affiliation]"
    return query

def run_search(query):
    try:
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
    except Exception as e:
        # print(e)
        return None

def create_dataframe(results):
    
    # Create empty dataframe
    dataframe = pd.DataFrame(columns=['emails'])

    # Loop over all pubmed articles in search results
    for i, paper in enumerate(results['PubmedArticle']):
        
        # Create empty dataframe for each paper
        single_paper_df = pd.DataFrame()

        try:
            # Get emails by regex across string of all author information
            string = str(results['PubmedArticle'][i]['MedlineCitation']['Article']['AuthorList'])
            emails = re.findall('\S+@[a-zA-Z]+.[a-zA-Z]+.[a-zA-Z]+.[a-zA-Z]+', string)
            
            if emails != []:
                emails_to_add = []

                for index, email in enumerate(emails):

                    # If email not already in list, add it
                    if emails[index] not in list(dataframe['emails']):
                        emails_to_add.append(emails[index])

                emails_series = pd.Series(emails_to_add, dtype=str)
                single_paper_df['emails'] = emails_series
                            
                # Add other paper information
                single_paper_df['paper_title'] = results['PubmedArticle'][i]['MedlineCitation']['Article']['ArticleTitle']
                single_paper_df['year_published'] = results['PubmedArticle'][i]['MedlineCitation']['Article']['ArticleDate'][0]['Year']
                
                dataframe = pd.concat([dataframe, single_paper_df])

        except Exception as e:
            # print(e)
            pass

    return dataframe

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

if __name__ == "__main__":
    query = create_search_term(args)
    results = run_search(query)
    if results:
        dataframe = create_dataframe(results)

        # Name file based on search term and date, and save
        file_name = f"./output/{query}_{str(date.today())}.csv"
        dataframe.to_csv(file_name)
        print(f"Search completed. {len(dataframe)} results found. File saved as: '{file_name}'")
    else:
        print("No results found.")
