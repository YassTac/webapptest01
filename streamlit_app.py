import streamlit as st
import pandas as pd
from Bio import Entrez

def pubmed_search(email, query, retmax):
    Entrez.email = email
    handle = Entrez.esearch(db="pubmed", term=query, retmax=retmax)
    record = Entrez.read(handle)
    ids = record["IdList"]

    articles = []
    for id in ids:
        handle = Entrez.efetch(db="pubmed", id=id, retmode="xml")
        record = Entrez.read(handle)
        article = {
            "title": record["MedlineCitation"]["Article"]["ArticleTitle"],
            "abstract": record["MedlineCitation"]["Article"]["Abstract"]["AbstractText"],
            "year": record["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"]["Year"],
            "journal": record["MedlineCitation"]["Article"]["Journal"]["Title"],
            "doi": record["PubmedData"]["ArticleIdList"][0]
        }
        articles.append(article)

    return pd.DataFrame(articles)

st.title("PubMed Search")

email = st.text_input("Enter your email address:")
query = st.text_input("Enter your search query:")
retmax = st.number_input("Enter the maximum number of articles to retrieve:", min_value=1, step=1)

if st.button("Search"):
    df = pubmed_search(email, query, retmax)
    st.write(df)
    st.download_button(
        label="Download as Excel",
        data=df.to_excel(index=False),
        file_name="pubmed_search_results.xlsx",
        mime="application/vnd.ms-excel",
    )
