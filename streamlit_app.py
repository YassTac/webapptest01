import streamlit as st
import pandas as pd
from Bio import Entrez
import io

def pubmed_search(email, query, retmax):
    Entrez.email = email
    handle = Entrez.esearch(db="pubmed", term=query, retmax=retmax)
    record = Entrez.read(handle)
    ids = record["IdList"]

    articles = []
    for id in ids:
        handle = Entrez.efetch(db="pubmed", id=id, retmode="xml")
        record = Entrez.read(handle)
        article = {}
        try:
            article["title"] = record["MedlineCitation"]["Article"]["ArticleTitle"]
        except KeyError:
            article["title"] = "N/A"
        try:
            article["abstract"] = "\n".join(record["MedlineCitation"]["Article"]["Abstract"]["AbstractText"])
        except KeyError:
            article["abstract"] = "N/A"
        try:
            article["year"] = record["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"]["Year"]
        except KeyError:
            article["year"] = "N/A"
        try:
            article["journal"] = record["MedlineCitation"]["Article"]["Journal"]["Title"]
        except KeyError:
            article["journal"] = "N/A"
        try:
            article["doi"] = record["PubmedData"]["ArticleIdList"][0]
        except (KeyError, IndexError):
            article["doi"] = "N/A"
        articles.append(article)

    return pd.DataFrame(articles)

st.title("PubMed Search")

email = st.text_input("Enter your email address:")
query = st.text_input("Enter your search query:")
retmax = st.number_input("Enter the maximum number of articles to retrieve:", min_value=1, step=1)

if st.button("Search"):
    df = pubmed_search(email, query, retmax)
    st.write(df)
    if not df.empty:
        excel_file = io.BytesIO()
        df.to_excel(excel_file, index=False)
        st.download_button(
            label="Download as Excel",
            data=excel_file.getvalue(),
            file_name="pubmed_search_results.xlsx",
            mime="application/vnd.ms-excel",
        )
    else:
        st.write("No results found.")
