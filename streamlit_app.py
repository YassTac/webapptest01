import streamlit as st
from Bio import Entrez
import pandas as pd
from io import BytesIO

# Configuration de l'email pour NCBI Entrez
Entrez.email = "votre_email@example.com"

# Fonction pour effectuer la requête NCBI
def search_ncbi(terms, retmax):
    try:
        handle = Entrez.esearch(db="pubmed", term=terms, retmax=retmax)
        record = Entrez.read(handle)
        handle.close()
        return record["IdList"]
    except Exception as e:
        st.error(f"Erreur lors de la recherche dans PubMed: {e}")
        return []

# Fonction pour récupérer les résumés des articles
def fetch_abstracts(id_list):
    ids = ",".join(id_list)
    handle = Entrez.efetch(db="pubmed", id=ids, rettype="abstract", retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    
    abstracts = []
    for record in records['PubmedArticle']:
        article = record['MedlineCitation']['Article']
        title = article.get('ArticleTitle', 'No title available')
        abstract = article.get('Abstract', {}).get('AbstractText', ['No abstract available'])[0]

        # Tentative améliorée pour récupérer l'année de publication
        year = 'Not available'
        if 'ArticleDate' in article and article['ArticleDate']:
            article_date = article['ArticleDate'][0]
            year = article_date.get('Year', 'Not available')
        elif 'JournalIssue' in article['Journal'] and 'PubDate' in article['Journal']['JournalIssue']:
            pub_date = article['Journal']['JournalIssue']['PubDate']
            year = pub_date.get('Year', 'Not available')

        journal = article['Journal'].get('Title', 'No journal title available')

        # Récupération du DOI
        doi = 'No DOI available'
        if 'ELocationID' in article:
            for elocation in article['ELocationID']:
                if elocation.attributes.get('EIdType') == 'doi':
                    doi = elocation.get('#text', 'No DOI available')
                    break
        
        abstracts.append({'Title': title, 'Abstract': abstract, 'Year': year, 'Journal': journal, 'DOI': doi})
        
    return abstracts


# Fonction pour convertir les données en Excel et générer un lien de téléchargement
def convert_to_excel(data):
    df = pd.DataFrame(data)
    output = BytesIO()
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        df.to_excel(writer, index=False, sheet_name='Articles')
    processed_data = output.getvalue()
    return processed_data

# Interface utilisateur Streamlit
st.title("Recherche d'articles scientifiques")
st.write("Entrez votre adresse email et les termes MeSH pour rechercher des articles scientifiques sur PubMed.")

# Champ pour entrer l'adresse email
email = st.text_input("Adresse email")

# Champ pour entrer les termes MeSH
mesh_terms = st.text_input("Termes MeSH")

# Sélection du nombre d'articles à récupérer
num_articles = st.number_input("Nombre d'articles à récupérer", min_value=1, max_value=1000, value=10)

if st.button("Rechercher"):
    if email and mesh_terms:
        # Mise à jour de l'email pour NCBI Entrez
        Entrez.email = email

        # Recherche des articles
        id_list = search_ncbi(mesh_terms, num_articles)
        if id_list:
            st.write(f"{len(id_list)} articles trouvés.")
            abstracts = fetch_abstracts(id_list)
            st.dataframe(abstracts)

            # Conversion des données en Excel
            excel_data = convert_to_excel(abstracts)
            st.download_button(label="Télécharger les résultats en Excel", data=excel_data, file_name='articles.xlsx', mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')
        else:
            st.write("Aucun article trouvé pour ces termes.")
    else:
        st.write("Veuillez entrer une adresse email et des termes MeSH.")

if st.button("Réinitialiser"):
    st.experimental_rerun()
