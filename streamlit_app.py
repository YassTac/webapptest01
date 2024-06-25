import streamlit as st
from Bio import Entrez

# Configuration de l'email pour NCBI Entrez
Entrez.email = "votre_email@example.com"

# Fonction pour effectuer la requête NCBI
def search_ncbi(terms):
    handle = Entrez.esearch(db="pubmed", term=terms, retmax=10)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]

# Fonction pour récupérer les résumés des articles
def fetch_abstracts(id_list):
    ids = ",".join(id_list)
    handle = Entrez.efetch(db="pubmed", id=ids, rettype="abstract", retmode="text")
    abstracts = handle.read()
    handle.close()
    return abstracts

# Interface utilisateur Streamlit
st.title("Recherche d'articles scientifiques")
st.write("Entrez votre adresse email et les termes MeSH pour rechercher des articles scientifiques sur PubMed.")

# Champ pour entrer l'adresse email
email = st.text_input("Adresse email")

# Champ pour entrer les termes MeSH
mesh_terms = st.text_input("Termes MeSH")

if st.button("Rechercher"):
    if email and mesh_terms:
        # Mise à jour de l'email pour NCBI Entrez
        Entrez.email = email

        # Recherche des articles
        id_list = search_ncbi(mesh_terms)
        if id_list:
            st.write(f"{len(id_list)} articles trouvés.")
            abstracts = fetch_abstracts(id_list)
            st.text_area("Résumés des articles", abstracts, height=300)
        else:
            st.write("Aucun article trouvé pour ces termes.")
    else:
        st.write("Veuillez entrer une adresse email et des termes MeSH.")

if st.button("Réinitialiser"):
    email = ""
    mesh_terms = ""
