import pandas as pd
from Bio import Entrez
import time
import os

def fetch_pubmed_articles(query, max_results=1000):
    """
    Busca artículos en PubMed y extrae el título y el abstract.
    """
    Entrez.email = "nijakes@gmail.com"  # Reemplaza con tu correo electrónico
    search_handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, sort="relevance")
    record = Entrez.read(search_handle)
    id_list = record["IdList"]

    if not id_list:
        print("No se encontraron artículos.")
        return None

    fetch_handle = Entrez.efetch(db="pubmed", id=id_list, rettype="fasta", retmode="xml")
    articles = Entrez.read(fetch_handle)
    
    data = []
    for article in articles['PubmedArticle']:
        try:
            # Extracción del título y el abstract
            title = article['MedlineCitation']['Article']['ArticleTitle']
            abstract_text = article['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
            data.append({"title": title, "abstract": abstract_text})
        except KeyError:
            continue
            
    return pd.DataFrame(data)

if __name__ == '__main__':
    # Reemplaza la consulta con tu búsqueda específica
    query = "nanoparticles AND drug delivery"
    df = fetch_pubmed_articles(query, max_results=100)

    if df is not None:
        # Define la ruta para guardar el archivo
        project_dir = os.path.join(os.path.dirname(__file__), '..', '..')
        data_dir = os.path.join(project_dir, 'data', 'raw')
        
        # Crea la carpeta si no existe
        os.makedirs(data_dir, exist_ok=True)
        
        # Guarda los datos en un archivo CSV
        file_path = os.path.join(data_dir, 'pubmed_articles.csv')
        df.to_csv(file_path, index=False)
        print(f"Datos extraídos y guardados en: {file_path}")
    else:
        print("La extracción de datos no pudo completarse.")