import pandas as pd
import os
import requests
from Bio import Entrez
import time

def fetch_pubmed_articles(query, max_results=1000):
    """Busca y extrae artículos de PubMed."""
    Entrez.email = "nijakes@gmail.com" # ¡IMPORTANTE! Reemplaza esto con tu correo
    print(f"Buscando en PubMed para la consulta: '{query}'")
    try:
        search_handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, sort="relevance")
        record = Entrez.read(search_handle)
        id_list = record["IdList"]

        if not id_list:
            print("No se encontraron artículos en PubMed.")
            return None

        fetch_handle = Entrez.efetch(db="pubmed", id=id_list, retmode="xml")
        articles = Entrez.read(fetch_handle)
        
        data = []
        for article in articles['PubmedArticle']:
            try:
                title = article['MedlineCitation']['Article']['ArticleTitle']
                abstract_text = article['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
                data.append({"title": title, "abstract": abstract_text})
            except KeyError:
                continue
        
        return pd.DataFrame(data)

    except Exception as e:
        print(f"Error al conectar con PubMed: {e}")
        return None


def fetch_scopus_articles(query, api_key, max_results=1000):
    """Busca y extrae artículos de Scopus."""
    # ¡IMPORTANTE! Necesitas una API key de Scopus para que esto funcione.
    url = "https://api.elsevier.com/content/search/scopus"
    headers = {
        "X-ELS-APIKey": api_key,
        "Accept": "application/json"
    }
    params = {
        "query": query,
        "count": max_results,
        "view": "COMPLETE"
    }
    
    print(f"Buscando en Scopus para la consulta: '{query}'")
    try:
        response = requests.get(url, headers=headers, params=params)
        response.raise_for_status()
        data = response.json()
        
        articles = data.get('search-results', {}).get('entry', [])
        
        article_list = []
        for article in articles:
            try:
                title = article.get('dc:title', '')
                # El abstract no siempre está disponible en la vista completa de la búsqueda
                abstract = article.get('prism:url', '') 
                article_list.append({"title": title, "abstract": abstract})
            except KeyError:
                continue

        return pd.DataFrame(article_list)
        
    except requests.exceptions.RequestException as e:
        print(f"Error al conectar con Scopus: {e}")
        return None
    

if __name__ == '__main__':
    # 1. Pide un ID de proyecto al usuario
    project_id = input("Introduce un ID para tu proyecto: ")

    # 2. Pide al usuario que elija la fuente de datos
    source = input("¿Quieres buscar en 'pubmed' o en 'scopus'? ")
    source = source.lower().strip() # Normaliza la entrada

    # Define la consulta de búsqueda
    query = "nanoparticles AND drug delivery"
    df = None

    if source == 'pubmed':
        df = fetch_pubmed_articles(query, max_results=100)
    elif source == 'scopus':
        scopus_api_key = input("ef8659e9fa072768e2e6a32e3d3b5f0b")
        df = fetch_scopus_articles(query, scopus_api_key, max_results=100)
    else:
        print("Fuente de datos no reconocida. Por favor, elige 'pubmed' o 'scopus'.")

    # Si se han obtenido datos, los guarda en una carpeta con el ID del proyecto
    if df is not None:
        # Define la ruta para la carpeta del proyecto
        project_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        output_dir = os.path.join(project_dir, 'data', 'raw', project_id)
        
        # Crea la carpeta si no existe
        os.makedirs(output_dir, exist_ok=True)
        
        # Guarda los datos en un archivo CSV
        file_path = os.path.join(output_dir, 'articles.csv')
        df.to_csv(file_path, index=False)
        print(f"Datos extraídos y guardados en: {file_path}")
    else:
        print("La extracción de datos no pudo completarse.")