import pandas as pd
import os
from Bio import Entrez
import time
import json

def fetch_pubmed_articles_by_year(query, year, max_results=100):
    """
    Busca artículos en PubMed para una consulta y un año específicos.
    Devuelve los 100 abstracts más relevantes y el total de artículos encontrados para ese año.
    """
    Entrez.email = "tu.email@ejemplo.com" # Reemplaza con tu correo electrónico
    search_term = f'("{query}")[Title/Abstract] AND {year}[pdat]'
    print(f"Buscando en PubMed para la consulta: '{search_term}'")

    try:
        # Obtener el número total de publicaciones para el año.
        search_handle = Entrez.esearch(db="pubmed", term=search_term, retmax=1)
        record = Entrez.read(search_handle)
        total_results = int(record["Count"])
        print(f"Se encontraron un total de {total_results} artículos para el año {year}.")

        # Obtener los 100 abstracts más relevantes.
        search_handle = Entrez.esearch(db="pubmed", term=search_term, retmax=max_results, sort="relevance")
        record = Entrez.read(search_handle)
        id_list = record["IdList"]

        if not id_list:
            print(f"No se encontraron abstracts para el año {year}.")
            return None, 0

        fetch_handle = Entrez.efetch(db="pubmed", id=id_list, retmode="xml")
        articles = Entrez.read(fetch_handle)

        data = []
        for article in articles['PubmedArticle']:
            try:
                pmid = article['MedlineCitation']['PMID']
                title = article['MedlineCitation']['Article']['ArticleTitle']
                abstract_text = article['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
                data.append({
                    "pmid": str(pmid),
                    "title": title,
                    "abstract": abstract_text
                })
            except KeyError:
                continue
        
        return pd.DataFrame(data), total_results

    except Exception as e:
        print(f"Error al conectar con PubMed: {e}")
        return None, 0

if __name__ == '__main__':
    # 1. Pide un ID de proyecto, término de búsqueda, y el rango de años al usuario
    project_id = input("Introduce un ID para tu proyecto: ")
    search_term = input("Introduce el término de búsqueda (ej. 'nanoparticle'): ")
    
    try:
        start_year = int(input("Introduce el año de inicio (ej. 2013): "))
        end_year = int(input("Introduce el año de fin (ej. 2024): "))
    except ValueError:
        print("Entrada de año no válida. Por favor, introduce un número.")
        exit()

    # Define la ruta base para el proyecto
    project_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data', 'raw', project_id)
    os.makedirs(project_dir, exist_ok=True)
    
    # Lista para almacenar el recuento de publicaciones por año
    yearly_counts = []

    for year in range(start_year, end_year + 1): # Se suma 1 para incluir el año de fin
        df, total_count = fetch_pubmed_articles_by_year(search_term, year, max_results=100)
        
        # Almacena el recuento para el año
        yearly_counts.append({'año': year, 'publicaciones': total_count})
        
        if df is not None and not df.empty:
            # Define la ruta para la carpeta de abstracts por año
            output_dir = os.path.join(project_dir, search_term, str(year))
            os.makedirs(output_dir, exist_ok=True)
            
            # Guarda cada artículo como un archivo JSON separado
            for _, row in df.iterrows():
                file_path = os.path.join(output_dir, f"{row['pmid']}.json")
                if not os.path.exists(file_path):
                    doc_dict = row.to_dict()
                    with open(file_path, 'w') as json_file:
                        json.dump(doc_dict, json_file, indent=4)
                    print(f"Artículo {row['pmid']} guardado en: {file_path}")
                else:
                    print(f"SKIP: El archivo {row['pmid']}.json ya existe.")
        else:
            print(f"No se encontraron abstracts para el año {year}.")

    # Guarda el recuento de publicaciones por año en un archivo CSV
    if yearly_counts:
        counts_df = pd.DataFrame(yearly_counts)
        csv_path = os.path.join(project_dir, f"publicaciones_por_año_{search_term}.csv")
        counts_df.to_csv(csv_path, index=False)
        print(f"\nRecuento de publicaciones guardado en: {csv_path}")