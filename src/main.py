import pandas as pd
import os
import spacy
import re
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from sklearn.feature_extraction.text import TfidfVectorizer, CountVectorizer
from sklearn.decomposition import LatentDirichletAllocation

# Carga el modelo de spaCy
try:
    nlp = spacy.load("en_core_web_sm")
except OSError:
    print("Descargando el modelo 'en_core_web_sm' de spaCy. Por favor, espera...")
    spacy.cli.download("en_core_web_sm")
    nlp = spacy.load("en_core_web_sm")


def preprocess_text(text):
    """Limpia y lematiza el texto de un abstract."""
    text = re.sub(r'[^a-zA-Z\s]', '', text, re.I|re.A)
    text = text.lower()
    doc = nlp(text)
    tokens = [token.lemma_ for token in doc if not token.is_stop and not token.is_punct and token.text.strip() != '']
    return " ".join(tokens)


def plot_top_words(model, feature_names, n_top_words, title, save_path):
    """Crea y guarda un gráfico de barras de los temas."""
    fig, axes = plt.subplots(1, model.n_components, figsize=(30, 10), sharex=True)
    axes = axes.flatten()
    for topic_idx, topic in enumerate(model.components_):
        ax = axes[topic_idx]
        top_features_ind = topic.argsort()[:-n_top_words - 1:-1]
        top_features = [feature_names[i] for i in top_features_ind]
        weights = [topic[i] for i in top_features_ind]
        
        ax.barh(top_features, weights, height=0.7)
        ax.set_title(f'Tema #{topic_idx + 1}', fontdict={'fontsize': 20})
        ax.invert_yaxis()
        ax.tick_params(axis='both', which='major', labelsize=15)
        for i in 'top right left'.split():
            ax.spines[i].set_visible(False)
    
    fig.suptitle(title, fontsize=25)
    plt.subplots_adjust(top=0.90, bottom=0.05, wspace=0.90)
    plt.savefig(save_path)
    plt.close(fig) # Cierra la figura para liberar memoria


def main():
    """Ejecuta el flujo completo de minería de literatura científica."""
    # Solicita el ID del proyecto
    project_id = input("Introduce un ID para tu proyecto: ")
    
    # Define las rutas del proyecto
    project_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    raw_data_path = os.path.join(project_dir, 'data', 'raw', 'pubmed_articles.csv')
    
    # Crea la carpeta de salida con el ID del proyecto
    output_dir = os.path.join(project_dir, 'data', 'processed', project_id)
    os.makedirs(output_dir, exist_ok=True)
    
    processed_data_path = os.path.join(output_dir, 'articles_with_topics.csv')
    topics_plot_path = os.path.join(output_dir, 'topics_bar_chart.png')
    co_occurrence_plot_path = os.path.join(output_dir, 'co_occurrence_graph.png')
    

    print("\n--- 1. Cargando y preprocesando datos ---")
    df = pd.read_csv(raw_data_path)
    df['clean_abstract'] = df['abstract'].apply(preprocess_text)
    print("¡Datos preprocesados con éxito!")

    print("\n--- 2. Realizando Topic Modeling ---")
    vectorizer = TfidfVectorizer(max_df=0.95, min_df=2, stop_words='english')
    data_vectorized = vectorizer.fit_transform(df['clean_abstract'])
    num_topics = 5
    lda = LatentDirichletAllocation(n_components=num_topics, random_state=42)
    lda.fit(data_vectorized)
    
    # Asignación de temas
    topic_probabilities = lda.transform(data_vectorized)
    most_likely_topic = np.argmax(topic_probabilities, axis=1)
    df['topic'] = most_likely_topic

    # Guarda el DataFrame
    df.to_csv(processed_data_path, index=False)
    print(f"¡Modelo de Topic Modeling entrenado y datos guardados en: {processed_data_path}!")

    # Guarda el gráfico de barras de los temas
    feature_names = vectorizer.get_feature_names_out()
    plot_top_words(lda, feature_names, 10, 'Temas en Abstracts de Nanopartículas', topics_plot_path)
    print(f"¡Gráfico de temas guardado en: {topics_plot_path}!")
    
    print("\n--- 3. Analizando Co-ocurrencia ---")
    filtered_df = df[(df['topic'] == 0) & (df['title'].str.contains('cancer', case=False))]
    
    if not filtered_df.empty:
        # Crea la matriz de co-ocurrencia
        co_vectorizer = CountVectorizer(stop_words='english', max_features=50)
        word_counts = co_vectorizer.fit_transform(filtered_df['clean_abstract'])
        co_occurrence_matrix = (word_counts.T * word_counts)
        co_occurrence_matrix.setdiag(0)

        # Crea y guarda el gráfico de red
        G = nx.from_scipy_sparse_array(co_occurrence_matrix)
        G = nx.relabel_nodes(G, {i: word for i, word in enumerate(co_vectorizer.get_feature_names_out())})
        edges_to_keep = [(u, v, d) for u, v, d in G.edges(data=True) if d['weight'] > 1]
        G_filtered = nx.Graph()
        G_filtered.add_edges_from(edges_to_keep)
        
        plt.figure(figsize=(15, 15))
        pos = nx.spring_layout(G_filtered, k=0.5, iterations=50)
        nx.draw_networkx_nodes(G_filtered, pos, node_size=200, node_color='skyblue')
        nx.draw_networkx_edges(G_filtered, pos, alpha=0.5)
        nx.draw_networkx_labels(G_filtered, pos, font_size=10, font_family="Arial")
        plt.title('Red de Co-ocurrencia para Artículos sobre Cáncer')
        plt.savefig(co_occurrence_plot_path)
        plt.close()
        
        print(f"¡Gráfico de co-ocurrencia guardado en: {co_occurrence_plot_path}!")
    else:
        print("No se encontraron artículos para el análisis de co-ocurrencia en el tema del cáncer.")


if __name__ == "__main__":
    main()