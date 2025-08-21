from setuptools import find_packages, setup

setup(
    name='src',
    packages=find_packages(),
    version='0.1.0',
    description='This project is a complete data science pipeline designed to extract valuable insights from scientific literature. Focusing on nanomedicine as a case study, it uses Natural Language Processing (NLP) techniques to mine research papers from sources like PubMed. The pipeline covers the entire process, from data acquisition and preprocessing to modeling (e.g., topic modeling and text classification) and the visualization of key trends.',
    author='Jakes Udabe',
    license='MIT',
)
