import pdfplumber
import fitz
import camelot
import os
import re
import pandas as pd

def extract_tables_from_pdf(pdf_path, output_folder):
    os.makedirs(output_folder, exist_ok=True)
    doc = fitz.open(pdf_path)

    for page_num in range(len(doc)):
        page = doc[page_num]
        text = page.get_text("text")

        # Table 1 detection (spans columns, complex structure)
        if page_num == 0:  # Table 1 is on the first page
            try:
                # Use camelot with more specific parameters, adjust if needed
                tables = camelot.read_pdf(pdf_path, pages=str(page_num+1), flavor='stream', table_areas=['50,150,550,700'])  # Adjust coordinates
                for i, table in enumerate(tables):
                    df = table.df
                    # Data Cleaning specific to table 1
                    df = df.iloc[1:] #remove first row
                    df = df.dropna(axis=1, how='all') #remove empty columns
                    df = df.rename(columns=df.iloc[0]).drop(df.index[0]).reset_index(drop=True) #set the first row as columns and drop it
                    df = df.drop(df.columns[[6]], axis=1) #drop the last column
                    df.columns = ['Parameter', 'Aspasomes', 'units', 'Aspasomes', 'units', 'Conventional liposomes', 'units'] #set the name of the columns
                    output_file = os.path.join(output_folder, f"table1_page{page_num+1}.csv")
                    df.to_csv(output_file, index=False)
                    print(f"Table 1 extracted to: {output_file}")
            except Exception as e:
                print(f"Error extracting Table 1: {e}")

        # Table 2 and 3 detection (full page width)
        elif page_num in [1, 2]:  # Tables 2 and 3 are on pages 2 and 3
            try:
                tables = camelot.read_pdf(pdf_path, pages=str(page_num+1), flavor='stream')
                for i, table in enumerate(tables):
                    df = table.df
                    # Data Cleaning specific to table 2 and 3
                    df = df.iloc[1:] #remove first row
                    df = df.dropna(axis=1, how='all') #remove empty columns
                    df = df.rename(columns=df.iloc[0]).drop(df.index[0]).reset_index(drop=True) #set the first row as columns and drop it
                    output_file = os.path.join(output_folder, f"table{i+2}_page{page_num+1}.csv") #Corrected file naming
                    df.to_csv(output_file, index=False)
                    print(f"Table {i+2} extracted to: {output_file}")
            except Exception as e:
                print(f"Error extracting Table {i+2}: {e}")

# Example usage:
pdf_path = r"C:\xxx\xxx.pdf"  # Use raw string for file path
output_folder = "extracted_tables"
extract_tables_from_pdf(pdf_path, output_folder)
print(f"Tables extracted to: {output_folder}")
