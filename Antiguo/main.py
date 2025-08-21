import fitz  # PyMuPDF
import pdfplumber
import pytesseract
from PIL import Image
import io
import os
import logging
from typing import Dict, List, Optional

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def extract_abstract(text: str) -> str:
    """
    Extracts the abstract from the first page text.
    Stops extraction when encountering "Keywords", "Key words", or "All rights reserved".
    """
    abstract_start = text.find("Abstract")
    if abstract_start == -1:
        return ""  # No abstract found

    # Define stopping points for the abstract
    stopping_phrases = ["Keywords", "Key words", "All rights reserved"]
    abstract_end = len(text)
    for phrase in stopping_phrases:
        phrase_pos = text.find(phrase, abstract_start)
        if phrase_pos != -1 and phrase_pos < abstract_end:
            abstract_end = phrase_pos

    abstract = text[abstract_start:abstract_end].strip()
    return abstract

def detect_column_layout(plumber_page) -> str:
    """
    Detects the column layout of a page (single-column or two-column).
    """
    words = plumber_page.extract_words()
    if not words:
        return "single"  # Default to single-column if no words are found

    # Calculate the horizontal midpoint of the page
    page_width = plumber_page.width
    midpoint = page_width / 2

    # Check if words are distributed across both halves of the page
    left_count = sum(1 for word in words if word["x0"] < midpoint)
    right_count = sum(1 for word in words if word["x0"] >= midpoint)

    if left_count > 0 and right_count > 0:
        return "two-column"
    else:
        return "single-column"

def extract_text_from_columns(plumber_page) -> str:
    """
    Extracts text from a page, handling both single-column and two-column layouts.
    """
    layout = detect_column_layout(plumber_page)
    page_width = plumber_page.width
    page_height = plumber_page.height

    if layout == "two-column":
        # Define bounding boxes for left and right columns
        left_bbox = (0, 0, page_width / 2, page_height)
        right_bbox = (page_width / 2, 0, page_width, page_height)

        # Extract text from left column first
        left_text = plumber_page.within_bbox(left_bbox).extract_text(x_tolerance=1, y_tolerance=1) or ""
        right_text = plumber_page.within_bbox(right_bbox).extract_text(x_tolerance=1, y_tolerance=1) or ""

        return left_text + "\n" + right_text
    else:
        # Single-column layout
        return plumber_page.extract_text(x_tolerance=1, y_tolerance=1) or ""

def extract_data_from_pdf(pdf_path: str, output_dir: str = "extracted_data") -> Optional[Dict[str, List[str]]]:
    """
    Extracts structured data (text and figures) from a scientific PDF.

    Args:
        pdf_path (str): The path to the PDF file.
        output_dir (str): The directory to save extracted data. Defaults to "extracted_data".

    Returns:
        dict: A dictionary containing the extracted data (text, figures). Returns None if errors.
    """
    try:
        os.makedirs(output_dir, exist_ok=True)

        extracted_data = {
            "abstract": "",
            "text": "",
            "figures": []
        }

        with fitz.open(pdf_path) as doc, pdfplumber.open(pdf_path) as pdf:
            for page_num, (page, plumber_page) in enumerate(zip(doc, pdf.pages)):
                # 1. Abstract Extraction (First Page Only)
                if page_num == 0:
                    first_page_text = plumber_page.extract_text(x_tolerance=1, y_tolerance=1) or ""
                    extracted_data["abstract"] = extract_abstract(first_page_text)

                # 2. Text Extraction (Handling Column Layouts)
                page_text = extract_text_from_columns(plumber_page)
                extracted_data["text"] += page_text + "\n"

                # 3. Figure Extraction
                images = page.get_images(full=True)
                for i, img in enumerate(images):
                    xref = img[0]
                    base_image = doc.extract_image(xref)
                    image_bytes = base_image["image"]
                    image_ext = base_image["ext"]
                    image = Image.open(io.BytesIO(image_bytes))

                    image_name = f"figure_page_{page_num+1}_figure_{i+1}.{image_ext}"
                    image_path = os.path.join(output_dir, image_name)
                    image.save(image_path)
                    extracted_data["figures"].append(image_name)

        # Save extracted text and abstract
        with open(os.path.join(output_dir, "extracted_text.txt"), "w", encoding="utf-8") as f:
            f.write(extracted_data["text"])
        with open(os.path.join(output_dir, "abstract.txt"), "w", encoding="utf-8") as f:
            f.write(extracted_data["abstract"])

        return extracted_data

    except Exception as e:
        logging.error(f"An error occurred: {e}")
        return None

# Example usage:
pdf_file_path = r"C:\Users\xxx\xxx.pdf"  # Use raw string for file path
extracted_data = extract_data_from_pdf(pdf_file_path)

if extracted_data:
    logging.info("Data extracted successfully!")
else:
    logging.error("Data extraction failed.")
