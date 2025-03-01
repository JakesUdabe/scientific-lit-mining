import pdfplumber
import os
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# Directories
PDF_DIR = "pdf_files"
TEXT_OUTPUT_DIR = "extracted_text"


def extract_text_from_pdf(pdf_path: str, output_dir: str):
    """
    Extracts text from a PDF and saves it as a .txt file.

    Args:
        pdf_path (str): Path to the PDF file.
        output_dir (str): Directory to save extracted text.
    """
    try:
        os.makedirs(output_dir, exist_ok=True)
        extracted_text = ""

        with pdfplumber.open(pdf_path) as pdf:
            for page in pdf.pages:
                page_text = page.extract_text(x_tolerance=1, y_tolerance=1) or ""
                extracted_text += page_text + "\n"

        text_filename = os.path.splitext(os.path.basename(pdf_path))[0] + ".txt"
        text_path = os.path.join(output_dir, text_filename)

        with open(text_path, "w", encoding="utf-8") as f:
            f.write(extracted_text)

        logging.info(f"Text extracted and saved to {text_path}")

    except Exception as e:
        logging.error(f"Error extracting text from {pdf_path}: {e}")


def process_pdfs(pdf_directory: str):
    """
    Processes all PDFs in the given directory to extract text.

    Args:
        pdf_directory (str): Directory containing PDF files.
    """
    if not os.path.exists(pdf_directory):
        logging.error(f"Directory '{pdf_directory}' not found.")
        return

    pdf_files = [f for f in os.listdir(pdf_directory) if f.lower().endswith(".pdf")]

    if not pdf_files:
        logging.warning(f"No PDF files found in '{pdf_directory}'.")
        return

    for pdf_file in pdf_files:
        pdf_path = os.path.join(pdf_directory, pdf_file)
        pdf_name = os.path.splitext(pdf_file)[0]
        text_subdir = os.path.join(TEXT_OUTPUT_DIR, pdf_name)

        logging.info(f"Processing {pdf_file} for text extraction...")
        extract_text_from_pdf(pdf_path, text_subdir)

    logging.info("Text extraction completed for all PDFs!")


if __name__ == "__main__":
    process_pdfs(PDF_DIR)
