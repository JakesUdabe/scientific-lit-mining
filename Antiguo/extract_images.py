import fitz  # PyMuPDF
import os
import io
import logging
from PIL import Image

# Set up logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# Directories
PDF_DIR = "pdf_files"
IMAGE_OUTPUT_DIR = "extracted_images"


def extract_images_from_pdf(pdf_path: str, output_dir: str):
    """
    Extracts images from a PDF and saves them in an organized directory.

    Args:
        pdf_path (str): Path to the PDF file.
        output_dir (str): Directory to save extracted images.
    """
    try:
        os.makedirs(output_dir, exist_ok=True)
        extracted_images = []

        with fitz.open(pdf_path) as doc:
            for page_num, page in enumerate(doc):
                images = page.get_images(full=True)
                for i, img in enumerate(images):
                    xref = img[0]
                    base_image = doc.extract_image(xref)
                    image_bytes = base_image["image"]
                    image_ext = base_image["ext"]
                    image = Image.open(io.BytesIO(image_bytes))

                    image_name = f"image_{page_num}_{i}.{image_ext}"
                    image_path = os.path.join(output_dir, image_name)
                    image.save(image_path)
                    extracted_images.append(image_name)

        logging.info(f"Extracted {len(extracted_images)} images from {pdf_path}")

    except Exception as e:
        logging.error(f"Error extracting images from {pdf_path}: {e}")


def process_pdfs(pdf_directory: str):
    """
    Processes all PDFs in the given directory to extract images.

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
        image_subdir = os.path.join(IMAGE_OUTPUT_DIR, pdf_name)

        logging.info(f"Processing {pdf_file} for image extraction...")
        extract_images_from_pdf(pdf_path, image_subdir)

    logging.info("Image extraction completed for all PDFs!")


if __name__ == "__main__":
    process_pdfs(PDF_DIR)
