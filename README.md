# 🧬 Bioinformatics Toolkit

A comprehensive web-based Bioinformatics Toolkit offering a suite of simulated tools for molecular biology analysis and database lookups. This project demonstrates a clear separation of concerns with a Flask-based Python backend and a beautifully designed, responsive HTML/Tailwind CSS frontend.

## ✨ Features

This toolkit provides the following bioinformatics functionalities (with simulated data for demonstration purposes):

* **Medication Interaction Checker:** Input two medication names to check for potential interactions, their severity, description, and recommendations.
* **BLAST Sequence Alignment:** Simulate aligning a DNA or protein sequence against a database to find homologous regions, complete with mock query coverage, E-value, and identity.
* **KEGG Database Lookup:** Retrieve simulated information for KEGG pathways, drugs, or orthologies using their respective IDs.
* **PCR Primer Design:** Generate simulated forward and reverse PCR primers for a given DNA template, including basic length, GC content, and approximate melting temperature (Tm) data.
* **DNA to RNA Transcription:** Transcribe a DNA sequence into its corresponding RNA sequence.
* **RNA to Protein Translation:** Translate an RNA sequence into its corresponding protein sequence, recognizing start and stop codons.
* **Complementary Sequence Generator:** Create the complementary strand for either a DNA or RNA sequence.
* **UniProt Lookup:** Fetch simulated protein information (name, species, function) using UniProt accession IDs.
* **ChEMBL Lookup:** Retrieve simulated chemical compound data (name, molecular weight, type) using ChEMBL IDs.

## 🛠️ Technologies Used

**Backend:**
* **Python 3.x:** The core programming language.
* **Flask:** A lightweight web framework for building the API.
* **Flask-CORS:** For handling Cross-Origin Resource Sharing, allowing the frontend to communicate with the backend.

**Frontend:**
* **HTML5:** Structure of the web page.
* **JavaScript (ES6+):** For client-side logic, interacting with the backend API, and dynamic content updates.
* **Tailwind CSS:** A utility-first CSS framework for rapid and responsive UI development, enabling the pastel and modern design.

## 🚀 Getting Started

Follow these instructions to get a copy of the project up and running on your local machine.

### Prerequisites

* Python 3.x installed on your system.
* `pip` (Python package installer) is usually included with Python.

### Installation

1.  **Clone the repository (or create the files manually):**
    If you're using Git, you would clone:
    ```bash
    git clone [https://github.com/YourUsername/Bioinformatics-Toolkit.git](https://github.com/YourUsername/Bioinformatics-Toolkit.git)
    cd Bioinformatics-Toolkit
    ```
    Otherwise, create two files in the same directory: `apptwo.py` and `pagetwo.html`.

2.  **Backend Setup (`apptwo.py`):**
    * **Create/Update `apptwo.py`**: Copy the content of the `apptwo.py` code provided previously into this file.
    * **Install Python dependencies**:
        ```bash
        pip install Flask Flask-Cors
        ```
    * **Run the Flask server**:
        Open your terminal or command prompt, navigate to the directory where `apptwo.py` is saved, and run:
        ```bash
        python apptwo.py
        ```
        You should see output similar to `* Running on http://127.0.0.1:5000`. Keep this terminal window open; the backend must be running for the frontend to work.

3.  **Frontend Setup (`pagetwo.html`):**
    * **Create/Update `pagetwo.html`**: Copy the content of the `pagetwo.html` code provided previously into this file.
    * **Open the HTML file**: Simply open `pagetwo.html` in your preferred web browser (e.g., Chrome, Firefox, Edge) by double-clicking it or using "File > Open".

## 💡 Usage

1.  Ensure both the **backend server (`apptwo.py`) is running** in a terminal and the **frontend page (`pagetwo.html`) is open** in your web browser.
2.  On the frontend, use the **sidebar navigation** to select a specific bioinformatics tool.
3.  **Enter the required input** (e.g., medication names, DNA sequence, UniProt ID) into the input fields.
4.  Click the **"Run" or "Get Info" button** relevant to the tool.
5.  The results will be displayed in the designated **results area** below the input fields, complete with detailed information and a clean design.

## 📁 Project Structure
## 📝 License

This project is open source and available under the MIT License. Feel free to modify and distribute.
