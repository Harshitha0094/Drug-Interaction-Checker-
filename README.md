# 🧬 Bio-Toolkit: Integrated Full-Stack Bioinformatics Suite

![Python](https://img.shields.io/badge/Python-3.9+-3776AB?style=for-the-badge&logo=python&logoColor=white)
![Flask](https://img.shields.io/badge/Flask-1.1.2+-000000?style=for-the-badge&logo=flask&logoColor=white)
![TailwindCSS](https://img.shields.io/badge/Tailwind_CSS-3.0+-38B2AC?style=for-the-badge&logo=tailwind-css&logoColor=white)

## 📖 Project Overview
The **Bio-Toolkit** is a high-performance web application designed to streamline common bioinformatics workflows and clinical decision-support tasks. Developed as a final-year project for **B.E. Bioinformatics Engineering**, it bridges the gap between raw genomic data and clinical application.

---

## 🛠️ Core Modules

### 1. 💊 Clinical Informatics
* **Drug-Drug Interaction (DDI) Checker:** A normalized lookup algorithm to cross-reference drug pairs.
* **Severity Classification:** Categorizes results into `Major`, `Moderate`, and `Minor` risks.
* **Database Integration:** Simulated mappings for **DrugBank** and **ChEMBL** standards.

### 2. 🧬 Genomic Sequence Engineering
* **Transcription Engine:** Automated DNA $\rightarrow$ RNA conversion.
* **Translation Engine:** A high-fidelity implementation of the **Standard Genetic Code** featuring:
    * Start codon (`AUG`) recognition.
    * Stop codon (`UAA`, `UAG`, `UGA`) termination logic.
* **PCR Primer Design:** Mathematical approach to designing oligonucleotide sequences based on Target Region and Melting Temperature ($T_m$).

### 3. 🧪 Systems Biology & Databases
* **UniProt:** Protein-level functional annotations and species data.
* **KEGG:** Mapping genomic data to metabolic and signaling pathways.
* **BLAST:** Sequence homology simulation for genomic research.

---

## 💻 Technical Stack

### **Backend (Python/Flask)**
* **RESTful API:** Decoupled endpoints for independent modular updates.
* **CORS Management:** Configured via `flask-cors` for secure frontend communication.
* **Bio-Algorithms:** Custom implementation of transcription/translation logic.

### **Frontend (JS/Tailwind CSS)**
* **SPA Architecture:** Single-page design for a seamless, "no-refresh" user experience.
* **Responsive UI:** Professional "Scientific Dashboard" aesthetic optimized for all screen sizes.

---

## 🚀 Installation & Setup

### 1. Install Dependencies
```bash
pip install flask flask-cors
