# 🧬 Bioinformatics Tools & Drug Interaction Checker

A full-stack web application that provides a suite of bioinformatics tools including a **Drug/Medication Interaction Checker**, BLAST search simulation, KEGG pathway lookup, primer design, sequence transcription/translation, complementary sequence generation, and database lookups for UniProt and ChEMBL.

---

## 📸 Preview

> A clean, responsive sidebar-based UI built with Tailwind CSS, backed by a Flask REST API.

---

## 🚀 Features

| Tool | Description |
|------|-------------|
| 💊 **Medication Interaction Checker** | Check interactions between two drugs with severity levels (Major / Moderate / Minor / None) |
| 🔬 **BLAST Tool** | Simulate sequence alignment and find significant hits |
| 🧪 **KEGG Pathway Lookup** | Retrieve KEGG pathway information by ID (e.g., `hsa00010`) |
| 🧫 **Primer Design** | Design forward and reverse primers from a DNA sequence and target region |
| 🔁 **Transcription Tool** | Convert a DNA sequence to an RNA sequence (T → U) |
| 🔤 **Translation Tool** | Translate an RNA sequence to a protein sequence using the standard genetic code |
| 🔗 **Complementary Sequence** | Generate the complementary strand of a DNA sequence |
| 🧬 **UniProt Lookup** | Retrieve protein information by UniProt ID (e.g., `P04637`) |
| 💉 **ChEMBL Lookup** | Retrieve compound information by ChEMBL ID (e.g., `CHEMBL12`) |

---

## 🛠️ Tech Stack

**Frontend**
- HTML5, JavaScript (Vanilla)
- [Tailwind CSS](https://tailwindcss.com/) (via CDN)
- Google Fonts (Inter)

**Backend**
- Python 3
- [Flask](https://flask.palletsprojects.com/)
- [Flask-CORS](https://flask-cors.readthedocs.io/)

---

## 📁 Project Structure

```
bioinformatics-tools/
│
├── app.py              # Flask backend with all API endpoints
├── index.html          # Frontend single-page application
└── README.md
```

---

## ⚙️ Installation & Setup

### 1. Clone the Repository

```bash
git clone https://github.com/your-username/bioinformatics-tools.git
cd bioinformatics-tools
```

### 2. Set Up the Backend

```bash
pip install flask flask-cors
```

### 3. Run the Flask Server

```bash
python app.py
```

The backend will start at `http://localhost:5000`.

### 4. Open the Frontend

Simply open `index.html` in your browser — no build step required.

> ⚠️ Make sure the Flask server is running before using the app.

---

## 🔌 API Endpoints

| Endpoint | Method | Request Body | Description |
|----------|--------|--------------|-------------|
| `/api/check-interaction` | POST | `{ "medication1": "aspirin", "medication2": "warfarin" }` | Check drug-drug interaction |
| `/api/blast` | POST | `{ "sequence": "ATGCGT..." }` | Simulate BLAST search |
| `/api/kegg` | POST | `{ "kegg_id": "hsa00010" }` | Get KEGG pathway info |
| `/api/primer-design` | POST | `{ "sequence": "...", "target_region": "10-50", "tm": "60" }` | Design primers |
| `/api/transcription` | POST | `{ "dna_sequence": "ATGC..." }` | Transcribe DNA to RNA |
| `/api/translation` | POST | `{ "rna_sequence": "AUGC..." }` | Translate RNA to protein |
| `/api/complementary` | POST | `{ "sequence": "ATGC..." }` | Get complementary DNA strand |
| `/api/uniprot` | POST | `{ "uniprot_id": "P04637" }` | Get UniProt protein info |
| `/api/chembl` | POST | `{ "chembl_id": "CHEMBL12" }` | Get ChEMBL compound info |

---

## 💊 Drug Interaction Severity Levels

| Severity | Meaning |
|----------|---------|
| 🔴 **Major** | Potentially life-threatening; avoid combination |
| 🟠 **Moderate** | May worsen condition or require close monitoring |
| 🟡 **Minor** | Limited clinical effect; usually manageable |
| 🟢 **None** | No significant interaction found |

### Example Interactions Supported

- `aspirin` + `warfarin` → **Major**: Increased bleeding risk
- `sildenafil` + `nitrates` → **Major**: Contraindicated (severe hypotension)
- `metronidazole` + `alcohol` → **Major**: Disulfiram-like reaction
- `ibuprofen` + `lisinopril` → **Moderate**: Increased kidney damage risk
- `amoxicillin` + `clarithromycin` → **Minor**: May increase stomach upset
- `metformin` + `atorvastatin` → **None**: No significant interaction

> **Note:** This is a simulated dataset for demonstration purposes only.

---

## 🗄️ Referenced Data Sources

- [DrugBank](https://www.drugbank.ca/)
- [PharmGKB](https://www.pharmgkb.org/)
- [SIDER](http://sideeffects.embl.de/)
- [ChEBI](https://www.ebi.ac.uk/chebi/)
- [FDA FAERS](https://www.fda.gov/drugs/questions-and-answers-fdas-adverse-event-reporting-system-faers)
- [DailyMed](https://dailymed.nlm.nih.gov/)
- [UniProt](https://www.uniprot.org/)
- [GenBank](https://www.ncbi.nlm.nih.gov/genbank/)
- [KEGG](https://www.genome.jp/kegg/)
- [ChEMBL](https://www.ebi.ac.uk/chembl/)

---

## 🧪 Usage Examples

### Check a Drug Interaction
1. Click **Medication Interaction** in the sidebar
2. Enter `aspirin` and `warfarin`
3. Click **Check Interaction**
4. Result: *"Major interaction: Increased risk of bleeding."*

### Transcribe a DNA Sequence
1. Click **Transcription Tool**
2. Enter `ATGCTTACG`
3. Result: `AUGCUUACG`

### Look Up a Protein
1. Click **UniProt**, enter `P04637`
2. Result: *Tumor protein p53, Homo sapiens*

---

## ⚠️ Disclaimer

This application is for **educational and demonstration purposes only**. It uses **simulated data** and should **not** be used for actual clinical or medical decision-making. Always consult a licensed healthcare professional or pharmacist.

---

## 📄 License

This project is open-source and available under the [MIT License](LICENSE).

---

## 🤝 Contributing

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/my-feature`)
3. Commit your changes (`git commit -m 'Add my feature'`)
4. Push to the branch (`git push origin feature/my-feature`)
5. Open a Pull Request

---

*Built with 🧬 for bioinformatics education and drug safety awareness.*
