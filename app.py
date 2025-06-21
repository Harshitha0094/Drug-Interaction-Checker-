# Import necessary libraries for Flask and CORS

from flask import Flask, request, jsonify
from flask_cors import CORS # type: ignore
import time
import random

# Initialize the Flask application
app = Flask(__name__)
# Enable CORS for all routes, allowing your frontend HTML to make requests
CORS(app)

# --- Backend API Endpoints for Bioinformatics Tools ---

@app.route('/api/check-interaction', methods=['POST'])
def check_interaction():
    """
    API endpoint to check medication interactions.
    Expects two medication names in the request JSON.
    """
    data = request.get_json()
    medication1 = data.get('medication1', '').strip().lower()
    medication2 = data.get('medication2', '').strip().lower()

    # Simulate a delay for demonstration purposes
    time.sleep(1.5)

    if not medication1 or not medication2:
        return jsonify({"result": "Please enter both medications.", "severity": "None"}), 400

    # Define simulated interactions. In a real app, this would query a database.
    interactions = {
        "aspirin-warfarin": {"result": "Major interaction: Increased risk of bleeding.", "severity": "Major"},
        "ibuprofen-lisinopril": {"result": "Moderate interaction: Increased risk of kidney damage.", "severity": "Moderate"},
        "amoxicillin-clarithromycin": {"result": "Minor interaction: May increase stomach upset.", "severity": "Minor"},
        "metformin-atorvastatin": {"result": "No significant interaction.", "severity": "None"},
        "losartan-ibuprofen": {"result": "Moderate interaction: Increased risk of kidney problems.", "severity": "Moderate"},
        "citalopram-tramadol": {"result": "Moderate interaction: Increased risk of serotonin syndrome.", "severity": "Moderate"},
        "gabapentin-morphine": {"result": "Moderate interaction: Increased risk of respiratory depression.", "severity": "Moderate"},
        "prednisone-insulin": {"result": "Moderate interaction: Prednisone can increase blood sugar levels, affecting insulin dosage.", "severity": "Moderate"},
        "digoxin-verapamil": {"result": "Moderate interaction: Verapamil can increase digoxin levels, leading to toxicity.", "severity": "Moderate"},
        "lithium-nsaids": {"result": "Moderate interaction: NSAIDs can increase lithium levels, increasing the risk of toxicity.", "severity": "Moderate"},
        "theophylline-caffeine": {"result": "Moderate interaction: Increased risk of CNS stimulation.", "severity": "Moderate"},
        "phenytoin-carbamazepine": {"result": "Variable: Can increase or decrease levels of either drug.", "severity": "Minor"},
        "isoniazid-rifampin": {"result": "Moderate interaction: Increased risk of liver damage.", "severity": "Moderate"},
        "enalapril-spironolactone": {"result": "Moderate interaction: Increased risk of hyperkalemia.", "severity": "Moderate"},
        "metoprolol-clonidine": {"result": "Moderate interaction: Risk of excessive bradycardia.", "severity": "Moderate"},
        "nifedipine-ranitidine": {"result": "Minor interaction: Ranitidine may increase nifedipine levels slightly.", "severity": "Minor"},
        "pioglitazone-gemfibrozil": {"result": "Moderate interaction: Increased risk of hypoglycemia and edema.", "severity": "Moderate"},
        "sildenafil-nitrates": {"result": "Contraindicated: Severe hypotension and cardiovascular events.", "severity": "Major"},
        "tamsulosin-warfarin": {"result": "Moderate interaction: Increased risk of bleeding.", "severity": "Moderate"},
        "venlafaxine-maois": {"result": "Contraindicated: Risk of serotonin syndrome.", "severity": "Major"},
        "zolpidem-alcohol": {"result": "Major interaction: Increased risk of CNS depression and respiratory depression.", "severity": "Major"},
        "atorvastatin-clarithromycin": {"result": "Moderate interaction: Increased risk of myopathy.", "severity": "Moderate"},
        "benazepril-hydrochlorothiazide": {"result": "Additive effect: Increased risk of hypotension.", "severity": "Minor"},
        "carvedilol-fluoxetine": {"result": "Moderate interaction: Increased risk of bradycardia.", "severity": "Moderate"},
        "doxycycline-antacids": {"result": "Moderate interaction: Antacids can decrease doxycycline absorption.", "severity": "Moderate"},
        "erythromycin-simvastatin": {"result": "Moderate interaction: Increased risk of myopathy.", "severity": "Moderate"},
        "furosemide-digoxin": {"result": "Moderate interaction: Increased risk of digoxin toxicity if potassium is depleted.", "severity": "Moderate"},
        "haloperidol-lithium": {"result": "Moderate interaction: Increased risk of neurotoxicity.", "severity": "Moderate"},
        "irbesartan-potassium supplements": {"result": "Moderate interaction: Increased risk of hyperkalemia.", "severity": "Moderate"},
        "ketoconazole-warfarin": {"result": "Moderate interaction: Increased risk of bleeding.", "severity": "Moderate"},
        "lamotrigine-valproate": {"result": "Moderate interaction: Increased risk of skin rash.", "severity": "Moderate"},
        "metronidazole-alcohol": {"result": "Disulfiram-like reaction: Nausea, vomiting, headache.", "severity": "Major"},
        "naltrexone-opioids": {"result": "Contraindicated: Opioid withdrawal.", "severity": "Major"},
        "olanzapine-benzodiazepines": {"result": "Moderate interaction: Increased risk of sedation and respiratory depression.", "severity": "Moderate"},
        "paroxetine-tamoxifen": {"result": "Moderate interaction: Paroxetine may reduce the effectiveness of tamoxifen.", "severity": "Moderate"},
        "quinidine-digoxin": {"result": "Moderate interaction: Quinidine increases digoxin levels.", "severity": "Moderate"},
        "rosiglitazone-insulin": {"result": "Moderate interaction: Increased risk of hypoglycemia.", "severity": "Moderate"},
        "spironolactone-trimethoprim": {"result": "Moderate interaction: Increased risk of hyperkalemia.", "severity": "Moderate"},
        "terazosin-sildenafil": {"result": "Moderate interaction: Increased risk of hypotension.", "severity": "Moderate"},
        "uricosurics-probenecid": {"result": "Therapeutic duplication: May not provide additional benefit.", "severity": "Minor"},
        "verapamil-beta-blockers": {"result": "Moderate interaction: Increased risk of bradycardia and heart block.", "severity": "Moderate"},
        "warfarin-aspirin": {"result": "Major interaction: Increased risk of bleeding.", "severity": "Major"},
        "xylometazoline-maois": {"result": "Moderate interaction: Risk of hypertensive crisis.", "severity": "Moderate"},
        "yohimbine-antidepressants": {"result": "Moderate interaction: May cause anxiety and increased blood pressure.", "severity": "Moderate"},
        "zafirlukast-warfarin": {"result": "Moderate interaction: Zafirlukast may increase warfarin levels.", "severity": "Moderate"},
        "alendronate-calcium supplements": {"result": "Moderate interaction: Calcium can decrease alendronate absorption.", "severity": "Moderate"},
        "bromocriptine-metoclopramide": {"result": "Moderate interaction: Metoclopramide can decrease the effectiveness of bromocriptine.", "severity": "Moderate"},
        "cyclosporine-ketoconazole": {"result": "Moderate interaction: Ketoconazole can increase cyclosporine levels.", "severity": "Moderate"},
        "diazepam-cimetidine": {"result": "Moderate interaction: Cimetidine can increase diazepam levels.", "severity": "Moderate"},
        "epinephrine-beta-blockers": {"result": "Moderate interaction: Beta-blockers can block the bronchodilating effects of epinephrine.", "severity": "Moderate"},
        "fluconazole-phenytoin": {"result": "Moderate interaction: Fluconazole can increase phenytoin levels.", "severity": "Moderate"},
        "gentamicin-vancomycin": {"result": "Additive nephrotoxicity: Increased risk of kidney damage.", "severity": "Moderate"},
        "hydralazine-nsaids": {"result": "Moderate interaction: NSAIDs can decrease the antihypertensive effect of hydralazine.", "severity": "Moderate"},
        "isoniazid-phenytoin": {"result": "Moderate interaction: Isoniazid can increase phenytoin levels.", "severity": "Moderate"},
        "ketorolac-anticoagulants": {"result": "Increased bleeding risk: Increased risk of bleeding.", "severity": "Major"},
        "levodopa-antipsychotics": {"result": "Antipsychotics can block dopamine receptors, reducing levodopa's effectiveness.", "severity": "Moderate"},
        "meperidine-maois": {"result": "Risk of serotonin syndrome: Risk of serotonin syndrome.", "severity": "Major"},
        "nifedipine-digoxin": {"result": "Moderate interaction: Nifedipine can increase digoxin levels.", "severity": "Moderate"},
        "omeprazole-diazepam": {"result": "Moderate interaction: Omeprazole can increase diazepam levels.", "severity": "Moderate"},
        "probenecid-penicillins": {"result": "Therapeutic benefit: Probenecid can increase penicillin levels.", "severity": "None"},
        "quinolone antibiotics-theophylline": {"result": "Increased theophylline levels: Increased risk of theophylline toxicity.", "severity": "Moderate"},
        "rifampin-oral contraceptives": {"result": "Decreased contraceptive effectiveness: Decreased contraceptive effectiveness.", "severity": "Moderate"},
        "sulfonamides-sulfonylureas": {"result": "Increased hypoglycemic effect: Increased risk of hypoglycemia.", "severity": "Moderate"},
        "tetracyclines-dairy products": {"result": "Decreased antibiotic absorption: Decreased antibiotic absorption.", "severity": "Minor"},
        "valproic acid-phenobarbital": {"result": "Increased phenobarbital levels: Increased risk of sedation.", "severity": "Moderate"},
        "zidovudine-probenecid": {"result": "Increased zidovudine levels: Increased risk of zidovudine toxicity.", "severity": "Moderate"}
    }

    # Normalize input for lookup
    key1 = f"{medication1}-{medication2}"
    key2 = f"{medication2}-{medication1}"

    interaction_result = interactions.get(key1) or interactions.get(key2)

    if interaction_result:
        return jsonify(interaction_result)
    else:
        return jsonify({"result": f"No significant interaction found between {medication1} and {medication2}.", "severity": "None"})

@app.route('/api/blast', methods=['POST'])
def run_blast():
    """
    API endpoint for the BLAST tool.
    Expects a sequence in the request JSON.
    """
    data = request.get_json()
    sequence = data.get('sequence', '').strip()

    time.sleep(2) # Simulate network delay and processing

    if not sequence:
        return jsonify({"result": "Please enter a sequence."}), 400
    if len(sequence) < 10:
        return jsonify({"result": "Sequence is too short for a meaningful BLAST."}), 400

    hits = random.randint(1, 5)
    return jsonify({"result": f"BLAST found {hits} significant hits for your sequence."})

@app.route('/api/kegg', methods=['POST'])
def get_kegg_info():
    """
    API endpoint for the KEGG tool.
    Expects a KEGG ID in the request JSON.
    """
    data = request.get_json()
    kegg_id = data.get('kegg_id', '').strip()

    time.sleep(1.5)

    if not kegg_id:
        return jsonify({"result": "Please enter a KEGG ID."}), 400
    if not kegg_id.startswith('hsa'):
        return jsonify({"result": "Invalid KEGG ID. Example: hsa00010"}), 400

    # More sophisticated KEGG simulation
    kegg_data = {
        'hsa00010': 'KEGG Pathway: Glycolysis/Gluconeogenesis - Central energy metabolism.',
        'hsa00020': 'KEGG Pathway: Citrate cycle (TCA cycle) - Aerobic respiration.',
        'hsa04110': 'KEGG Pathway: Cell cycle - Regulation of cell division.',
        'hsa05200': 'KEGG Pathway: Pathways in cancer - Various molecular pathways associated with cancer.'
    }
    result = kegg_data.get(kegg_id, f"KEGG ID: {kegg_id} - Information not found (Simulated).")
    return jsonify({"result": result})

@app.route('/api/primer-design', methods=['POST'])
def design_primers():
    """
    API endpoint for the Primer Design tool.
    Expects sequence, target_region, and tm in the request JSON.
    """
    data = request.get_json()
    sequence = data.get('sequence', '').strip().upper()
    target_region = data.get('target_region', '').strip()
    tm = data.get('tm', '').strip()

    time.sleep(2)

    if not sequence or not target_region or not tm:
        return jsonify({"result": "Please enter sequence, region, and Tm."}), 400

    try:
        start_str, end_str = target_region.split('-')
        start = int(start_str) - 1 # Adjust to 0-indexed
        end = int(end_str) # Keep as 1-indexed for end slice, Python slice is exclusive
        melting_temp = float(tm)
    except (ValueError, IndexError):
        return jsonify({"result": "Invalid target region or melting temperature."}), 400

    if start < 0 or end > len(sequence) or start >= end:
        return jsonify({"result": "Target region out of bounds or invalid."}), 400

    # Basic primer design logic (highly simplified)
    primer_length = 20
    if (start + primer_length) > len(sequence) or (end - primer_length) < 0:
        return jsonify({"result": "Sequence too short for primers in specified region."}), 400

    forward_primer = sequence[start : start + primer_length]
    # To get a reverse primer for amplification, we need the reverse complement of the end section
    # This is a very simplified example and doesn't account for optimal primer placement/properties
    reverse_complement_seq = ''.join([{'A':'T', 'T':'A', 'C':'G', 'G':'C'}.get(base, '?') for base in sequence[::-1]])
    reverse_primer = reverse_complement_seq[len(sequence) - end : len(sequence) - end + primer_length]


    return jsonify({
        "result": "Primers designed successfully (Simulated).",
        "forward_primer": forward_primer,
        "reverse_primer": reverse_primer,
        "tm_used": f"{melting_temp}°C"
    })

@app.route('/api/transcription', methods=['POST'])
def transcribe_sequence():
    """
    API endpoint for the Transcription tool.
    Expects a DNA sequence in the request JSON.
    """
    data = request.get_json()
    dna_sequence = data.get('dna_sequence', '').strip().upper()

    time.sleep(1)

    if not dna_sequence:
        return jsonify({"result": "Please enter a DNA sequence."}), 400

    # Transcription: T -> U
    rna_sequence = dna_sequence.replace('T', 'U')
    return jsonify({"result": "Transcription successful.", "rna_sequence": rna_sequence})

@app.route('/api/translation', methods=['POST'])
def translate_sequence():
    """
    API endpoint for the Translation tool.
    Expects an RNA sequence in the request JSON.
    """
    data = request.get_json()
    rna_sequence = data.get('rna_sequence', '').strip().upper()

    time.sleep(1.5)

    if not rna_sequence:
        return jsonify({"result": "Please enter an RNA sequence."}), 400

    # Codon table for translation (standard genetic code)
    codon_table = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M', # AUG is start
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',

        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',

        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*', # Stop codons
        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',

        'UGU': 'C', 'UGC': 'C', 'UGA': '*', # Stop codon
        'UGG': 'W',
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }

    protein_sequence = []
    # Start translation from the first 'AUG' if found, otherwise from beginning
    start_index = rna_sequence.find('AUG')
    if start_index == -1: # No start codon, translate from beginning
        start_index = 0

    for i in range(start_index, len(rna_sequence) - 2, 3):
        codon = rna_sequence[i:i+3]
        amino_acid = codon_table.get(codon, '?')
        if amino_acid == '*':
            break # Stop translation at stop codon
        protein_sequence.append(amino_acid)

    if not protein_sequence:
        return jsonify({"result": "Could not translate. Check RNA sequence for valid codons or start codon.", "protein_sequence": ""})

    return jsonify({"result": "Translation successful.", "protein_sequence": "".join(protein_sequence)})

@app.route('/api/complementary', methods=['POST'])
def get_complement():
    """
    API endpoint for the Complementary tool.
    Expects a DNA sequence in the request JSON.
    """
    data = request.get_json()
    sequence = data.get('sequence', '').strip().upper()

    time.sleep(1)

    if not sequence:
        return jsonify({"result": "Please enter a DNA sequence."}), 400

    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    complement = ''.join([complement_map.get(base, '?') for base in sequence])
    return jsonify({"result": "Complementary sequence generated.", "complement_sequence": complement})

@app.route('/api/uniprot', methods=['POST'])
def get_uniprot_info():
    """
    API endpoint for the UniProt tool.
    Expects a UniProt ID in the request JSON.
    """
    data = request.get_json()
    uniprot_id = data.get('uniprot_id', '').strip().upper()

    time.sleep(1.5)

    if not uniprot_id:
        return jsonify({"result": "Please enter a UniProt ID."}), 400

    uniprot_data = {
        'P04637': 'UniProt ID: P04637, Protein Name: Tumor protein p53, Species: Homo sapiens, Function: Tumor suppressor.',
        'P00533': 'UniProt ID: P00533, Protein Name: Epidermal growth factor receptor (EGFR), Species: Homo sapiens, Function: Cell growth, differentiation, and survival.',
        'Q13153': 'UniProt ID: Q13153, Protein Name: PTEN, Species: Homo sapiens, Function: Tumor suppressor, cell signaling.',
    }
    result = uniprot_data.get(uniprot_id, f"UniProt ID: {uniprot_id} - Information not found (Simulated).")
    return jsonify({"result": result})

@app.route('/api/chembl', methods=['POST'])
def get_chembl_info():
    """
    API endpoint for the ChEMBL tool.
    Expects a ChEMBL ID in the request JSON.
    """
    data = request.get_json()
    chembl_id = data.get('chembl_id', '').strip().upper()

    time.sleep(1.5)

    if not chembl_id:
        return jsonify({"result": "Please enter a ChEMBL ID."}), 400

    chembl_data = {
        'CHEMBL12': 'ChEMBL ID: CHEMBL12, Compound Name: Aspirin, Molecular Weight: 180.16 g/mol, Type: Small molecule, ATC: N02BA01.',
        'CHEMBL1766': 'ChEMBL ID: CHEMBL1766, Compound Name: Metformin, Molecular Weight: 129.16 g/mol, Type: Small molecule, ATC: A10BA02.',
        'CHEMBL108': 'ChEMBL ID: CHEMBL108, Compound Name: Insulin, Molecular Weight: ~5808 Da, Type: Protein, ATC: A10AB01.',
    }
    result = chembl_data.get(chembl_id, f"ChEMBL ID: {chembl_id} - Information not found (Simulated).")
    return jsonify({"result": result})

# Run the Flask app
if __name__ == '__main__':
    # You can change the port if 5000 is in use
    app.run(debug=True, port=5000)
