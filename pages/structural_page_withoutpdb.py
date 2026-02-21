import streamlit as st
import pandas as pd
import requests
import json
import time
import base64
import gzip
import io
import py3Dmol
from stmol import showmol

# -------------------------------------------------
# PAGE CONFIG & CSS: UNIFIED DARK THEME
# -------------------------------------------------
st.set_page_config(page_title="Structural Summary", layout="wide")

st.markdown("""
<style>
    /* Unified Dark Theme */
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;600;700&display=swap');
    
    .block-container {
        padding-top: 2.5rem !important; 
        padding-bottom: 1rem !important;
        font-family: 'Inter', sans-serif !important;
    }
    
    header {
        visibility: hidden;
    }

    /* Main Background */
    .stApp {
        background-color: #0F0F0F;
        color: #EAEAEA;
    }
    
    .stMarkdown, .stMarkdown p, .stMarkdown span, h1, h2, h3 {
        background-color: transparent !important;
        color: #EAEAEA;
        font-family: 'Inter', sans-serif !important;
    }

    h1 { font-weight: 700 !important; font-size: 2.2rem !important; margin-bottom: 0.5rem !important; color: #FFFFFF !important; }
    h2 { font-weight: 700 !important; font-size: 1.8rem !important; margin-bottom: 0.5rem !important; color: #4DA6FF !important; }

    /* Cards & Containers */
    .section-box {
        background-color: #1A1A1A;
        border: 1px solid #333;
        padding: 24px;
        border-radius: 12px;
        margin-bottom: 20px;
    }

    .section-title {
        font-size: 20px;
        font-weight: 600;
        color: #4DA6FF;
        display: block;
    }

    /* Button Styling */
    .stButton>button {
        background-color: #1A1A1A;
        color: #EAEAEA;
        border: 1px solid #333;
        border-radius: 8px;
        padding: 0.5rem 1rem;
        font-weight: 600;
        transition: all 0.3s ease;
    }
    .stButton>button:hover {
        border-color: #4DA6FF;
        color: #4DA6FF;
    }

    code {
        background-color: #262626 !important;
        color: #4DA6FF !important;
        padding: 2px 6px !important;
        border-radius: 4px !important;
    }
</style>
""", unsafe_allow_html=True)

# -----------------------------
# Navigation & Session Check
# -----------------------------
if st.button("⬅ Go To Dashboard"):
    st.switch_page("dashboard.py")

if "selected_gene" not in st.session_state:
    st.warning("No gene selected. Please search from dashboard first.")
    st.stop()

gene = st.session_state["selected_gene"]
st.title(f"Structural Homology: {gene}")

# -----------------------------
# Data Loading
# -----------------------------
@st.cache_data
def load_structural_data():
    df = pd.read_csv("data/alphafill_report.csv")
    return df

df_struct = load_structural_data()
gene_row = df_struct[df_struct["Gene"].str.upper() == gene.upper()]

if gene_row.empty:
    st.error("No structural data available for this gene.")
    st.stop()

row = gene_row.iloc[0]

# -----------------------------
# SECTION 1: STRUCTURAL DETAILS
# -----------------------------
# Opening the container
st.markdown('<div class="section-box"><div class="section-title">Structural Details (Alphafold and Alphafill)</div>', unsafe_allow_html=True)

# Putting columns inside the container
col1, col2 = st.columns(2)
with col1:
    st.markdown(f"**UniProt ID:** `{row['UniProt']}`")
    st.markdown(f"**Template PDB:** `{row['Template_PDB']}`")
    st.markdown(f"**Sequence Identity:** {row['Identity']}")
with col2:
    st.markdown(f"**Best Ligand:** {row['Best_Ligand']}")
    st.markdown(f"**RMSD:** {row['RMSD']}")
    st.markdown(f"**Clash Count:** {row['Clash_Count']}")

# Closing the container
st.markdown('</div>', unsafe_allow_html=True)

# -----------------------------
# SECTION 2: QUALITY ASSESSMENT
# -----------------------------
st.markdown('<div class="section-box"><div class="section-title">Homologous template identification</div>', unsafe_allow_html=True)
if row["Best_Ligand"] == "No Structure":
    st.warning("No experimentally resolved structure available for this gene.")
else:
    identity_val = float(row["Identity"].replace("%",""))
    if identity_val >= 40:
        st.success("✅ High structural similarity template identified.")
    elif identity_val >= 30:
        st.info("ℹ️ Moderate structural similarity template identified.")
    else:
        st.warning("⚠️ Low structural similarity – interpret with caution.")
st.markdown('</div>', unsafe_allow_html=True)

# -----------------------------
# SECTION 3: DOGSITE ANALYSIS
# -----------------------------
st.markdown('<div class="section-box"><div class="section-title">DoGSite Pocket Analysis</div>', unsafe_allow_html=True)

pdb_id = row['Template_PDB']

if pd.notna(pdb_id) and str(pdb_id).strip() and str(pdb_id).lower() != "nan":
    st.markdown(f"<p style='color: #999;'>Analyzing PDB Template: <strong>{pdb_id}</strong></p>", unsafe_allow_html=True)
    
    if st.button("▶ Run DoGSite Analysis"):
        with st.spinner("Communicating with ProteinsPlus API..."):
            API_BASE_URL = "https://proteins.plus/api/dogsite_rest"
            HEADERS = {"Accept": "application/json", "Content-Type": "application/json"}
            JOB_PARAMETERS = {
                "dogsite": {
                    "pdbCode": pdb_id,
                    "analysisDetail": "0",
                    "bindingSitePredictionGranularity": "1",
                    "ligand": "", "chain": ""
                }
            }

            try:
                post_response = requests.post(API_BASE_URL, headers=HEADERS, data=json.dumps(JOB_PARAMETERS), timeout=30)
                post_response.raise_for_status()
                job_location = post_response.json().get("location")
                
                # Polling
                result_data = None
                attempts = 0
                while attempts < 60:
                    get_response = requests.get(job_location, headers=HEADERS, timeout=30)
                    if get_response.status_code == 200:
                        result_data = get_response.json()
                        break
                    elif get_response.status_code == 202:
                        time.sleep(10)
                        attempts += 1
                    else:
                        st.error("API error during polling.")
                        break
                
                if result_data:
                    pockets = result_data.get("pockets", [])
                    if not pockets:
                        st.warning("No pockets detected.")
                    else:
                        best_pocket_url = pockets[0]
                        
                        # Fetch PDB Structure
                        PDB_URL = f"https://alphafold.ebi.ac.uk/files/AF-{pdb_id}-F1-model_v6.pdb"
                        pdb_res = requests.get(PDB_URL, timeout=30)
                        if pdb_res.status_code != 200:
                            pdb_res = requests.get(f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=30)
                        
                        # Fetch and decompress pocket
                        pocket_res = requests.get(best_pocket_url, timeout=30)
                        with gzip.open(io.BytesIO(pocket_res.content), "rb") as f:
                            pocket_binary = f.read()
                        pocket_b64 = base64.b64encode(pocket_binary).decode("utf-8")

                        # Visualization inside the card
                        st.markdown("---")
                        st.markdown('<div class="section-title">🧬 3D Pocket Visualization</div>', unsafe_allow_html=True)
                        view = py3Dmol.view(width=800, height=600)
                        view.addModel(pdb_res.text, "pdb")
                        view.setStyle({"model": 0}, {"cartoon": {"color": "lightgray"}})
                        view.addVolumetricData(pocket_b64, "ccp4", {"isoval": 0.5, "color": "red", "opacity": 0.8})
                        view.zoomTo()
                        showmol(view, height=600, width=800)
                        st.caption("Red surface indicates predicted druggable cavity.")

            except Exception as e:
                st.error(f"Analysis failed: {str(e)}")
else:
    st.info(f"ℹ️ No valid PDB ID available for {gene}.")

st.markdown('</div>', unsafe_allow_html=True)