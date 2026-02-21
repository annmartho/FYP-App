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


#with col_button:
if st.button("Go To Dashboard"):
    st.switch_page("dashboard.py")
# -------------------------------------------------
# PAGE CONFIG
# -------------------------------------------------
st.set_page_config(page_title="Structural Druggability", layout="wide")

# -------------------------------------------------
# CSS: UNIFIED DARK THEME
# -------------------------------------------------
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

    h1 { font-weight: 700 !important; font-size: 2.2rem !important; margin-bottom: 0.5rem !important; }
    h2 { font-weight: 700 !important; font-size: 1.8rem !important; margin-bottom: 0.5rem !important; }
    h3 { font-weight: 600 !important; font-size: 1.2rem !important; }

    /* Cards & Containers */
    .metric-card, .section-box {
        background-color: #1A1A1A;
        border: 1px solid #333;
        padding: 24px;
        border-radius: 12px;
        transition: all 0.3s ease;
        box-shadow: 0 4px 12px rgba(0, 0, 0, 0.3);
    }
    
    .metric-card:hover {
        border-color: #444;
        background-color: #252525;
        transform: translateY(-2px);
    }

    .section-label {
        font-size: 11px;
        text-transform: uppercase;
        letter-spacing: 1.5px;
        color: #999;
        margin-bottom: 8px;
        font-weight: 600;
    }

    .metric-value {
        font-size: 28px;
        font-weight: 700;
        color: #FFFFFF;
    }

    /* Status Badges */
    .status-container {
        display: flex;
        align-items: center;
        gap: 15px;
        margin-bottom: 2rem;
    }

    .status-pill {
        padding: 6px 16px;
        border-radius: 20px;
        font-size: 12px;
        font-weight: 700;
        text-transform: uppercase;
        letter-spacing: 1px;
    }

    /* Status Colors */
    .high-priority { background-color: #7A1F1F; color: #FFC1C1; }
    .mid-priority { background-color: #6B5E1A; color: #FFF3C1; }
    .low-priority { background-color: #333333; color: #AAAAAA; }

    /* Dividers */
    hr {
        margin: 2em 0;
        border: 0;
        border-top: 1px solid #333;
    }
    
    .stExpander {
        border: 1px solid #333 !important;
        background-color: #1A1A1A !important;
        border-radius: 12px !important;
    }
    
    code {
        background-color: #262626 !important;
        color: #999 !important;
        border: 1px solid #333 !important;
        border-radius: 4px !important;
        padding: 2px 6px !important;
    }

    .section-title {
        font-size: 20px;
        font-weight: 600;
        color: #4DA6FF;
        margin-bottom: 16px;
        display: block;
    }
</style>
""", unsafe_allow_html=True)



# -------------------------------------------------
# DATA LOADING 
# -------------------------------------------------
@st.cache_data
def load_data():
    try:
        pdb_map = pd.read_csv("data/best_pdbs_corrected_final.csv")
        struct_df = pd.read_csv("data/Structural_validation_with_PDB.csv")
        
        pdb_map.columns = [c.strip() for c in pdb_map.columns]
        struct_df.columns = [c.strip() for c in struct_df.columns]
        
        merged = struct_df.merge(pdb_map, left_on="Gene", right_on="gene", how="left")
        
        merged["Druggability_Evidence"] = merged["Druggability_Evidence"].fillna("Low Confidence").astype(str)
        merged["Known_Ligands"] = merged["Known_Ligands"].fillna("None").astype(str)
        merged["Binding_Affinity"] = merged["Binding_Affinity"].fillna("No Data").astype(str)
        merged["all_pdbs"] = merged["all_pdbs"].fillna("").astype(str)
        
        return merged
    except:
        return pd.DataFrame({
            "Gene": ["EGFR", "BRCA1"],
            "Druggability_Evidence": ["Very High", "Low Confidence"],
            "PDB_ID": ["1M17", "1JNX"],
            "uniprot": ["P00533", "P38398"],
            "Binding_Affinity": ["0.5 nM", "No Data"],
            "Known_Ligands": ["Erlotinib, Gefitinib", "None"],
            "all_pdbs": ["1M17;2J6M;3W2S", "1JNX;1LZE"]
        })

df = load_data()

# -------------------------------------------------
# MAIN DASHBOARD TOP (Search Section)
# -------------------------------------------------
st.markdown("<h1>Structural Druggability Validations</h1>", unsafe_allow_html=True)
st.markdown("<h2>Proteins with experimentally determined protein structure</h2>", unsafe_allow_html=True)

# Check if gene is selected from another page (via session state)
if "selected_gene" in st.session_state:
    gene = st.session_state["selected_gene"]
    st.markdown(f"<p style='color: #8ba3cb; font-size: 14px;'>Gene: <strong>{gene}</strong></p>", unsafe_allow_html=True)
else:
    # Fallback: show selectbox if no gene is selected
    gene = st.selectbox("Search / Select Target Gene:", sorted(df["Gene"].unique()))

st.markdown("<hr style='margin-top: 1em; margin-bottom: 2em;'>", unsafe_allow_html=True)

# -------------------------------------------------
# DASHBOARD CONTENT
# -------------------------------------------------
if gene:
    # perform a case-insensitive, stripped match and guard against empty results
    filtered = df[df["Gene"].astype(str).str.strip().str.upper() == str(gene).strip().upper()]

    if filtered.empty:
        sample_genes = ", ".join(sorted(df["Gene"].dropna().astype(str).unique()[:20]))
        st.error(f"No structural data found for '{gene}'.\nSamples: {sample_genes}...")
        st.stop()

    row = filtered.iloc[0]

    # Clean Header for the selected profile
    st.markdown(f"## {gene} Profile")
    
    # Status Badge
    evidence = str(row["Druggability_Evidence"]).lower()
    if "very high" in evidence or "high" in evidence:
        status_class = "high-priority"
    elif "medium" in evidence:
        status_class = "mid-priority"
    else:
        status_class = "low-priority"

    st.markdown(f"""
        <div class="status-container">
            <span class="status-pill {status_class}">{row['Druggability_Evidence']}</span>
        </div>
    """, unsafe_allow_html=True)

    # Metrics
    c1, c2, c3 = st.columns(3)
    
    metrics = [
        ("Primary PDB", row["PDB_ID"]),
        ("UniProt Accession", row["uniprot"]),
        ("Binding Affinity", row["Binding_Affinity"])
    ]

    for col, (label, val) in zip([c1, c2, c3], metrics):
        col.markdown(f"""
            <div class="metric-card">
                <div class="section-label">{label}</div>
                <div class="metric-value">{val}</div>
            </div>
        """, unsafe_allow_html=True)

    st.markdown("<hr>", unsafe_allow_html=True)

    # Lower Section
    col_l, col_r = st.columns([1, 1])

    with col_l:
        st.markdown("### Co-Crystallized Ligands")
        ligands = str(row["Known_Ligands"]).strip()
        if ligands.lower() not in ["none", "nan", ""]:
            for l in ligands.split(","):
                st.markdown(f"🔹 {l.strip()}")
        else:
            st.markdown("*No ligands mapped to this structure.*")

    with col_r:
        st.markdown("### Structural Library")
        with st.expander("Explore all PDB entries"):
            all_structs = str(row["all_pdbs"]).strip()
            if all_structs:
                for p in all_structs.split(";"):
                    st.code(p.strip())
            else:
                st.write("*No additional structures available.*")

    st.markdown("<hr>", unsafe_allow_html=True)

    # -------------------------------------------------
    # DOGSITE ANALYSIS SECTION
    # -------------------------------------------------
    st.markdown("<h2>🔬 DoGSite Pocket Analysis</h2>", unsafe_allow_html=True)
    
    pdb_id = row["PDB_ID"]
    
    # Only proceed if we have a valid PDB ID
    if pd.notna(pdb_id) and str(pdb_id).strip() and str(pdb_id).lower() != "nan":
        col_dogsite_left, col_dogsite_right = st.columns([1, 1])
        
        with col_dogsite_left:
            st.markdown(f"<p style='color: #8ba3cb; font-size: 14px;'>Analyzing PDB: <strong>{pdb_id}</strong></p>", unsafe_allow_html=True)
            
            run_analysis = st.button("▶ Run DoGSite Analysis")
        
        if run_analysis:
            with st.spinner("Submitting job to ProteinsPlus..."):
                API_BASE_URL = "https://proteins.plus/api/dogsite_rest"
                HEADERS = {
                    "Accept": "application/json",
                    "Content-Type": "application/json"
                }

                JOB_PARAMETERS = {
                    "dogsite": {
                        "pdbCode": pdb_id,
                        "analysisDetail": "0",
                        "bindingSitePredictionGranularity": "1",
                        "ligand": "",
                        "chain": ""
                    }
                }

                try:
                    post_response = requests.post(
                        API_BASE_URL,
                        headers=HEADERS,
                        data=json.dumps(JOB_PARAMETERS),
                        timeout=30
                    )
                    post_response.raise_for_status()
                    response_data = post_response.json()
                    job_location = response_data.get("location")
                except Exception as e:
                    st.error(f"Failed to submit job: {str(e)}")
                    st.stop()

            st.success("Job Submitted!")

            # Poll for Results
            result_data = None

            with st.spinner("Running pocket detection..."):
                attempts = 0
                max_attempts = 60  # 10 minutes with 10-second intervals
                
                while attempts < max_attempts:
                    try:
                        get_response = requests.get(job_location, headers=HEADERS, timeout=30)

                        if get_response.status_code == 200:
                            result_data = get_response.json()
                            break
                        elif get_response.status_code == 202:
                            time.sleep(10)
                            attempts += 1
                        else:
                            st.error(f"API error: {get_response.status_code}")
                            st.stop()
                    except Exception as e:
                        st.error(f"Error polling results: {str(e)}")
                        st.stop()
                
                if attempts >= max_attempts:
                    st.error("Analysis timed out. Please try again.")
                    st.stop()

            st.success("DoGSite Analysis Completed!")

            # Display Pocket Info
            pockets = result_data.get("pockets", [])

            if not pockets:
                st.warning("No pockets detected for this structure.")
            else:
                best_pocket_url = pockets[0]

                # Fetch AlphaFold or PDB structure
                try:
                    PDB_URL = f"https://alphafold.ebi.ac.uk/files/AF-{pdb_id}-F1-model_v6.pdb"
                    pdb_response = requests.get(PDB_URL, timeout=30)

                    if pdb_response.status_code != 200:
                        # Fallback to RCSB PDB
                        PDB_URL = f"https://files.rcsb.org/download/{pdb_id}.pdb"
                        pdb_response = requests.get(PDB_URL, timeout=30)

                    if pdb_response.status_code != 200:
                        st.error("Could not fetch structure data.")
                        st.stop()

                    pdb_data = pdb_response.text
                except Exception as e:
                    st.error(f"Error fetching structure: {str(e)}")
                    st.stop()

                # Download + Decompress CCP4
                try:
                    pocket_response = requests.get(best_pocket_url, timeout=30)
                    compressed_data = io.BytesIO(pocket_response.content)

                    with gzip.open(compressed_data, "rb") as f:
                        pocket_binary = f.read()

                    pocket_b64 = base64.b64encode(pocket_binary).decode("utf-8")
                except Exception as e:
                    st.error(f"Error processing pocket data: {str(e)}")
                    st.stop()

                # 3D Visualization
                st.subheader("3D Pocket Visualization")

                view = py3Dmol.view(width=800, height=600)

                # Protein
                view.addModel(pdb_data, "pdb")
                view.setStyle({"model": 0}, {"cartoon": {"color": "lightgray"}})

                # Pocket surface
                view.addVolumetricData(
                    pocket_b64,
                    "ccp4",
                    {"isoval": 0.5, "color": "red", "opacity": 0.8}
                )

                view.zoomTo()

                showmol(view, height=600, width=800)

                st.caption("Predicted Pocket (DoGSite Scorer) – Red surface indicates druggable cavity region.")
    else:
        st.info(f"ℹ️ No valid PDB ID available for {gene}. Cannot run DoGSite analysis.")