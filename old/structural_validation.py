#run pip3 install py3Dmol stmol 
#python3 -m pip install ipython-genutils ipywidgets
import streamlit as st
import requests
import json
import time
import base64
import gzip
import io
import py3Dmol
from stmol import showmol

st.set_page_config(page_title="Structural Validation")

st.title("🧬 Structural Validation – DoGSite Scorer")

if "selected_gene" not in st.session_state:
    st.warning("No gene selected. Please search from dashboard first.")
    st.stop()

gene = st.session_state["selected_gene"]
st.subheader(f"Gene: {gene}")

# -----------------------------
# User Input: PDB ID
# -----------------------------

pdb_id = st.text_input("Enter UniProt / AlphaFold ID", value="Q9Y216")

if st.button("Run DoGSite Analysis"):

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

        post_response = requests.post(
            API_BASE_URL,
            headers=HEADERS,
            data=json.dumps(JOB_PARAMETERS)
        )

        post_response.raise_for_status()
        response_data = post_response.json()
        job_location = response_data.get("location")

    st.success("Job Submitted!")

    # -----------------------------
    # Poll for Results
    # -----------------------------

    result_data = None

    with st.spinner("Running pocket detection..."):

        while True:
            get_response = requests.get(job_location, headers=HEADERS)

            if get_response.status_code == 200:
                result_data = get_response.json()
                break

            elif get_response.status_code == 202:
                time.sleep(10)
            else:
                st.error("API error while polling.")
                st.stop()

    st.success("DoGSite Analysis Completed!")

    # -----------------------------
    # Display Pocket Info
    # -----------------------------

    pockets = result_data.get("pockets", [])

    if not pockets:
        st.warning("No pockets detected.")
        st.stop()


    # Use first pocket for visualization
    best_pocket_url = pockets[0]

    # -----------------------------
    # Fetch AlphaFold PDB
    # -----------------------------

    PDB_URL = f"https://alphafold.ebi.ac.uk/files/AF-{pdb_id}-F1-model_v6.pdb"
    pdb_response = requests.get(PDB_URL)

    if pdb_response.status_code != 200:
        st.error("Could not fetch AlphaFold structure.")
        st.stop()

    pdb_data = pdb_response.text

    # -----------------------------
    # Download + Decompress CCP4
    # -----------------------------

    pocket_response = requests.get(best_pocket_url)
    compressed_data = io.BytesIO(pocket_response.content)

    with gzip.open(compressed_data, "rb") as f:
        pocket_binary = f.read()

    pocket_b64 = base64.b64encode(pocket_binary).decode("utf-8")

    # -----------------------------
    # 3D Visualization
    # -----------------------------

    st.subheader("🧬 3D Pocket Visualization")

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
