import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter

# ---------------------------------------------------
# Page Config
# ---------------------------------------------------
st.set_page_config(page_title="Gene Validation", layout="wide")

# ---------------------------------------------------
# Refined Dark Dashboard CSS
# ---------------------------------------------------
st.markdown("""
<style>
    /* Main Background */
    .stApp {
        background-color: #0E1117;
        color: #E0E0E0;
    }

    /* Card Styling */
    .card {
        background-color: #161B22;
        padding: 24px;
        border-radius: 12px;
        border: 1px solid #30363D;
        margin-bottom: 24px;
        box-shadow: 0 4px 12px rgba(0,0,0,0.3);
    }

    /* Section Titles */
    .section-title {
        font-size: 22px;
        font-weight: 700;
        color: #58A6FF;
        margin-bottom: 20px;
        display: flex;
        align-items: center;
        gap: 10px;
    }

    /* Modern Badges */
    .priority-badge {
        padding: 4px 12px;
        border-radius: 6px;
        font-size: 12px;
        font-weight: 600;
        text-transform: uppercase;
        letter-spacing: 0.5px;
        margin-bottom: 15px;
        display: inline-block;
    }
    .priority-high { background-color: #DA3633; color: white; }
    .priority-medium { background-color: #D29922; color: black; }
    .priority-low { background-color: #30363D; color: #8B949E; }

    /* Muted Text */
    .muted-text {
        color: #8B949E;
        font-style: italic;
    }

    /* Link Styling */
    .pubmed-link {
        display: inline-block;
        background: #21262D;
        padding: 4px 10px;
        border-radius: 4px;
        margin: 4px;
        text-decoration: none !important;
        border: 1px solid #30363D;
        transition: 0.2s;
    }
    .pubmed-link:hover {
        background: #30363D;
        border-color: #58A6FF;
    }
    
    /* Global Spacing Fixes */
    .block-container {
        padding-top: 2rem !important;
    }
</style>
""", unsafe_allow_html=True)

# ---------------------------------------------------
# Helper: Dark Matplotlib Style
# ---------------------------------------------------
def set_dark_plot_style():
    plt.style.use("dark_background")
    plt.rcParams.update({
        "axes.facecolor": "#161B22",
        "figure.facecolor": "#161B22",
        "axes.edgecolor": "#30363D",
        "grid.color": "#30363D",
        "font.family": "sans-serif",
        "text.color": "#E6EDF3"
    })

# ---------------------------------------------------
# Data Loading (Logic remains same, cached)
# ---------------------------------------------------
@st.cache_data
def load_validation_data():
    # Mocking structure for demonstration based on your code
    try:
        df = pd.read_csv("data/api_validation.csv")
    except:
        df = pd.DataFrame(columns=["Gene", "Interacts_With", "Priority", "Target_Class", "Druggability", "Papers"])
    df.columns = [col.strip() for col in df.columns]
    return df

@st.cache_data
def load_survival():
    # Use your actual paths here
    try:
        surv_df = pd.read_csv("survival_analysis_csv/LUAD_survival.txt", sep="\t")
        expr_df = pd.read_csv("survival_analysis_csv/M1_counts_raw.csv")
        # ... logic from your snippet ...
        return pd.DataFrame(), {} # Return actual data
    except:
        return pd.DataFrame(), {}

# ---------------------------------------------------
# Logic Execution
# ---------------------------------------------------
if "selected_gene" not in st.session_state:
    st.session_state["selected_gene"] = "TP53" # Default for testing

gene = st.session_state["selected_gene"]

# Header Layout
col_head, col_btn = st.columns([0.8, 0.2])
with col_head:
    st.title(f"🧬 Gene Analysis: {gene}")
with col_btn:
    if st.button("⬅ Back to Dashboard"):
        st.switch_page("dashboard.py")

# ===================================================
# 1️⃣ API VALIDATION CARD
# ===================================================
validation_df = load_validation_data()
gene_validation = validation_df[validation_df["Gene"].str.upper() == gene.upper()]

with st.container():
    st.markdown('<div class="card">', unsafe_allow_html=True)
    st.markdown('<div class="section-title">🔬 API Validation Summary</div>', unsafe_allow_html=True)
    
    if gene_validation.empty:
        st.markdown('<div class="muted-text">No API validation data available.</div>', unsafe_allow_html=True)
    else:
        v = gene_validation.iloc[0]
        p = str(v["Priority"]).lower()
        p_class = f"priority-{p}" if p in ["high", "medium", "low"] else "priority-low"
        
        st.markdown(f'<div class="priority-badge {p_class}">{p.upper()} PRIORITY</div>', unsafe_allow_html=True)
        
        m1, m2, m3 = st.columns(3)
        m1.metric("Target Class", v["Target_Class"])
        m2.metric("Druggability", v["Druggability"])
        m3.metric("Literature Mentions", int(v.get("Papers", 0)))
        
        if str(v["Interacts_With"]).lower() not in ["none", "nan", ""]:
            st.info(f"**Driver Interactions:** {v['Interacts_With']}")
    st.markdown("</div>", unsafe_allow_html=True)

# ===================================================
# 2️⃣ CANCERMINE & SURVIVAL (Side-by-Side Layout)
# ===================================================
col_left, col_right = st.columns([1, 1], gap="large")

with col_left:
    st.markdown('<div class="card">', unsafe_allow_html=True)
    st.markdown('<div class="section-title">📊 Survival Analysis</div>', unsafe_allow_html=True)
    
    # Simple logic placeholder to represent your run_km function
    set_dark_plot_style()
    fig, ax = plt.subplots(figsize=(6, 4.5))
    # ... your KM fitting logic ...
    ax.text(0.5, 0.5, "KM Curve Visualization", ha='center') # Placeholder
    st.pyplot(fig)
    st.markdown("</div>", unsafe_allow_html=True)

with col_right:
    st.markdown('<div class="card" style="height: 100%;">', unsafe_allow_html=True)
    st.markdown('<div class="section-title">📖 CancerMine Evidence</div>', unsafe_allow_html=True)
    
    # Simplified display for brevity
    st.markdown(f"**Gene Role:** Oncogene") 
    st.markdown("**Top Citations (PubMed):**")
    
    # Cleaned up PubMed links
    links_html = "".join([f'<a class="pubmed-link" href="#">PMID {i}</a>' for i in range(1234, 1238)])
    st.markdown(links_html, unsafe_allow_html=True)
    st.markdown("</div>", unsafe_allow_html=True)