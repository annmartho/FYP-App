import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter

# ---------------------------------------------------
# PAGE CONFIG & CSS: UNIFIED DARK THEME
# ---------------------------------------------------
st.set_page_config(page_title="Validation", layout="wide")

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

    /* Cards & Containers */
    .stElementContainer div[data-testid="stMarkdownContainer"] > div.section-box {
        background-color: #1A1A1A;
        padding: 24px;
        border-radius: 12px;
        border: 1px solid #333;
        margin-bottom: 20px;
    }

    .section-title {
        font-size: 20px;
        font-weight: 600;
        color: #4DA6FF;
        display: block;
    }

    .priority-badge {
        padding: 6px 14px;
        border-radius: 20px;
        font-size: 14px;
        font-weight: bold;
        display: inline-block;
        margin-bottom: 15px;
    }
    .priority-high { background-color: #7A1F1F; color: #FFC1C1; }
    .priority-medium { background-color: #6B5E1A; color: #FFF3C1; }
    .priority-low { background-color: #333333; color: #AAAAAA; }

    .pubmed-link {
        color: #4DA6FF !important;
        text-decoration: none;
        padding: 4px 8px;
        background: #262626;
        border-radius: 4px;
        margin-right: 8px;
        font-size: 13px;
    }

    /* Dividers */
    hr {
        margin: 2em 0;
        border: 0;
        border-top: 1px solid #333;
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

    .stExpander {
        border: 1px solid #333 !important;
        background-color: #1A1A1A !important;
        border-radius: 12px !important;
    }
</style>
""", unsafe_allow_html=True)

# ---------------------------------------------------
# Data Loading (Optimized for Memory)
# ---------------------------------------------------
@st.cache_data
def load_validation_data():
    """Load fast, small datasets"""
    val_df = pd.read_csv("data/api_validation.csv")
    val_df.columns = [col.strip() for col in val_df.columns]
    
    imp_genes = pd.read_csv("survival_analysis_csv/imp_genes.csv")
    hugo_to_ens = dict(zip(imp_genes["HUGO_Symbol"], imp_genes["Ensembl_ID"]))
    
    collated = pd.read_csv("cancermine_csv/collated_filtered.csv")
    sentences = pd.read_csv("cancermine_csv/sentences_filtered.csv")
    collated["gene_normalized"] = collated["gene_normalized"].str.upper()
    
    return val_df, hugo_to_ens, collated, sentences

@st.cache_data
def load_expression_data():
    """Load and cache expression data once"""
    expr_df = pd.read_csv("survival_analysis_csv\M1_counts_raw.csv", low_memory=False)
    expr_df = expr_df.rename(columns={expr_df.columns[0]: "Gene"}).set_index("Gene")
    return expr_df

@st.cache_data
def load_survival_data():
    """Load and cache survival data once"""
    surv_df = pd.read_csv("survival_analysis_csv/LUAD_survival.txt", sep="\t")
    surv_df["sample"] = surv_df["sample"].str.slice(0, 12)
    return surv_df

def get_survival_data_for_gene(gene_name, expr_df, surv_df):
    """Extract survival data for a specific gene from cached data"""
    if gene_name not in expr_df.index:
        return pd.DataFrame()
    
    gene_expr = expr_df.loc[[gene_name]].T.reset_index()
    gene_expr.columns = ["sample", gene_name]
    gene_expr["sample"] = gene_expr["sample"].astype(str)
    
    # Convert to log scale
    gene_expr[gene_name] = np.log2(gene_expr[gene_name].astype(float) + 1)
    
    # Merge with survival data
    merged = surv_df.merge(gene_expr, on="sample", how="inner")
    return merged

# Load fast data immediately
val_df, hugo_to_ens, collated_df, sentences_df = load_validation_data()

# Gene Check
if "selected_gene" not in st.session_state:
    st.warning("Please select a gene from the dashboard.")
    st.stop()

gene = st.session_state["selected_gene"]
no_pdb = ['MTMR7',
 'ADAMTS16',
 'CIT',
 'LTBP4',
 'MMRN1',
 'PAX9',
 'SOX5',
 'STXBP6',
 'PDZD2',
 'ITIH5',
 'VWA3B',
 'DENND2A',
 'ATP10B',
 'MEOX2',
 'RGS22',
 'STXBP5L',
 'GDF10',
 'C1orf162',
 'PLCL1',
 'SYNPO2',
 'RP1',
 'CLEC14A',
 'TBX4']

# ---------------------------------------------------
# UI Layout
# ---------------------------------------------------
col_title, col_button = st.columns([4, 1])

with col_title:
    st.markdown(f"## Target Validation: {gene}")

with col_button:
    if st.button("📊 Structural Data"):
        if gene in no_pdb:
            st.switch_page("pages/structural_page_withoutpdb.py")
        else:
            st.switch_page("pages/structural_page_withpdb.py")

# 1. API Validation Section
with st.container():
    # We wrap the content in a single markdown block to avoid the "blob" fragments
    gene_val = val_df[val_df["Gene"].str.upper() == gene.upper()]
    
    if gene_val.empty:
        st.info(f"No repository validation data available for {gene}.")
    else:
        v = gene_val.iloc[0]
        p = str(v["Priority"]).lower()
        
        # Start Section Box
        st.markdown(f"""
        <div class="section-box">
            <div class="section-title">🔬 Repository Validation Details</div>
            <div class="priority-badge priority-{p}">{p.upper()} PRIORITY</div>
        </div>
        """, unsafe_allow_html=True)
        
        c1, c2, c3 = st.columns(3)
        c1.metric("Target Class", v["Target_Class"])
        c2.metric("Druggability", v["Druggability"])
        c3.metric("Literature Mentions", int(v["Papers"]))

# 2. Side-by-Side Analysis
col_left, col_right = st.columns(2, gap="medium")

with col_left:
    st.markdown("""
    <div class="section-box" style="background-color: #1A1A1A; border: 1px solid #333; padding: 12px; border-radius: 12px; margin-bottom: 10px;">
        <div class="section-title">Survival Analysis - Kaplan Meier Plot</div>
    """, unsafe_allow_html=True)
    
    if gene in hugo_to_ens:
        ens_id = hugo_to_ens[gene]
        
        # Load cached data and extract for this gene
        expr_df = load_expression_data()
        surv_df = load_survival_data()
        merged = get_survival_data_for_gene(ens_id, expr_df, surv_df)
        
        if not merged.empty and ens_id in merged.columns:
            df_surv = merged.dropna(subset=["OS", "OS.time", ens_id])
            if not df_surv.empty:
                median_val = df_surv[ens_id].median()
                df_surv["group"] = np.where(df_surv[ens_id] >= median_val, "High", "Low")
                
                fig, ax = plt.subplots(figsize=(5, 4))
                plt.style.use("dark_background")
                ax.set_facecolor("#1A1A1A")
                fig.patch.set_facecolor("#0F0F0F")
                
                kmf = KaplanMeierFitter()
                for label in ["High", "Low"]:
                    mask = df_surv["group"] == label
                    if mask.sum() > 0:
                        kmf.fit(df_surv[mask]["OS.time"], df_surv[mask]["OS"], label=f"{label} Expr")
                        kmf.plot(ax=ax)
                
                st.pyplot(fig)
            else:
                st.write("Insufficient data for survival analysis.")
        else:
            st.write("Expression data missing.")
    else:
        st.write("Gene ID not found in survival dataset.")
    
    st.markdown("</div>", unsafe_allow_html=True)

with col_right:
    st.markdown("""
    <div class="section-box" style="background-color: #1A1A1A; border: 1px solid #333; padding:10px; border-radius: 12px; margin-bottom: 10px;">
        <div class="section-title">CancerMine Evidence</div>
    """, unsafe_allow_html=True)
    
    cm_data = collated_df[collated_df["gene_normalized"] == gene.upper()]
    
    if cm_data.empty:
        st.write("No literature evidence found.")
    else:
        for _, row in cm_data.iterrows():
            st.markdown(f"**{row['cancer_normalized']}** — *{row['role']}*")
            st.progress(min(int(row["citation_count"])/500, 1.0))
            
            # Filter PMIDs
            pmids = sentences_df[(sentences_df["gene_normalized"] == gene.upper()) & 
                                 (sentences_df["cancer_normalized"] == row["cancer_normalized"])]["pmid"].unique()[:5]
            
            links_html = "".join([f'<a class="pubmed-link" href="https://pubmed.ncbi.nlm.nih.gov/{int(p)}/" target="_blank">PMID:{int(p)}</a>' for p in pmids])
            st.markdown(links_html, unsafe_allow_html=True)
            st.write("") # Spacer
    
    st.markdown("</div>", unsafe_allow_html=True)

if st.button("⬅ Back to Dashboard"):
    st.switch_page("dashboard.py")