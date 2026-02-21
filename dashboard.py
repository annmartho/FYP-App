import streamlit as st
import pandas as pd

# -------------------------------------------------
# PAGE CONFIG 
# -------------------------------------------------
st.set_page_config(page_title="Hybrid QNN Druggability Dashboard", layout="wide")

# -------------------------------------------------
# CSS: CLEAN MODERN TECH + MIDNIGHT BLUE
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

    /* Cards & Containers */
    .metric-card, .section-box {
        background-color: #1A1A1A;
        border: 1px solid #333;
        padding: 24px;
        border-radius: 12px;
        transition: all 0.3s ease;
        box-shadow: 0 4px 12px rgba(0, 0, 0, 0.3);
        margin-bottom: 1rem;
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
        margin-bottom: 1.5rem;
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
    .low-priority { background-color: #333333; color: #AAAAAA; }
    .priority-high { background-color: #7A1F1F; color: #FFC1C1; }
    .priority-medium { background-color: #6B5E1A; color: #FFF3C1; }
    .priority-low { background-color: #333333; color: #AAAAAA; }

    /* Dividers */
    hr {
        margin: 2em 0;
        border: 0;
        border-top: 1px solid #333;
    }
    
    /* Input field styling */
    .stSelectbox div[data-baseweb="select"] {
        background-color: #1A1A1A !important;
        color: #EAEAEA !important;
        border: 1px solid #333 !important;
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

    /* Expanders & Code */
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

    .pubmed-link {
        color: #4DA6FF !important;
        text-decoration: none;
        padding: 4px 8px;
        background: #262626;
        border-radius: 4px;
        margin-right: 8px;
        font-size: 13px;
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
        return pd.read_csv("data/QNN_RESULTS.csv")
    except:
        # Mock data for demonstration if file isn't found
        return pd.DataFrame({
            "HUGO_ID": ["EGFR", "BRCA1", "TP53", "KRAS", "BRAF"],
            "Rank": [1, 45, 102, 5, 2],
            "Score": [0.9876, 0.4532, 0.1234, 0.8845, 0.9543],
            "Predicted": ["Druggable", "Non-Druggable", "Non-Druggable", "Druggable", "Druggable"],
            "Ground_Truth": ["Druggable", "Non-Druggable", "Non-Druggable", "Druggable", "Druggable"]
        })

@st.cache_data
def load_feature_vector_data():
    """Load feature vector data for all genes"""
    feature_df = pd.read_csv("data/FeatureVectorAll.csv", index_col=0)
    return feature_df

df = load_data()

# -------------------------------------------------
# MAIN DASHBOARD TOP 
# -------------------------------------------------
st.markdown("<h1>🧬 Hybrid QNN Gene Druggability Dashboard 🧬</h1>", unsafe_allow_html=True)

# Simplified to just the selectbox (which supports typing to search automatically)
all_genes = sorted(df["HUGO_ID"].astype(str).unique())
gene_input = st.selectbox("Search / Select Target Gene (HUGO ID):", all_genes)

st.markdown("<hr style='margin-top: 0.5em; margin-bottom: 2em;'>", unsafe_allow_html=True)

# -------------------------------------------------
# DASHBOARD CONTENT
# -------------------------------------------------
if gene_input:
    result = df[df["HUGO_ID"] == gene_input]

    if not result.empty:
        rank = result["Rank"].values[0]
        score = float(result["Score"].values[0])
        predicted = str(result["Predicted"].values[0])
        ground_truth = str(result["Ground_Truth"].values[0])

        col_title, col_button = st.columns([4, 1])
        
        with col_title:
            st.markdown(f"## {gene_input} Analysis")
        
        with col_button:
            def _is_unlabeled(val):
                if pd.isna(val):
                    return True
                s = str(val).strip().lower()
                if s == "":
                    return True
                # treat values containing 'unlabeled' or starting with 'n/a' as unlabeled
                if "unlabeled" in s or s.startswith("n/a"):
                    return True
                return s in ("none", "unknown", "na", "nan", "-")

            if _is_unlabeled(ground_truth):
                if st.button("🔬 Validation"):
                    st.switch_page("pages/validation.py")
            else:
                st.markdown(
                    "<div style='text-align:center; padding-top:12px;'><span class='status-pill' style='background-color: #1A4D2E; color: #A8E6D3;'>Established Drug Target</span></div>",
                    unsafe_allow_html=True,
                )
        
        # Pill colors logic (Green for Druggable, Red/Gray for Non-Druggable)
        pred_class = "high-priority" if "Druggable" in predicted and "Non" not in predicted else "low-priority"
        truth_class = "high-priority" if "Druggable" in ground_truth and "Non" not in ground_truth else "low-priority"

        st.markdown(f"""
            <div class="status-container">
                <span style="color: #8ba3cb; font-size: 13px; font-weight: 600; text-transform: uppercase; margin-left: 20px;">Ground Truth:</span> 
                <span class="status-pill {truth_class}">{ground_truth}</span>
            </div>
        """, unsafe_allow_html=True)

        # Metrics Row: show Target Rank only if unlabeled
        is_unlabeled = _is_unlabeled(ground_truth)
        num_cols = 3 if is_unlabeled else 2
        cols = st.columns(num_cols)
        
        cols[0].markdown(f"""
            <div class="metric-card">
                <div class="section-label">QNN Rank</div>
                <div class="metric-value">#{rank}</div>
            </div>
        """, unsafe_allow_html=True)

        cols[1].markdown(f"""
            <div class="metric-card">
                <div class="section-label">QNN Druggability Score</div>
                <div class="metric-value">{score:.4f}</div>
            </div>
        """, unsafe_allow_html=True)

        # Load and display Target Rank only if unlabeled
        if is_unlabeled:
            final_rank_display = "N/A"
            try:
                pri_df = pd.read_csv("data/serialized_gene_priorities.csv")
                if "Gene" in pri_df.columns and "Final Rank" in pri_df.columns:
                    match = pri_df[pri_df["Gene"].astype(str).str.strip().str.upper() == gene_input.strip().upper()]
                    if not match.empty:
                        fr = match.iloc[0]["Final Rank"]
                        final_rank_display = f"#{int(fr)}" if pd.notna(fr) else "N/A"
            except Exception:
                final_rank_display = "N/A"

            cols[2].markdown(f"""
                <div class="metric-card">
                    <div class="section-label">Target Rank</div>
                    <div class="metric-value">{final_rank_display}</div>
                </div>
            """, unsafe_allow_html=True)

        st.markdown("<br>", unsafe_allow_html=True)

        # Gene Feature Vector Section: render metrics inside the styled container
        feature_df = load_feature_vector_data()
        gene_features = feature_df[feature_df["Hugo_id"].str.upper() == gene_input.upper()]

        if gene_features.empty:
            st.info(f"No feature vector data available for {gene_input}.")
        else:
            feature_row = gene_features.iloc[0]

            feature_items = [
                ("Mutation Score", "Mutation_Value"),
                ("Proteomics Score", "Proteomics_Value"),
                ("Methylation Score", "Methylation_Value"),
                ("Expression Score", "Expression_Value")
            ]

            # Build an HTML grid of feature metrics without leading indentation
            item_htmls = []
            for label, col_name in feature_items:
                if col_name in feature_row.index:
                    value = feature_row[col_name]
                    val_str = f"{value:.4f}" if isinstance(value, (int, float)) else str(value)
                else:
                    val_str = "N/A"

                item_htmls.append(
                    f'<div style="flex:1; min-width:150px; margin:8px 8px 8px 0;">'
                    + f'<div style="font-size:12px; text-transform:uppercase; color:#999; font-weight:600;">{label}</div>'
                    + f'<div style="font-size:20px; font-weight:700; color:#FFFFFF; margin-top:6px;">{val_str}</div>'
                    + '</div>'
                )

            grid_html = '<div style="display:flex; gap:16px; flex-wrap:wrap; align-items:flex-start;">' + ''.join(item_htmls) + '</div>'

            container_html = (
                '<div style="background-color: #1A1A1A; border: 1px solid #333; padding: 24px; border-radius: 12px; margin-bottom: 1rem;">'
                + '<div style="font-size: 20px; font-weight: 600; color: #4DA6FF; margin-bottom: 16px;">Gene Feature Vector</div>'
                + grid_html
                + '</div>'
            )

            st.markdown(container_html, unsafe_allow_html=True)

        st.markdown("<br>", unsafe_allow_html=True)

        # Save to session state
        st.session_state["selected_gene"] = gene_input
        st.session_state["rank"] = rank
        st.session_state["score"] = score
            
st.markdown("<hr>", unsafe_allow_html=True)

# Full Table display
if st.checkbox("Show Full QNN Ranking Table"):
    st.dataframe(df, use_container_width=True)