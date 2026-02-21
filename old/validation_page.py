import streamlit as st
import pandas as pd

# ---------------------------------------
# Page Setup
# ---------------------------------------
st.set_page_config(layout="wide")

st.markdown("""
<style>
.big-font {
    font-size:22px !important;
    font-weight:600;
}
.badge-high {color:white;background:#1f3b57;padding:6px 12px;border-radius:20px;}
.badge-medium {color:white;background:#4a6c8c;padding:6px 12px;border-radius:20px;}
.badge-low {color:white;background:#9e9e9e;padding:6px 12px;border-radius:20px;}
</style>
""", unsafe_allow_html=True)

st.title("🧬 Gene Validation Explorer")
st.caption("Automated multi-source validation summary")

# ---------------------------------------
# Load and Clean CSV
# ---------------------------------------
@st.cache_data
def load_validation_data():
    df = pd.read_csv("data/api_validation.csv")

    # Strip column whitespace
    df.columns = [col.strip() for col in df.columns]

    # Clean and normalize fields
    df["Interacts_With"] = df["Interacts_With"].fillna("None").astype(str)
    df["Priority"] = df["Priority"].fillna("Low").astype(str)
    df["Target_Class"] = df["Target_Class"].fillna("Unclassified").astype(str)
    df["Druggability"] = df["Druggability"].fillna("Low").astype(str)
    df["Papers"] = df["Papers"].fillna(0)

    return df

df = load_validation_data()

# ---------------------------------------
# Gene Selection
# ---------------------------------------
gene = st.selectbox(
    "Select Gene",
    sorted(df["Gene"].unique())
)

# ---------------------------------------
# Display Results
# ---------------------------------------
if gene:

    gene_data = df[df["Gene"] == gene].iloc[0]

    st.divider()

    # ================================
    # Priority Badge
    # ================================
    priority = gene_data["Priority"]

    if priority.lower() == "high":
        badge_class = "badge-high"
    elif priority.lower() == "medium":
        badge_class = "badge-medium"
    else:
        badge_class = "badge-low"

    st.markdown(
        f"<div class='big-font'>Priority Status: "
        f"<span class='{badge_class}'>{priority}</span></div>",
        unsafe_allow_html=True
    )

    st.divider()

    # ================================
    # Core Metrics Row
    # ================================
    col1, col2, col3 = st.columns(3)

    col1.metric("Target Class", gene_data["Target_Class"])
    col2.metric("Druggability", gene_data["Druggability"])
    col3.metric("Literature Mentions", int(gene_data["Papers"]))

    st.divider()

    # ================================
    # Interaction Section
    # ================================
    st.subheader("Driver Interaction Network")

    interactions = str(gene_data["Interacts_With"]).strip()

    if interactions.lower() not in ["none", "nan", ""]:
        genes = interactions.split(", ")

        st.write(f"Interacts with {len(genes)} reference driver(s):")

        for g in genes:
            st.markdown(f"- **{g}**")

    else:
        st.info("No high-confidence interactions with reference drivers detected.")

    st.divider()

    # ================================
    # Summary Block
    # ================================
    st.subheader("Validation Summary")

    summary_text = f"""
    The gene **{gene}** is classified as **{priority} priority** based on:

    • Target class: {gene_data['Target_Class']}  
    • Druggability assessment: {gene_data['Druggability']}  
    • Literature evidence count: {int(gene_data['Papers'])}
    """

    if interactions.lower() not in ["none", "nan", ""]:
        summary_text += "\n• Demonstrates mechanistic linkage to known LUAD driver genes."

    st.markdown(summary_text)
