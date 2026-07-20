# utils.py

import streamlit as st
import pandas as pd
from pathlib import Path

def load_css(css_file_path):
    """
    Safely resolves the path to the CSS file and injects its contents
    into the Streamlit app.
    """
    css_path = Path(css_file_path)
    if css_path.exists():
        with open(css_path) as f:
            st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)
    else:
        st.warning(f"CSS file not found at {css_file_path}")

@st.cache_data
def load_data(filepath="proteinatlas_search.tsv"):
    """
    Loads the Human Protein Atlas dataset from a TSV file and caches the result.
    """
    return pd.read_csv(filepath, sep="\t")

def detect_columns(df):
    """
    Scans the dataframe columns and automatically categorizes them into:
    - normal_cols: tissue-specific RNA expression columns [nTPM]
    - cell_cols: cancer cell line expression columns [nTPM]
    - metadata: a dictionary containing column names for gene, description,
      chromosome, protein class, disease involvement, and evidence.
    """
    expression_cols = [c for c in df.columns if "[nTPM]" in c]
    normal_cols = []
    cell_cols = []
    
    for col in expression_cols:
        col_lower = col.lower()
        if "tissue" in col_lower:
            normal_cols.append(col)
        else:
            cell_cols.append(col)

    metadata = {
        "gene": None,
        "description": None,
        "protein_class": None,
        "disease": None,
        "chromosome": None,
        "evidence": None
    }
    
    for c in df.columns:
        lc = c.lower()
        if lc == "gene":
            metadata["gene"] = c
        elif "description" in lc:
            metadata["description"] = c
        elif "protein" in lc and "class" in lc:
            metadata["protein_class"] = c
        elif "disease" in lc:
            metadata["disease"] = c
        elif "chromosome" in lc:
            metadata["chromosome"] = c
        elif "evidence" in lc:
            metadata["evidence"] = c
            
    return normal_cols, cell_cols, metadata

def filter_data(df, normal_col, expression_range, protein_col=None, selected_proteins=None, disease_col=None, selected_diseases=None):
    """
    Filters the dataset based on protein class, disease involvement,
    and expression range.
    """
    filtered = df.copy()
    
    if protein_col and selected_proteins:
        filtered = filtered[filtered[protein_col].isin(selected_proteins)]
        
    if disease_col and selected_diseases:
        filtered = filtered[filtered[disease_col].isin(selected_diseases)]
        
    filtered = filtered[
        (filtered[normal_col] >= expression_range[0]) &
        (filtered[normal_col] <= expression_range[1])
    ]
    
    return filtered
