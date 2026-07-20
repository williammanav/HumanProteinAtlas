import streamlit as st
import pandas as pd
import numpy as np
from pathlib import Path
import plotly.express as px
import plotly.graph_objects as go

# ----------------------------
# Page Config
# ----------------------------
st.set_page_config(
    page_title="Human Protein Atlas Dashboard",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)


# ----------------------------
# Helper & Utility Functions
# ----------------------------

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


# ----------------------------
# Chart Functions
# ----------------------------

def gene_expression_bar(df, gene_col, gene, normal_col, cell_lines):
    """
    Bar chart comparing one gene across
    normal tissue and selected cancer cell lines.
    """
    row = df[df[gene_col] == gene]

    if row.empty:
        return None

    samples = ["Normal Tissue"]
    values = [row.iloc[0][normal_col]]

    for col in cell_lines:
        samples.append(col.replace("[nTPM]", ""))
        values.append(row.iloc[0][col])

    plot_df = pd.DataFrame({
        "Sample": samples,
        "Expression": values
    })

    fig = px.bar(
        plot_df,
        x="Sample",
        y="Expression",
        color="Sample",
        text_auto=".2f",
        title=f"{gene} Expression"
    )

    fig.update_layout(
        template="plotly_dark",
        height=420,
        showlegend=False,
        xaxis_title="",
        yaxis_title="Expression (nTPM)"
    )

    return fig


def top_genes_chart(df, expression_col, gene_col, top_n=10):
    top = (
        df[[gene_col, expression_col]]
        .sort_values(expression_col, ascending=False)
        .head(top_n)
    )

    fig = px.bar(
        top,
        x=expression_col,
        y=gene_col,
        orientation="h",
        color=expression_col,
        text_auto=".1f"
    )

    fig.update_layout(
        template="plotly_dark",
        height=420,
        yaxis=dict(autorange="reversed"),
        xaxis_title="Expression (nTPM)",
        yaxis_title=""
    )

    return fig


def expression_distribution(df, expression_col):
    fig = px.histogram(
        df,
        x=expression_col,
        nbins=40,
        opacity=0.85
    )

    fig.update_layout(
        template="plotly_dark",
        height=420,
        xaxis_title="Expression (nTPM)",
        yaxis_title="Number of Genes"
    )

    return fig


def expression_heatmap(df, gene_col, cell_lines, top_n=20):
    heat_df = (
        df[[gene_col] + cell_lines]
        .head(top_n)
        .set_index(gene_col)
    )

    fig = px.imshow(
        heat_df,
        labels=dict(
            x="Cell Line",
            y="Gene",
            color="nTPM"
        ),
        aspect="auto"
    )

    fig.update_layout(
        template="plotly_dark",
        height=650
    )

    return fig


def correlation_heatmap(df, cell_lines):
    corr = df[cell_lines].corr()

    fig = px.imshow(
        corr,
        text_auto=".2f",
        color_continuous_scale="Viridis"
    )

    fig.update_layout(
        template="plotly_dark",
        height=500
    )

    return fig


def fold_change_plot(df, gene_col, normal_col, cancer_col):
    plot_df = df[[gene_col, normal_col, cancer_col]].copy()

    plot_df["FoldChange"] = (
        plot_df[cancer_col] /
        (plot_df[normal_col] + 0.01)
    )

    fig = px.scatter(
        plot_df,
        x=normal_col,
        y=cancer_col,
        hover_name=gene_col,
        color="FoldChange"
    )

    fig.update_layout(
        template="plotly_dark",
        height=500,
        xaxis_title="Normal Tissue",
        yaxis_title="Cancer Cell Line"
    )

    fig.add_shape(
        type="line",
        x0=0,
        y0=0,
        x1=plot_df[normal_col].max(),
        y1=plot_df[normal_col].max(),
        line=dict(
            dash="dash",
            color="white"
        )
    )

    return fig


def cellline_boxplot(df, cell_lines):
    box_df = df[cell_lines].melt(
        var_name="Cell Line",
        value_name="Expression"
    )

    fig = px.box(
        box_df,
        x="Cell Line",
        y="Expression",
        color="Cell Line"
    )

    fig.update_layout(
        template="plotly_dark",
        height=500,
        showlegend=False
    )

    return fig


def protein_class_pie(df, protein_col):
    counts = (
        df[protein_col]
        .value_counts()
        .head(8)
        .reset_index()
    )

    counts.columns = ["Protein Class", "Count"]

    fig = px.pie(
        counts,
        values="Count",
        names="Protein Class",
        hole=0.45
    )

    fig.update_layout(
        template="plotly_dark",
        height=450
    )

    return fig


def normal_vs_cancer(df, gene_col, normal_col, cancer_col):
    fig = px.scatter(
        df,
        x=normal_col,
        y=cancer_col,
        hover_name=gene_col,
        opacity=0.75
    )

    fig.update_layout(
        template="plotly_dark",
        height=500,
        xaxis_title="Normal Tissue",
        yaxis_title="Cancer Cell Line"
    )

    return fig


# ----------------------------
# Load CSS
# ----------------------------
load_css("style.css")

# ----------------------------
# Load Dataset
# ----------------------------
df = load_data("proteinatlas_search.tsv")

# ----------------------------
# Automatically Detect Columns
# ----------------------------
normal_cols, cell_cols, metadata = detect_columns(df)

gene_col = metadata["gene"]
desc_col = metadata["description"]
protein_col = metadata["protein_class"]
disease_col = metadata["disease"]
chromosome_col = metadata["chromosome"]
evidence_col = metadata["evidence"]

if gene_col is None:
    st.error("Gene column not found.")
    st.stop()
if len(normal_cols) == 0:
    st.error("No normal tissue column found.")
    st.stop()
if len(cell_cols) == 0:
    st.error("No cell line columns found.")
    st.stop()

# ----------------------------
# Sidebar
# ----------------------------
st.sidebar.title("🧬 HPA Dashboard")
st.sidebar.caption("Gene Expression Explorer")

st.sidebar.divider()

st.sidebar.subheader("Filters")

# ----------------------------
# Gene Search
# ----------------------------
gene = st.sidebar.selectbox(
    "Search Gene",
    sorted(df[gene_col].dropna().unique())
)

# ----------------------------
# Tissue
# ----------------------------
if len(normal_cols):
    normal_col = st.sidebar.selectbox(
        "Normal Tissue",
        normal_cols
    )
else:
    st.sidebar.error("No normal tissue column found.")
    st.stop()

# ----------------------------
# Cell Lines
# ----------------------------
selected_cell_lines = st.sidebar.multiselect(
    "Cancer Cell Lines",
    cell_cols,
    default=cell_cols[:5]
)

# ----------------------------
# Protein Class
# ----------------------------
if protein_col:
    protein = st.sidebar.multiselect(
        "Protein Class",
        sorted(df[protein_col].dropna().unique())
    )
else:
    protein = []

# ----------------------------
# Disease
# ----------------------------
if disease_col:
    disease = st.sidebar.multiselect(
        "Disease",
        sorted(df[disease_col].dropna().unique())
    )
else:
    disease = []

# ----------------------------
# Expression Slider
# ----------------------------
all_values = []
cols = [normal_col] + selected_cell_lines

for c in cols:
    if c in df.columns:
        all_values.extend(df[c].fillna(0).tolist())

max_exp = max(1, int(np.nanmax(all_values)))

expression_range = st.sidebar.slider(
    "Expression Range (nTPM)",
    0,
    max_exp,
    (0, max_exp)
)

st.sidebar.divider()

st.sidebar.info(
    """
Human Protein Atlas (HPA)

RNA expression in normal tissues and cancer cell lines.

Explore gene expression patterns interactively.
"""
)

# ----------------------------
# Filter Data
# ----------------------------
filtered = filter_data(
    df=df,
    normal_col=normal_col,
    expression_range=expression_range,
    protein_col=protein_col,
    selected_proteins=protein,
    disease_col=disease_col,
    selected_diseases=disease
)

# ----------------------------
# Header
# ----------------------------
st.title("🧬 Human Protein Atlas Gene Expression Dashboard")

st.markdown(
"""
Explore and compare transcriptomic expression levels between
normal human tissues and cancer cell lines.
"""
)

st.divider()

# ----------------------------
# KPI Calculations
# ----------------------------
total_genes = len(filtered)
num_cell_lines = len(cell_cols)
selected_lines = len(selected_cell_lines)

if filtered.empty:
    avg_expression = 0
    high_expression = 0
else:
    avg_expression = round(filtered[normal_col].mean(), 2)
    high_expression = (filtered[normal_col] > 100).sum()

# ----------------------------
# KPI Cards
# ----------------------------
c1, c2, c3, c4, c5 = st.columns(5)

with c1:
    st.metric(
        "Total Genes",
        f"{total_genes:,}"
    )

with c2:
    st.metric(
        "Cancer Cell Lines",
        num_cell_lines
    )

with c3:
    st.metric(
        "Average Expression",
        f"{avg_expression} nTPM"
    )

with c4:
    st.metric(
        "Genes >100 nTPM",
        high_expression
    )

with c5:
    st.metric(
        "Selected Cell Lines",
        selected_lines
    )

st.divider()

# ==========================================================
# Layout
# ==========================================================
row1_col1, row1_col2, row1_col3 = st.columns([2, 1, 1])

with row1_col1:
    st.subheader("Gene Expression Comparison")
    if not filtered.empty:
        fig_bar = gene_expression_bar(df, gene_col, gene, normal_col, selected_cell_lines)
        if fig_bar is not None:
            st.plotly_chart(fig_bar, use_container_width=True)
        else:
            st.warning("No data found for selected gene.")
    else:
        st.warning("No filtered data available.")

with row1_col2:
    st.subheader("Top Expressed Genes")
    if not filtered.empty:
        fig_top = top_genes_chart(filtered, normal_col, gene_col)
        st.plotly_chart(fig_top, use_container_width=True)
    else:
        st.warning("No filtered data available.")

with row1_col3:
    st.subheader("Expression Distribution")
    if not filtered.empty:
        fig_dist = expression_distribution(filtered, normal_col)
        st.plotly_chart(fig_dist, use_container_width=True)
    else:
        st.warning("No filtered data available.")

st.divider()

row2_col1, row2_col2, row2_col3 = st.columns([2, 2, 1])

with row2_col1:
    tab_heat, tab_box = st.tabs(["Expression Heatmap", "Cell Line Boxplot"])
    
    with tab_heat:
        if not filtered.empty:
            if selected_cell_lines:
                fig_heat = expression_heatmap(filtered, gene_col, selected_cell_lines)
                st.plotly_chart(fig_heat, use_container_width=True)
            else:
                st.warning("Select at least one cancer cell line.")
        else:
            st.warning("No filtered data available.")
            
    with tab_box:
        if not filtered.empty:
            if selected_cell_lines:
                fig_box = cellline_boxplot(filtered, selected_cell_lines)
                st.plotly_chart(fig_box, use_container_width=True)
            else:
                st.warning("Select at least one cancer cell line.")
        else:
            st.warning("No filtered data available.")

with row2_col2:
    st.subheader("Fold Change")
    if not filtered.empty:
        if selected_cell_lines:
            fc_cell_line = st.selectbox(
                "Compare Normal Tissue with Cell Line:",
                selected_cell_lines,
                key="fold_change_cell_line"
            )
            fig_fc = fold_change_plot(filtered, gene_col, normal_col, fc_cell_line)
            st.plotly_chart(fig_fc, use_container_width=True)
        else:
            st.warning("Please select at least one cancer cell line.")
    else:
        st.warning("No filtered data available.")

with row2_col3:
    tab_corr, tab_pie = st.tabs(["Correlation", "Protein Classes"])
    
    with tab_corr:
        if not filtered.empty:
            if len(selected_cell_lines) >= 2:
                fig_corr = correlation_heatmap(filtered, selected_cell_lines)
                st.plotly_chart(fig_corr, use_container_width=True)
            else:
                st.warning("Select at least 2 cancer cell lines to show correlation.")
        else:
            st.warning("No filtered data available.")
            
    with tab_pie:
        if not filtered.empty:
            if protein_col:
                fig_pie = protein_class_pie(filtered, protein_col)
                st.plotly_chart(fig_pie, use_container_width=True)
            else:
                st.info("Protein class column not available.")
        else:
            st.warning("No filtered data available.")

st.divider()

st.subheader("Gene Information")

gene_row = filtered[
    filtered[gene_col] == gene
]

if not gene_row.empty:
    info = gene_row.iloc[0]
    col1, col2, col3 = st.columns(3)

    with col1:
        st.write("**Gene**")
        st.write(info[gene_col])

        if chromosome_col and chromosome_col in info:
            st.write("**Chromosome**")
            st.write(info[chromosome_col])

    with col2:
        if protein_col and protein_col in info:
            st.write("**Protein Class**")
            st.write(info[protein_col])

        if evidence_col and evidence_col in info:
            st.write("**Evidence**")
            st.write(info[evidence_col])

    with col3:
        if disease_col and disease_col in info:
            st.write("**Disease**")
            st.write(info[disease_col])

        st.write("**Normal Expression**")
        st.write(round(info[normal_col], 2))

st.divider()

st.subheader("Gene Expression Data")

table_cols = [gene_col]
for c in [
    desc_col,
    chromosome_col,
    protein_col,
    disease_col,
    evidence_col,
]:
    if c:
        table_cols.append(c)
table_cols.append(normal_col)
table_cols.extend(selected_cell_lines)
table_cols = [c for c in table_cols if c in filtered.columns]

st.dataframe(
    filtered[table_cols],
    use_container_width=True,
    height=450
)

csv = filtered.to_csv(index=False).encode()

st.download_button(
    "⬇ Download Filtered Dataset",
    csv,
    "Filtered_HPA.csv",
    "text/csv"
)
