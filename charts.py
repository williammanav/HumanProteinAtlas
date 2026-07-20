# charts.py

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go


# -------------------------------------------------
# Gene Expression Comparison
# -------------------------------------------------
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


# -------------------------------------------------
# Top Expressed Genes
# -------------------------------------------------
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


# -------------------------------------------------
# Expression Distribution
# -------------------------------------------------
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


# -------------------------------------------------
# Heatmap
# -------------------------------------------------
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


# -------------------------------------------------
# Correlation Heatmap
# -------------------------------------------------
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


# -------------------------------------------------
# Fold Change Scatter
# -------------------------------------------------
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


# -------------------------------------------------
# Boxplot
# -------------------------------------------------
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


# -------------------------------------------------
# Pie Chart
# -------------------------------------------------
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


# -------------------------------------------------
# Scatter Plot
# -------------------------------------------------
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