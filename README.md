# Human Protein Atlas (HPA) Gene Expression Dashboard

🧬 An interactive, high-performance Streamlit dashboard designed to explore, compare, and analyze RNA transcriptomic expression levels (nTPM) between normal human tissues and cancer cell lines.

---

## Features

- **Interactive Filtering**: Filter by gene name, normal tissue type, cancer cell line list, protein classes, disease involvement, and expression ranges.
- **Gene Info Lookup**: Dynamically details specific gene chromosome mappings, evidence levels, disease relevance, and normal expression levels.
- **Rich Data Visualizations**:
  - **Gene Expression Comparison**: Side-by-side bar chart comparison of the selected gene.
  - **Top Expressed Genes**: Horizontal ranking of highly expressed genes under the selected tissue criteria.
  - **Expression Distribution**: Histogram of gene expression frequencies.
  - **Heatmap & Boxplot Tabs**: Dynamic toggle between cell-line heatmaps and cell-line expression range boxplots.
  - **Fold Change Comparison**: Dynamic comparison scatter plot showing fold changes relative to a chosen cell line.
  - **Correlation & Classification Tabs**: Dynamic correlation heatmap comparing selected cell lines paired with a protein class distribution pie chart.
- **Data Export**: Inspect the filtered dataset in a responsive data table and download the selection directly as a CSV.
- **Premium Aesthetics**: Customized dark theme styling (`style.css`) featuring custom component borders, hover micro-animations, custom scrollbars, and colored indicator metrics.

---

## File Structure

```text
├── main.py                 # Streamlit application layout, tabs, and page setup
├── utils.py                # Data loading (cached), CSS injector, column detector, and filters
├── charts.py               # Custom Plotly figure generation functions (dark template)
├── style.css               # Premium CSS stylesheets and component selectors
├── requirements.txt        # Application package dependencies
└── proteinatlas_search.tsv # Human Protein Atlas TSV search dataset
```

---

## Installation & Running

### 1. Install Dependencies
Make sure you have python installed, then install all project requirements:
```bash
pip install -r requirements.txt
```

### 2. Launch the Application
Run the Streamlit server from your project workspace directory:
```bash
streamlit run main.py
```

The application will automatically launch in your default web browser at `http://localhost:8501`.
