from flask import Flask, render_template, request
import pandas as pd
import os

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

import umap.umap_ as umap

import gseapy as gp

import matplotlib

from matplotlib.colors import LinearSegmentedColormap

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from werkzeug.utils import secure_filename

app = Flask(__name__)

UPLOAD_FOLDER = "uploads"
CHART_FOLDER = "static/charts"

os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(CHART_FOLDER, exist_ok=True)

app.config["UPLOAD_FOLDER"] = UPLOAD_FOLDER

latest_df = None

latest_results = {
    "uploaded_filename": "No dataset uploaded",
    "num_rows": 0,
    "num_genes": 0,
    "num_conditions": 0,
    "avg_expression": 0,
    "table": None,
    "chart": None,
    "condition_chart": None,
    "top_degs": [],
    "expression_bars": [],
    "heatmap_cells": [],
    "condition_bars": [],
    "top_variable_genes": [],
    "enriched_pathways": [],
    "enrichment_chart": None,
    "qc_metrics": {
        "rows": 0,
        "missing_removed": 0,
        "genes": 0,
        "conditions": 0,
        "completeness": 0,
    },
}

def find_column(df, possible_names):
    for col in df.columns:
        if col in possible_names:
            return col
    return None


def run_analysis(df, gene_col, expr_col, cond_col, filename):
    global latest_results, latest_df

    original_rows = len(df)

    df[expr_col] = pd.to_numeric(df[expr_col], errors="coerce")
    df = df.dropna(subset=[expr_col])

    missing_removed = original_rows - len(df)

    df = df.rename(
        columns={
            gene_col: "gene",
            expr_col: "expression_level",
            cond_col: "condition",
        }
    )

    latest_df = df.copy()

    # ---- BASIC METRICS ----
    num_rows = len(df)
    num_genes = df["gene"].nunique()
    num_conditions = df["condition"].nunique()
    avg_expression = round(df["expression_level"].mean(), 2)

    total_cells = df.shape[0] * df.shape[1]
    missing_cells = df.isna().sum().sum()
    completeness = (
        round((1 - (missing_cells / total_cells)) * 100, 1) if total_cells else 0
    )

    table = df.head(10).to_html(classes="data", index=False)

    # ---- BAR CHART DATA ----
    gene_avg = (
        df.groupby("gene")["expression_level"].mean().sort_values(ascending=False)
    )

    top_expression = gene_avg.head(8)
    max_expression = top_expression.max() if not top_expression.empty else 1

    expression_bars = [
        {
            "gene": gene,
            "value": round(val, 2),
            "height": round((val / max_expression) * 100, 0) if max_expression else 0,
        }
        for gene, val in top_expression.items()
    ]

    # ---- MATPLOTLIB CHARTS (for other pages) ----
    plt.style.use('dark_background')

    chart_path = os.path.join(CHART_FOLDER, "chart.png")

    plt.figure(figsize=(7, 5), facecolor='#0B1E2D')
    gene_avg.plot(kind="bar", color="#3D8BFF")

    plt.title("Average Expression by Gene", color="white")
    plt.xlabel("Gene", color="white")
    plt.ylabel("Expression Level", color="white")

    plt.xticks(color="white")
    plt.yticks(color="white")

    plt.tight_layout()
    plt.savefig(chart_path)
    plt.close()

    condition_avg = (
        df.groupby(["gene", "condition"])["expression_level"].mean().unstack()
    )

    # ---- UMAP PLOT ----
    umap_chart_path = os.path.join(CHART_FOLDER, "umap_chart.png")
    umap_chart = None

    umap_data = condition_avg.fillna(0).T

    if umap_data.shape[0] >= 3 and umap_data.shape[1] >= 2:
        scaled_umap_data = StandardScaler().fit_transform(umap_data)

        reducer = umap.UMAP(
            n_components=2, random_state=42, n_neighbors=min(5, umap_data.shape[0] - 1)
        )

        umap_result = reducer.fit_transform(scaled_umap_data)

        plt.figure(figsize=(7, 5), facecolor="#0B1E2D")
        ax = plt.gca()
        ax.set_facecolor("#071421")

        plt.scatter(
            umap_result[:, 0],
            umap_result[:, 1],
            s=90,
            color="#7C3AED",
            edgecolors="#7DD3FC",
            linewidths=1.2,
        )

        for i, condition in enumerate(umap_data.index):
            plt.text(
                umap_result[i, 0],
                umap_result[i, 1],
                f"  {condition}",
                color="white",
                fontsize=9,
            )

        plt.title("UMAP Projection", color="white")
        plt.xlabel("UMAP 1", color="white")
        plt.ylabel("UMAP 2", color="white")

        plt.xticks(color="white")
        plt.yticks(color="white")

        plt.grid(alpha=0.15)
        plt.tight_layout()
        plt.savefig(umap_chart_path, transparent=False)
        plt.close()

        umap_chart = "charts/umap_chart.png"

    # ---- PCA PLOT ----
    pca_chart_path = os.path.join(CHART_FOLDER, "pca_chart.png")

    pca_data = condition_avg.fillna(0).T

    pca_chart = None

    if pca_data.shape[0] >= 2 and pca_data.shape[1] >= 2:
        scaled_data = StandardScaler().fit_transform(pca_data)

        pca = PCA(n_components=2)
        pca_result = pca.fit_transform(scaled_data)

        plt.figure(figsize=(7, 5), facecolor="#0B1E2D")
        ax = plt.gca()
        ax.set_facecolor("#071421")

        plt.scatter(
            pca_result[:, 0],
            pca_result[:, 1],
            s=90,
            color="#00C896",
            edgecolors="#7DD3FC",
            linewidths=1.2,
        )

        for i, condition in enumerate(pca_data.index):
            plt.text(
                pca_result[i, 0],
                pca_result[i, 1],
                f"  {condition}",
                color="white",
                fontsize=9,
            )

        plt.title("PCA Projection", color="white")
        plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)", color="white")
        pc2_var = (
        pca.explained_variance_ratio_[1] * 100
        if len(pca.explained_variance_ratio_) > 1
        else 0
        )

        plt.ylabel(
            f"PC2 ({pc2_var:.1f}%)",
            color="white"
        )

        plt.xticks(color="white")
        if pc2_var < 0.01:
            plt.yticks([])
        else:
            plt.yticks(color="white")

        plt.grid(alpha=0.15)
        plt.ticklabel_format(style="plain")
        plt.axhline(0, color="#123056", linewidth=1)
        plt.tight_layout()
        plt.savefig(pca_chart_path, transparent=False)
        plt.close()

        pca_chart = "charts/pca_chart.png"
    # ---- REAL HEATMAP IMAGE ----
    heatmap_chart_path = os.path.join(CHART_FOLDER, "heatmap_chart.png")

    # Select most variable genes
    gene_variance = condition_avg.var(axis=1)

    # ---- TOP VARIABLE GENES ----
    top_variable_genes = [
        {
            "gene": gene,
            "variance": round(value, 2),
            "width": round(
                (value / gene_variance.max()) * 100,
                0
            ) if gene_variance.max() else 0,
        }
        for gene, value in (
            gene_variance
            .sort_values(ascending=False)
            .head(5)
            .items()
        )
    ]

    top_genes = (
        gene_variance
        .sort_values(ascending=False)
        .head(25)
        .index
    )

    heatmap_data = condition_avg.loc[top_genes].fillna(0)

    # Log normalization
    import numpy as np
    heatmap_data = np.log1p(heatmap_data)

    plt.figure(figsize=(8, 5), facecolor='#0B1E2D')

    genescope_cmap = LinearSegmentedColormap.from_list(
        "genescope",
        [
            "#071421",  # deep background
            "#123056",  # navy
            "#2563EB",  # biotech blue
            "#00C896",  # teal
            "#7DD3FC",  # light cyan
        ],
    )

    plt.imshow(
        heatmap_data,
        aspect='auto',
        cmap=genescope_cmap
    )
    # Blend chart with UI
    plt.gca().set_facecolor('#071421')
    plt.box(False)

    plt.colorbar(label='Expression')

    plt.xticks(
        range(len(heatmap_data.columns)),
        heatmap_data.columns,
        rotation=20,
        color='white'
    )

    plt.yticks(
        range(len(heatmap_data.index)),
        heatmap_data.index,
        color='white'
    )

    plt.title(
        "Gene × Condition Heatmap",
        color='white'
    )

    plt.tight_layout()
    plt.savefig(
        heatmap_chart_path,
        transparent=False
    )

    plt.close()

    # Custom grouped bar data for Visualizations page
    condition_bars = []

    condition_avg_for_ui = condition_avg.fillna(0)
    max_condition_value = (
        condition_avg_for_ui.max().max() if not condition_avg_for_ui.empty else 1)

    for gene, row in condition_avg_for_ui.head(8).iterrows():
        condition_bars.append(
            {
                "gene": gene,
                "conditions": [
                    {
                        "name": condition,
                        "value": round(value, 2),
                        "height": (
                            round((value / max_condition_value) * 100, 0)
                            if max_condition_value
                            else 0
                        ),
                    }
                    for condition, value in row.items()
                ],
            })

    condition_chart_path = os.path.join(CHART_FOLDER, "condition_chart.png")

    plt.figure(figsize=(7, 5), facecolor='#0B1E2D')
    condition_avg.plot(kind="bar", color=["#3D8BFF", "#00C896"])

    plt.title("Expression by Gene and Condition", color="white")
    plt.xlabel("Gene", color="white")
    plt.ylabel("Expression Level", color="white")

    plt.xticks(color="white")
    plt.yticks(color="white")

    plt.tight_layout()
    plt.savefig(condition_chart_path)
    plt.close()

    # ---- HEATMAP (FIXED VERSION) ----
    top_heatmap_genes = gene_avg.head(8).index.tolist()
    heatmap_df = df[df["gene"].isin(top_heatmap_genes)].head(48)

    max_val = heatmap_df["expression_level"].max() if not heatmap_df.empty else 1
    min_val = heatmap_df["expression_level"].min() if not heatmap_df.empty else 0

    range_val = max_val - min_val if max_val != min_val else 1

    heatmap_cells = [
        {
            "gene": row["gene"],
            "condition": row["condition"],
            "value": round(row["expression_level"], 2),
            "intensity": round(
                ((row["expression_level"] - min_val) / range_val) * 100,
                0,
            ),
        }
        for _, row in heatmap_df.iterrows()
    ]

    # ---- DEGs ----
    deg = condition_avg.fillna(0)
    condition_names = list(condition_avg.columns)

    if len(condition_names) >= 2:
        c1, c2 = condition_names[:2]

        deg["diff"] = deg[c2] - deg[c1]
        deg["abs_diff"] = deg["diff"].abs()

        top_df = deg.sort_values("abs_diff", ascending=False).head(5)
        max_diff = top_df["abs_diff"].max() if not top_df.empty else 1

        top_degs = [
            {
                "gene": gene,
                "value": round(row["diff"], 2),
                "width": round((row["abs_diff"] / max_diff) * 100, 0),
                "direction": "up" if row["diff"] > 0 else "down",
            }
            for gene, row in top_df.iterrows()
        ]
    else:
        top_degs = []

    # ---- PATHWAY ENRICHMENT ----
    enriched_pathways = []
    enrichment_chart = None

    gene_list = [gene["gene"] for gene in top_degs]

    if len(gene_list) >= 3:
        try:
            enr = gp.enrichr(
                gene_list=gene_list,
                gene_sets="KEGG_2021_Human",
                organism="Human",
                outdir=None,
            )

            if enr.results is not None and not enr.results.empty:
                enrichment_df = enr.results.head(8)

                enriched_pathways = [
                    {
                        "pathway": row["Term"],
                        "p_value": round(row["P-value"], 5),
                        "score": round(row["Combined Score"], 2),
                        "width": (
                            round(
                                (
                                    row["Combined Score"]
                                    / enrichment_df["Combined Score"].max()
                                )
                                * 100,
                                0,
                            )
                            if enrichment_df["Combined Score"].max()
                            else 0
                        ),
                    }
                    for _, row in enrichment_df.iterrows()
                ]

                enrichment_chart_path = os.path.join(CHART_FOLDER, "enrichment_chart.png")

                plt.figure(figsize=(8, 5), facecolor="#0B1E2D")
                ax = plt.gca()
                ax.set_facecolor("#071421")

                plot_df = enrichment_df.sort_values("Combined Score", ascending=True)

                plt.barh(plot_df["Term"], plot_df["Combined Score"], color="#00C896")

                plt.title("Top Enriched Pathways", color="white")
                plt.xlabel("Combined Score", color="white")
                plt.ylabel("Pathway", color="white")

                plt.xticks(color="white")
                plt.yticks(color="white", fontsize=8)

                plt.grid(axis="x", alpha=0.15)
                plt.tight_layout()
                plt.savefig(enrichment_chart_path, transparent=False)
                plt.close()

                enrichment_chart = "charts/enrichment_chart.png"

        except Exception as e:
            print("Pathway enrichment failed:", e)

    latest_results = {
        "uploaded_filename": filename,
        "num_rows": num_rows,
        "num_genes": num_genes,
        "num_conditions": num_conditions,
        "avg_expression": avg_expression,
        "table": table,
        "chart": "charts/chart.png",
        "condition_chart": "charts/condition_chart.png",
        "heatmap_chart": "charts/heatmap_chart.png",
        "pca_chart": pca_chart,
        "umap_chart": umap_chart,
        "top_degs": top_degs,
        "enriched_pathways": enriched_pathways,
        "enrichment_chart": enrichment_chart,
        "expression_bars": expression_bars,
        "heatmap_cells": heatmap_cells,
        "condition_bars": condition_bars,
        "top_variable_genes": top_variable_genes,
        "qc_metrics": {
            "rows": num_rows,
            "missing_removed": missing_removed,
            "genes": num_genes,
            "conditions": num_conditions,
            "completeness": completeness,
        },
    }

    return render_template("index.html", **latest_results)


def prepare_uploaded_dataframe(df):
    """
    Supports two formats:

    1. Long format:
       gene, expression, condition

    2. Wide format:
       sample_id, condition, BRCA1, EGFR, TP53, ...
    """

    df.columns = df.columns.str.strip().str.lower()

    # Existing long-format support
    gene_col = find_column(df, ["gene", "gene_name", "symbol"])
    expr_col = find_column(
        df, ["expression", "expression_level", "tpm", "value", "count"]
    )
    cond_col = find_column(
        df, ["condition", "group", "status", "sample_type", "class", "label"]
    )

    if gene_col and expr_col and cond_col:
        return df, gene_col, expr_col, cond_col

    # Wide-format support
    possible_id_cols = ["sample", "sample_id", "patient", "patient_id", "id"]
    possible_condition_cols = [
        "condition",
        "group",
        "status",
        "sample_type",
        "class",
        "label",
    ]

    id_col = find_column(df, possible_id_cols)
    cond_col = find_column(df, possible_condition_cols)

    if cond_col:
        metadata_cols = [col for col in [id_col, cond_col] if col]

        gene_columns = [col for col in df.columns if col not in metadata_cols]

        long_df = df.melt(
            id_vars=[cond_col],
            value_vars=gene_columns,
            var_name="gene",
            value_name="expression",
        )

        return long_df, "gene", "expression", cond_col

    return None, None, None, None


@app.route("/")
def index():
    return render_template("index.html", **latest_results)


@app.route("/upload", methods=["POST"])
def upload():
    file = request.files.get("file")

    if not file or file.filename == "":
        return "No file uploaded."

    filename = secure_filename(file.filename)
    filepath = os.path.join(app.config["UPLOAD_FOLDER"], filename)
    file.save(filepath)

    df = pd.read_csv(filepath)

    prepared_df, gene_col, expr_col, cond_col = prepare_uploaded_dataframe(df)

    if prepared_df is None:
        return """
        Missing required columns.

        Supported formats:

        Long format:
        gene, expression, condition

        Wide format:
        sample_id, condition, BRCA1, EGFR, TP53, ...
        """

    return run_analysis(prepared_df, gene_col, expr_col, cond_col, filename)


@app.route("/analysis")
def analysis():
    return render_template("analysis.html", **latest_results)


@app.route("/visualizations")
def visualizations():
    return render_template("visualizations.html", **latest_results)
@app.route("/reports")
def reports():
    return render_template("reports.html", **latest_results)

@app.route("/download-report")
def download_report():
    top_genes_text = ""

    for gene in latest_results.get("top_degs", []):
        top_genes_text += f"- {gene['gene']}: {gene['value']}\n"

    report = f"""
GeneScope Analysis Report

Dataset: {latest_results.get('uploaded_filename')}
Rows/Samples: {latest_results.get('num_rows')}
Genes Detected: {latest_results.get('num_genes')}
Conditions: {latest_results.get('num_conditions')}
Average Expression: {latest_results.get('avg_expression')}

Top Differentially Expressed Genes:
{top_genes_text if top_genes_text else "No DEG results available."}

Quality Summary:
Rows Processed: {latest_results.get('qc_metrics', {}).get('rows')}
Missing Values Removed: {latest_results.get('qc_metrics', {}).get('missing_removed')}
Data Completeness: {latest_results.get('qc_metrics', {}).get('completeness')}%

Summary:
GeneScope analyzed the uploaded expression dataset and identified genes with the largest expression differences between conditions.
"""

    return (
        report,
        200,
        {
            "Content-Type": "text/plain",
            "Content-Disposition": "attachment; filename=genescope_report.txt",
        },
    )


@app.route("/datasets")
def datasets():
    return render_template("datasets.html", **latest_results)


@app.route("/gene-explorer", methods=["GET", "POST"])
def gene_explorer():
    global latest_df

    selected_gene = None
    gene_table = None
    gene_summary = None

    if request.method == "POST":
        gene_query = request.form.get("gene_query", "").strip()

        if latest_df is not None and gene_query:
            df_filtered = latest_df[
                latest_df["gene"].astype(str).str.lower() == gene_query.lower()
            ]

            if not df_filtered.empty:
                selected_gene = gene_query.upper()
                gene_table = df_filtered.head(20).to_html(classes="data", index=False)

                cond_avg = (
                    df_filtered.groupby("condition")["expression_level"]
                    .mean()
                    .sort_values(ascending=False)
                )

                max_val = cond_avg.max() if not cond_avg.empty else 1

                gene_summary = [
                    {
                        "condition": condition,
                        "value": round(value, 2),
                        "width": round((value / max_val) * 100, 0) if max_val else 0,
                    }
                    for condition, value in cond_avg.items()
                ]

    return render_template(
        "gene_explorer.html",
        selected_gene=selected_gene,
        gene_table=gene_table,
        gene_summary=gene_summary,
        **latest_results
    )


@app.route("/pca-umap")
def pca_umap():
    return render_template("pca_umap.html", **latest_results)


@app.route("/pathways")
def pathways():
    return render_template("pathways.html", **latest_results)


@app.route("/raw-data")
def raw_data():
    raw_table = None

    if latest_df is not None:
        raw_table = latest_df.to_html(classes="data", index=False)

    return render_template("raw_data.html", raw_table=raw_table, **latest_results)


@app.route("/load-tcga-demo")
def load_tcga_demo():
    data_path = os.path.join(app.config["UPLOAD_FOLDER"], "data.csv")
    labels_path = os.path.join(app.config["UPLOAD_FOLDER"], "labels.csv")

    if not os.path.exists(data_path) or not os.path.exists(labels_path):
        return "Please place data.csv and labels.csv inside the uploads folder."

    data_df = pd.read_csv(data_path)
    labels_df = pd.read_csv(labels_path)

    data_df = data_df.rename(columns={"Unnamed: 0": "sample_id"})
    labels_df = labels_df.rename(
        columns={"Unnamed: 0": "sample_id", "Class": "condition"}
    )

    merged_df = data_df.merge(labels_df, on="sample_id")

    gene_cols = [
        col for col in merged_df.columns if col not in ["sample_id", "condition"]
    ]

    # Pick top 100 most variable genes so the app does not lag
    top_genes = (
        merged_df[gene_cols].var().sort_values(ascending=False).head(100).index.tolist()
    )

    reduced_df = merged_df[["sample_id", "condition"] + top_genes]

    long_df = reduced_df.melt(
        id_vars=["sample_id", "condition"],
        value_vars=top_genes,
        var_name="gene",
        value_name="expression",
    )

    return run_analysis(
        long_df, "gene", "expression", "condition", "TCGA Cancer RNA-Seq Demo"
    )


if __name__ == "__main__":
    app.run(debug=True)
