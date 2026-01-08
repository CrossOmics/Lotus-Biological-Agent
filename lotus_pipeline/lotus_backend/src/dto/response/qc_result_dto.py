from dataclasses import dataclass

@dataclass
class QCResultDTO:
    """
    Data Transfer Object for QC Calculation Results.
    Contains the dataset ID and relative paths to all generated assets.
    """
    dataset_id: str
    metrics_json_path: str      # Path to the interactive JSON for UI
    violin_plot_path: str       # Path to the static Violin PDF
    scatter_mt_path: str        # Path to the static Scatter (MT vs Counts) PDF
    scatter_genes_path: str     # Path to the static Scatter (Genes vs Counts) PDF