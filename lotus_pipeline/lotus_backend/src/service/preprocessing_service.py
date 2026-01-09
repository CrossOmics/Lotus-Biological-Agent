import matplotlib.pyplot as plt
from typing import Optional, Dict
from fastapi import Depends
from infrastructure.filesystem.constants.filesystem_constants import QC_SUBSPACE
from lotus.preprocessing import calculate_qc_metrics
from lotus.visualization import violin, scatter
from infrastructure.database.dao.dataset_dao import DatasetDAO
from infrastructure.filesystem.storage import AssetStorage
from dto.response.qc_result_dto import QCResultDTO

from loguru import logger

from util.path_utils import get_project_relative


class PreprocessingService:
    def __init__(self, dataset_dao: DatasetDAO = Depends(), storage: AssetStorage = Depends()):
        self.dataset_dao = dataset_dao
        self.storage = storage

    def qc_calculation(
            self,
            project_id: str,
            dataset_id: str,
            organism: str,
            custom_prefix: Optional[Dict] = None
    ) -> QCResultDTO:  # Updated return type hint
        """
        Performs QC metric calculations and generates visualization assets.

        Steps:
        1. Load AnnData.
        2. Define gene prefixes based on organism (Human/Mouse).
        3. Calculate metrics via Scanpy.
        4. Save interactive JSON statistics for UI.
        5. Save static PDF plots for reports.
        """

        # Load Dataset
        dataset = self.dataset_dao.get_dataset_by_business_id(dataset_id)
        if not dataset:
            raise ValueError(f"Dataset {dataset_id} not found.")

        # Load raw data from storage
        adata = self.storage.load_anndata(dataset.dataset_path)

        # 2. Define Gene Prefixes (Logic for mt/MT detection)
        # Defaults
        mt_prefix = "MT-"
        ribo_prefix = ("RPS", "RPL")
        hb_pattern = "^HB[^(P)]"

        # Organism Logic
        if organism.lower() == "mouse":
            mt_prefix = "mt-"
            ribo_prefix = ("Rps", "Rpl")
            hb_pattern = "^Hb[^(p)]"
        elif organism.lower() == "human":
            pass  # Use defaults

        # Override with custom prefixes if provided
        if custom_prefix:
            mt_prefix = custom_prefix.get("mt", mt_prefix)
            ribo_prefix = tuple(custom_prefix.get("ribo", ribo_prefix))
            hb_pattern = custom_prefix.get("hb", hb_pattern)

        # 3. Annotate Genes & Calculate Metrics
        logger.info(f"Calculating QC with prefixes: MT={mt_prefix}, Ribo={ribo_prefix}")

        # Boolean indexing for gene groups
        adata.var['mt'] = adata.var_names.str.startswith(mt_prefix)
        adata.var['ribo'] = adata.var_names.str.startswith(ribo_prefix)
        # Using regex for hemoglobin
        adata.var['hb'] = adata.var_names.str.contains(hb_pattern, regex=True)

        # Run QC calculation
        qc_vars = ['mt', 'ribo', 'hb']

        calculate_qc_metrics(
            adata,
            qc_vars=qc_vars,
            percent_top=None,
            log1p=False,
            inplace=True
        )

        # Generate & Save UI Data (JSON)
        # Extract columns required for frontend Violin/Scatter plots
        obs_cols = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo', 'pct_counts_hb']
        # Filter columns that actually exist (in case hb is missing)
        existing_cols = [col for col in obs_cols if col in adata.obs.columns]

        qc_df = adata.obs[existing_cols].copy()
        qc_df['barcode'] = qc_df.index.astype(str)

        # Convert to list of dicts for JSON serialization
        qc_json = qc_df.to_dict(orient='records')

        # Save JSON to workspace
        json_path = self.storage.save_file(
            content=qc_json,
            relative_key=get_project_relative(project_id, QC_SUBSPACE, 'qc_metrics.json'),
            project_id=project_id
        )

        # 5. Generate & Save Static Plots (PDF)
        # Set show=False to prevent display, capture fig via plt.gcf()

        # 5.1 Violin Plot
        violin(
            adata,
            ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
            jitter=0.4,
            multi_panel=True,
            show=False
        )
        fig_violin = plt.gcf()
        violin_path = self.storage.save_file(
            fig_violin,
            get_project_relative(project_id, QC_SUBSPACE, 'qc_violin_overview.pdf')
        )
        plt.close(fig_violin)

        # 5.2 Scatter: Counts vs MT
        scatter(
            adata,
            x='total_counts',
            y='pct_counts_mt',
            show=False
        )
        fig_mt = plt.gcf()
        scatter_mt_path = self.storage.save_file(
            fig_mt,
            get_project_relative(project_id, QC_SUBSPACE, "qc_scatter_mt.pdf")
        )
        plt.close(fig_mt)

        # 5.3 Scatter: Counts vs Genes
        scatter(
            adata,
            x='total_counts',
            y='n_genes_by_counts',
            show=False
        )
        fig_genes = plt.gcf()
        scatter_genes_path = self.storage.save_file(
            fig_genes,
            get_project_relative(project_id, QC_SUBSPACE, "qc_scatter_genes.pdf")
        )
        plt.close(fig_genes)

        # Skip saving intermediate AnnData state

        # Return DTO with all paths
        return QCResultDTO(
            dataset_id=dataset_id,
            metrics_json_path=json_path,
            violin_plot_path=violin_path,
            scatter_mt_path=scatter_mt_path,
            scatter_genes_path=scatter_genes_path
        )


    