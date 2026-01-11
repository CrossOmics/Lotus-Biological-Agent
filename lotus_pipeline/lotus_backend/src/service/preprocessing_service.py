import uuid
from datetime import datetime

import matplotlib.pyplot as plt
from typing import Optional, Dict

from anndata import AnnData
from fastapi import Depends

from dto.request.run_hvg_request import RunHVGRequest
from dto.response.hvg_result_dto import HVGResultDTO
from infrastructure.database.dao.analysis_snapshots_dao import AnalysisSnapshotsDAO
from infrastructure.filesystem.constants.filesystem_constants import QC_SUBSPACE, SNAPSHOTS_SUBSPACE, CACHE_SUBSPACE
import lotus
from lotus.preprocessing import calculate_qc_metrics, filter_cells, filter_genes, normalize_total, log1p, scale
from lotus.visualization import violin, scatter
from infrastructure.database.dao.project_meta_dao import ProjectMetaDAO

from infrastructure.database.dao.dataset_dao import DatasetDAO
from infrastructure.filesystem.storage import AssetStorage
from dto.response.qc_result_dto import QCResultDTO
from dto.request.filter_qc_request import FilterQCRequest
from dto.response.filter_qc_response import FilterQCResponse
from loguru import logger

from util.path_utils import get_project_relative


class PreprocessingService:
    def __init__(self, dataset_dao: DatasetDAO = Depends(),
                 snapshot_dao: AnalysisSnapshotsDAO = Depends(),
                 project_dao: ProjectMetaDAO = Depends(),
                 storage: AssetStorage = Depends()):
        self.dataset_dao = dataset_dao
        self.snapshot_dao = snapshot_dao
        self.project_dao = project_dao
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

        cache_file_name = f"{dataset_id}_qc_metrics.h5ad"
        cache_key = get_project_relative(project_id, CACHE_SUBSPACE, cache_file_name)

        self._compute_qc_metrics_inplace(adata, organism, custom_prefix, cache_key)

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

    def apply_filter(self, request: FilterQCRequest) -> FilterQCResponse:
        """
        Applies QC filtering logic and persists the result.

        Strategy:
        1. Load Raw Data (Clean Slate).
        2. Attempt to merge QC metrics from the Lightweight Cache (Obs/Var only).
        3. If cache miss/invalid, fallback to re-calculating metrics.
        4. Apply filters.
        5. Save Snapshot.
        """

        # 1. Validation & Loading Raw Data
        dataset = self.dataset_dao.get_dataset_by_business_id(request.dataset_id)
        if not dataset:
            raise ValueError(f"Dataset {request.dataset_id} not found.")

        try:
            adata = self.storage.load_anndata(dataset.dataset_path)
        except FileNotFoundError:
            raise FileNotFoundError(f"Physical file missing for dataset: {dataset.dataset_path}")
        except Exception as e:
            raise RuntimeError(f"Failed to load source file: {str(e)}")

        # 2. [OPTIMIZED] Try Loading Lightweight Cache & Merge
        cache_file_name = f"{request.dataset_id}_qc_metrics.h5ad"
        cache_key = get_project_relative(request.project_id, CACHE_SUBSPACE, cache_file_name)
        # Use the generic merge method
        is_merged = self.storage.load_and_merge_anndata(adata, cache_key)
        # Verify if critical metrics exist (Business logic validation)
        metrics_ready = is_merged and ('total_counts' in adata.obs)

        # Fallback: Recalculate Metrics (If cache failed or is invalid)
        if not metrics_ready:
            logger.info("Fallback: Recalculating QC metrics...")
            # Retrieve Organism Info
            project = None
            try:
                project = self.project_dao.get_project_by_id(request.project_id)
            except:
                pass

            organism = project.organism if project else "Human"
            self._compute_qc_metrics_inplace(adata, organism, cache_key)

        # 4. Apply Filtering
        initial_cells = adata.n_obs
        logger.info(f"Applying filters. Initial cells: {initial_cells}")

        if request.min_genes > 0:
            filter_cells(adata, min_genes=request.min_genes)
        if request.min_cells > 0:
            filter_genes(adata, min_cells=request.min_cells)

        # Safe filtering now that metrics are guaranteed
        if request.max_counts is not None and 'total_counts' in adata.obs:
            adata = adata[adata.obs['total_counts'] <= request.max_counts, :]

        if request.pct_mt_max is not None and 'pct_counts_mt' in adata.obs:
            adata = adata[adata.obs['pct_counts_mt'] <= request.pct_mt_max, :]

        if request.pct_hb_max is not None and 'pct_counts_hb' in adata.obs:
            adata = adata[adata.obs['pct_counts_hb'] <= request.pct_hb_max, :]

        logger.info(f"Filtering complete. Remaining cells: {adata.n_obs}")

        # 5. Save Snapshot
        timestamp_str = datetime.now().strftime("%Y%m%d_%H%M%S")
        unique_suffix = uuid.uuid4().hex[:4]
        snapshot_id = f"snap_{timestamp_str}_{unique_suffix}"

        file_name = f"node_root_qc_{snapshot_id}.h5ad"
        relative_key = get_project_relative(request.project_id, SNAPSHOTS_SUBSPACE, file_name)

        try:
            saved_path = self.storage.save_anndata(adata, relative_key)
        except Exception as e:
            raise RuntimeError(f"Failed to save snapshot to storage: {str(e)}")

        # 6. Create DB Record
        params = request.model_dump()

        db_record = self.snapshot_dao.create_snapshot(
            dataset_id=dataset.dataset_id,
            snapshot_id=snapshot_id,
            branch_name="QC Filtered",
            snapshot_path=saved_path,
            params_json=params,
            user_notes=f"Filtered from {initial_cells} to {adata.n_obs} cells."
        )

        if not db_record:
            raise RuntimeError("Database integrity error: Failed to create snapshot record.")

        return FilterQCResponse(
            snapshot_id=snapshot_id,
            snapshot_path=saved_path,
            n_obs_remaining=adata.n_obs,
            n_vars_remaining=adata.n_vars
        )

    # Internal Methods
    def _compute_qc_metrics_inplace(
            self,
            adata: AnnData,
            organism: str,
            custom_prefix: Optional[Dict] = None,
            cache_key: str = ''
    ) -> None:
        """
        Internal helper to determine prefixes based on organism and run scanpy.pp.calculate_qc_metrics.
        Modifies the AnnData object in-place.
        """
        # cache QC Calculation score
        # Record the column names before calculation
        obs_cols_before = set(adata.obs.columns)
        var_cols_before = set(adata.var.columns)

        # 1. Define Defaults
        mt_prefix = "MT-"
        ribo_prefix = ("RPS", "RPL")
        hb_pattern = "^HB[^(P)]"

        # 2. Organism Logic
        if organism and organism.lower() == "mouse":
            mt_prefix = "mt-"
            ribo_prefix = ("Rps", "Rpl")
            hb_pattern = "^Hb[^(p)]"

        # 3. Custom Overrides
        if custom_prefix:
            mt_prefix = custom_prefix.get("mt", mt_prefix)
            ribo_prefix = tuple(custom_prefix.get("ribo", ribo_prefix))
            hb_pattern = custom_prefix.get("hb", hb_pattern)

        logger.info(f"Computing QC metrics (Organism: {organism}). Prefixes: MT={mt_prefix}")

        # 4. Annotate Genes
        adata.var['mt'] = adata.var_names.str.startswith(mt_prefix)
        adata.var['ribo'] = adata.var_names.str.startswith(ribo_prefix)
        adata.var['hb'] = adata.var_names.str.contains(hb_pattern, regex=True)

        # 5. Run Scanpy Calculation
        qc_vars = ['mt', 'ribo', 'hb']
        calculate_qc_metrics(
            adata,
            qc_vars=qc_vars,
            percent_top=None,
            log1p=False,
            inplace=True
        )

        # Newly added columns = current columns - previous columns
        obs_new_cols = list(set(adata.obs.columns) - obs_cols_before)
        var_new_cols = list(set(adata.var.columns) - var_cols_before)

        #  Minimal storage Create an AnnData object containing only incremental data
        self.storage.save_incremental_anndata(
            adata_source=adata,
            obs_cols=obs_new_cols,
            var_cols=var_new_cols,
            relative_key=cache_key
        )

    def apply_hvg(self, request: RunHVGRequest) -> HVGResultDTO:
        """
        Executes Feature Selection & Scaling (Normalization -> HVG -> Scale).

        The flow follows the standard Scanpy "Golden Order":
        1. Normalize & Log1p.
        2. Identify Highly Variable Genes (HVG).
        3. Backup 'raw' state (normalized, full gene set) for DE analysis.
        4. Subset to HVGs and Scale (Z-score) for PCA/Clustering.
        """

        # Resolve Source Data Path
        # Logic: If source_snapshot_id is provided, use it.
        # Otherwise, find the latest "QC Filtered" snapshot from DB.

        source_path = None
        hint_message = "Please run QC calculation and QC filtering before proceeding to the Highly Variable Genes step."
        if request.source_snapshot_id:
            snapshot = self.snapshot_dao.get_snapshot_by_business_id(request.source_snapshot_id)
            if not snapshot:
                raise ValueError(
                    f"Source snapshot {request.source_snapshot_id} is missing. {hint_message}")
            source_path = snapshot.snapshot_path
        else:
            # Auto-find latest QC snapshot
            latest_qc = self.snapshot_dao.get_latest_snapshot(
                request.dataset_id, branch_name="QC Filtered"
            )
            if not latest_qc:
                raise ValueError(f"No QC snapshot found. Please provide source_snapshot_id. {hint_message}")
            source_path = latest_qc.snapshot_path
            # Update the source snapshot id
            request.source_snapshot_id = latest_qc.snapshot_id

        # Load Data
        try:
            adata = self.storage.load_anndata(source_path)
        except Exception as e:
            raise RuntimeError(f"Failed to load source file: {e}")

        # Normalization & Log1p
        logger.info("Step A: Normalizing and Log-transforming...")

        # Safety check: Ensure raw counts are preserved in layers if not already
        if 'counts' not in adata.layers:
            adata.layers['counts'] = adata.X.copy()

        # Normalize total counts per cell (default 1e4)
        normalize_total(adata, target_sum=request.target_sum)

        # Logarithmize the data (log(x+1))
        log1p(adata)

        # Identify Highly Variable Genes (HVG)
        logger.info(f"Step B: Identifying top {request.n_top_genes} HVGs...")

        lotus.preprocessing.highly_variable_genes(
            adata,
            n_top_genes=request.n_top_genes,
            flavor=request.flavor,
            subset=False
        )

        # Visualization (Dispersion Plot)
        # Generate the plot
        lotus.visualization.highly_variable_genes(adata, show=False)
        fig = plt.gcf()

        # Save plot to storage
        plot_filename = "hvg_dispersion.pdf"
        hvg_plot_path = self.storage.save_file(
            fig,
            get_project_relative(request.project_id, QC_SUBSPACE, plot_filename)
        )
        plt.close(fig)

        # Backup Raw
        # We freeze the "Normalized + Full Gene" state into .raw
        # This is required for future Differential Expression (DE) analysis
        logger.info("Step D: Backing up full normalized data to .raw...")
        adata.raw = adata

        # Step E: Subset & Scale
        logger.info("Step E: Subsetting to HVGs and Scaling...")

        # 1. Physical subset: Keep only HVGs in X
        adata = adata[:, adata.var.highly_variable]

        # 2. Scale: Z-score normalization (make mean=0, std=1)
        # Note: This introduces negative numbers
        scale(adata, max_value=10)

        # Persistence (Snapshot)
        timestamp_str = datetime.now().strftime("%Y%m%d_%H%M%S")
        unique_suffix = uuid.uuid4().hex[:4]
        snapshot_id = f"snap_{timestamp_str}_{unique_suffix}"

        # Save the processed AnnData (PCA Ready)
        file_name = f"node_hvg_snap_{snapshot_id}.h5ad"
        relative_key = get_project_relative(request.project_id, SNAPSHOTS_SUBSPACE, file_name)

        try:
            saved_path = self.storage.save_anndata(adata, relative_key)
        except Exception as e:
            raise RuntimeError(f"Failed to save snapshot: {e}")

        # Create DB Record
        n_genes = int(sum(adata.var['highly_variable'])) if 'highly_variable' in adata.var else adata.n_vars

        # TODO: parent_snapshot_id=request.source_snapshot_id,  # Link lineage and change the snapshot schema
        self.snapshot_dao.create_snapshot(
            dataset_id=request.dataset_id,
            snapshot_id=snapshot_id,
            branch_name="HVG Selected",
            snapshot_path=saved_path,
            parent_snapshot_id=request.source_snapshot_id,
            params_json=request.model_dump(),
            user_notes=f"Selected {n_genes} HVGs via {request.flavor}."
        )

        return HVGResultDTO(
            snapshot_id=snapshot_id,
            snapshot_path=saved_path,
            hvg_plot_path=hvg_plot_path,
            n_genes_found=adata.n_vars,  # Should be equal to n_top_genes
            msg="HVG selection and scaling complete."
        )
