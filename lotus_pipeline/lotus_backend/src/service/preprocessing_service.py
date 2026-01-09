import uuid
from datetime import datetime

import matplotlib.pyplot as plt
from typing import Optional, Dict

from anndata import AnnData
from fastapi import Depends

from infrastructure.database.dao.analysis_snapshots_dao import AnalysisSnapshotsDAO
from infrastructure.filesystem.constants.filesystem_constants import QC_SUBSPACE, SNAPSHOTS_SUBSPACE, CACHE_SUBSPACE
from lotus.preprocessing import calculate_qc_metrics, filter_cells, filter_genes
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

        # cache QC Calculation score
        # Record the column names before calculation
        obs_cols_before = set(adata.obs.columns)
        var_cols_before = set(adata.var.columns)

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
        logger.debug(f"Incremental columns to cache (Obs): {obs_new_cols}")

        #  Minimal storage Create an AnnData object containing only incremental data
        try:
            adata_cache = AnnData(
                X=None,  # Do not store X
                obs=adata.obs[obs_new_cols].copy(),  # Store only new obs columns
                var=adata.var[var_new_cols].copy()  # Store only new var columns
            )

            # Save cache
            cache_file_name = f"{dataset_id}_qc_metrics.h5ad"
            cache_relative_key = get_project_relative(
                project_id, CACHE_SUBSPACE, cache_file_name
            )
            self.storage.save_anndata(adata_cache, cache_relative_key)

        except Exception as e:
            logger.warning(f"Failed to cache incremental metrics: {e}")

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
        """

        # 1. Validation & Loading Raw Data
        dataset = self.dataset_dao.get_dataset_by_business_id(request.dataset_id)
        if not dataset:
            raise ValueError(f"Dataset {request.dataset_id} not found.")

        try:
            # Load Raw Data
            adata = self.storage.load_anndata(dataset.dataset_path)
        except FileNotFoundError:
            raise FileNotFoundError(f"Physical file missing for dataset: {dataset.dataset_path}")
        except Exception as e:
            raise RuntimeError(f"Failed to load source file: {str(e)}")

        # 2. Try Loading Lightweight Cache & Merge
        cache_file_name = f"{request.dataset_id}_qc_metrics.h5ad"
        cache_key = get_project_relative(request.project_id, CACHE_SUBSPACE, cache_file_name)

        metrics_loaded = False

        try:
            logger.info(f"Attempting to load QC Metrics Cache: {cache_key}")
            adata_cache = self.storage.load_anndata(cache_key)

            if not adata.obs_names.equals(adata_cache.obs_names):
                logger.warning("Cache index mismatch. Ignoring cache.")
            else:
                new_obs_cols = adata_cache.obs.columns.difference(adata.obs.columns)
                new_var_cols = adata_cache.var.columns.difference(adata.var.columns)

                if len(new_obs_cols) > 0:
                    adata.obs = list(
                        adata.obs.join(adata_cache.obs[new_obs_cols]))  # Join is safer/faster for aligned index
                    for col in new_obs_cols:
                        adata.obs[col] = adata_cache.obs[col]

                if len(new_var_cols) > 0:
                    for col in new_var_cols:
                        adata.var[col] = adata_cache.var[col]

                # Check validity
                if 'total_counts' in adata.obs:
                    metrics_loaded = True
                    logger.success("Lightweight Cache Hit & Merged successfully.")
                else:
                    logger.warning("Cache loaded but missing 'total_counts' column.")

        except FileNotFoundError:
            logger.info("Cache Miss. Proceeding to fallback.")
        except Exception as e:
            logger.error(f"Error merging cache: {e}. Proceeding to fallback.")

        # 3. Fallback: Recalculate Metrics
        if not metrics_loaded:
            logger.info("Fallback: Recalculating QC metrics...")

            try:
                project = self.project_dao.get_project_by_project_id(request.project_id)
            except AttributeError:
                project = None

            organism = project.organism if project else "Human"
            logger.info(f"Using organism: {organism}")

            # Define Prefixes
            mt_prefix = "mt-" if organism.lower() == "mouse" else "MT-"
            hb_pattern = "^Hb[^(p)]" if organism.lower() == "mouse" else "^HB[^(P)]"

            # Annotate
            adata.var['mt'] = adata.var_names.str.startswith(mt_prefix)
            adata.var['hb'] = adata.var_names.str.contains(hb_pattern, regex=True)

            # Calculate
            calculate_qc_metrics(
                adata,
                qc_vars=['mt', 'hb'],
                percent_top=None,
                log1p=False,
                inplace=True
            )

            # Self-Repair: Save the recalculation to cache for next time
            try:
                logger.info(f"Caching recalculated metrics to: {cache_key}")
                adata_cache_new = AnnData(
                    X=None,
                    obs=adata.obs.copy(),
                    var=adata.var.copy()
                )
                self.storage.save_anndata(adata_cache_new, cache_key)
            except Exception as e:
                logger.warning(f"Failed to save fallback cache: {e}")

        # 4. Apply Filtering (Business as usual)
        initial_cells = adata.n_obs

        if request.min_genes > 0:
            filter_cells(adata, min_genes=request.min_genes)
        if request.min_cells > 0:
            filter_genes(adata, min_cells=request.min_cells)
        if request.max_counts:
            # Safe access now
            if 'total_counts' in adata.obs:
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

        # 6. DB Record
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