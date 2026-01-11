import uuid

import matplotlib.pyplot as plt
from fastapi import Depends

from dto.request.run_neighbors_request import RunNeighborsRequest
from dto.request.run_pca_request import RunPCARequest
from dto.response.neighbor_result_dto import NeighborsResultDTO
from dto.response.pca_result_dto import PCAResultDTO

from infrastructure.database.dao.analysis_snapshots_dao import AnalysisSnapshotsDAO
from infrastructure.database.dao.dataset_dao import DatasetDAO
from infrastructure.filesystem.storage import AssetStorage
from infrastructure.filesystem.constants.filesystem_constants import SNAPSHOTS_SUBSPACE, DIMRED_SUBSPACE
from lotus import preprocessing, visualization

from util.path_utils import get_project_relative
from util.snapshot_utils import resolve_source_snapshot, generate_snapshot_id
from loguru import logger


class DimReductionService:
    def __init__(
            self,
            dataset_dao: DatasetDAO = Depends(),
            snapshot_dao: AnalysisSnapshotsDAO = Depends(),
            storage: AssetStorage = Depends()
    ):
        self.dataset_dao = dataset_dao
        self.snapshot_dao = snapshot_dao
        self.storage = storage

    def run_pca(self, request: RunPCARequest) -> PCAResultDTO:
        """
        Executes Principal Component Analysis (PCA).

        Logic:
        1. Resolve Source: Use provided ID or find latest 'HVG Selected' snapshot.
        2. Compute PCA
        3. Extract Metrics: Variance ratio for Elbow Plot.
        4. Visualize: PC1 vs PC2 Scatter.
        5. Persist: Save snapshot and update DB.
        """

        # Resolve Source Snapshot (HVG Result)
        source_snapshot = resolve_source_snapshot(
            request.dataset_id,
            request.source_snapshot_id,
            target_branch="HVG Selected",
            snapshot_dao=self.snapshot_dao
        )
        # Load Data
        adata = self.storage.load_anndata(source_snapshot.snapshot_path)
        # dynamic n_comps adjustment
        n_samples = adata.n_obs
        n_features = adata.n_vars
        min_dim = min(n_samples, n_features)

        actual_n_comps = request.n_comps

        # ARPACK requires n_comps < min_dim; auto‑reduce for very small datasets
        if request.n_comps >= min_dim:
            new_comps = max(1, min_dim - 1)
            logger.info(
                f"Requested {request.n_comps} PCs, but data dims are ({n_samples}, {n_features}). "
                f"Auto-adjusting n_comps to {new_comps}."
            )
            actual_n_comps = new_comps

        # Too few observations/genes to run PCA
        if actual_n_comps < 1:
            raise ValueError(
                f"Dataset too small for PCA (dims: {adata.shape}). Check upstream QC filters."
            )

        # Execute Scanpy PCA
        # 'arpack' is efficient for finding a few components from large matrices
        preprocessing.pca(
            adata,
            n_comps=actual_n_comps,
            svd_solver=request.svd_solver
        )

        # Generate Outputs & Visualizations

        # Extract Data for Frontend Interactive Plot
        variance_ratios = adata.uns['pca']['variance_ratio'].tolist()

        # Static Plot 1: Elbow Plot (Variance Ratio)
        visualization.pca_variance_ratio(adata, log=True, show=False)
        fig_var = plt.gcf()
        variance_plot_path = self.storage.save_file(
            fig_var,
            get_project_relative(request.project_id, DIMRED_SUBSPACE, f"pca_variance_{uuid.uuid4().hex[:6]}.pdf")
        )
        plt.close(fig_var)

        # Static Plot 2: PCA Scatter (PC1 vs PC2)
        # We color by 'total_counts' if available to check for depth bias, otherwise just shape
        color_key = 'total_counts' if 'total_counts' in adata.obs else None

        visualization.pca(adata, color=color_key, show=False)
        fig_scatter = plt.gcf()

        scatter_filename = f"pca_scatter_{uuid.uuid4().hex[:6]}.pdf"
        scatter_path = self.storage.save_file(
            fig_scatter,
            get_project_relative(request.project_id, DIMRED_SUBSPACE, scatter_filename)
        )
        plt.close(fig_scatter)

        # Persist Result
        # Save snapshot
        snapshot_id = generate_snapshot_id()
        file_name = f"node_pca_{snapshot_id}.h5ad"
        relative_key = get_project_relative(request.project_id, SNAPSHOTS_SUBSPACE, file_name)

        saved_path = self.storage.save_anndata(adata, relative_key)

        # Update Database
        self.snapshot_dao.create_snapshot(
            dataset_id=request.dataset_id,
            snapshot_id=snapshot_id,
            branch_name="PCA Computed",
            snapshot_path=saved_path,
            parent_snapshot_id=source_snapshot.snapshot_id,  # Link lineage
            params_json=request.model_dump(),
            thumbnail_json={"pca_scatter": scatter_path, "variance_plot_path": variance_plot_path},
            user_notes=f"PCA computed with {request.n_comps} components."
        )

        return PCAResultDTO(
            variance_plot_data=variance_ratios,
            pca_scatter_path=scatter_path,
            snapshot_path=saved_path,
            snapshot_id=snapshot_id
        )

    def build_neighborhood_graph(self, request: RunNeighborsRequest) -> NeighborsResultDTO:
        """
        Constructs the neighborhood graph based on PCA results.

        Logic:
        1. Resolve Source: Use provided ID or find latest 'PCA Computed' snapshot.
        2. Compute Neighbors: sc.pp.neighbors using user-selected n_pcs.
        3. Persist: Save snapshot and update DB.
        """

        # Resolve Source Snapshot (PCA Result)
        source_snapshot = resolve_source_snapshot(
            request.dataset_id,
            request.source_snapshot_id,
            target_branch="PCA Computed",
            snapshot_dao=self.snapshot_dao
        )

        # Load Data
        adata = self.storage.load_anndata(source_snapshot.snapshot_path)

        n_obs = adata.n_obs
        actual_n_neighbors = request.n_neighbors

        if actual_n_neighbors >= n_obs:
            # Ensure there are at least 2 neighbors
            # unless the total number of cells < 3, then clustering is impossible
            new_neighbors = max(2, n_obs - 1)

            logger.warning(
                f"Requested {request.n_neighbors} neighbors, but only have {n_obs} cells. "
                f"Auto-adjusting n_neighbors to {new_neighbors}."
            )
            actual_n_neighbors = new_neighbors

        # The number of cells is too small
        if n_obs < 3:
            raise ValueError(
                f"Too few cells ({n_obs}) to build a neighborhood graph. Please check upstream QC/Filter parameters.")

        # Validate PCA dimension
        if "X_pca" not in adata.obsm:
            raise ValueError("PCA results not found. Please run PCA before neighbors.")

        available_pcs = adata.obsm["X_pca"].shape[1]
        if request.n_pcs > available_pcs:
            logger.warning(
                f"Requested {request.n_pcs} PCs, but only {available_pcs} PCs available. "
                f"Auto-adjusting n_pcs to {available_pcs}."
            )
            actual_n_pcs = available_pcs
        else:
            actual_n_pcs = request.n_pcs

        # Compute Neighbors
        preprocessing.neighbors(
            adata,
            n_neighbors=actual_n_neighbors,
            n_pcs=actual_n_pcs
        )

        # Persist Result
        snapshot_id = generate_snapshot_id()
        file_name = f"node_graph_{snapshot_id}.h5ad"
        relative_key = get_project_relative(request.project_id, SNAPSHOTS_SUBSPACE, file_name)

        saved_path = self.storage.save_anndata(adata, relative_key)

        # Helper to count connections for debugging info
        # 'connectivities' is a sparse matrix
        n_conns = adata.obsp['connectivities'].nnz if 'connectivities' in adata.obsp else 0

        # Update Database
        self.snapshot_dao.create_snapshot(
            dataset_id=request.dataset_id,
            snapshot_id=snapshot_id,
            branch_name="Neighborhood Graph",
            snapshot_path=saved_path,
            parent_snapshot_id=source_snapshot.snapshot_id,
            params_json=request.model_dump(),
            user_notes=f"Neighbors graph built (n_pcs={request.n_pcs}, k={request.n_neighbors})."
        )

        return NeighborsResultDTO(
            snapshot_path=saved_path,
            snapshot_id=snapshot_id,
            n_connectivities=n_conns
        )
