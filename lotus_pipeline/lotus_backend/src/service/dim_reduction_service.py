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
            target_branch="HVG Selected"
        )

        # Load Data
        adata = self.storage.load_anndata(source_snapshot.snapshot_path)

        # Execute Scanpy PCA
        # 'arpack' is efficient for finding a few components from large matrices
        preprocessing.pca(
            adata,
            n_comps=request.n_comps,
            svd_solver=request.svd_solver
        )

        # Generate Outputs

        # A. Elbow Plot Data (Variance Ratio)
        # This list allows the frontend to render an interactive line chart
        variance_ratios = adata.uns['pca']['variance_ratio'].tolist()

        # B. Static Scatter Plot (PC1 vs PC2)
        visualization.pca_scatter(adata, show=False)
        fig = plt.gcf()

        scatter_filename = f"pca_scatter_{uuid.uuid4().hex[:6]}.pdf"
        scatter_path = self.storage.save_file(
            fig,
            get_project_relative(request.project_id, DIMRED_SUBSPACE, scatter_filename)
        )
        plt.close(fig)

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
            thumbnail_json={"pca_scatter": scatter_path},
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
            target_branch="PCA Computed"
        )

        # Load Data
        adata = self.storage.load_anndata(source_snapshot.snapshot_path)

        # Compute Neighbors
        # n_pcs is the critical decision point from the Elbow Plot
        preprocessing.neighbors(
            adata,
            n_neighbors=request.n_neighbors,
            n_pcs=request.n_pcs
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
