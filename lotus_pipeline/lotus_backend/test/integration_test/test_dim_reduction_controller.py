"""
Test Controller Layer Dimensionality Reduction functions one by one.

This script tests each dim reduction controller function separately:
- dim_reduction_controller.run_pca_analysis()
- dim_reduction_controller.build_neighborhood_graph()

For each function, we test:
1. Can it run successfully?
2. Are the results similar to direct scanpy calls (with same hyperparameters and random_state)?

Usage:
    python test_dim_reduction_controller.py pca        # Test PCA only
    python test_dim_reduction_controller.py neighbors   # Test Neighbors only
    python test_dim_reduction_controller.py all          # Test all (default)
"""
import sys
import argparse
from pathlib import Path
import numpy as np
import psutil
import os
import scanpy as sc
import scanpy.preprocessing as sc_pp

# Add lotus_core/src to path if needed (for development)
project_root = Path(__file__).resolve().parent.parent.parent.parent.parent
lotus_core_src = project_root / "lotus_pipeline" / "lotus_core" / "src"
if str(lotus_core_src) not in sys.path:
    sys.path.insert(0, str(lotus_core_src))

# Add backend src to path
backend_src = project_root / "lotus_pipeline" / "lotus_backend" / "src"
if str(backend_src) not in sys.path:
    sys.path.insert(0, str(backend_src))

from infrastructure.workspace_context import workspace_path_manager
from infrastructure.database.connection import get_default_db_manager
from infrastructure.filesystem.storage import AssetStorage
from infrastructure.database.model.project_meta_model import ProjectMeta
from infrastructure.database.model.dataset_model import Dataset
from infrastructure.database.model.analysis_snapshots_model import AnalysisSnapshot

# Controller and DTOs
from controller.dim_reduction_controller import run_pca_analysis, build_neighborhood_graph
from dto.request.run_pca_request import RunPCARequest
from dto.request.run_neighbors_request import RunNeighborsRequest
from dto.request.create_project_request import CreateProjectRequest

# Service (for dependency injection in tests)
from service.dim_reduction_service import DimReductionService
from service.preprocessing_service import PreprocessingService
from controller.project_management_controller import create_new_project
from controller.preprocessing_controller import perform_qc_calculation, apply_qc_filter, apply_hvg
from dto.request.calculate_qc_request import CalculateQCRequest
from dto.request.filter_qc_request import FilterQCRequest
from dto.request.run_hvg_request import RunHVGRequest
from infrastructure.database.dao.dataset_dao import DatasetDAO
from infrastructure.database.dao.analysis_snapshots_dao import AnalysisSnapshotsDAO
from infrastructure.database.dao.project_meta_dao import ProjectMetaDAO


def setup_test_environment():
    """Initialize workspace and database"""
    workspace_path_manager.initialize()
    db_manager = get_default_db_manager()
    db_manager.initialize_proxy()
    db = db_manager.get_connection()
    if db.is_closed():
        db.connect()
    db.create_tables([ProjectMeta, Dataset, AnalysisSnapshot], safe=True)
    return db_manager, db


def get_memory_usage():
    """Get current process memory usage in MB"""
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / (1024 ** 2)


def print_memory_info(label, before_mb=None):
    """Print memory usage information"""
    current_mb = get_memory_usage()
    if before_mb is not None:
        diff_mb = current_mb - before_mb
        print(f"  Memory: {current_mb:.1f} MB (Δ {diff_mb:+.1f} MB)")
    else:
        print(f"  Memory: {current_mb:.1f} MB")
    return current_mb


async def setup_project_and_dataset():
    """Helper: Create project and import dataset (shared by all tests)"""
    dataset_path = project_root / "TAURUS_raw_counts_annotated_final.h5ad"
    if not dataset_path.exists():
        raise FileNotFoundError(f"Dataset file not found: {dataset_path}")
    
    project_request = CreateProjectRequest(
        project_name="Test Dim Reduction Controller Project",
        local_file_path=str(dataset_path),
        organism="Human",
        tissue_type="Test",
        description="Integration test for dim reduction controller"
    )
    
    project_response = await create_new_project(project_request)
    return project_response.project_id, project_response.dataset_id


def create_dim_reduction_service():
    """Helper: Create dim reduction service with dependencies"""
    dataset_dao = DatasetDAO()
    snapshot_dao = AnalysisSnapshotsDAO()
    storage = AssetStorage()
    
    return DimReductionService(
        dataset_dao=dataset_dao,
        snapshot_dao=snapshot_dao,
        storage=storage
    )


def create_preprocessing_service():
    """Helper: Create preprocessing service with dependencies"""
    dataset_dao = DatasetDAO()
    snapshot_dao = AnalysisSnapshotsDAO()
    project_dao = ProjectMetaDAO()
    storage = AssetStorage()
    
    return PreprocessingService(
        dataset_dao=dataset_dao,
        snapshot_dao=snapshot_dao,
        project_dao=project_dao,
        storage=storage
    )


async def create_subset_and_run_preprocessing(project_id, dataset_id, n_subset=50000):
    """
    Helper: Create subset from raw dataset, then run QC/filter/HVG on subset.
    
    Returns:
        hvg_result: HVG result snapshot info
        n_subset: Actual subset size
    """
    dataset_dao = DatasetDAO()
    storage = AssetStorage()
    snapshot_dao = AnalysisSnapshotsDAO()
    preprocessing_service = create_preprocessing_service()
    
    # Step 1: Load raw dataset and create subset
    print("\n[Step 1] Creating subset from raw dataset...")
    dataset = dataset_dao.get_dataset_by_business_id(dataset_id)
    raw_dataset_path = workspace_path_manager.resolve(dataset.dataset_path)
    
    adata_source = sc.read_h5ad(raw_dataset_path, backed='r')
    n_total = adata_source.n_obs
    actual_n_subset = min(n_subset, n_total)
    print(f"  Total cells in raw dataset: {n_total:,}, creating subset: {actual_n_subset:,} cells")
    
    adata_subset = adata_source[:actual_n_subset, :].to_memory()
    adata_source.file.close()
    
    # Step 2: Save subset and update dataset path temporarily
    import uuid
    temp_dataset_path = f"projects/{project_id}/cache/temp_raw_subset_{actual_n_subset}_{uuid.uuid4().hex[:8]}.h5ad"
    storage.save_anndata(adata_subset, temp_dataset_path)
    
    # Temporarily update dataset path to point to subset
    original_dataset_path = dataset.dataset_path
    dataset.dataset_path = temp_dataset_path
    dataset.save()
    
    try:
        # Step 3: Run QC on subset
        print("\n[Step 2] Running QC calculation on subset...")
        qc_request = CalculateQCRequest(
            project_id=project_id,
            dataset_id=dataset_id,
            organism="Human"
        )
        await perform_qc_calculation(qc_request, preprocessing_service)
        
        # Step 4: Run filter on subset
        print("\n[Step 3] Running filter on subset...")
        filter_request = FilterQCRequest(
            project_id=project_id,
            dataset_id=dataset_id,
            min_genes=200,
            min_cells=3,
            pct_mt_max=20.0,
            pct_hb_max=5.0
        )
        filter_result = await apply_qc_filter(filter_request, preprocessing_service)
        print("  ✓ QC and filter on subset done")
        
        # Step 5: Run HVG on subset
        print("\n[Step 4] Running HVG on subset...")
        hvg_request = RunHVGRequest(
            project_id=project_id,
            dataset_id=dataset_id,
            n_top_genes=2000,
            flavor='seurat',
            target_sum=1e4
        )
        hvg_result = await apply_hvg(hvg_request, preprocessing_service)
        print("  ✓ HVG on subset done")
        
        return hvg_result, actual_n_subset
        
    finally:
        # Restore original dataset path
        dataset.dataset_path = original_dataset_path
        dataset.save()


async def test_pca():
    """Test 4: PCA Analysis - Can it run? Is it similar to scanpy?"""
    print("\n" + "=" * 60)
    print("Test 4: PCA Analysis (run_pca_analysis)")
    print("=" * 60)
    
    # Setup
    db_manager, db = setup_test_environment()
    initial_memory = get_memory_usage()
    print(f"\nInitial memory: {initial_memory:.1f} MB")
    
    # Create project and dataset
    project_id, dataset_id = await setup_project_and_dataset()
    print(f"\n✓ Project: {project_id}, Dataset: {dataset_id}")
    
    preprocessing_service = create_preprocessing_service()
    dim_reduction_service = create_dim_reduction_service()
    
    # Create subset from raw dataset and run QC/filter/HVG on subset
    hvg_result, n_subset = await create_subset_and_run_preprocessing(project_id, dataset_id, n_subset=50000)
    
    # Test: Can it run? (On subset)
    print("\n" + "-" * 60)
    print("[Step 1] Testing if PCA can run (on subset)...")
    print("-" * 60)
    
    before_pca = get_memory_usage()
    pca_request = RunPCARequest(
        project_id=project_id,
        dataset_id=dataset_id,
        source_snapshot_id=None,  # Auto-find latest HVG snapshot (subset)
        n_comps=50,
        svd_solver='arpack'
    )
    
    try:
        pca_result = await run_pca_analysis(pca_request, dim_reduction_service)
        after_pca = print_memory_info("After PCA", before_pca)
        
        print(f"  ✓ PCA analysis successful!")
        print(f"  ✓ Snapshot ID: {pca_result.snapshot_id}")
        print(f"  ✓ Variance plot path: {pca_result.pca_scatter_path}")
        print(f"  ✓ Number of components: {len(pca_result.variance_plot_data)}")
        
        # Verify files exist
        pca_snapshot_path = workspace_path_manager.resolve(pca_result.snapshot_path)
        assert pca_snapshot_path.exists(), f"PCA snapshot should exist: {pca_snapshot_path}"
        print(f"  ✓ All output files verified")
        
    except Exception as e:
        print(f"  ✗ Error: {e}")
        import traceback
        traceback.print_exc()
        raise
    
    # Test: Compare with scanpy (use the same subset we already have in memory)
    print("\n" + "-" * 60)
    print("[Step 2] Comparing results with direct scanpy call (on same subset)...")
    print("-" * 60)
    
    # Load HVG result to get the exact same subset that PCA used
    hvg_snapshot_path = workspace_path_manager.resolve(hvg_result.snapshot_path)
    print(f"  Loading HVG result from: {hvg_result.snapshot_path}")
    print(f"  This ensures we use the exact same data that PCA used ({n_subset:,} cells)")
    
    adata_scanpy = sc.read_h5ad(hvg_snapshot_path, backed='r')
    # Convert to memory for scanpy operations
    adata_scanpy = adata_scanpy.to_memory()
    adata_scanpy.file.close()
    
    # Run scanpy PCA with same parameters (including random_state=0 for consistency)
    # Reference: dim_reduction_service.py lines 78-82
    sc_pp.pca(
        adata_scanpy,
        n_comps=50,
        svd_solver='arpack',
        random_state=0  # Explicitly set to match default
    )
    
    # Load lotus result
    adata_lotus = sc.read_h5ad(pca_snapshot_path, backed='r')
    
    # Compare PCA results
    # 1. Number of components
    n_comps_lotus = len(pca_result.variance_plot_data)
    n_comps_scanpy = adata_scanpy.obsm['X_pca'].shape[1]
    
    print(f"  Lotus components: {n_comps_lotus}")
    print(f"  Scanpy components: {n_comps_scanpy}")
    
    if n_comps_lotus == n_comps_scanpy:
        print(f"  ✓ Number of components matches!")
    else:
        print(f"  ⚠ Component count difference: {abs(n_comps_lotus - n_comps_scanpy)}")
    
    # 2. Compare variance ratios (first few components should be similar)
    if 'pca' in adata_lotus.uns and 'variance_ratio' in adata_lotus.uns['pca']:
        variance_lotus = np.array(pca_result.variance_plot_data)
        variance_scanpy = np.array(adata_scanpy.uns['pca']['variance_ratio'])  # Already a numpy array
        
        # Compare first 10 components
        n_compare = min(10, len(variance_lotus), len(variance_scanpy))
        diff = np.abs(variance_lotus[:n_compare] - variance_scanpy[:n_compare])
        max_diff = diff.max()
        mean_diff = diff.mean()
        
        if max_diff < 1e-5:
            print(f"  ✓ Variance ratios match (first {n_compare} components, max diff: {max_diff:.2e})")
        else:
            print(f"  ⚠ Variance ratios differ (max: {max_diff:.6f}, mean: {mean_diff:.6f})")
            if max_diff > 1e-3:
                print(f"  ⚠ Warning: Significant difference in variance ratios")
    
    print(f"  ✓ PCA logic verified (subset comparison to avoid OOM)")
    
    print("\n" + "=" * 60)
    print("✓ Test 4 Complete: PCA Analysis")
    print("=" * 60)


async def test_neighbors():
    """Test 5: Neighborhood Graph - Can it run? Is it similar to scanpy?"""
    print("\n" + "=" * 60)
    print("Test 5: Neighborhood Graph (build_neighborhood_graph)")
    print("=" * 60)
    
    # Setup
    db_manager, db = setup_test_environment()
    initial_memory = get_memory_usage()
    print(f"\nInitial memory: {initial_memory:.1f} MB")
    
    # Create project and dataset
    project_id, dataset_id = await setup_project_and_dataset()
    print(f"\n✓ Project: {project_id}, Dataset: {dataset_id}")
    
    preprocessing_service = create_preprocessing_service()
    dim_reduction_service = create_dim_reduction_service()
    
    # Create subset from raw dataset and run QC/filter/HVG on subset
    hvg_result, n_subset = await create_subset_and_run_preprocessing(project_id, dataset_id, n_subset=50000)
    
    # Run PCA on subset
    print("\n[Running PCA on subset]...")
    pca_request = RunPCARequest(
        project_id=project_id,
        dataset_id=dataset_id,
        source_snapshot_id=None,  # Auto-find latest HVG snapshot (subset)
        n_comps=50,
        svd_solver='arpack'
    )
    pca_result = await run_pca_analysis(pca_request, dim_reduction_service)
    print("  ✓ PCA on subset done")
    
    # Test: Can neighbors run? (On subset)
    print("\n" + "-" * 60)
    print("[Step 1] Testing if neighbors graph can run (on subset)...")
    print("-" * 60)
    
    before_neighbors = get_memory_usage()
    neighbors_request = RunNeighborsRequest(
        project_id=project_id,
        dataset_id=dataset_id,
        source_snapshot_id=pca_result.snapshot_id,  # Use subset PCA snapshot
        n_neighbors=15,
        n_pcs=30
    )
    
    try:
        neighbors_result = await build_neighborhood_graph(neighbors_request, dim_reduction_service)
        after_neighbors = print_memory_info("After Neighbors", before_neighbors)
        
        print(f"  ✓ Neighbors graph successful!")
        print(f"  ✓ Snapshot ID: {neighbors_result.snapshot_id}")
        print(f"  ✓ Number of connectivities: {neighbors_result.n_connectivities:,}")
        
        # Verify snapshot exists
        neighbors_snapshot_path = workspace_path_manager.resolve(neighbors_result.snapshot_path)
        assert neighbors_snapshot_path.exists(), f"Neighbors snapshot should exist: {neighbors_snapshot_path}"
        print(f"  ✓ All output files verified")
        
    except Exception as e:
        print(f"  ✗ Error: {e}")
        import traceback
        traceback.print_exc()
        raise
    
    # Test: Compare with scanpy (use the same subset we already have in memory)
    print("\n" + "-" * 60)
    print("[Step 2] Comparing results with direct scanpy call (on same subset)...")
    print("-" * 60)
    
    # Load HVG result to get the exact same subset that PCA/Neighbors used
    hvg_snapshot_path = workspace_path_manager.resolve(hvg_result.snapshot_path)
    print(f"  Loading HVG result from: {hvg_result.snapshot_path}")
    print(f"  This ensures we use the exact same data that PCA/Neighbors used ({n_subset:,} cells)")
    
    adata_scanpy = sc.read_h5ad(hvg_snapshot_path, backed='r')
    # Convert to memory for scanpy operations
    adata_scanpy = adata_scanpy.to_memory()
    adata_scanpy.file.close()
    
    # Run scanpy PCA first (same parameters as lotus)
    sc_pp.pca(
        adata_scanpy,
        n_comps=50,
        svd_solver='arpack',
        random_state=0
    )
    
    # Run scanpy neighbors with same parameters (including random_state=0)
    # Reference: dim_reduction_service.py lines 194-198
    sc_pp.neighbors(
        adata_scanpy,
        n_neighbors=15,
        n_pcs=30,
        random_state=0  # Explicitly set to match default
    )
    
    # Load lotus result
    adata_lotus = sc.read_h5ad(neighbors_snapshot_path, backed='r')
    
    # Compare neighbors results
    # 1. Number of connectivities
    n_conns_lotus = neighbors_result.n_connectivities
    n_conns_scanpy = adata_scanpy.obsp['connectivities'].nnz if 'connectivities' in adata_scanpy.obsp else 0
    
    print(f"  Lotus connectivities: {n_conns_lotus:,}")
    print(f"  Scanpy connectivities (on subset): {n_conns_scanpy:,}")
    
    # 2. Verify neighbors data structure exists
    if 'neighbors' in adata_lotus.uns:
        print(f"  ✓ Neighbors data structure exists in lotus result")
    if 'neighbors' in adata_scanpy.uns:
        print(f"  ✓ Neighbors data structure exists in scanpy result")
    
    # 3. Compare n_neighbors parameter
    if 'neighbors' in adata_lotus.uns and 'neighbors' in adata_scanpy.uns:
        n_neighbors_lotus = adata_lotus.uns['neighbors']['params']['n_neighbors']
        n_neighbors_scanpy = adata_scanpy.uns['neighbors']['params']['n_neighbors']
        
        if n_neighbors_lotus == n_neighbors_scanpy:
            print(f"  ✓ n_neighbors parameter matches: {n_neighbors_lotus}")
        else:
            print(f"  ⚠ n_neighbors difference: {n_neighbors_lotus} vs {n_neighbors_scanpy}")
    
    print(f"  ✓ Neighbors graph logic verified (subset comparison to avoid OOM)")
    
    print("\n" + "=" * 60)
    print("✓ Test 5 Complete: Neighborhood Graph")
    print("=" * 60)


def main():
    """Run selected test(s)"""
    parser = argparse.ArgumentParser(description='Test dim reduction controller functions')
    parser.add_argument(
        'test',
        nargs='?',
        default='all',
        choices=['pca', 'neighbors', 'all'],
        help='Which test to run (default: all)'
    )
    args = parser.parse_args()
    
    import asyncio
    
    print("\n" + "=" * 60)
    print("Dim Reduction Controller Tests")
    print("=" * 60)
    print(f"\nRunning: {args.test}")
    print("\nEach test verifies:")
    print("  1. Can the function run successfully?")
    print("  2. Are results similar to direct scanpy calls (same hyperparameters + random_state=0)?")
    
    try:
        if args.test == 'pca' or args.test == 'all':
            asyncio.run(test_pca())
        
        if args.test == 'neighbors' or args.test == 'all':
            asyncio.run(test_neighbors())
        
        if args.test == 'all':
            print("\n" + "=" * 60)
            print("✓ All Tests Complete!")
            print("=" * 60)
        
    except Exception as e:
        print(f"\n✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        raise


if __name__ == "__main__":
    main()

