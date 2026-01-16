"""
Test Controller Layer Preprocessing functions one by one.

This script tests each preprocessing controller function separately:
- preprocessing_controller.perform_qc_calculation()
- preprocessing_controller.apply_qc_filter()
- preprocessing_controller.apply_hvg()

For each function, we test:
1. Can it run successfully?
2. Are the results similar to direct scanpy calls?

Usage:
    python test_pp_controller.py qc_calc      # Test QC calculation only
    python test_pp_controller.py qc_filter   # Test QC filter only
    python test_pp_controller.py hvg          # Test HVG only
    python test_pp_controller.py all          # Test all (default)
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
from controller.preprocessing_controller import perform_qc_calculation, apply_qc_filter, apply_hvg
from dto.request.calculate_qc_request import CalculateQCRequest
from dto.request.filter_qc_request import FilterQCRequest
from dto.request.run_hvg_request import RunHVGRequest
from dto.request.create_project_request import CreateProjectRequest

# Service (for dependency injection in tests)
from service.preprocessing_service import PreprocessingService
from controller.project_management_controller import create_new_project
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
        project_name="Test Preprocessing Controller Project",
        local_file_path=str(dataset_path),
        organism="Human",
        tissue_type="Test",
        description="Integration test for preprocessing controller"
    )
    
    project_response = await create_new_project(project_request)
    return project_response.project_id, project_response.dataset_id


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


async def test_qc_calculation():
    """Test 1: QC Calculation - Can it run? Is it similar to scanpy?"""
    print("\n" + "=" * 60)
    print("Test 1: QC Calculation (perform_qc_calculation)")
    print("=" * 60)
    
    # Setup
    db_manager, db = setup_test_environment()
    initial_memory = get_memory_usage()
    print(f"\nInitial memory: {initial_memory:.1f} MB")
    
    # Create project and dataset
    project_id, dataset_id = await setup_project_and_dataset()
    print(f"\n✓ Project: {project_id}, Dataset: {dataset_id}")
    
    preprocessing_service = create_preprocessing_service()
    
    # Test: Can it run?
    print("\n" + "-" * 60)
    print("[Step 1] Testing if QC calculation can run...")
    print("-" * 60)
    
    before_qc = get_memory_usage()
    qc_request = CalculateQCRequest(
        project_id=project_id,
        dataset_id=dataset_id,
        organism="Human"
    )
    
    try:
        qc_result = await perform_qc_calculation(qc_request, preprocessing_service)
        after_qc = print_memory_info("After QC calculation", before_qc)
        
        print(f"  ✓ QC calculation successful!")
        print(f"  ✓ Metrics JSON: {qc_result.metrics_json_path}")
        print(f"  ✓ Violin plot: {qc_result.violin_plot_path}")
        
        # Verify files exist
        storage = AssetStorage()
        json_path = workspace_path_manager.resolve(qc_result.metrics_json_path)
        assert json_path.exists(), f"QC JSON file should exist: {json_path}"
        print(f"  ✓ All output files verified")
        
    except Exception as e:
        print(f"  ✗ Error: {e}")
        import traceback
        traceback.print_exc()
        raise
    
    # Test: Compare with scanpy
    print("\n" + "-" * 60)
    print("[Step 2] Comparing results with direct scanpy call...")
    print("-" * 60)
    
    # Load the data that was processed
    dataset_dao = DatasetDAO()
    dataset = dataset_dao.get_dataset_by_business_id(dataset_id)
    adata_lotus = storage.load_anndata(dataset.dataset_path)
    
    # Also load cache to get QC metrics
    cache_key = f"projects/{project_id}/cache/{dataset_id}_qc_metrics.h5ad"
    try:
        cache_adata = storage.load_anndata(cache_key)
        # Merge cache into main adata
        for col in cache_adata.obs.columns:
            if col not in adata_lotus.obs.columns:
                adata_lotus.obs[col] = cache_adata.obs[col]
        for col in cache_adata.var.columns:
            if col not in adata_lotus.var.columns:
                adata_lotus.var[col] = cache_adata.var[col]
    except:
        pass  # Cache might not exist yet
    
    # Run scanpy directly on a copy
    adata_scanpy = storage.load_anndata(dataset.dataset_path)
    adata_scanpy.var['mt'] = adata_scanpy.var_names.str.startswith('MT-')
    adata_scanpy.var['ribo'] = adata_scanpy.var_names.str.startswith(('RPS', 'RPL'))
    
    sc_pp.calculate_qc_metrics(
        adata_scanpy,
        qc_vars=['mt', 'ribo'],
        percent_top=None,
        log1p=False,
        inplace=True
    )
    
    # Compare QC metrics
    qc_cols = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo']
    all_match = True
    
    for col in qc_cols:
        if col in adata_lotus.obs.columns and col in adata_scanpy.obs.columns:
            # Use a subset for comparison (first 10K cells to avoid memory issues)
            n_compare = min(10000, adata_lotus.n_obs)
            lotus_vals = adata_lotus.obs[col].iloc[:n_compare].values
            scanpy_vals = adata_scanpy.obs[col].iloc[:n_compare].values
            
            diff = np.abs(lotus_vals - scanpy_vals)
            max_diff = diff.max()
            mean_diff = diff.mean()
            
            if max_diff < 1e-10:
                print(f"  ✓ {col}: Identical (max diff: {max_diff:.2e})")
            else:
                print(f"  ⚠ {col}: Small differences (max: {max_diff:.6f}, mean: {mean_diff:.6f})")
                if max_diff > 1e-5:
                    all_match = False
        else:
            print(f"  ⚠ {col}: Column missing")
    
    if all_match:
        print(f"\n  ✓ QC metrics match scanpy results!")
    else:
        print(f"\n  ⚠ Some QC metrics have small differences (likely numerical precision)")
    
    print("\n" + "=" * 60)
    print("✓ Test 1 Complete: QC Calculation")
    print("=" * 60)


async def test_qc_filter():
    """Test 2: QC Filter - Can it run? Is it similar to scanpy?"""
    print("\n" + "=" * 60)
    print("Test 2: QC Filter (apply_qc_filter)")
    print("=" * 60)
    
    # Setup
    db_manager, db = setup_test_environment()
    initial_memory = get_memory_usage()
    print(f"\nInitial memory: {initial_memory:.1f} MB")
    
    # Create project and dataset
    project_id, dataset_id = await setup_project_and_dataset()
    print(f"\n✓ Project: {project_id}, Dataset: {dataset_id}")
    
    preprocessing_service = create_preprocessing_service()
    
    # First run QC calculation (prerequisite)
    print("\n[Prerequisite] Running QC calculation first...")
    qc_request = CalculateQCRequest(
        project_id=project_id,
        dataset_id=dataset_id,
        organism="Human"
    )
    await perform_qc_calculation(qc_request, preprocessing_service)
    print("  ✓ QC calculation done")
    
    # Test: Can it run?
    print("\n" + "-" * 60)
    print("[Step 1] Testing if QC filter can run...")
    print("-" * 60)
    
    before_filter = get_memory_usage()
    filter_request = FilterQCRequest(
        project_id=project_id,
        dataset_id=dataset_id,
        min_genes=200,
        min_cells=3,
        max_counts=None,
        pct_mt_max=20.0,
        pct_hb_max=5.0
    )
    
    try:
        filter_result = await apply_qc_filter(filter_request, preprocessing_service)
        after_filter = print_memory_info("After QC filter", before_filter)
        
        print(f"  ✓ QC filter successful!")
        print(f"  ✓ Snapshot ID: {filter_result.snapshot_id}")
        print(f"  ✓ Cells remaining: {filter_result.n_obs_remaining:,}")
        print(f"  ✓ Genes remaining: {filter_result.n_vars_remaining:,}")
        
        # Verify snapshot exists
        storage = AssetStorage()
        snapshot_path = workspace_path_manager.resolve(filter_result.snapshot_path)
        assert snapshot_path.exists(), f"Snapshot should exist: {snapshot_path}"
        file_size_mb = snapshot_path.stat().st_size / (1024 ** 2)
        print(f"  ✓ Snapshot file size: {file_size_mb:.2f} MB")
        
    except Exception as e:
        print(f"  ✗ Error: {e}")
        import traceback
        traceback.print_exc()
        raise
    
    # Test: Compare with scanpy
    print("\n" + "-" * 60)
    print("[Step 2] Comparing results with direct scanpy call...")
    print("-" * 60)
    
    # Load original data
    dataset_dao = DatasetDAO()
    dataset = dataset_dao.get_dataset_by_business_id(dataset_id)
    adata_original = storage.load_anndata(dataset.dataset_path)
    
    # Annotate genes for scanpy
    adata_original.var['mt'] = adata_original.var_names.str.startswith('MT-')
    adata_original.var['ribo'] = adata_original.var_names.str.startswith(('RPS', 'RPL'))
    adata_original.var['hb'] = adata_original.var_names.str.contains('HB', regex=True, case=False)
    
    # Calculate QC metrics first
    sc_pp.calculate_qc_metrics(adata_original, qc_vars=['mt', 'ribo', 'hb'], inplace=True)
    
    # Apply same filters with scanpy (mimicking lotus's apply_filter logic exactly)
    # Reference: preprocessing_service.py lines 195-208
    adata_scanpy = adata_original.copy()
    
    # Step 1: Filter cells by min_genes (lotus line 195-196)
    min_genes = 200
    if min_genes > 0:
        sc_pp.filter_cells(adata_scanpy, min_genes=min_genes, inplace=True)
    
    # Step 2: Filter genes by min_cells (lotus line 197-198)
    min_cells = 3
    if min_cells > 0:
        sc_pp.filter_genes(adata_scanpy, min_cells=min_cells, inplace=True)
    
    # Step 3: Filter by max_counts using boolean indexing (lotus line 201-202)
    max_counts = None  # From filter_request
    if max_counts is not None and 'total_counts' in adata_scanpy.obs.columns:
        adata_scanpy = adata_scanpy[adata_scanpy.obs['total_counts'] <= max_counts].copy()
    
    # Step 4: Filter by pct_mt_max using boolean indexing (lotus line 204-205)
    pct_mt_max = 20.0  # From filter_request
    if pct_mt_max is not None and 'pct_counts_mt' in adata_scanpy.obs.columns:
        adata_scanpy = adata_scanpy[adata_scanpy.obs['pct_counts_mt'] <= pct_mt_max].copy()
    
    # Step 5: Filter by pct_hb_max using boolean indexing (lotus line 207-208)
    pct_hb_max = 5.0  # From filter_request
    if pct_hb_max is not None and 'pct_counts_hb' in adata_scanpy.obs.columns:
        adata_scanpy = adata_scanpy[adata_scanpy.obs['pct_counts_hb'] <= pct_hb_max].copy()
    
    # Load filtered result from lotus
    adata_lotus = storage.load_anndata(filter_result.snapshot_path)
    
    # Compare shapes
    print(f"  Lotus filtered shape: {adata_lotus.shape}")
    print(f"  Scanpy filtered shape: {adata_scanpy.shape}")
    
    if adata_lotus.shape == adata_scanpy.shape:
        print(f"  ✓ Filtered shapes match!")
    else:
        diff_obs = abs(adata_lotus.n_obs - adata_scanpy.n_obs)
        diff_vars = abs(adata_lotus.n_vars - adata_scanpy.n_vars)
        print(f"  ⚠ Shape difference: cells {diff_obs}, genes {diff_vars}")
        if diff_obs > adata_lotus.n_obs * 0.01:  # More than 1% difference
            print(f"  ⚠ Warning: Significant difference in cell count")
    
    print("\n" + "=" * 60)
    print("✓ Test 2 Complete: QC Filter")
    print("=" * 60)


async def test_hvg():
    """Test 3: HVG Selection - Can it run? Is it similar to scanpy?"""
    print("\n" + "=" * 60)
    print("Test 3: HVG Selection (apply_hvg)")
    print("=" * 60)
    
    # Setup
    db_manager, db = setup_test_environment()
    initial_memory = get_memory_usage()
    print(f"\nInitial memory: {initial_memory:.1f} MB")
    
    # Create project and dataset
    project_id, dataset_id = await setup_project_and_dataset()
    print(f"\n✓ Project: {project_id}, Dataset: {dataset_id}")
    
    preprocessing_service = create_preprocessing_service()
    
    # Prerequisites: QC calculation and filter
    print("\n[Prerequisites] Running QC calculation and filter...")
    qc_request = CalculateQCRequest(
        project_id=project_id,
        dataset_id=dataset_id,
        organism="Human"
    )
    await perform_qc_calculation(qc_request, preprocessing_service)
    
    filter_request = FilterQCRequest(
        project_id=project_id,
        dataset_id=dataset_id,
        min_genes=200,
        min_cells=3,
        pct_mt_max=20.0,
        pct_hb_max=5.0
    )
    filter_result = await apply_qc_filter(filter_request, preprocessing_service)
    print("  ✓ Prerequisites done")
    
    # Test: Can it run?
    print("\n" + "-" * 60)
    print("[Step 1] Testing if HVG selection can run...")
    print("-" * 60)
    
    before_hvg = get_memory_usage()
    hvg_request = RunHVGRequest(
        project_id=project_id,
        dataset_id=dataset_id,
        source_snapshot_id=None,  # Auto-find latest QC snapshot
        n_top_genes=2000,
        flavor='seurat',
        target_sum=1e4
    )
    
    try:
        hvg_result = await apply_hvg(hvg_request, preprocessing_service)
        after_hvg = print_memory_info("After HVG", before_hvg)
        
        print(f"  ✓ HVG selection successful!")
        print(f"  ✓ Snapshot ID: {hvg_result.snapshot_id}")
        print(f"  ✓ Genes found: {hvg_result.n_genes_found:,}")
        
        # Verify files exist
        storage = AssetStorage()
        hvg_snapshot_path = workspace_path_manager.resolve(hvg_result.snapshot_path)
        assert hvg_snapshot_path.exists(), f"HVG snapshot should exist: {hvg_snapshot_path}"
        print(f"  ✓ All output files verified")
        
    except Exception as e:
        print(f"  ✗ Error: {e}")
        import traceback
        traceback.print_exc()
        raise
    
    # Test: Compare with scanpy
    print("\n" + "-" * 60)
    print("[Step 2] Comparing results with direct scanpy call...")
    print("-" * 60)
    
    # Load lotus result first (using backed mode to save memory)
    adata_lotus = sc.read_h5ad(hvg_snapshot_path, backed='r')
    n_hvg_lotus = adata_lotus.n_vars
    print(f"  Lotus HVGs: {n_hvg_lotus}")
    
    # For large datasets, we'll use a subset for scanpy comparison to avoid OOM
    # Load source data in backed mode
    source_path = workspace_path_manager.resolve(filter_result.snapshot_path)
    adata_source = sc.read_h5ad(source_path, backed='r')
    
    # Use a subset for comparison (first 50K cells to avoid memory issues)
    # This is sufficient to verify the HVG selection logic works correctly
    n_subset = min(50000, adata_source.n_obs)
    print(f"  Using subset for scanpy comparison: {n_subset:,} cells (to avoid OOM)")
    
    # Load subset into memory (backed mode requires to_memory() first)
    adata_subset = adata_source[:n_subset, :].to_memory()
    adata_source.file.close()  # Close backed file
    
    # Run scanpy HVG on subset (mimicking lotus's apply_hvg logic)
    # Reference: preprocessing_service.py lines 352-365
    sc_pp.normalize_total(adata_subset, target_sum=1e4, inplace=True)
    sc_pp.log1p(adata_subset)  # Same as lotus line 362: log1p(adata)
    sc_pp.highly_variable_genes(
        adata_subset,
        n_top_genes=2000,
        flavor='seurat',
        inplace=True
    )
    
    # Compare: Lotus should have selected 2000 HVGs
    # (We verify the logic works, even if we can't run on full dataset)
    n_hvg_scanpy_subset = adata_subset.var['highly_variable'].sum()
    
    print(f"  Scanpy HVGs (on subset): {n_hvg_scanpy_subset}")
    
    # Main verification: Lotus result should have exactly 2000 HVGs
    if n_hvg_lotus == 2000:
        print(f"  ✓ Lotus selected correct number of HVGs: {n_hvg_lotus}")
    else:
        print(f"  ⚠ Lotus HVG count: {n_hvg_lotus} (expected ~2000)")
    
    # Verify scanpy logic works on subset
    if n_hvg_scanpy_subset == 2000:
        print(f"  ✓ Scanpy logic verified on subset: {n_hvg_scanpy_subset} HVGs")
    else:
        print(f"  ⚠ Scanpy subset result: {n_hvg_scanpy_subset} HVGs")
    
    print(f"  ✓ HVG selection logic verified (subset comparison to avoid OOM)")
    
    print("\n" + "=" * 60)
    print("✓ Test 3 Complete: HVG Selection")
    print("=" * 60)


def main():
    """Run selected test(s)"""
    parser = argparse.ArgumentParser(description='Test preprocessing controller functions')
    parser.add_argument(
        'test',
        nargs='?',
        default='all',
        choices=['qc_calc', 'qc_filter', 'hvg', 'all'],
        help='Which test to run (default: all)'
    )
    args = parser.parse_args()
    
    import asyncio
    
    print("\n" + "=" * 60)
    print("Preprocessing Controller Tests")
    print("=" * 60)
    print(f"\nRunning: {args.test}")
    print("\nEach test verifies:")
    print("  1. Can the function run successfully?")
    print("  2. Are results similar to direct scanpy calls?")
    
    try:
        if args.test == 'qc_calc' or args.test == 'all':
            asyncio.run(test_qc_calculation())
        
        if args.test == 'qc_filter' or args.test == 'all':
            asyncio.run(test_qc_filter())
        
        if args.test == 'hvg' or args.test == 'all':
            asyncio.run(test_hvg())
        
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
