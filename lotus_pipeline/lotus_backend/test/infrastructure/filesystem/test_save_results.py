import pytest
import scanpy as sc
from pathlib import Path
from infrastructure.filesystem.storage import AssetStorage
from infrastructure.workspace_context import workspace_path_manager
from lotus.preprocessing import neighbors
from lotus.clustering import leiden
from lotus.io import read_h5ad

CURRENT_FILE = Path(__file__).resolve()

# Manually define the user data root path
PROJECT_ROOT = CURRENT_FILE.parents[4]
TEST_DATA_ROOT = PROJECT_ROOT / "data_root"


def test_scanpy_clustering_persistence():
    """
    Verifies that Scanpy clustering results can be correctly persisted 
    to a specific relative path using AssetStorage.
    """

    # Step 1: Prepare Data (Mock Data)
    print("\n[Step 1] Generating dummy single-cell data...")
    # Generate random blobs: 100 cells x 50 genes
    adata = sc.datasets.blobs(n_variables=50, n_observations=100)

    # Step 2: Run Clustering Algorithm (In-Memory)
    print("[Step 2] Running Leiden clustering...")
    # Compute neighbors map
    neighbors(adata, n_neighbors=10, n_pcs=10)

    # Run Leiden clustering
    leiden(adata, resolution=0.5)

    # Verify clustering results exist in memory
    assert 'leiden' in adata.obs, "Leiden key missing in adata.obs after clustering"
    print(f"         Clusters found: {adata.obs['leiden'].unique().tolist()}")

    # Step 3: Persist to Storage (File I/O)
    storage = AssetStorage()

    # The expected relative path key required by the business logic
    relative_key = "cluster/test/test01.h5ad"

    print(f"[Step 3] Saving to workspace key: '{relative_key}'")
    # This method should automatically handle absolute path resolution
    saved_path_key = storage.save_anndata(adata, relative_key)

    # Step 4: Verification
    # 4.1 Verify the returned key matches the input key
    assert saved_path_key == "cluster/test/test01.h5ad"

    # 4.2 Verify the physical file exists on disk
    expected_abs_path = TEST_DATA_ROOT / "cluster" / "test" / "test01.h5ad"
    print(f"         Checking physical path: {expected_abs_path}")
    assert expected_abs_path.exists(), "Physical .h5ad file was not created!"

    # 4.3 Verify file size (ensure it's not empty)
    file_size = expected_abs_path.stat().st_size
    print(f"         File size: {file_size / 1024:.2f} KB")
    assert file_size > 0, "Created file is empty!"

    # 4.4 Verify Data Integrity (Reloading)
    print("[Step 4] Reloading data to verify integrity...")
    # Read the file back from the disk
    loaded_adata = read_h5ad(expected_abs_path)

    # Check if the 'leiden' clustering result is preserved
    assert 'leiden' in loaded_adata.obs, "Persisted file lost the 'leiden' clustering info!"
    assert loaded_adata.n_obs == 100, "Cell count mismatch!"

    print("TEST PASSED: Data successfully clustered, saved, and reloaded.")


def test_storage_load_capability():
    """
    Integration Test: Load Flow
    Verifies that AssetStorage can correctly load an existing .h5ad file
    using its relative key.
    """

    # Step 1: Setup - Create a file to load
    # We create a dummy file first so we have something to test the 'load' function with.
    print("Step 1: Pre-creating a dummy file for loading test...")
    dummy_adata = sc.datasets.blobs(n_variables=20, n_observations=50)

    # We use AssetStorage to save it first (trusting Test 1 passed),
    # or we could manually write it to disk. Let's use storage for consistency.
    storage = AssetStorage()
    target_load_key = "cluster/test/test_load.h5ad"
    storage.save_anndata(dummy_adata, target_load_key)

    # Step 2: Execute Load
    print(f"Step 2: Attempting to load key: '{target_load_key}'")
    # This is the core function we are testing
    loaded_adata = storage.load_anndata(target_load_key)

    # Step 3: Verification
    print("Step 3: Verifying loaded data integrity...")

    # 3.1 Check object type
    assert loaded_adata is not None, "Loaded object is None!"

    # 3.2 Check dimensions (Shape)
    # Original was (50, 20), loaded should be (50, 20)
    print(f"         Original shape: {dummy_adata.shape}")
    print(f"         Loaded shape:   {loaded_adata.shape}")
    assert loaded_adata.n_obs == 50, "Cell count mismatch in loaded data"
    assert loaded_adata.n_vars == 20, "Gene count mismatch in loaded data"

    # 3.3 Check if it matches the original data structure
    # (Optional: Check if .obs indices match)
    assert dummy_adata.obs_names[0] == loaded_adata.obs_names[0], "Index mismatch"

    print("TEST 2 PASSED: Load functionality works.")
