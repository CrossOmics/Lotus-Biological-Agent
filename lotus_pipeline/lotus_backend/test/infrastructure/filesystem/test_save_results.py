import pytest
import scanpy as sc
from pathlib import Path
import sys
from infrastructure.filesystem.storage import AssetStorage
from infrastructure.workspace_context import workspace_path_manager
from lotus.preprocessing import neighbors
from lotus.clustering import leiden
from lotus.io import read_h5ad

CURRENT_FILE = Path(__file__).resolve()

# Manually define the user data root path
PROJECT_ROOT = CURRENT_FILE.parents[4]
TEST_DATA_ROOT = PROJECT_ROOT / "data_root"


@pytest.fixture(scope="module")
def setup_workspace():
    """
    [Fixture] Setup and Teardown for the test environment.
    1. Ensure 'data_root' exists.
    2. Initialize the Singleton Context Manager.
    3. Cleanup test files after execution.
    """
    print(f"\n[Setup] Initializing workspace at: {TEST_DATA_ROOT}")

    # 1. Ensure the physical workspace directory exists
    if not TEST_DATA_ROOT.exists():
        TEST_DATA_ROOT.mkdir(parents=True, exist_ok=True)

    # 2. Initialize the singleton (Simulating App Startup)
    # This sets the global root path for the AssetStorage to use
    workspace_path_manager.initialize(str(TEST_DATA_ROOT))


def test_scanpy_clustering_persistence(setup_workspace):
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

    # The expected relative path key required by your business logic
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

    print("✅ TEST PASSED: Data successfully clustered, saved, and reloaded.")

