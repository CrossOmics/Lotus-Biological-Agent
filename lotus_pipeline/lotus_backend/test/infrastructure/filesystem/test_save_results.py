import pytest
import scanpy as sc
from infrastructure.filesystem.storage import AssetStorage
from infrastructure.workspace_context import workspace_path_manager
from lotus.preprocessing import neighbors
from lotus.clustering import leiden
from lotus.io import read_h5ad


@pytest.fixture(scope="function")
def temp_workspace(tmp_path):
    """
    Fixture: Provide isolated temporary workspace for each test.
    Automatically cleans up after test completion.
    """
    # Initialize workspace with temp directory
    workspace_path_manager.reset()
    workspace_path_manager.initialize(str(tmp_path))

    yield tmp_path

    # Cleanup: Reset workspace context
    workspace_path_manager.reset()


@pytest.fixture
def mock_adata():
    """
    Fixture: Generate mock single-cell data for testing.
    Returns: AnnData with 100 cells x 50 genes.
    """
    return sc.datasets.blobs(n_variables=50, n_observations=100)


def test_scanpy_clustering_persistence(temp_workspace, mock_adata):
    """
    CI/CD Test: Verify Scanpy clustering results persist correctly via AssetStorage.

    Pipeline:
        1. Run Leiden clustering (in-memory)
        2. Save to storage with relative key
        3. Verify file creation and data integrity
        4. Cleanup automatically via fixture
    """

    # Step 1: Execute clustering pipeline
    neighbors(mock_adata, n_neighbors=10, n_pcs=10)
    leiden(mock_adata, resolution=0.5)

    assert 'leiden' in mock_adata.obs, "Clustering failed: 'leiden' key missing"

    # Step 2: Persist to storage
    storage = AssetStorage()
    relative_key = "cluster/test/test01.h5ad"
    saved_key = storage.save_anndata(mock_adata, relative_key)

    # Step 3: Assertions
    assert saved_key == relative_key, "Returned key mismatch"

    # Verify physical file exists
    expected_path = temp_workspace / relative_key
    assert expected_path.exists(), f"File not created at: {expected_path}"
    assert expected_path.stat().st_size > 0, "File is empty"

    # Step 4: Verify data integrity via reload
    reloaded_adata = read_h5ad(expected_path)
    assert 'leiden' in reloaded_adata.obs, "Clustering info lost after persistence"
    assert reloaded_adata.n_obs == 100, "Cell count mismatch"
    assert reloaded_adata.n_vars == 50, "Gene count mismatch"

    # Cleanup handled by temp_workspace fixture


def test_storage_load_capability(temp_workspace):
    """
    CI/CD Test: Verify AssetStorage load functionality.

    Pipeline:
        1. Create and save test data
        2. Load via relative key
        3. Verify data integrity
        4. Cleanup automatically
    """

    # Step 1: Setup - Create test data
    dummy_adata = sc.datasets.blobs(n_variables=20, n_observations=50)

    storage = AssetStorage()
    target_key = "cluster/test/test_load.h5ad"
    storage.save_anndata(dummy_adata, target_key)

    # Step 2: Execute load operation
    loaded_adata = storage.load_anndata(target_key)

    # Step 3: Assertions
    assert loaded_adata is not None, "Load returned None"
    assert loaded_adata.n_obs == 50, "Cell count mismatch"
    assert loaded_adata.n_vars == 20, "Gene count mismatch"
    assert loaded_adata.obs_names[0] == dummy_adata.obs_names[0], "Index mismatch"

    # Cleanup handled by temp_workspace fixture


def test_storage_overwrites_existing_file(temp_workspace, mock_adata):
    """
    CI/CD Test: Verify storage overwrites existing files correctly.

    Pipeline:
        1. Save initial data
        2. Modify data and save to same key
        3. Verify file was overwritten (not duplicated)
        4. Cleanup automatically
    """

    storage = AssetStorage()
    key = "cluster/test/overwrite.h5ad"

    # First save
    storage.save_anndata(mock_adata, key)
    first_size = (temp_workspace / key).stat().st_size

    # Modify data
    mock_adata.obs['new_column'] = 'test_value'

    # Second save (overwrite)
    storage.save_anndata(mock_adata, key)
    second_size = (temp_workspace / key).stat().st_size

    # Verify overwrite occurred (size should differ)
    assert second_size != first_size, "File was not overwritten"

    # Verify new data exists
    reloaded = storage.load_anndata(key)
    assert 'new_column' in reloaded.obs, "Overwritten file missing new data"


def test_storage_handles_nested_paths(temp_workspace, mock_adata):
    """
    CI/CD Test: Verify storage creates nested directories automatically.

    Pipeline:
        1. Save to deeply nested path
        2. Verify all intermediate directories created
        3. Verify file is loadable
        4. Cleanup automatically
    """

    storage = AssetStorage()
    nested_key = "level1/level2/level3/deep_file.h5ad"

    storage.save_anndata(mock_adata, nested_key)

    # Verify nested directories exist
    expected_path = temp_workspace / nested_key
    assert expected_path.exists(), "Nested file not created"
    assert expected_path.parent.exists(), "Parent directories not created"

    # Verify loadable
    loaded = storage.load_anndata(nested_key)
    assert loaded.n_obs == mock_adata.n_obs


def test_storage_load_nonexistent_file(temp_workspace):
    """
    CI/CD Test: Verify storage raises error for missing files.

    Expected: FileNotFoundError when loading non-existent key.
    """

    storage = AssetStorage()

    with pytest.raises(FileNotFoundError):
        storage.load_anndata("nonexistent/path.h5ad")


def test_storage_save_invalid_key(temp_workspace, mock_adata):
    """
    CI/CD Test: Verify storage rejects invalid relative keys.

    Expected: ValueError for absolute paths or path traversal attempts.
    """

    storage = AssetStorage()

    # Test absolute path rejection
    with pytest.raises(ValueError, match="Security violation"):
        storage.save_anndata(mock_adata, "/absolute/path.h5ad")

    # Test path traversal rejection
    with pytest.raises(ValueError, match="Security violation"):
        storage.save_anndata(mock_adata, "../escape/path.h5ad")


@pytest.mark.parametrize("n_obs,n_vars", [
    (10, 5),
    (100, 50),
    (1000, 100),
])
def test_storage_handles_various_sizes(temp_workspace, n_obs, n_vars):
    """
    CI/CD Test: Verify storage handles various data sizes.

    Parameterized test for different AnnData dimensions.
    """

    adata = sc.datasets.blobs(n_variables=n_vars, n_observations=n_obs)
    storage = AssetStorage()
    key = f"size_test/{n_obs}x{n_vars}.h5ad"

    storage.save_anndata(adata, key)
    loaded = storage.load_anndata(key)

    assert loaded.n_obs == n_obs
    assert loaded.n_vars == n_vars


# Integration Test Suite
def test_full_pipeline_integration(temp_workspace):
    """
    CI/CD Integration Test: Complete analysis pipeline.

    Simulates real workflow:
        1. Generate data
        2. Preprocess (neighbors)
        3. Cluster (Leiden)
        4. Save results
        5. Load and verify
        6. Auto-cleanup
    """

    # Generate data
    adata = sc.datasets.blobs(n_variables=50, n_observations=100)

    # Preprocessing
    neighbors(adata, n_neighbors=10, n_pcs=10)

    # Clustering
    leiden(adata, resolution=0.5)

    # Persist
    storage = AssetStorage()
    key = "integration/pipeline_result.h5ad"
    storage.save_anndata(adata, key)

    # Reload and verify complete pipeline
    result = storage.load_anndata(key)

    assert 'leiden' in result.obs, "Pipeline result incomplete"
    assert 'neighbors' in result.uns, "Neighbor graph lost"
    assert result.n_obs == 100
    assert result.n_vars == 50

    # Success: All operations completed and cleaned up automatically
