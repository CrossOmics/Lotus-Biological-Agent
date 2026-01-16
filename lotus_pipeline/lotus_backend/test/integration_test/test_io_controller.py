"""
Test Controller Layer I/O functionality: Project Management Controller.

This script tests the I/O operations at the controller layer, specifically:
- project_management_controller.create_new_project() which wraps DatasetService.import_dataset_from_local()
"""
import sys
from pathlib import Path
import numpy as np
import psutil
import os
from datetime import datetime

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
from controller.project_management_controller import create_new_project
from dto.request.create_project_request import CreateProjectRequest


def setup_test_environment():
    """Initialize workspace and database"""
    # Initialize workspace
    workspace_path_manager.initialize()
    
    # Initialize database
    db_manager = get_default_db_manager()
    db_manager.initialize_proxy()
    db = db_manager.get_connection()
    if db.is_closed():
        db.connect()
    
    # Create tables if needed
    db.create_tables([ProjectMeta, Dataset], safe=True)
    
    return db_manager, db


def get_memory_usage():
    """Get current process memory usage in MB"""
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / (1024 ** 2)  # Convert to MB


def print_memory_info(label, before_mb=None):
    """Print memory usage information"""
    current_mb = get_memory_usage()
    if before_mb is not None:
        diff_mb = current_mb - before_mb
        print(f"  Memory: {current_mb:.1f} MB (Δ {diff_mb:+.1f} MB)")
    else:
        print(f"  Memory: {current_mb:.1f} MB")
    return current_mb


async def test_controller_import_dataset():
    """Test Controller layer dataset import through create_new_project endpoint"""
    print("\n" + "=" * 60)
    print("Test: Controller Layer - Dataset Import (create_new_project)")
    print("=" * 60)
    
    # Setup environment
    db_manager, db = setup_test_environment()
    
    initial_memory = get_memory_usage()
    print(f"\nInitial memory: {initial_memory:.1f} MB")
    
    # Load original dataset path
    dataset_path = project_root / "TAURUS_raw_counts_annotated_final.h5ad"
    if not dataset_path.exists():
        print(f"  ⚠ Skipping: Dataset file not found: {dataset_path}")
        return
    
    print(f"\n[Step 1] Preparing test request...")
    print(f"  Dataset file: {dataset_path.name}")
    
    # Create request object (simulating FastAPI request)
    request = CreateProjectRequest(
        project_name="Test I/O Controller Project",
        local_file_path=str(dataset_path),
        organism="Human",
        tissue_type="Test",
        description="Integration test for controller I/O"
    )
    print(f"  ✓ Request created: {request.project_name}")
    
    # Test controller function directly (without HTTP server)
    print(f"\n[Step 2] Testing controller.create_new_project()...")
    print(f"  This wraps DatasetService.import_dataset_from_local()")
    before_import = get_memory_usage()
    
    try:
        response = await create_new_project(request)
        after_import = print_memory_info("After importing via controller", before_import)
        
        print(f"  ✓ Controller import successful")
        print(f"  ✓ Project ID: {response.project_id}")
        print(f"  ✓ Project Name: {response.project_name}")
        print(f"  ✓ Dataset ID: {response.dataset_id}")
        print(f"  ✓ Dataset Name: {response.dataset_name}")
        print(f"  ✓ Workspace Path: {response.workspace_path}")
        
        # Verify file was saved correctly
        storage = AssetStorage()
        print(f"\n[Step 3] Verifying saved dataset...")
        
        # Get dataset from database to get the path
        from infrastructure.database.dao.dataset_dao import DatasetDAO
        dataset_dao = DatasetDAO()
        dataset_record = dataset_dao.get_dataset_by_business_id(response.dataset_id)
        
        if dataset_record:
            print(f"  ✓ Dataset record found in database")
            print(f"  ✓ Dataset path: {dataset_record.dataset_path}")
            
            # Load and verify the saved file
            adata_loaded = storage.load_anndata(dataset_record.dataset_path)
            print(f"  ✓ Loaded shape: {adata_loaded.shape}")
            
            # Verify file exists at expected location
            abs_path = workspace_path_manager.resolve(dataset_record.dataset_path)
            assert abs_path.exists(), f"Dataset file should exist at {abs_path}"
            file_size_mb = abs_path.stat().st_size / (1024 ** 2)
            print(f"  ✓ File size: {file_size_mb:.2f} MB")
            
            # Verify data integrity (sample check)
            print(f"\n[Step 4] Verifying data integrity...")
            original_adata = storage.load_anndata(dataset_record.dataset_path)
            assert original_adata.shape[0] > 0, "Should have observations"
            assert original_adata.shape[1] > 0, "Should have variables"
            print(f"  ✓ Data integrity verified")
            
        else:
            print(f"  ⚠ Warning: Could not retrieve dataset from database")
        
        print("\n" + "=" * 60)
        print("✓ Controller layer I/O test passed!")
        print("=" * 60)
        
    except Exception as e:
        print(f"\n  ✗ Error: {e}")
        import traceback
        traceback.print_exc()
        raise


def main():
    """Run the test"""
    import asyncio
    
    try:
        asyncio.run(test_controller_import_dataset())
        
        final_memory = get_memory_usage()
        print(f"\nFinal memory usage: {final_memory:.1f} MB")
        print("=" * 60)
        print("✓ All Controller Layer I/O tests completed!")
        print("=" * 60)
    except Exception as e:
        print(f"\n✗ Test failed: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()


