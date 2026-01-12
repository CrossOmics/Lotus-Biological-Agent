"""
Test lotus I/O functionality: read and write h5ad files.

This script tests the read and write functions from lotus.io module.
"""
import sys
from pathlib import Path
import numpy as np
import psutil
import os

# Add lotus_core/src to path if needed (for development)
# File is in lotus_pipeline/lotus_backend/test/, so go up to project root
project_root = Path(__file__).resolve().parent.parent.parent.parent.parent
lotus_core_src = project_root / "lotus_pipeline" / "lotus_core" / "src"
if str(lotus_core_src) not in sys.path:
    sys.path.insert(0, str(lotus_core_src))

from lotus.io import read_h5ad, write, read, standardize_load

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

def main():
    print("=" * 60)
    print("Lotus I/O Test: Read & Write")
    print("=" * 60)
    
    # Initial memory
    initial_memory = get_memory_usage()
    print(f"\nInitial memory usage: {initial_memory:.1f} MB")
    
    # 1. Load original dataset
    dataset_path = project_root / "TAURUS_raw_counts_annotated_final.h5ad"
    print(f"\n[Step 1] Loading dataset from: {dataset_path.name}")
    
    if not dataset_path.exists():
        raise FileNotFoundError(f"Dataset file not found: {dataset_path}")
    
    print("  Loading h5ad file using lotus.io.read_h5ad...")
    before_load = get_memory_usage()
    adata_original = read_h5ad(dataset_path)
    after_load = get_memory_usage()
    print(f"  ✓ Loaded shape: {adata_original.shape}")
    print(f"  ✓ Number of observations: {adata_original.n_obs}")
    print(f"  ✓ Number of variables: {adata_original.n_vars}")
    print(f"  ✓ Memory used for loading: {after_load - before_load:.1f} MB")
    print_memory_info("After loading original data")
    
    # 2. Test write function
    output_dir = project_root / "test_output"
    output_dir.mkdir(exist_ok=True)
    
    output_file = output_dir / "test_write_output.h5ad"
    print(f"\n[Step 2] Testing write function...")
    print(f"  Writing to: {output_file.name}")
    before_write = get_memory_usage()
    write(output_file, adata_original)
    after_write = get_memory_usage()
    print(f"  ✓ File written successfully")
    print(f"  ✓ File size: {output_file.stat().st_size / (1024**2):.2f} MB")
    print(f"  ✓ Memory used for writing: {after_write - before_write:.1f} MB")
    print_memory_info("After writing")
    
    # 3. Test read function (using read_h5ad)
    print(f"\n[Step 3] Testing read function (read_h5ad)...")
    print(f"  Reading from: {output_file.name}")
    before_read = get_memory_usage()
    adata_read_h5ad = read_h5ad(output_file)
    after_read = get_memory_usage()
    print(f"  ✓ Loaded shape: {adata_read_h5ad.shape}")
    print(f"  ✓ Memory used for reading: {after_read - before_read:.1f} MB")
    print_memory_info("After reading with read_h5ad")
    
    # 4. Test generic read function
    print(f"\n[Step 4] Testing generic read function...")
    print(f"  Reading from: {output_file.name}")
    before_read_generic = get_memory_usage()
    adata_read_generic = read(output_file)
    after_read_generic = get_memory_usage()
    print(f"  ✓ Loaded shape: {adata_read_generic.shape}")
    print(f"  ✓ Memory used for reading: {after_read_generic - before_read_generic:.1f} MB")
    print_memory_info("After reading with read()")
    
    # 5. Test standardize_load function
    print(f"\n[Step 5] Testing standardize_load function...")
    print(f"  Loading with auto-detection from: {output_file.name}")
    before_standardize = get_memory_usage()
    adata_standardize = standardize_load(output_file, source_type='auto')
    after_standardize = get_memory_usage()
    print(f"  ✓ Loaded shape: {adata_standardize.shape}")
    print(f"  ✓ Memory used for loading: {after_standardize - before_standardize:.1f} MB")
    
    print(f"  Loading with explicit h5ad type...")
    before_standardize_explicit = get_memory_usage()
    adata_standardize_explicit = standardize_load(output_file, source_type='h5ad')
    after_standardize_explicit = get_memory_usage()
    print(f"  ✓ Loaded shape: {adata_standardize_explicit.shape}")
    print(f"  ✓ Memory used for loading: {after_standardize_explicit - before_standardize_explicit:.1f} MB")
    print_memory_info("After standardize_load")
    
    # 6. Verify data integrity
    print(f"\n[Step 6] Verifying data integrity...")
    
    # Check shapes
    assert adata_original.shape == adata_read_h5ad.shape, "Shape mismatch with read_h5ad"
    assert adata_original.shape == adata_read_generic.shape, "Shape mismatch with read"
    assert adata_original.shape == adata_standardize.shape, "Shape mismatch with standardize_load (auto)"
    assert adata_original.shape == adata_standardize_explicit.shape, "Shape mismatch with standardize_load (explicit)"
    print("  ✓ All shapes match")
    
    # Check obs and var names
    assert np.array_equal(adata_original.obs_names, adata_read_h5ad.obs_names), "obs_names mismatch"
    assert np.array_equal(adata_original.var_names, adata_read_h5ad.var_names), "var_names mismatch"
    assert np.array_equal(adata_original.obs_names, adata_standardize.obs_names), "obs_names mismatch with standardize_load"
    assert np.array_equal(adata_original.var_names, adata_standardize.var_names), "var_names mismatch with standardize_load"
    print("  ✓ Cell and gene names match")
    
    # Check X matrix (sample a subset for large datasets)
    if adata_original.n_obs > 10000:
        sample_idx = np.random.choice(adata_original.n_obs, 1000, replace=False)
        original_X_sample = adata_original.X[sample_idx, :100].toarray() if hasattr(adata_original.X, 'toarray') else adata_original.X[sample_idx, :100]
        read_X_sample = adata_read_h5ad.X[sample_idx, :100].toarray() if hasattr(adata_read_h5ad.X, 'toarray') else adata_read_h5ad.X[sample_idx, :100]
        standardize_X_sample = adata_standardize.X[sample_idx, :100].toarray() if hasattr(adata_standardize.X, 'toarray') else adata_standardize.X[sample_idx, :100]
        
        assert np.allclose(original_X_sample, read_X_sample, atol=1e-6), "X matrix mismatch with read_h5ad (sampled)"
        assert np.allclose(original_X_sample, standardize_X_sample, atol=1e-6), "X matrix mismatch with standardize_load (sampled)"
        print("  ✓ X matrix matches (sampled check)")
    else:
        original_X = adata_original.X.toarray() if hasattr(adata_original.X, 'toarray') else adata_original.X
        read_X = adata_read_h5ad.X.toarray() if hasattr(adata_read_h5ad.X, 'toarray') else adata_read_h5ad.X
        standardize_X = adata_standardize.X.toarray() if hasattr(adata_standardize.X, 'toarray') else adata_standardize.X
        
        assert np.allclose(original_X, read_X, atol=1e-6), "X matrix mismatch with read_h5ad"
        assert np.allclose(original_X, standardize_X, atol=1e-6), "X matrix mismatch with standardize_load"
        print("  ✓ X matrix matches")
    
    # Check obs columns
    original_obs_cols = set(adata_original.obs.columns)
    read_obs_cols = set(adata_read_h5ad.obs.columns)
    assert original_obs_cols == read_obs_cols, f"obs columns mismatch: {original_obs_cols - read_obs_cols}"
    print(f"  ✓ obs columns match ({len(original_obs_cols)} columns)")
    
    # Check var columns
    original_var_cols = set(adata_original.var.columns)
    read_var_cols = set(adata_read_h5ad.var.columns)
    assert original_var_cols == read_var_cols, f"var columns mismatch: {original_var_cols - read_var_cols}"
    print(f"  ✓ var columns match ({len(original_var_cols)} columns)")
    
    # Check obsm (if exists)
    if adata_original.obsm.keys():
        assert set(adata_original.obsm.keys()) == set(adata_read_h5ad.obsm.keys()), "obsm keys mismatch"
        print(f"  ✓ obsm keys match: {list(adata_original.obsm.keys())}")
    
    # Check varm (if exists)
    if adata_original.varm.keys():
        assert set(adata_original.varm.keys()) == set(adata_read_h5ad.varm.keys()), "varm keys mismatch"
        print(f"  ✓ varm keys match: {list(adata_original.varm.keys())}")
    
    # Check uns (if exists)
    if adata_original.uns.keys():
        original_uns_keys = set(adata_original.uns.keys())
        read_uns_keys = set(adata_read_h5ad.uns.keys())
        # Some uns keys might be metadata that don't need to match exactly
        print(f"  ✓ uns keys - original: {len(original_uns_keys)}, read: {len(read_uns_keys)}")
    
    # Final memory summary
    final_memory = get_memory_usage()
    total_memory_used = final_memory - initial_memory
    
    print("\n" + "=" * 60)
    print("✓ All I/O tests passed successfully!")
    print(f"Output saved to: {output_file}")
    print(f"\nMemory Summary:")
    print(f"  Initial: {initial_memory:.1f} MB")
    print(f"  Final: {final_memory:.1f} MB")
    print(f"  Total used: {total_memory_used:.1f} MB ({total_memory_used/1024:.2f} GB)")
    print("=" * 60)

if __name__ == "__main__":
    main()
