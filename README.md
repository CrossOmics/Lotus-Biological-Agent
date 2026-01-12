# Lotus Biological Agent

A biological analysis engine for single-cell RNA sequencing data.

## Integration Tests

### I/O Functionality Test

The integration test (`lotus_pipeline/lotus_backend/test/integration_test/test_ibd_io.py`) validates the core I/O functionality of the lotus pipeline.

#### Tested Functions

The test covers the following lotus I/O APIs:

1. **`read_h5ad()`** - Read h5ad format files
2. **`write()`** - Write AnnData objects to h5ad files
3. **`read()`** - Generic read function that auto-detects file format
4. **`standardize_load()`** - Standardized loading with format auto-detection and data validation

#### Test Results

**Dataset:** TAURUS dataset (987,743 cells × 33,075 genes)

**Test Coverage:**
- ✅ Data loading and saving
- ✅ Data integrity verification (shapes, names, X matrix, obs/var columns)
- ✅ Memory usage monitoring

**Performance Metrics:**
- **Memory Usage:**
  - Initial: ~200 MB
  - Loading dataset: +1,152 MB (total ~1.35 GB)
  - Final memory: ~1.17 GB
  - **Memory efficiency:** Loading 987K cells × 33K genes uses only ~1.15 GB (sparse matrix format saves ~99% memory compared to dense matrix)

- **File Size:**
  - Original file: 12 GB (uncompressed)
  - Saved file: 3.6 GB (gzip compressed)
  - **Compression ratio:** ~70% size reduction

**Key Findings:**
- All I/O functions work correctly
- Data integrity is preserved after save/load operations
- Sparse matrix format provides significant memory savings for single-cell data
- Gzip compression effectively reduces file size for sparse datasets

#### Running the Test

```bash
cd lotus_pipeline/lotus_backend/test/integration_test
python test_ibd_io.py
```

**Requirements:**
- Dataset: `TAURUS_raw_counts_annotated_final.h5ad` in project root
- Environment: `lotus_agent` conda environment with lotus package installed
