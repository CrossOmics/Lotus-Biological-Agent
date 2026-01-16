# Lotus Biological Agent

A biological analysis engine for single-cell RNA sequencing data.

## Integration Tests

Integration tests validate the Controller layer functionality, ensuring that lotus preprocessing functions produce identical results to direct scanpy calls.

### Test Files

1. **`test_io_controller.py`** - Tests I/O operations at the Controller layer
2. **`test_pp_controller.py`** - Tests Preprocessing operations at the Controller layer (can run individually)

### Test Coverage

#### Controller Layer Tests

| Controller | Function | Test File | Status |
|-----------|---------|-----------|--------|
| Project Management | `create_new_project()` | `test_io_controller.py` | ✅ Passed |
| Preprocessing | `perform_qc_calculation()` | `test_pp_controller.py` | ✅ Passed |
| Preprocessing | `apply_qc_filter()` | `test_pp_controller.py` | ✅ Passed |
| Preprocessing | `apply_hvg()` | `test_pp_controller.py` | ✅ Passed |

**Coverage: 4/6 controller functions (67%)**

### Running Tests

```bash
# Test individual preprocessing functions
python test_pp_controller.py qc_calc      # Test QC calculation only
python test_pp_controller.py qc_filter   # Test QC filter only
python test_pp_controller.py hvg          # Test HVG only
python test_pp_controller.py all          # Test all (default)

# Test I/O controller
python test_io_controller.py
```

### Test Results

**Dataset:** TAURUS dataset (987,743 cells × 33,075 genes)

#### Test 1: QC Calculation (`perform_qc_calculation`)
- ✅ Function runs successfully
- ✅ Generates QC metrics JSON and visualization PDFs
- ✅ **Results match scanpy:** All QC metrics (n_genes_by_counts, total_counts, pct_counts_mt, pct_counts_ribo) are identical (max diff: 0.00e+00)

#### Test 2: QC Filter (`apply_qc_filter`)
- ✅ Function runs successfully
- ✅ Filters cells: 987,743 → 657,111 cells
- ✅ Filters genes: 33,075 → 31,113 genes
- ✅ **Results match scanpy:** Filtered shapes are identical (657,111 × 31,113)

#### Test 3: HVG Selection (`apply_hvg`)
- ✅ Function runs successfully
- ✅ Selects 2,000 highly variable genes
- ✅ Generates HVG dispersion plot
- ✅ **Results match scanpy:** HVG count matches (2,000 genes)

#### I/O Controller Test (`create_new_project`)
- ✅ Dataset import successful
- ✅ File storage and database persistence verified
- ✅ Data integrity maintained

### Performance Metrics

**Memory Usage:**
- Initial: ~370 MB
- QC Calculation: +130-770 MB (peak ~1.1 GB)
- QC Filter: ~460 MB (after filtering, memory reduced)
- HVG Selection: +1,300 MB (peak ~1.8 GB, includes scaling which densifies sparse matrix)
- **Memory efficiency:** Sparse matrix format saves ~99% memory compared to dense matrix

**File Sizes:**
- Original dataset: 12 GB (uncompressed)
- Saved files: 3.6 GB (gzip compressed)
- QC filtered snapshot: 2.27 GB
- **Compression ratio:** ~70% size reduction

### Test Strategy

Each test verifies:
1. **Functionality:** Can the function run successfully?
2. **Correctness:** Are results similar to direct scanpy calls?

Tests use the same dataset and compare lotus results with scanpy results to ensure the abstraction layer doesn't introduce bugs.
