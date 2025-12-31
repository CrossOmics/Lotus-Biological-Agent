from pathlib import Path
import numpy as np
from scipy import sparse
import anndata as ad
import pandas as pd
import pytest

from service import dataset_service
from infrastructure.workspace_context import workspace_path_manager


def generate_mock_file(data_dir: Path, file_type: str) -> Path:
    """
    Create small biological-like dataset in various formats.

    :param data_dir: Path to store test files
    :param file_type: 'csv', 'tsv', 'mtx', 'h5ad', 'xlsx'
    :return: Path to the generated file
    """
    data_dir.mkdir(exist_ok=True, parents=True)

    genes = ["GeneA", "GeneB", "GeneC"]
    cells = ["Cell1", "Cell2", "Cell3"]

    # Dense matrix for CSV/TSV/Excel
    matrix = np.array([[1, 0, 2], [0, 3, 0], [1, 1, 1]])
    df = pd.DataFrame(matrix, index=genes, columns=cells)
    df.index.name = "gene"

    file_path = data_dir / f"mock_data.{file_type}"

    if file_type == "csv":
        df.to_csv(file_path)
    elif file_type == "tsv":
        df.to_csv(file_path, sep="\t")
    elif file_type == "xlsx":
        file_path = data_dir / "mock_data.xlsx"
        df = df.astype(float)  # ensure numeric
        df.to_excel(file_path, index=True)
    elif file_type == "mtx":
        # For MTX, use scipy.sparse
        sparse_matrix = sparse.csr_matrix(matrix)
        from scipy.io import mmwrite
        mmwrite(file_path, sparse_matrix)
    elif file_type == "h5ad":
        adata = ad.AnnData(X=matrix, obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=genes))
        adata.write(file_path)
    else:
        raise ValueError(f"Unsupported file type: {file_type}")

    return file_path


path = Path(__file__).resolve().parent.parent
TEST_DATA_DIR = path / "test_data"


@pytest.mark.parametrize("file_type", ["csv", "tsv", "mtx", "h5ad", "xlsx"])
def test_import_various_formats(file_type):
    file_path = generate_mock_file(TEST_DATA_DIR, file_type)

    service = dataset_service.DatasetService()
    meta = service.import_dataset_from_local(str(file_path), project_id="test_id_1", dataset_name=f"Test_{file_type}")

    print(f"[Test {file_type}] Metadata:", meta)

    # Basic assertions
    assert meta["n_obs"] == 3
    assert meta["n_vars"] == 3
    saved_file = workspace_path_manager.root / meta["relative_path"]
    assert saved_file.exists()
    assert saved_file.stat().st_size > 0
