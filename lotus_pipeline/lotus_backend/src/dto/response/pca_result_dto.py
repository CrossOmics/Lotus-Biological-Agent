from typing import List

from pydantic import BaseModel


class PCAResultDTO(BaseModel):
    variance_plot_data: List[float]
    pca_scatter_path: str
    snapshot_path: str
    snapshot_id: str
