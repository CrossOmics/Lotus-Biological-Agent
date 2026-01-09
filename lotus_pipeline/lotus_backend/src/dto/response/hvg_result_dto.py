from pydantic import BaseModel


class HVGResultDTO(BaseModel):
    snapshot_id: str
    snapshot_path: str
    hvg_plot_path: str  # Path to 'hvg_dispersion.pdf' or '.png'
    n_genes_found: int
    msg: str = "HVG selection and scaling complete."
