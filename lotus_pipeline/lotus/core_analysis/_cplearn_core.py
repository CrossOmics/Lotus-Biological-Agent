from typing import Literal, Any, Mapping
import numpy as np
import pandas as pd
from anndata import AnnData
from cplearn.corespect import CorespectModel
from cplearn.corespect.config import CoreSpectConfig


def corespect(
        adata: AnnData,
        *,
        use_rep: str = "X_pca",
        key_added: str = "corespect",
        # --- CoreSpect Config Parameters ---
        q: int = 20,
        r: int = 10,
        core_frac: float = 0.2,
        densify: bool = False,
        granularity: float = 0.5,
        resolution: float = 0.5,
        # --- Run Parameters ---
        fine_grained: bool = True,
        propagate: bool = True,
        copy: bool = False,
) -> AnnData | None:
    """
    Cluster cells using the CoreSpect algorithm.

    This algorithm identifies stable 'core' cells and propagates labels to the rest.

    Parameters:
        adata: The annotated data matrix.
        use_rep: The representation to use (e.g., 'X_pca', 'X_umap').
        key_added: Key under which to add the cluster labels.
        q: Determines neighborhood size for the underlying q-NN graph.
        r: Neighborhood radius parameter for ascending random walk with FlowRank.
        core_frac: Fraction of points in the top-layer (core size).
        densify: Whether to densify different parts of the data to reduce fragmentation.
        granularity: Higher granularity finds more local cores but implies weaker clusters.
        resolution: Resolution for clustering with Leiden inside CoreSpect.
        fine_grained: Whether to run the fine-grained core refinement step.
        propagate: Whether to propagate labels from cores to all cells.
        copy: Whether to return a copy of the AnnData object or modify it in-place.

    Returns:
        If copy=True: Returns a new AnnData object with the results.
        If copy=False: Returns None and updates the input `adata` in-place.

        Updates:
            - adata.obs[key_added]: The cluster labels.
            - adata.obs[key_added + '_is_core']: Boolean indicating if a cell belongs to the stable core.
            - adata.uns[key_added]: Dictionary containing parameters and layer info.
    """

    # 0. Handle In-Place vs Copy Mode
    # If copy=True: Create a deep copy of the input AnnData object.
    # This ensures the original object remains unchanged. We will perform
    # the analysis on this new instance and return it at the end.
    if copy:
        adata = adata.copy()

    # 1. Prepare Data
    if use_rep not in adata.obsm:
        raise KeyError(f"Could not find representation '{use_rep}' in adata.obsm")

    X = adata.obsm[use_rep]

    # 2. Configure Model
    # Wrap parameters into the config object as required by CoreSpect
    cfg = CoreSpectConfig(
        q=q,
        r=r,
        core_frac=core_frac,
        densify=densify,
        granularity=granularity,
        resolution=resolution
    ).configure()

    # 3. Initialize and Run Model
    # unpack() unpacks the config into separate sub-configs (flowrank_cfg, etc.)
    model = CorespectModel(X, **cfg.unpack())

    # Run the core logic to perform clustering and core identification
    model.run(fine_grained=fine_grained, propagate=propagate)

    # 4. Write Results back to AnnData
    # Note: 'adata' here refers to either the deep copy (if copy=True)
    # or the original reference (if copy=False).

    # A. Write Cluster Labels
    # Convert to categorical string for consistent plotting compatibility
    labels = model.labels_.astype(str)
    # Handle potentially unassigned points (-1) if propagation was False
    labels[labels == "-1"] = "Unassigned"
    adata.obs[key_added] = pd.Categorical(labels)

    # B. Write Core Identity (Boolean mask)
    # model.layers_[0] contains indices of the core cells
    is_core = np.zeros(adata.n_obs, dtype=bool)
    if model.layers_ and len(model.layers_) > 0:
        core_indices = model.layers_[0]
        is_core[core_indices] = True

    core_key = f"{key_added}_is_core"
    adata.obs[core_key] = is_core

    # C. Write Parameters and Layers to uns (Metadata)
    adata.uns[key_added] = {
        'params': {
            'q': q, 'r': r, 'core_frac': core_frac,
            'densify': densify, 'resolution': resolution,
            'fine_grained': fine_grained, 'propagate': propagate
        },
        # Storing layers can be useful for advanced visualization (Coremap)
        # We store them as a list of arrays
        'layers_indices': model.layers_
    }

    # 5. Return Logic
    # If copy=True: Return the modified new instance.
    # If copy=False: Return None (standard Scanpy in-place behavior).
    return adata if copy else None
