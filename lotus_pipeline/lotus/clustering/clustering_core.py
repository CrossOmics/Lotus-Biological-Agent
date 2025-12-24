from typing import Literal, Sequence, Mapping, Any, Tuple
import numpy as np
import scanpy as sc
from anndata import AnnData
from scipy.sparse import csr_matrix, csc_matrix


def leiden(
        adata: AnnData,
        resolution: float = 1,
        *,
        restrict_to: Tuple[str, Sequence[str]] | None = None,
        random_state: int | np.random.RandomState | None = 0,
        key_added: str = 'leiden',
        adjacency: csr_matrix | csc_matrix | None = None,
        directed: bool | None = None,
        use_weights: bool = True,
        n_iterations: int = -1,
        partition_type: Any | None = None,
        neighbors_key: str | None = None,
        obsp: str | None = None,
        copy: bool = False,
        flavor: Literal['leidenalg', 'igraph'] = 'leidenalg',
        **clustering_args,
) -> AnnData | None:
    """
    Cluster cells into subgroups using the Leiden algorithm.

    This requires having run `neighbors()` or `bbknn()` first.

    Parameters:
        adata: The annotated data matrix.
        resolution: A parameter value controlling the coarseness of the clustering.
                    Higher values lead to more clusters.
        restrict_to: Restrict the clustering to the categories within the key for sample annotation.
        random_state: Change the initialization of the optimization.
        key_added: adata.obs key under which to add the cluster labels.
        adjacency: Sparse adjacency matrix of the graph, defaults to neighbors connectivities.
        directed: Whether to treat the graph as directed or undirected.
        use_weights: If True, edge weights from the graph are used in the computation.
        n_iterations: How many iterations of the Leiden clustering algorithm to perform.
        partition_type: Type of partition to use.
        neighbors_key: Use neighbors connectivities as adjacency.
        obsp: Use .obsp[obsp] as adjacency.
        copy: Whether to copy adata or modify it inplace.
        flavor: Which package’s implementation to use ('leidenalg', 'igraph').
        **clustering_args: Any further arguments to pass to find_partition().

    Returns:
        Returns None if copy=False, else returns an AnnData object.
        Sets `adata.obs[key_added]` (cluster labels).
    """
    return sc.tl.leiden(
        adata,
        resolution=resolution,
        restrict_to=restrict_to,
        random_state=random_state,
        key_added=key_added,
        adjacency=adjacency,
        directed=directed,
        use_weights=use_weights,
        n_iterations=n_iterations,
        partition_type=partition_type,
        neighbors_key=neighbors_key,
        obsp=obsp,
        copy=copy,
        flavor=flavor,
        **clustering_args
    )


def louvain(
        adata: AnnData,
        resolution: float | None = None,
        *,
        random_state: int | np.random.RandomState | None = 0,
        restrict_to: Tuple[str, Sequence[str]] | None = None,
        key_added: str = 'louvain',
        adjacency: csr_matrix | csc_matrix | None = None,
        flavor: Literal['vtraag', 'igraph', 'rapids'] = 'vtraag',
        directed: bool = True,
        use_weights: bool = False,
        partition_type: Any | None = None,
        partition_kwargs: Mapping[str, Any] = {},
        neighbors_key: str | None = None,
        obsp: str | None = None,
        copy: bool = False,
) -> AnnData | None:
    """
    Cluster cells into subgroups using the Louvain algorithm.

    This requires having run `neighbors()` or `bbknn()` first.

    Parameters:
        adata: The annotated data matrix.
        resolution: Resolution parameter (higher means more clusters). Defaults to 1.0.
        random_state: Change the initialization of the optimization.
        restrict_to: Restrict the clustering to the categories within the key for sample annotation.
        key_added: Key under which to add the cluster labels.
        adjacency: Sparse adjacency matrix of the graph.
        flavor: Choose between packages for computing the clustering ('vtraag', 'igraph', 'rapids').
        directed: Interpret the adjacency matrix as directed graph?
        use_weights: Use weights from knn graph.
        partition_type: Type of partition to use.
        partition_kwargs: Key word arguments to pass to partitioning.
        neighbors_key: Use neighbors connectivities as adjacency.
        obsp: Use .obsp[obsp] as adjacency.
        copy: Copy adata or modify it inplace.

    Returns:
        Returns None if copy=False, else returns an AnnData object.
        Sets `adata.obs[key_added]` (cluster labels).
    """
    return sc.tl.louvain(
        adata,
        resolution=resolution,
        random_state=random_state,
        restrict_to=restrict_to,
        key_added=key_added,
        adjacency=adjacency,
        flavor=flavor,
        directed=directed,
        use_weights=use_weights,
        partition_type=partition_type,
        partition_kwargs=partition_kwargs,
        neighbors_key=neighbors_key,
        obsp=obsp,
        copy=copy
    )

def corespect_clustering():

