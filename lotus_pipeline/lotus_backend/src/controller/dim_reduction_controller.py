from fastapi import APIRouter, HTTPException, status, Depends
from typing import Annotated

from dto.request.run_neighbors_request import RunNeighborsRequest
from dto.request.run_pca_request import RunPCARequest
from dto.response.neighbor_result_dto import NeighborsResultDTO
from service.dim_reduction_service import DimReductionService

from dto.response.pca_result_dto import PCAResultDTO
from loguru import logger

router = APIRouter(prefix="/api/v1/dim-reduction", tags=["DimReduction"])
# injected instance
ServiceDep = Annotated[DimReductionService, Depends()]


@router.post(
    "/pca",
    response_model=PCAResultDTO,
    status_code=status.HTTP_201_CREATED,
    summary="Run PCA Analysis",
    description="Performs Principal Component Analysis on highly variable genes and returns Elbow Plot data."
)
async def run_pca_analysis(
        request: RunPCARequest,
        service: ServiceDep
):
    """
    Executes PCA, extracts variance ratios for the Elbow Plot,
    and saves the resulting snapshot.
    """
    try:
        logger.info(f"Starting PCA for project: {request.project_id}, dataset: {request.dataset_id}")
        result = service.run_pca(request)
        return result

    except ValueError as e:
        # Handles missing snapshots or invalid dataset IDs
        logger.error(f"PCA validation error: {str(e)}")
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))
    except Exception as e:
        # Handles unexpected calculation or I/O errors
        logger.exception("Unexpected error during PCA execution")
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=str(e))


@router.post(
    "/neighbors",
    response_model=NeighborsResultDTO,
    status_code=status.HTTP_201_CREATED,
    summary="Build Neighborhood Graph",
    description="Computes cell-cell neighbor relations based on a specified number of PCs."
)
async def build_neighborhood_graph(
        request: RunNeighborsRequest,
        service: ServiceDep
):
    """
    Constructs the k-nearest neighbor graph required for UMAP and Clustering.
    User provides 'n_pcs' based on the Elbow Plot results.
    """
    try:
        logger.info(f"Building neighbors graph for project: {request.project_id} using {request.n_pcs} PCs")
        result = service.build_neighborhood_graph(request)
        return result

    except ValueError as e:
        # Handles cases where PCA snapshot is missing
        logger.error(f"Neighbors graph validation error: {str(e)}")
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))
    except Exception as e:
        logger.exception("Unexpected error during neighbors graph construction")
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=str(e))
