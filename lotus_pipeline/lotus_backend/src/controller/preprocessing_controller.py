from fastapi import APIRouter, HTTPException, status, Depends
from typing import Annotated

from dto.request.calculate_qc_request import CalculateQCRequest
from dto.request.filter_qc_request import FilterQCRequest
from dto.response.filter_qc_response import FilterQCResponse
from dto.response.qc_result_dto import QCResultDTO
from service.preprocessing_service import PreprocessingService

router = APIRouter(prefix="/api/v1/preprocessing", tags=["Preprocessing"])
# injected instance
ServiceDep = Annotated[PreprocessingService, Depends()]


@router.post("/qc/calculate", response_model=QCResultDTO, status_code=status.HTTP_200_OK)
async def perform_qc_calculation(
        request: CalculateQCRequest,
        service: ServiceDep
):
    """
    Calculate QC Metrics.
    """
    try:
        result = service.qc_calculation(
            project_id=request.project_id,
            dataset_id=request.dataset_id,
            organism=request.organism,
            custom_prefix=request.custom_prefixes
        )
        return result

    except ValueError as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=str(e))


@router.post("/qc/filter", response_model=FilterQCResponse,status_code=status.HTTP_201_CREATED)
async def apply_qc_filter(request: FilterQCRequest,service: ServiceDep):
    """
    Apply QC Filtering (Step 2).

    This endpoint:
    1. Accepts filtering thresholds (min_genes, mt_pct, etc.).
    2. Physically filters the AnnData object.
    3. Creates a new Snapshot (H5AD file + DB Record).
    """
    try:
        # Delegate business logic to Service
        result = service.apply_filter(request)
        return result

    except ValueError as e:
        # Case: Dataset ID provided in request does not exist in DB
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=str(e)
        )

    except FileNotFoundError as e:
        # Case: DB record exists, but physical .h5ad file is missing
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=str(e)
        )

    except RuntimeError as e:
        # Case: Database write failed or Storage write failed
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=str(e)
        )

    except Exception as e:
        # Catch-all for code bugs or unhandled libraries
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"System Error during filtering: {str(e)}"
        )
