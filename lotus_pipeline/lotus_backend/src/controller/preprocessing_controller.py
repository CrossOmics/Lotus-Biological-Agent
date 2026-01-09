from fastapi import APIRouter, HTTPException, status, Depends
from typing import Annotated

from dto.request.calculate_qc_request import CalculateQCRequest
from dto.request.filter_qc_request import FilterQCRequest
from dto.request.run_hvg_request import RunHVGRequest
from dto.response.filter_qc_response import FilterQCResponse
from dto.response.hvg_result_dto import HVGResultDTO
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


@router.post("/qc/filter", response_model=FilterQCResponse, status_code=status.HTTP_201_CREATED)
async def apply_qc_filter(request: FilterQCRequest, service: ServiceDep):
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


@router.post("/hvg", response_model=HVGResultDTO, status_code=status.HTTP_201_CREATED,
             summary="Feature Selection & Scaling",
             description="Performs Normalization -> Log1p -> HVG Identification -> Raw Backup -> Scaling."
             )
async def apply_hvg(request: RunHVGRequest, service: ServiceDep
                    ):
    """
    SFeature Selection (HVG) & Scaling.
    """
    try:
        # Normal business logic to Service
        return service.apply_hvg(request)

    except ValueError as e:
        error_msg = str(e).lower()
        # Map ValueError to 404 if resource is missing, else 400
        if "not found" in error_msg:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"DATASET_NOT_FOUND: {str(e)}"
            )
        else:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail=f"INVALID_PARAMS: {str(e)}"
            )

    except FileNotFoundError as e:
        # Physical file missing despite DB record existence
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"DATASET_NOT_FOUND: Physical file missing - {str(e)}"
        )

    except RuntimeError as e:
        # Scanpy calculation crashes (OOM, Data integrity issues)
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"CALCULATION_ERROR: {str(e)}"
        )

    except Exception as e:
        # Catch-all for unexpected system errors
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"SYSTEM_ERROR: {str(e)}"
        )
