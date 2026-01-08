from fastapi import APIRouter, HTTPException, status, Depends
from typing import Annotated

from dto.request.calculate_qc_request import CalculateQCRequest
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
