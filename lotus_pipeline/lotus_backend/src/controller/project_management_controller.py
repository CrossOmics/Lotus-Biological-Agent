from fastapi import APIRouter, HTTPException, status

from dto.request.create_project_request import CreateProjectRequest
from dto.response.project_response import ProjectResponse
# Service Imports
from service.dataset_service import DatasetService
from infrastructure.filesystem.constants.filesystem_constants import USER_PROJECT_ROOT
from service.project_management_service import ProjectManagementService

router = APIRouter(prefix="/api/v1/projects/import", tags=["Project Management"])


@router.post("/", response_model=ProjectResponse, status_code=status.HTTP_201_CREATED)
async def create_new_project(request: CreateProjectRequest):
    """
    Creates a new project and ingests the initial raw dataset.

    Workflow:
    1. Create Project Metadata in DB.
    2. Trigger Dataset Service to ingest raw data (Copy/Convert/Save).
    3. Link Dataset to Project in DB.
    """

    # Instantiate Services
    # (In a larger app, use Dependency Injection via Depends())
    project_service = ProjectManagementService()
    dataset_service = DatasetService()

    try:
        # Create Project (Business Logic)
        # This service creates the 'project_meta' record and project folder structure
        new_project = project_service.create_project(
            project_name=request.project_name,
            organism=request.organism,
            tissue_type=request.tissue_type,
            description=request.description
        )

        if not new_project:
            raise HTTPException(
                status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                detail="Failed to create project metadata."
            )

        # Import Dataset (IO + DB)
        # We pass the ProjectMeta object (new_project) to link the FK automatically
        new_dataset = dataset_service.import_dataset_from_local(
            local_file_path=request.local_file_path,
            project=new_project,
            dataset_name=None  # Use original filename by default
        )

        if not new_dataset:
            # In a real app, consider rolling back (deleting) the project here
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail="Failed to import dataset. Check file path or format."
            )

        # 4. Return Success Response
        return ProjectResponse(
            project_id=new_project.project_id,
            project_name=new_project.project_name,
            dataset_id=new_dataset.dataset_id,
            dataset_name=new_dataset.dataset_name,
            workspace_path=new_project.project_path
        )

    except FileNotFoundError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        # Catch-all for unexpected errors
        raise HTTPException(status_code=500, detail=f"System Error: {str(e)}")
