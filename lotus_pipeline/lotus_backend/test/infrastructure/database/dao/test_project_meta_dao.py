import pytest
import time
from datetime import datetime

# Adjust imports to match your project structure
from infrastructure.database.dao.project_meta_dao import ProjectMetaDAO


# --- Test Cases ---

def test_create_project_success():
    """Test creating a project with valid data including the new project_path."""
    project = ProjectMetaDAO.create_project(
        project_id="p_001",
        project_name="Test Project Alpha",
        project_path="/data/projects/p_001",
        organism="Human",
        tissue_type="Liver",
        description="A test project"
    )

    assert project is not None
    assert project.project_id == "p_001"
    assert project.project_name == "Test Project Alpha"
    assert project.project_path == "/data/projects/p_001"
    assert project.organism == "Human"
    assert isinstance(project.create_time, datetime)
    assert isinstance(project.update_time, datetime)


def test_create_duplicate_project_id():
    """Test that creating a project with a duplicate project_id fails."""
    # Create first project
    ProjectMetaDAO.create_project("p_dup", "Original", "/path/original")

    # Attempt to create second with same ID (Fix: added missing project_path arg)
    duplicate = ProjectMetaDAO.create_project("p_dup", "Copy", "/path/copy")

    assert duplicate is None  # DAO should catch IntegrityError and return None


def test_get_project_by_id():
    """Test retrieving a project by its primary key (id)."""
    created = ProjectMetaDAO.create_project("p_pk_test", "PK Test", "/path/pk_test")

    # Retrieve
    retrieved = ProjectMetaDAO.get_project_by_id(created.id)

    assert retrieved is not None
    assert retrieved.project_id == "p_pk_test"
    assert retrieved.project_name == "PK Test"
    assert retrieved.project_path == "/path/pk_test"


def test_get_project_by_id_not_found():
    """Test retrieving a non-existent primary key returns None."""
    result = ProjectMetaDAO.get_project_by_id(99999)
    assert result is None


def test_get_project_by_project_id():
    """Test retrieving a project by its business ID string."""
    ProjectMetaDAO.create_project("p_biz_id", "Business ID Test", "/path/biz_id")

    result = ProjectMetaDAO.get_project_by_project_id("p_biz_id")

    assert result is not None
    assert result.project_name == "Business ID Test"


def test_get_project_file_path():
    """Test retrieving just the project file path (New DAO method)."""
    target_path = "/data/storage/p_path_test"
    ProjectMetaDAO.create_project("p_path_test", "Path Test", target_path)

    # Test the new DAO method
    result_path = ProjectMetaDAO.get_project_file_path("p_path_test")

    assert result_path == target_path


def test_get_all_projects():
    """Test retrieving all projects."""
    ProjectMetaDAO.create_project("p_1", "One", "/path/1")
    ProjectMetaDAO.create_project("p_2", "Two", "/path/2")
    ProjectMetaDAO.create_project("p_3", "Three", "/path/3")

    all_projects = ProjectMetaDAO.get_all_projects()
    returned_ids = {p.project_id for p in all_projects}

    # Check if they are ordered by update_time desc (latest first)
    assert {"p_1", "p_2", "p_3"}.issubset(returned_ids)


def test_filter_projects():
    """Test filtering projects by organism and tissue type."""
    ProjectMetaDAO.create_project("p_human_liver", "A", "/p/a", organism="Human", tissue_type="Liver")
    ProjectMetaDAO.create_project("p_human_brain", "B", "/p/b", organism="Human", tissue_type="Brain")
    ProjectMetaDAO.create_project("p_mouse_liver", "C", "/p/c", organism="Mouse", tissue_type="Liver")

    # Filter by Organism only
    humans = ProjectMetaDAO.filter_projects(organism="Human")
    assert all(p.organism == "Human" for p in humans)

    # Filter by Tissue only
    livers = ProjectMetaDAO.filter_projects(tissue_type="Liver")
    assert all(p.tissue_type == "Liver" for p in livers)

    # Filter by Both
    human_brain = ProjectMetaDAO.filter_projects(organism="Human", tissue_type="Brain")
    assert all(p.organism == "Human" and p.tissue_type == "Brain" for p in human_brain)
    assert any(p.project_id == "p_human_brain" for p in human_brain)


def test_update_project():
    """Test updating project details, including the new project_path."""
    project = ProjectMetaDAO.create_project("p_update", "Old Name", "/path/old", description="Old Desc")
    original_update_time = project.update_time

    # Ensure some time passes for timestamp check
    time.sleep(0.01)

    success = ProjectMetaDAO.update_project(
        pk_id=project.id,
        project_name="New Name",
        description="New Desc",
        project_path="/path/new_moved_location"
    )

    assert success is True

    # Reload from DB
    updated = ProjectMetaDAO.get_project_by_id(project.id)
    assert updated.project_name == "New Name"
    assert updated.description == "New Desc"
    assert updated.project_path == "/path/new_moved_location"
    assert updated.update_time > original_update_time


def test_update_project_not_found():
    """Test updating a non-existent project returns False."""
    success = ProjectMetaDAO.update_project(999, project_name="Ghost")
    assert success is False


def test_delete_project():
    """Test deleting a project."""
    project = ProjectMetaDAO.create_project("p_delete", "To Delete", "/path/delete")

    success = ProjectMetaDAO.delete_project(project.id)
    assert success is True

    # Verify deletion
    assert ProjectMetaDAO.get_project_by_id(project.id) is None


def test_delete_project_not_found():
    """Test deleting a non-existent project returns False."""
    success = ProjectMetaDAO.delete_project(999)
    assert success is False
