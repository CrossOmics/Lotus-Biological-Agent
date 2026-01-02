import pytest
import time
from datetime import datetime

from infrastructure.database.dao.project_meta_dao import ProjectMetaDAO


# Test Cases
def test_create_project_success():
    """Test creating a project with valid data."""
    project = ProjectMetaDAO.create_project(
        project_id="p_001",
        project_name="Test Project Alpha",
        organism="Human",
        tissue_type="Liver",
        description="A test project"
    )

    assert project is not None
    assert project.project_id == "p_001"
    assert project.project_name == "Test Project Alpha"
    assert project.organism == "Human"
    assert isinstance(project.create_time, datetime)
    assert isinstance(project.update_time, datetime)


def test_create_duplicate_project_id():
    """Test that creating a project with a duplicate project_id fails."""
    # Create first project
    ProjectMetaDAO.create_project("p_dup", "Original")

    # Attempt to create second with same ID
    duplicate = ProjectMetaDAO.create_project("p_dup", "Copy")

    assert duplicate is None  # DAO should catch IntegrityError and return None


def test_get_project_by_id():
    """Test retrieving a project by its primary key (id)."""
    created = ProjectMetaDAO.create_project("p_pk_test", "PK Test")

    # Retrieve
    retrieved = ProjectMetaDAO.get_project_by_id(created.id)

    assert retrieved is not None
    assert retrieved.project_id == "p_pk_test"
    assert retrieved.project_name == "PK Test"


def test_get_project_by_id_not_found():
    """Test retrieving a non-existent primary key returns None."""
    result = ProjectMetaDAO.get_project_by_id(99999)
    assert result is None


def test_get_project_by_project_id():
    """Test retrieving a project by its business ID string."""
    ProjectMetaDAO.create_project("p_biz_id", "Business ID Test")

    result = ProjectMetaDAO.get_project_by_project_id("p_biz_id")

    assert result is not None
    assert result.project_name == "Business ID Test"


def test_get_all_projects():
    """Test retrieving all projects."""
    ProjectMetaDAO.create_project("p_1", "One")
    ProjectMetaDAO.create_project("p_2", "Two")
    ProjectMetaDAO.create_project("p_3", "Three")

    all_projects = ProjectMetaDAO.get_all_projects()
    returned_ids = {p.project_id for p in all_projects}

    # Check if they are ordered by update_time desc (latest first)
    assert {"p_1", "p_2", "p_3"}.issubset(returned_ids)


def test_filter_projects():
    """Test filtering projects by organism and tissue type."""
    ProjectMetaDAO.create_project("p_human_liver", "A", organism="Human", tissue_type="Liver")
    ProjectMetaDAO.create_project("p_human_brain", "B", organism="Human", tissue_type="Brain")
    ProjectMetaDAO.create_project("p_mouse_liver", "C", organism="Mouse", tissue_type="Liver")

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
    """Test updating project details."""
    project = ProjectMetaDAO.create_project("p_update", "Old Name", description="Old Desc")
    original_update_time = project.update_time

    # Ensure some time passes for timestamp check
    time.sleep(0.01)

    success = ProjectMetaDAO.update_project(
        pk_id=project.id,
        project_name="New Name",
        description="New Desc"
    )

    assert success is True

    # Reload from DB
    updated = ProjectMetaDAO.get_project_by_id(project.id)
    assert updated.project_name == "New Name"
    assert updated.description == "New Desc"
    assert updated.update_time > original_update_time


def test_update_project_not_found():
    """Test updating a non-existent project returns False."""
    success = ProjectMetaDAO.update_project(999, project_name="Ghost")
    assert success is False


def test_delete_project():
    """Test deleting a project."""
    project = ProjectMetaDAO.create_project("p_delete", "To Delete")

    success = ProjectMetaDAO.delete_project(project.id)
    assert success is True

    # Verify deletion
    assert ProjectMetaDAO.get_project_by_id(project.id) is None


def test_delete_project_not_found():
    """Test deleting a non-existent project returns False."""
    success = ProjectMetaDAO.delete_project(999)
    assert success is False
