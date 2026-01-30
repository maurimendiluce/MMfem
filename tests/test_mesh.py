import pytest
from MMfem.mesh import (
    rectangle_mesh,
    l_shaped_mesh,
    L_mesh,
    get_boundary_labels,
    boundary_string,
    get_dirichlet_boundaries,
)


# -----------------------------------------------------------------------------
# Rectangle mesh
# -----------------------------------------------------------------------------

def test_rectangle_mesh_boundaries_default():
    mesh = rectangle_mesh(0, 1, 0, 1, h=0.2)
    boundaries = set(mesh.GetBoundaries())

    assert boundaries == {"bottom", "right", "top", "left"}


def test_rectangle_mesh_custom_boundaries():
    mesh = rectangle_mesh(
        0, 1, 0, 1, h=0.2,
        boundary_names={"bottom": "inlet", "top": "outlet"}
    )
    boundaries = set(mesh.GetBoundaries())

    assert boundaries == {"inlet", "right", "outlet", "left"}


def test_rectangle_mesh_invalid_boundary_key():
    with pytest.raises(ValueError):
        rectangle_mesh(
            0, 1, 0, 1, h=0.2,
            boundary_names={"foo": "bar"}
        )


# -----------------------------------------------------------------------------
# L-shaped mesh
# -----------------------------------------------------------------------------

def test_l_shaped_mesh_boundaries():
    mesh = l_shaped_mesh(h=0.3)
    boundaries = set(mesh.GetBoundaries())

    assert boundaries == {"side", "flow_in", "flow_out"}


def test_L_mesh_backward_compatibility():
    mesh1 = l_shaped_mesh(h=0.3)
    mesh2 = L_mesh(h=0.3)

    assert set(mesh1.GetBoundaries()) == set(mesh2.GetBoundaries())


# -----------------------------------------------------------------------------
# Boundary utilities
# -----------------------------------------------------------------------------

def test_get_boundary_labels_returns_list():
    mesh = rectangle_mesh(0, 1, 0, 1, h=0.3)
    labels = get_boundary_labels(mesh)

    assert isinstance(labels, list)
    assert set(labels) == {"floor", "right", "top", "left"}


def test_boundary_string_format():
    mesh = rectangle_mesh(0, 1, 0, 1, h=0.3)
    labels_str = boundary_string(mesh)

    for label in ["floor", "right", "top", "left"]:
        assert label in labels_str
    assert "|" in labels_str


def test_get_dirichlet_boundaries_single_neumann():
    mesh = rectangle_mesh(0, 1, 0, 1, h=0.3)
    dirichlet = get_dirichlet_boundaries(mesh, "top")

    assert "top" not in dirichlet
    for label in ["floor", "right", "left"]:
        assert label in dirichlet


def test_get_dirichlet_boundaries_multiple_neumann():
    mesh = rectangle_mesh(0, 1, 0, 1, h=0.3)
    dirichlet = get_dirichlet_boundaries(mesh, ["top", "right"])

    assert "top" not in dirichlet
    assert "right" not in dirichlet
    assert "floor" in dirichlet
    assert "left" in dirichlet


def test_get_dirichlet_boundaries_invalid_label():
    mesh = rectangle_mesh(0, 1, 0, 1, h=0.3)

    with pytest.raises(ValueError):
        get_dirichlet_boundaries(mesh, "nonexistent")
