import pytest

from module_1_assignment.physics import compute_kinetic_energy


def test_compute_kinetic_energy_standard() -> None:
    assert compute_kinetic_energy(10.0, 5.0) == pytest.approx(125.0)


def test_compute_kinetic_energy_zero_velocity() -> None:
    assert compute_kinetic_energy(100.0, 0.0) == pytest.approx(0.0)


def test_compute_kinetic_energy_negative_velocity() -> None:
    assert compute_kinetic_energy(2.0, -3.0) == pytest.approx(9.0)


def test_compute_kinetic_energy_negative_mass_rejected() -> None:
    with pytest.raises(ValueError, match="non-negative"):
        compute_kinetic_energy(-1.0, 3.0)
