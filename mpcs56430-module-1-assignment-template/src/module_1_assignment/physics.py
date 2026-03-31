"""Physics utilities for simple kinematic simulations."""


def compute_kinetic_energy(mass: float, velocity: float) -> float:
    """Calculate kinetic energy (J): 0.5 * mass * velocity^2.

    Args:
        mass: Mass in kilograms (must be non-negative).
        velocity: Velocity in meters per second.

    Returns:
        Kinetic energy in joules.
    """
    if mass < 0:
        raise ValueError("mass must be non-negative")
    return 0.5 * mass * (velocity ** 2)
