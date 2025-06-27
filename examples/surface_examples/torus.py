from os import minor
import matplotlib.pyplot as plt
import numpy as np

from distmeshpy import distmeshsurface
from distmeshpy.utils import huniform


def mesh_torus(
    major_radius: float = 0.8,
    minor_radius: float = 0.2,
):
    def fd(p):
        # (R - sqrt(x^2 + y^2))^2 + z^2 = r^2
        return (
            np.sum(p**2, axis=1) + major_radius**2 - minor_radius**2
        ) - 2 * major_radius * np.sqrt(p[:, 0] ** 2 + p[:, 1] ** 2)

    tol = 0.1
    bbox = (
        (-major_radius - minor_radius - tol, major_radius + minor_radius + tol),
        (-major_radius - minor_radius - tol, major_radius + minor_radius + tol),
        (-minor_radius - tol, minor_radius + tol),
    )
    return distmeshsurface(fd, huniform, 0.1, bbox)


if __name__ == "__main__":
    points, triangles = mesh_torus()

    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    ax.plot_trisurf(points[:, 0], points[:, 1], triangles, points[:, 2])
    ax.set_aspect("equal")
    plt.show()
