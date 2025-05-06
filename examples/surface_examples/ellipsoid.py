import matplotlib.pyplot as plt
import numpy as np

from pydistmesh import distmeshsurface
from pydistmesh.utils import huniform


def mesh_ellipsoid(
    axis_1: float = 2.0,
    axis_2: float = 1.5,
    axis_3: float = 1.0,
):
    def fd(p):
        return (
            p[:, 0] ** 2 / axis_1**2
            + p[:, 1] ** 2 / axis_2**2
            + p[:, 2] ** 2 / axis_3**2
            - 1.0
        )

    tol = 0.1
    bbox = (
        (-axis_1 - tol, axis_1 + tol),
        (-axis_2 - tol, axis_2 + tol),
        (-axis_3 - tol, axis_3 + tol),
    )
    return distmeshsurface(fd, huniform, 0.2, bbox)


if __name__ == "__main__":
    points, triangles = mesh_ellipsoid()

    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    ax.plot_trisurf(points[:, 0], points[:, 1], triangles, points[:, 2])
    ax.set_aspect("equal")
    plt.show()
