import numpy as np
import matplotlib.pyplot as plt

from pydistmesh import distmesh2d
from pydistmesh.utils import dcircle, drectangle, ddiff


def mesh_hollow_rectangle():
    def fd(p):
        return ddiff(drectangle(p, -1, 1, -1, 1), dcircle(p, 0, 0, 0.5))

    def fh(p):
        return 0.05 + 0.3 * dcircle(p, 0, 0, 0.5)

    return distmesh2d(
        fd,
        fh,
        0.05,
        ((-1.0, 1.0), (-1.0, 1.0)),
        np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]], dtype=np.float64),
    )


if __name__ == "__main__":
    points, triangles = mesh_hollow_rectangle()

    fig, ax = plt.subplots()
    ax.triplot(points[:, 0], points[:, 1], triangles)
    ax.set_aspect("equal")
    plt.show()
