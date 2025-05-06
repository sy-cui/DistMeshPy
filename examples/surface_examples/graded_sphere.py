import matplotlib.pyplot as plt

from pydistmesh import distmeshsurface
from pydistmesh.utils import dsphere


def mesh_graded_sphere():
    def fh(p):
        return 0.05 + 0.5 * dsphere(p, 0.0, 0.0, 1.0, 0.0)

    bbox = (
        (-1.1, 1.1),
        (-1.1, 1.1),
        (-1.1, 1.1),
    )
    return distmeshsurface(lambda x: dsphere(x, 0, 0, 0, 1), fh, 0.15, bbox)


if __name__ == "__main__":
    points, triangles = mesh_graded_sphere()

    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    ax.plot_trisurf(points[:, 0], points[:, 1], triangles, points[:, 2])
    ax.set_aspect("equal")
    plt.show()
