import matplotlib.pyplot as plt

from pydistmesh import distmeshsurface
from pydistmesh.utils import huniform


def mesh_egg():
    a = 1
    b = a / 1.5
    w = 0.2 * a

    def fd(p):
        return (
            p[:, 2] ** 2 / a**2
            + (p[:, 0] ** 2 + p[:, 1] ** 2)
            * (a**2 + 2 * w * p[:, 2] + w**2)
            / (b * a) ** 2
            - 1
        )

    bbox = (
        (-b - 0.1, b + 0.1),
        (-b - 0.1, b + 0.1),
        (-a - 0.1, a + 0.1),
    )
    return distmeshsurface(fd, huniform, 0.1, bbox)


if __name__ == "__main__":
    points, triangles = mesh_egg()

    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    ax.plot_trisurf(points[:, 0], points[:, 1], triangles, points[:, 2])
    ax.set_aspect("equal")
    plt.show()
