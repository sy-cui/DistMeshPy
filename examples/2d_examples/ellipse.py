import matplotlib.pyplot as plt

from pydistmesh import distmesh2d
from pydistmesh.utils import huniform


def mesh_ellipse(
    semi_major_axis: float = 2.0,
    semi_minor_axis: float = 1.0,
):
    def func(p):
        return (
            p[:, 0] ** 2 / semi_major_axis**2 + p[:, 1] ** 2 / semi_minor_axis**2 - 1.0
        )

    return distmesh2d(
        func,
        huniform,
        0.1,
        ((-semi_major_axis, semi_major_axis), (-semi_minor_axis, semi_minor_axis)),
    )


if __name__ == "__main__":
    points, triangles = mesh_ellipse()

    fig, ax = plt.subplots()
    ax.triplot(points[:, 0], points[:, 1], triangles)
    ax.set_aspect("equal")
    plt.show()
