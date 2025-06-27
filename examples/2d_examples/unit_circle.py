import matplotlib.pyplot as plt

from distmeshpy import distmesh2d
from distmeshpy.utils import dcircle, huniform


def mesh_unit_circle():
    return distmesh2d(
        lambda p: dcircle(p, 0.0, 0.0, 1.0),  # Distance function for a unit circle
        huniform,  # Uniform size function
        0.2,  # Initial edge size
        ((-1.1, 1.1), (-1.1, 1.1)),  # Bounding box
    )


if __name__ == "__main__":
    points, triangles = mesh_unit_circle()

    fig, ax = plt.subplots()
    ax.triplot(points[:, 0], points[:, 1], triangles)
    ax.set_aspect("equal")
    plt.show()
