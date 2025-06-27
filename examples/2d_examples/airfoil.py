import matplotlib.pyplot as plt
import numpy as np

from distmeshpy import distmesh2d
from distmeshpy.utils import dcircle, ddiff


def mesh_airfoil():
    """Generate mesh for NACA0012 airfoil based on the example from MATLAB."""
    hlead = 0.01
    htrail = 0.04
    hmax = 2
    circx = 2
    circr = 4
    a = 0.12 / 0.2 * np.array([0.2969, -0.1260, -0.3516, 0.2843, -0.1036])
    arr = np.array([a[4], a[3], a[2], a[1], 0.0])

    def fd(pts, **kwargs):
        return ddiff(
            dcircle(pts, circx, 0, circr),
            (np.abs(pts[:, 1]) - np.polyval(arr, pts[:, 0])) ** 2
            - a[0] ** 2 * pts[:, 0],
        )

    def fh(pts, **kwargs):
        return np.minimum(
            np.minimum(
                hlead + 0.3 * dcircle(pts, 0, 0, 0),
                htrail + 0.3 * dcircle(pts, 1, 0, 0),
            ),
            hmax,
        )

    fixx = 1 - htrail * np.cumsum(1.3 ** np.arange(0, 5))
    fixy = a[0] * np.sqrt(fixx) + np.polyval(arr, fixx)
    fixx = fixx.reshape(-1, 1)
    fixy = fixy.reshape(-1, 1)
    fix = np.array(
        [
            [circx - circr, 0.0],
            [circx + circr, 0.0],
            [circx, -circr],
            [circx, circr],
            [0.0, 0.0],
            [1.0, 0.0],
        ]
    )
    fix = np.vstack((fix, np.block([[fixx, fixy], [fixx, -fixy]])))
    box = ((circx - circr, circx + circr), (-circr, circr))
    h0 = min([hlead, htrail, hmax])

    return distmesh2d(fd, fh, h0, box, fix)


if __name__ == "__main__":
    points, triangles = mesh_airfoil()

    fig, ax = plt.subplots()
    ax.triplot(points[:, 0], points[:, 1], triangles)
    ax.set_aspect("equal")
    plt.show()
