import matplotlib.pyplot as plt
from distmeshpy import distmesh3d
from distmeshpy.utils import dsphere, huniform
from distmeshpy.core.tetrahedron_utils import tetrahedron_find_surface

if __name__ == "__main__":
    pos, tet = distmesh3d(
        dist_func=lambda x: dsphere(x, 0.0, 0.0, 0.0, 1.0),
        size_func=huniform,
        edge_size=0.14,
        bounding_box=((-1.1, 1.1), (-1.1, 1.1), (-1.1, 1.1)),
    )
    
    # Toggle this block to show only half of the mesh
    import numpy as np    
    avg_pos = np.mean(pos[tet], axis=1)
    tet = tet[np.where(avg_pos[:, 0] > 0)[0], :]

    # Matplotlib does have good support for plotting tetrahedral mesh.
    # Instead we find the surface triangles and plot them.
    surface_tri = tetrahedron_find_surface(pos, tet)

    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    ax.plot_trisurf(pos[:, 0], pos[:, 1], surface_tri, pos[:, 2])
    ax.set_box_aspect([1, 1, 1])
    plt.show()
