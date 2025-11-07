from __future__ import annotations
import pytest
import numpy as np
import shutil
import os

from frida.mesh import GmshReader

# def test_gmsh_square_regular():
#     # The Python API is entirely defined in the `gmsh.py' module (which contains the
#     # full documentation of all the functions in the API):
#     import gmsh
#     import sys

#     # Before using any functions in the Python API, Gmsh must be initialized:
#     gmsh.initialize()

#     # Next we add a new model named "t1" (if gmsh.model.add() is not called a new
#     # unnamed model will be created on the fly, if necessary):
#     gmsh.model.add("t1")

#     # The Python API provides direct access to each supported geometry (CAD)
#     # kernel. The built-in kernel is used in this first tutorial: the corresponding
#     # API functions have the `gmsh.model.geo' prefix.

#     # The first type of `elementary entity' in Gmsh is a `Point'. To create a point
#     # with the built-in CAD kernel, the Python API function is
#     # gmsh.model.geo.addPoint():
#     # - the first 3 arguments are the point coordinates (x, y, z)
#     # - the next (optional) argument is the target mesh size close to the point
#     # - the last (optional) argument is the point tag (a stricly positive integer
#     #   that uniquely identifies the point)
#     lc = 1e-2
#     gmsh.model.geo.addPoint(0, 0, 0, lc, 1)

#     # Note that in addition to the default ``camelCase'' function name `addPoint',
#     # the Python API also defines a ``snake case'' alias, i.e. `add_point'. You can
#     # use either interchangeably; all the tutorials are written using the camelCase
#     # convention for consistency.

#     # The distribution of the mesh element sizes will be obtained by interpolation
#     # of these mesh sizes throughout the geometry. Another method to specify mesh
#     # sizes is to use general mesh size Fields (see `t10.py'). A particular case is
#     # the use of a background mesh (see `t7.py').
#     #
#     # If no target mesh size of provided, a default uniform coarse size will be used
#     # for the model, based on the overall model size.
#     #
#     # We can then define some additional points. All points should have different
#     # tags:
#     gmsh.model.geo.addPoint(.1, 0, 0, lc, 2)
#     gmsh.model.geo.addPoint(.1, .3, 0, lc, 3)

#     # If the tag is not provided explicitly, a new tag is automatically created, and
#     # returned by the function:
#     p4 = gmsh.model.geo.addPoint(0, .3, 0, lc)

#     # Curves are Gmsh's second type of elementery entities, and, amongst curves,
#     # straight lines are the simplest. The API to create straight line segments with
#     # the built-in kernel follows the same conventions: the first 2 arguments are
#     # point tags (the start and end points of the line), and the last (optional one)
#     # is the line tag.
#     #
#     # In the commands below, for example, the line 1 starts at point 1 and ends at
#     # point 2.
#     #
#     # Note that curve tags are separate from point tags - hence we can reuse tag `1'
#     # for our first curve. And as a general rule, elementary entity tags in Gmsh
#     # have to be unique per geometrical dimension.
#     gmsh.model.geo.addLine(1, 2, 1)
#     gmsh.model.geo.addLine(3, 2, 2)
#     gmsh.model.geo.addLine(3, p4, 3)
#     gmsh.model.geo.addLine(4, 1, p4)

#     # The third elementary entity is the surface. In order to define a simple
#     # rectangular surface from the four curves defined above, a curve loop has first
#     # to be defined. A curve loop is defined by an ordered list of connected curves,
#     # a sign being associated with each curve (depending on the orientation of the
#     # curve to form a loop). The API function to create curve loops takes a list
#     # of integers as first argument, and the curve loop tag (which must be unique
#     # amongst curve loops) as the second (optional) argument:
#     gmsh.model.geo.addCurveLoop([4, 1, -2, 3], 1)

#     # We can then define the surface as a list of curve loops (only one here,
#     # representing the external contour, since there are no holes--see `t4.py' for
#     # an example of a surface with a hole):
#     gmsh.model.geo.addPlaneSurface([1], 1)

#     # Before they can be meshed (and, more generally, before they can be used by API
#     # functions outside of the built-in CAD kernel functions), the CAD entities must
#     # be synchronized with the Gmsh model, which will create the relevant Gmsh data
#     # structures. This is achieved by the gmsh.model.geo.synchronize() API call for
#     # the built-in CAD kernel. Synchronizations can be called at any time, but they
#     # involve a non trivial amount of processing; so while you could synchronize the
#     # internal CAD data after every CAD command, it is usually better to minimize
#     # the number of synchronization points.
#     gmsh.model.geo.synchronize()

#     # At this level, Gmsh knows everything to display the rectangular surface 1 and
#     # to mesh it. An optional step is needed if we want to group elementary
#     # geometrical entities into more meaningful groups, e.g. to define some
#     # mathematical ("domain", "boundary"), functional ("left wing", "fuselage") or
#     # material ("steel", "carbon") properties.
#     #
#     # Such groups are called "Physical Groups" in Gmsh. By default, if physical
#     # groups are defined, Gmsh will export in output files only mesh elements that
#     # belong to at least one physical group. (To force Gmsh to save all elements,
#     # whether they belong to physical groups or not, set the `Mesh.SaveAll' option
#     # to 1.) Physical groups are also identified by tags, i.e. stricly positive
#     # integers, that should be unique per dimension (0D, 1D, 2D or 3D). Physical
#     # groups can also be given names.
#     #
#     # Here we define a physical curve that groups the left, bottom and right curves
#     # in a single group (with prescribed tag 5); and a physical surface with name
#     # "My surface" (with an automatic tag) containing the geometrical surface 1:
#     gmsh.model.addPhysicalGroup(1, [1, 2, 4], 5)
#     gmsh.model.addPhysicalGroup(2, [1], name="My surface")

#     # We can then generate a 2D mesh...
#     gmsh.model.mesh.generate(2)

#     # ... and save it to disk
#     gmsh.write("t1.msh")

#     # Remember that by default, if physical groups are defined, Gmsh will export in
#     # the output mesh file only those elements that belong to at least one physical
#     # group. To force Gmsh to save all elements, you can use
#     #
#     # gmsh.option.setNumber("Mesh.SaveAll", 1)

#     # By default, Gmsh saves meshes in the latest version of the Gmsh mesh file
#     # format (the `MSH' format). You can save meshes in other mesh formats by
#     # specifying a filename with a different extension. For example
#     #
#     #   gmsh.write("t1.unv")
#     #
#     # will save the mesh in the UNV format. You can also save the mesh in older
#     # versions of the MSH format: simply set
#     #
#     #   gmsh.option.setNumber("Mesh.MshFileVersion", x)
#     #
#     # for any version number `x'. As an alternative, you can also not specify the
#     # format explicitly, and just choose a filename with the `.msh2' or `.msh4'
#     # extension.

#     # To visualize the model we can run the graphical user interface with
#     # `gmsh.fltk.run()'. Here we run it only if "-nopopup" is not provided in the
#     # command line arguments:
#     if '-nopopup' not in sys.argv:
#         gmsh.fltk.run()

#     # Note that starting with Gmsh 3.0, models can be built using other geometry
#     # kernels than the default "built-in" kernel. To use the OpenCASCADE CAD kernel
#     # instead of the built-in kernel, you should use the functions with the
#     # `gmsh.model.occ' prefix.
#     #
#     # Different CAD kernels have different features. With OpenCASCADE, instead of
#     # defining the surface by successively defining 4 points, 4 curves and 1 curve
#     # loop, one can define the rectangular surface directly with
#     #
#     # gmsh.model.occ.addRectangle(.2, 0, 0, .1, .3)
#     #
#     # After synchronization with the Gmsh model with
#     #
#     # gmsh.model.occ.synchronize()
#     #
#     # the underlying curves and points could be accessed with
#     # gmsh.model.getBoundary().
#     #
#     # See e.g. `t16.py', `t18.py', `t19.py' or `t20.py' for complete examples based
#     # on OpenCASCADE, and `examples/api' for more.

#     # This should be called when you are done using the Gmsh Python API:
#     gmsh.finalize()




@pytest.mark.parametrize("order", [1, 2])
def test_gmsh_from_geo_file(order: int):
    path_geo = "tests/fixtures/test_circle.geo"
    
    command = f"gmsh -2 -order {order} {path_geo}"
    path_msh = f"{path_geo.replace(".geo", ".msh")}"
    if os.path.exists(path_msh):
        os.remove(path_msh)
    
    os.system(command)
    
    # os.system(f"gmsh -2 -order {order} -format m {path_geo}")
    assert os.path.exists(path_msh)
    
    reader = GmshReader()
    mesh = reader.read_planar_mesh(path=path_msh)
    
    if os.path.exists(path_msh):
        os.remove(path_msh)
        
    
    # import pyvista as pv
    # import meshio
    
    # points = mesh.nodes
    
    # strips = np.hstack((3+np.zeros((mesh.elements.shape[0],1), dtype=np.uint8), mesh.elements)).reshape(-1,1).squeeze().tolist()
    # mesh_pv = pv.PolyData(points, strips=strips)
    
    
    # # Create an array of points
    # points = np.array(
    #     [
    #         [1.0, 0.0, 0.0],
    #         [0.0, 0.0, 0.0],
    #         [1.0, 1.0, 0.0],
    #         [0.0, 1.0, 0.0],
    #         [1.0, 2.0, 0.0],
    #         [0.0, 2.0, 0.0],
    #         [1.0, 3.0, 0.0],
    #         [0.0, 3.0, 0.0],
    #     ],
    # )
    # strips = np.array([8, 0, 1, 2, 3, 4, 5, 6, 7])
    
    # pl = pv.Plotter()
    # pl.add_mesh(mesh_pv, show_edges=True)
    # pl.add_point_labels(mesh_pv.points, range(mesh_pv.n_points))
    # pl.camera_position = 'yx'
    # pl.camera.zoom(1.2)
    # pl.show()
    # pl.close()
    
    
    
    # mesh = meshio.read(path_msh)

    # print(mesh.points)        # N x 3 numpy array of node coordinates
    # print(mesh.cells_dict)    # Dictionary of element types and connectivity

    # # strips = mesh.
    
    # pl = pv.Plotter()
    # pl.add_mesh(mesh, show_edges=True)
    # pl.add_point_labels(mesh.points, range(mesh.n_points))
    # pl.camera_position = 'yx'
    # pl.camera.zoom(1.2)
    # pl.show()
    
    
    # mesh = pyvista.Plane(i_resolution=10, j_resolution=10)
    # mesh.cell_data["colors"] = [
    #     [1.0, 0, 0],
    #     [0, 1, 0],
    #     [0, 0, 1],
    #     [1, 1, 0]
    # ]

    # plotter = pyvista.Plotter()
    # plotter.add_mesh(mesh, scalars="colors", rgb=True)
    # plotter.view_xy()
    # plotter.show()
    
    
    