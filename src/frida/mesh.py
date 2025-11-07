from __future__ import annotations
from abc import ABC, abstractmethod

import numpy as np

# import gmsh
import meshio


class PlanarMesh(ABC):
    N_elements: int
    N_nodes: int
    N_edges: int

    @classmethod
    @abstractmethod
    def from_meshio(cls, mesh_io): ...

    @abstractmethod
    def compute_edges(self): ...

    def __repr__(self, prefix: str = "") -> str:
        s = prefix + "\n"
        s += f" N_elements: {self.N_elements}\n"
        s += f" N_nodes: {self.N_nodes}\n"
        s += f" N_edges: {self.N_edges}\n"
        return s


class TriangularMesh(PlanarMesh):
    def __init__(self, nodes, elements):
        self.nodes = nodes
        self.N_nodes = nodes.shape[0]

        self.elements = elements
        assert elements.shape[1] in (3, 6)
        self.N_elements = elements.shape[0]
        self.order = 1 if elements.shape[1] == 3 else 2

        self.edges = self.compute_edges()
        self.N_edges = self.edges.shape[0]

    def compute_edges(self):
        edges_l = []
        for element in self.elements:
            edges_l.extend([[element[i], element[i + 1]] for i in range(2)])

        edges = np.array(edges_l)
        edges = np.sort(edges, axis=1)
        edges = np.unique(edges, axis=0)
        return edges

    @classmethod
    def from_meshio(cls, mesh_io) -> TriangularMesh:
        nodes = mesh_io.points
        cells = mesh_io.cells_dict

        err_msg = f"either triangle3 or triangle6 must be provided"
        assert "triangle6" in cells.keys() or "triangle" in cells.keys(), err_msg
        elements = (
            cells["triangle"] if "triangle" in cells.keys() else cells["triangle6"]
        )

        return cls(nodes=nodes, elements=elements)

    def __repr__(self):
        return super().__repr__(self.__class__.__name__)


# class TetraMesh(PlanarMesh):
#     def __init__(self, nodes, elements): ...

#     def compute_edges(self): ...

#     @classmethod
#     def from_meshio(cls, mesh_io) -> TetraMesh: ...

#     def __repr__(self):
#         return super().__repr__(self.__class__.__name__)


class GmshReader:
    def __init__(self):
        # gmsh.initialize()
        ...

    def read_planar_mesh(self, path: str) -> PlanarMesh:
        mesh_io = meshio.read(path)
        mesh = TriangularMesh.from_meshio(mesh_io)
        print(f"Loaded {len(mesh.nodes)} nodes")
        print(f"Number of element: {len(mesh.elements)}")
        return mesh

    # def finalize(self):
    #     gmsh.finalize()
