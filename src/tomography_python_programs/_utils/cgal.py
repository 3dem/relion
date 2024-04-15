import numpy as np

from typing import Tuple

from CGAL.CGAL_Polyhedron_3 import Polyhedron_modifier, Polyhedron_3
from CGAL.CGAL_Kernel import Point_3
from CGAL.CGAL_Polygon_mesh_processing import isotropic_remeshing, compute_vertex_normals, remove_isolated_vertices

def _vec2tuple(vec):
    return vec.x(), vec.y(), vec.z()

class Polyhedron(Polyhedron_3):
    def __init__(self, vertices: np.array, indices: np.array):
        super().__init__()
        m = Polyhedron_modifier()

        m.begin_surface(len(vertices), len(indices))
        for vertex in vertices:
            m.add_vertex(Point_3(*vertex))

        for facet in indices:
            m.begin_facet()
            for i in facet:
                    m.add_vertex_to_facet(int(i))
            m.end_facet()

        self.delegate(m)
        m.clear()

    def isotropic_remeshing(self, edge_length: float):
        isotropic_remeshing(self.facets(), edge_length, self, 1)

    def remove_isolated_vertices(self):
        remove_isolated_vertices(self)

    def vertices_array_shape(self) -> Tuple[int, int]:
        return self.size_of_vertices(), 3
    
    def indices_array_shape(self) -> Tuple[int, int]:
        return self.size_of_facets(), 3

    def compute_vertex_normals(self) -> np.array:
        out = np.zeros(self.vertices_array_shape(), dtype=float, order="C")
        l = list()
        compute_vertex_normals(self, l)
        out[:] = list(map(_vec2tuple, l))
        return out

    def vertices_array(self) -> np.array:
        out = np.zeros(self.vertices_array_shape(), dtype=float)
        for i, v in enumerate(self.vertices()):
            p = v.point()
            out[i] = p.x(), p.y(), p.z()
        return out

    def indices_array(self) -> np.array:
        out = np.zeros(self.indices_array_shape(), dtype=int)
        raw_vertices = list(self.vertices())
        for i, f in enumerate(self.facets()):
            halfedge = f.halfedge()
            out[i, 0] = raw_vertices.index(halfedge.vertex())
            out[i, 1] = raw_vertices.index(halfedge.next().vertex())
            out[i, 2] = raw_vertices.index(halfedge.next().next().vertex())
        return out