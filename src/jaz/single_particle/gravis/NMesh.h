#ifndef __LIBGRAVIS_NMESH_H__
#define __LIBGRAVIS_NMESH_H__
/******************************************************************************
**        Title: NMesh.h
**  Description: N-Sided Mesh representation.
**
**       Author: Jean-Sebastien Pierrard, 2005
**               Brian Amberg
**               Computer Science Department, University Basel (CH)
**
******************************************************************************/

#include <vector>
#include "tArray.h"
#include "tRGBA.h"
#include "t2Vector.h"
#include "t3Vector.h"
#include "tImage.h"
#include "NTuple.h"
#include "Mesh.h"

namespace gravis
{

  /*! \brief N-Mesh data structure.
   *
   * The below is more or less correct, but we may now have double
   * precision values and faces of arbitrary (but equal) length. (Usefull
   * for multicube)
   *
   * A Mesh contains vertex, normal, texture coordinate (uvw) and material information.
   * For the three types of primitives (triangle, line, point) there are index arrays
   * referencing above information. For example for lines, lvi indexes into
   * vertex, and lti into texture coordinates. The vertices and colors
   * for the 4th lines in the mesh are then vertex[lvi[3][0]], vertex[lvi[3][1]],
   * color[lci[3][0]] and color[lci[3][1]].
   *
   * fvi.size(), lvi.size() and pvi.size() implicitly specify how many triangles, lines
   * and points there are in the mesh. All other index arrays must either be of the
   * same length as the corresponding vertex index array, or of length 0.
   *
   * How is missing information handled? If for example no normals are assigned to
   * any triangles, fni.size() would be zero. If normals are assigned for some triangles,
   * but not for others, the fni-tuples for the respective triangles must have entries
   * of -1 (which is the 'invalid index').
   */
  template <class T, size_t N>
  class NMesh
  {
      typedef t3Vector<T> Vector;
      typedef NTuple<int, N> Tuple;
    public:

      tArray<Vector> vertex; /*!< \brief  Vertex array. */
      tArray<Vector> normal; /*!< \brief  Normal array. */
      tArray<Vector> texcrd; /*!< \brief  Texture coordinate array. */
      tArray<fRGBA> color; /*!< \brief  Color array. */
      std::vector<Material> material; /*!< \brief  Material array. */

      tArray<Tuple> fvi; /*!< \brief  Face vertex indices. */
      tArray<Tuple> fni; /*!< \brief  Face normal indices. */
      tArray<Tuple> fti; /*!< \brief  Face texcrd indices. */
      tArray<Tuple> fci; /*!< \brief  Face color indices. */
      tArray<int> fmi;    /*!< \brief  Face material indices. */

      tArray<Tuple2> lvi; /*!< \brief  Line vertex indices. */
      tArray<Tuple2> lti; /*!< \brief  Line texcrd indices. */
      tArray<Tuple2> lci; /*!< \brief  Line texcrd indices. */

      tArray<int> pvi; /*!< \brief  Point vertex indices. */
      tArray<int> pci; /*!< \brief  Point color indices. */

      tArray<Tuple> adjacent; /*!< \brief Adjacency list. See generateAdjacencyList(). */

      // Create a deep copy of the mesh
      void clone(Mesh& out) const
      {
        out.vertex   = vertex.clone();
        out.normal   = normal.clone();
        out.texcrd   = texcrd.clone();
        out.color    = color.clone();
        //out.material = material.save_clone();
        out.material = material;

        out.fvi = fvi.clone();
        out.fni = fni.clone();
        out.fti = fti.clone();
        out.fci = fci.clone();
        out.fmi = fmi.clone();

        out.lvi = lvi.clone();
        out.lti = lti.clone();
        out.lci = lci.clone();

        out.pvi = pvi.clone();
        out.pci = pci.clone();

        out.adjacent = adjacent.clone();
      }
  };

  typedef NMesh<double, 3> d3Mesh;
  typedef NMesh<double, 4> d4Mesh;
  typedef NMesh<float, 3> f3Mesh;
  typedef NMesh<float, 4> f4Mesh;

  // This should work nicely, as we have binary compatibility
  const Mesh& convert(const NMesh<float, 3> &in)
  {
    return *reinterpret_cast<Mesh*>(&in);
  }

} // namespace gravis

#endif
