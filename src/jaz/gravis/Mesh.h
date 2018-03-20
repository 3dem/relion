#ifndef __LIBGRAVIS_MESH_H__
#define __LIBGRAVIS_MESH_H__
/******************************************************************************
**        Title: Mesh.h
**  Description: Mesh representation.
**
**       Author: Jean-Sebastien Pierrard, 2005
**               Computer Science Department, University Basel (CH)
**
******************************************************************************/

#include "tArray.h"
#include "tRGBA.h"
#include "t2Vector.h"
#include "t3Vector.h"
#include "tImage.h"
#include "Tuple.h"
#include <vector>

namespace gravis
{

  class Material
  {
    public:
      Material (std::string n="")
        : name(n), ambient(0.1,1.0), diffuse(0.5,1.0), specular(1.0,1.0),
          shininess(10.0), hasTexture(false), textureName(), hasEnvMap(false), envMapName(), hasNormalMap(false), normalMapName() {}

      std::string name;

      fRGBA ambient;
      fRGBA diffuse;
      fRGBA specular;

      float shininess; /*!< \brief Phong exponent. */

      bool hasTexture; /*!< \brief whether a (diffuse) texture is defined for this material. */
      std::string textureName; /*!< \brief Filename of the (diffuse) texture. */

      bool hasEnvMap;
      std::string envMapName;

      bool hasNormalMap;
      std::string normalMapName;
  };

  /*! \brief Mesh data structure.
  *
  * A Mesh contains vertex, normal, texture coordinate (uvw) and material information.
  * For the three types of primitives (triangle, line, point) there are index arrays
  * referencing above information. For example for lines, lvi indexes into
  * vertex, and lti into texture coordinates. The vertices and colors
  * for the 4th lines in the mesh are then vertex[lvi[3][0]], vertex[lvi[3][1]],
  * color[lci[3][0]] and color[lci[3][1]].
  *
  * tvi.size(), lvi.size() and pvi.size() implicitly specify how many triangles, lines
  * and points there are in the mesh. All other index arrays must either be of the
  * same length as the corresponding vertex index array, or of length 0.
  *
  * How is missing information handled? If for example no normals are assigned to
  * any triangles, tni.size() would be zero. If normals are assigned for some triangles,
  * but not for others, the tni-tuples for the respective triangles must have entries
  * of -1 (which is the 'invalid index').
  */
  class Mesh
  {
    public:

      tArray<f3Vector> vertex; /*!< \brief  Vertex array. */
      tArray<f3Vector> normal; /*!< \brief  Normal array. */
      tArray<f3Vector> texcrd; /*!< \brief  Texture coordinate array. */
      tArray<fRGBA> color; /*!< \brief  Color array. */
      std::vector<Material> material; /*!< \brief  Material array. */

      tArray<Tuple3> tvi; /*!< \brief  Triangle vertex indices. */
      tArray<Tuple3> tni; /*!< \brief  Triangle normal indices. */
      tArray<Tuple3> tti; /*!< \brief  Triangle texcrd indices. */
      tArray<Tuple3> tci; /*!< \brief  Triangle color indices. */
      tArray<int> tmi;    /*!< \brief  Triangle material indices. */

      tArray<Tuple2> lvi; /*!< \brief  Line vertex indices. */
      tArray<Tuple2> lti; /*!< \brief  Line texcrd indices. */
      tArray<Tuple2> lci; /*!< \brief  Line texcrd indices. */

      tArray<int> pvi; /*!< \brief  Point vertex indices. */
      tArray<int> pci; /*!< \brief  Point color indices. */

      tArray<Tuple3> adjacent; /*!< \brief Adjacency list. See generateAdjacencyList(). */

      Mesh() :
        vertex(),
        normal(),
        texcrd(),
        color(),
        material(),
        tvi(),
        tni(),
        tti(),
        tci(),
        tmi(),
        lvi(),
        lti(),
        lci(),
        pvi(),
        pci(),
        adjacent() {}

      // Create a deep copy of the mesh
      void clone(Mesh& out) const
      {
        out.vertex   = vertex.clone();
        out.normal   = normal.clone();
        out.texcrd   = texcrd.clone();
        out.color    = color.clone();
        //out.material = material.save_clone();
        out.material = material;

        out.tvi = tvi.clone();
        out.tni = tni.clone();
        out.tti = tti.clone();
        out.tci = tci.clone();
        out.tmi = tmi.clone();

        out.lvi = lvi.clone();
        out.lti = lti.clone();
        out.lci = lci.clone();

        out.pvi = pvi.clone();
        out.pci = pci.clone();

        out.adjacent = adjacent.clone();
      }

      void generateNormals()
      {
        const int numFaces = tvi.size();
        tni.setSize(numFaces);
        normal.setSize(numFaces);

        for (int i = 0; i < numFaces; i++)
        {
          f3Vector a = (vertex[tvi[i][1]] - vertex[tvi[i][0]]);
          f3Vector b = (vertex[tvi[i][2]] - vertex[tvi[i][0]]);
          normal[i] = cross(a, b).normalize();
          tni[i] = Tuple3(i, i, i);
        }
      }

      void generatePerVertexNormals(unsigned int propagations=0)
      {
        tArray<int> ncount;
        f3Vector norm;
        const int numFaces = tvi.size();
        tni.setSize(numFaces);

        normal.setSize(vertex.size());
        ncount.setSize(vertex.size());
        for (unsigned int i = 0; i < ncount.size(); i++) ncount[i] = 0;
        for (unsigned int i = 0; i < normal.size(); i++) normal[i] = f3Vector(0.0f,0.0f,0.0f);


        for (int i = 0; i < numFaces; i++)
        {
          if(tvi[i].c0 < 0 || tvi[i].c1 < 0 || tvi[i].c2 < 0) continue;
          f3Vector a = (vertex[tvi[i][1]] - vertex[tvi[i][0]]);
          f3Vector b = (vertex[tvi[i][2]] - vertex[tvi[i][0]]);
          norm = cross(a, b).normalize();
          tni[i] = tvi[i];
          normal[tvi[i][0]] += norm;
          normal[tvi[i][1]] += norm;
          normal[tvi[i][2]] += norm;
          ncount[tvi[i][0]]++;
          ncount[tvi[i][1]]++;
          ncount[tvi[i][2]]++;
        }
        for (unsigned int i = 0; i < normal.size(); i++)
        {
          if(ncount[i] != 0)
            normal[i] /= ncount[i];
          normal[i] = normal[i].normalize();
        }

        tArray<f3Vector> nnormal;
        nnormal.setSize(ncount.size());
        for(unsigned int j=0; j<propagations; ++j)
        {
          for (unsigned int i = 0; i < ncount.size(); i++)
          {
            nnormal[i] = normal[i]*2.0;
            ncount[i]  = 2;
          }
          for (int i = 0; i < numFaces; i++)
          {
            if(tvi[i].c0 < 0 || tvi[i].c1 < 0 || tvi[i].c2 < 0) continue;
            for (int k = 0; k < 3; k++)
            {
              f3Vector norm = normal[tvi[i][k]] + normal[tvi[i][(k+2)%3]];
              nnormal[tvi[i][(k+1)%3]] += norm;
              ncount[tvi[i][(k+1)%3]]  += 2;
            }
          }
          for (unsigned int i = 0; i < normal.size(); i++)
          {
            if(ncount[i] != 0)
              nnormal[i] /= ncount[i];
            normal[i] = nnormal[i].normalize();
          }
        }
      }

      class Node
      {
        public:
          int count;
          Tuple2 faces[200];
          Node() : count(0) {}
          void addFace(const Tuple2& t)
          {
            if (count == 200) GRAVIS_THROW2(Exception, "Node in mesh has cardinality greater than 200!");
            faces[count++] = t;
          }
      };

      /*! \brief Generate the adjacency list.
      *
      * The adjacency list (adjacent) contains entries for each triangle.
      * Each entry specifies the adjacent triangle for each edge.
      *
      * The complexity of the algorithm is linear in the number of faces.
      *
      * \throw gravis::Exception if any vertex has a cardinality greater 20
      */
      void generateAdjacencyList()
      {

        const int numFaces = tvi.size();
        const int numVert = vertex.size();

        adjacent.setSize(numFaces);

        std::vector<Node> nodeFaces(numVert);
        for (int i = 0; i < numFaces; i++)
        {
          for (int j = 0; j < 3; j++)
          {
            nodeFaces[tvi[i][j]].addFace(Tuple2(i, j));
          }
        }

        // foreach face
        for (int f = 0; f < numFaces; f++)
        {
          Tuple3& ft = tvi[f];
          Tuple3& at = adjacent[f];
          // foreach edge
          for (int e = 0; e < 3; e++)
          {
            // already found adjacent face for this edge?
            if (at[e] >= 0) continue;
            // vertices for this edge
            int v1 = ft[e];
            int v2 = ft[(e+1)%3];
            // faces using these vertices
            Node& node1 = nodeFaces[v1];
            Node& node2 = nodeFaces[v2];
            for (int i = 0; i < node1.count; i++)
            {
              int f1 = node1.faces[i][0];
              if (f1 == f) continue; // self
              for (int j = 0; j < node2.count; j++)
              {
                if (f1 == node2.faces[j][0])
                {
                  adjacent[f][e] = f1;
                  adjacent[f1][node2.faces[j][1]] = f;
                }
              }
            }
          }
        }
      }
  };

} // namespace gravis

#endif
