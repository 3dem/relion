#ifndef __LIBGRAVIS_T_MESH_H__
#define __LIBGRAVIS_T_MESH_H__
/******************************************************************************
**        Title: tMesh.h
**  Description: Templated mesh representation using std:;vector.
**
**       Author: Jean-Sebastien Pierrard, 2005
**               Brian Amberg, 2006
**               Computer Science Department, University Basel (CH)
**
******************************************************************************/

#define ATT_PURE __attribute__ ((pure))

#include "tRGBA.h"
#include "t2Vector.h"
#include "t3Vector.h"
#include "tImage.h"
#include "Tuple.h"
#include <vector>
#include <limits>
#include "Mesh.h"

#include <boost/filesystem.hpp>

namespace gravis
{


  template <class T>
  class tMaterial
  {
    public:

      /**
       * Helper class with lazy loading images with an associated filename.
       * It may also contain no image, for a textureless mesh.
       **/
      class ImageFile
      {
        private:
          mutable bool loaded;
          mutable tImage< tRGBA<T> > image;
          boost::filesystem::path filename;

        public:

          ImageFile() : loaded(false), image(), filename() {};
          ImageFile(const std::string& fn) : loaded(false), image(), filename(fn) {};

          inline const bool& isLoaded() const
          {
            return loaded;
          }
          /**
           * Load the image into memory, if it is not yet loaded
           **/
          void load() const
          {
            if (!loaded)
            {
              reload();
            }
          };
          /**
           * Load the image into memory, even if it is already loaded
           **/
          void reload() const
          {
            if (isSet())
            {
              loaded = true;
              image.read( getFilenameNative() );
            }
          };
          /**
           * Load the given image
           **/
          void load(const std::string& filename)
          {
            setFilename(filename);
            reload();
          };
          /**
           * Change the name of the image file
           *
           * The empty filename resets the image.
           *
           * Convenience function, assuming the filename is in native format.
           **/
          void setFilename(const char* filename)
          {
            setFilename(boost::filesystem::path(filename));
          };
          /**
           * Change the name of the image file
           *
           * The empty filename resets the image.
           *
           * Convenience function, assuming the filename is in native format.
           **/
          void setFilename(const std::string& filename)
          {
            setFilename(boost::filesystem::path(filename));
          };
          /**
           * Change the name of the image file
           *
           * The empty filename resets the image.
           **/
          void setFilename(const boost::filesystem::path& filename)
          {
            this->filename = filename;
          };
          /**
           * Return the image filename in native format.
           **/
          const std::string getFilenameNative() const
          {
            return filename.string();
          }
          /**
           * Return the image filename
           **/
          const boost::filesystem::path& getFilename() const
          {
            return filename;
          }
          /**
           * Delete the texture
           **/
          void reset()
          {
            filename = "";
            image.resize(0,0);
          }
          /**
           * Do we represent the NULL image
           **/
          ATT_PURE bool isSet() const
          {
            return filename != "";
          }
          /**
           * Associate a texture from a tImage
           **/
          void set(const std::string& filename, const tImage<tRGBA<T> > &image)
          {
            loaded = true;
            this->filename = filename;
            this->image = image;
          }
          /**
           * Access the image.
           * There seems to be a problem with changing the image.
           **/
          tImage< tRGBA<T> > &getImage()
          {
            if (!loaded) load();
            return image;
          }
          /**
           * Access the image
           **/
          const tImage< tRGBA<T> > &getImage() const
          {
            if (!loaded) load();
            return image;
          }
          /**
           * Set the image
           **/
          void setImage(const tImage<tRGBA<T> > &img)
          {
            image=img;
          };
      };

    public:
      tMaterial(std::string n="")
        : name(n), ambient(T(0.1),T(1.0)), diffuse(T(0.9),T(1.0)), specular(T(0.6),T(1.0)),
          shininess(T(25.0)), texture(), envMap(), normalMap() {}

      tMaterial(const Material& o) :
        name(o.name),
        ambient(o.ambient),
        diffuse(o.diffuse),
        specular(o.specular),
        shininess(o.shininess),
        texture(o.textureName),
        envMap(o.envMapName),
        normalMap(o.normalMapName)
      {
      }

      tMaterial(const tMaterial& o) :
        name(o.name),
        ambient(o.ambient),
        diffuse(o.diffuse),
        specular(o.specular),
        shininess(o.shininess),
        texture(o.texture),
        envMap(o.envMap),
        normalMap(o.normalMap)
      {
      }

      tMaterial& operator=(const tMaterial& o)
      {
        name= o.name;
        ambient=o.ambient;
        diffuse=o.diffuse;
        specular=o.specular;
        shininess=o.shininess;
        texture=o.texture;
        envMap=o.envMap;
        normalMap=o.normalMap;
        return *this;
      }

      std::string name;

      tRGBA<T> ambient;
      tRGBA<T> diffuse;
      tRGBA<T> specular;

      T shininess; /*!< \brief Phong exponent. */

      ImageFile texture;
      ImageFile envMap;
      ImageFile normalMap;
      
      /// convert this brian-material to a gravis material
      Material getGravisMaterial() const
      {
        Material gravisMaterial;
        gravisMaterial.name = name;
        
        gravisMaterial.ambient = ambient;
        gravisMaterial.diffuse = diffuse;
        gravisMaterial.specular = specular;
        
        gravisMaterial.shininess = shininess;
        
        gravisMaterial.hasTexture  = texture.isSet();
        gravisMaterial.textureName = texture.getFilenameNative();
        gravisMaterial.hasEnvMap  = envMap.isSet();
        gravisMaterial.envMapName = envMap.getFilenameNative();
        gravisMaterial.hasNormalMap  = normalMap.isSet();
        gravisMaterial.normalMapName = normalMap.getFilenameNative();

        return gravisMaterial;
      }
      
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
   * of -1 (which is the 'invalid index' pointing to the default entry in the
   * corresponding defaultVectors).
   */
  template <class T>
  class tMesh
  {
    public:

      std::vector< t3Vector<T>  > vertex; /*!< \brief  Vertex array. */
      std::vector< t3Vector<T>  > normal; /*!< \brief  Normal array. */
      std::vector< t3Vector<T>  > texcrd; /*!< \brief  Texture coordinate array. */
      std::vector< tRGBA<T>     > color; /*!< \brief  Color array. */
      std::vector< tMaterial<T> > material; /*!< \brief  Material array. */

      std::vector<Tuple3> tvi; /*!< \brief  Triangle vertex indices. */
      std::vector<Tuple3> tni; /*!< \brief  Triangle normal indices. */
      std::vector<Tuple3> tti; /*!< \brief  Triangle texcrd indices. */
      std::vector<Tuple3> tci; /*!< \brief  Triangle color indices. */
      std::vector<int   > tmi;    /*!< \brief  Triangle material indices. */

      std::vector<Tuple2> lvi; /*!< \brief  Line vertex indices. */
      std::vector<Tuple2> lti; /*!< \brief  Line texcrd indices. */
      std::vector<Tuple2> lci; /*!< \brief  Line colour indices. */

      std::vector<int> pvi; /*!< \brief  Point vertex indices. */
      std::vector<int> pci; /*!< \brief  Point color indices. */

      std::vector<Tuple3> adjacent; /*!< \brief Adjacency list. See generateAdjacencyList(). */

      tMesh() :
        vertex(0,t3Vector<T>(std::numeric_limits<T>::quiet_NaN())), normal(), texcrd(), color(0, tRGBA<T>(1)), material(),
        tvi(0, Tuple3(-1,-1,-1)),
        tni(0, Tuple3(-1,-1,-1)),
        tti(0, Tuple3(-1,-1,-1)),
        tci(0, Tuple3(-1,-1,-1)),
        tmi(0, -1),
        lvi(0, Tuple2(-1,-1)),
        lti(0, Tuple2(-1,-1)),
        lci(0, Tuple2(-1,-1)),
        pvi(0, -1),
        pci(0, -1),
        adjacent(0,Tuple3(-1,-1,-1))
      {
      }

      tMesh(const tMesh& o) :
        vertex(o.vertex), normal(o.normal), texcrd(o.texcrd), color(o.color), material(o.material),
        tvi(o.tvi), tni(o.tni), tti(o.tti), tci(o.tci), tmi(o.tmi),
        lvi(o.lvi), lti(o.lti), lci(o.lci),
        pvi(o.pvi), pci(o.pci),
        adjacent(o.adjacent)
      {
      }

      tMesh(const Mesh& o) :
        vertex(o.vertex), normal(o.normal), texcrd(o.texcrd), color(o.color), material(o.material.size()),
        tvi(o.tvi), tni(o.tni), tti(o.tti), tci(o.tci), tmi(o.tmi),
        lvi(o.lvi), lti(o.lti), lci(o.lci),
        pvi(o.pvi), pci(o.pci),
        adjacent(o.adjacent)
      {
        for (size_t i=0; i<o.material.size(); ++i)
          material[i] = tMaterial<T>(o.material[i]);
      }

      /**
       * Exception safe swap operator
       **/
      void swap(tMesh& o)
      {
        vertex.swap(o.vertex);
        normal.swap(o.normal);
        texcrd.swap(o.texcrd);
        color.swap(o.color);
        material.swap(o.material);
        tvi.swap(o.tvi);
        tni.swap(o.tni);
        tti.swap(o.tti);
        tci.swap(o.tci);
        tmi.swap(o.tmi);
        lvi.swap(o.lvi);
        lti.swap(o.lti);
        lci.swap(o.lci);
        pvi.swap(o.pvi);
        pci.swap(o.pci);
        adjacent.swap(o.adjacent);
      }

      /**
       * Exception safe assignment operator
       **/
      tMesh& operator=(const tMesh& o)
      {
        tMesh tmp(o);
        tmp.swap(*this);
        return *this;
      }

      /// Generate a gravis mesh from this brian mesh
      Mesh getGravisMesh() const
      {
        Mesh gravisMesh;
        
        gravisMesh.vertex = vertex;
        gravisMesh.normal = normal;
        gravisMesh.texcrd = texcrd;
        gravisMesh.color  = color;
        
        gravisMesh.tvi    = tvi;
        gravisMesh.tni    = tni;
        gravisMesh.tti    = tti;
        gravisMesh.tci    = tci;
        gravisMesh.tmi    = tmi;
        
        gravisMesh.lvi    = lvi;
        gravisMesh.lti    = lti;
        gravisMesh.lci    = lci;
        
        gravisMesh.pvi    = pvi;
        gravisMesh.pci    = pci;
        
        gravisMesh.adjacent = adjacent;
        
        gravisMesh.material.resize( material.size() );
        for (size_t i=0; i<material.size(); ++i)
          gravisMesh.material[i] = material[i].getGravisMaterial();
        
        return gravisMesh;
      }
      
      void generateNormals()
      {
        const int numFaces = int(tvi.size());
        tni.resize(numFaces);
        normal.resize(numFaces);

        for (int i = 0; i < numFaces; i++)
        {
          t3Vector<T> a = (vertex[tvi[i][1]] - vertex[tvi[i][0]]);
          t3Vector<T> b = (vertex[tvi[i][2]] - vertex[tvi[i][0]]);
          normal[i] = cross(a, b).normalize();
          tni[i] = Tuple3(i, i, i);
        }
      }

      void generatePerVertexNormals()
      {
        std::vector<int> ncount;
        t3Vector<T> norm;
        const int numFaces = int(tvi.size());
        tni.resize(numFaces);

        normal.resize(vertex.size());
        ncount.resize(vertex.size());
        for (unsigned int i = 0; i < ncount.size(); i++)
        {
          ncount[i] = 0;
          normal[i] = t3Vector<T>(T(0));
        }


        for (int i = 0; i < numFaces; i++)
        {
          t3Vector<T> a = (vertex[tvi[i][1]] - vertex[tvi[i][0]]);
          t3Vector<T> b = (vertex[tvi[i][2]] - vertex[tvi[i][0]]);
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
            normal[i] /= T(ncount[i]);
          normal[i] = normal[i].normalize();
        }
      }

      class Node
      {
        public:
          int count;
          Tuple2 faces[20];
          Node() : count(0) {}
          void addFace(const Tuple2& t)
          {
            if (count == 20) GRAVIS_THROW2(Exception, "Node in mesh has cardinality greater than 20!");
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

        adjacent.resize(numFaces);

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

  /// See tMesh::swap().
  template<typename T>
  inline void
  swap(tMesh<T>& __x, tMesh<T>& __y)
  {
    __x.swap(__y);
  }

  typedef tMaterial<float>  fMaterial;
  typedef tMaterial<double> dMaterial;
  typedef tMesh<float>  fMesh;
  typedef tMesh<double> dMesh;

} // namespace gravis

#endif
