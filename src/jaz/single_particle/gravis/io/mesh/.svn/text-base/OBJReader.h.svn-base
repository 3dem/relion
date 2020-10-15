/******************************************************************************
**        Title: OBJReader.h
**  Description: Class to import a .obj file info a Mesh.
**
**       Author: Jean-Sebastien Pierrard, 2005
**               Computer Science Department, University Basel (CH)
**
******************************************************************************/
#ifndef _OBJREADER_H_
#define _OBJREADER_H_

#include <vector>
#include <string>

#include "../../Mesh.h"

namespace gravis
{

  /*! \brief Class for reading Wavefront OBJ files.
   *
   * Features:
   * - Currently only reads polygonal information (no points, no lines,
   *   and especially no NURBS, bezier patches etc.)
   * - The old maplib and usemap directives are ignored; texture maps must be
   *   specified in materials (with map_Kd)
   * - map_Ks, map_d, map_refl and map_bump are crrently not read
   * - all grouping stuff is ignored (and consequently smoothing groups as
   *   a way to implicitly specify normals)
   * - illum is ignored (it's a rendering setting)
   * - many other things in materials are ignored (Td, Ni, sharpness...)
   */
  class OBJReader
  {
    public:

      struct Texcoord
      {
        Texcoord():u(0.0),v(0.0),w(0.0) {}
        Texcoord(float u,float v,float w):u(u),v(v),w(w) {}
        float u, v, w;
      };

      struct Vertex
      {
        Vertex () : vidx(-1), nidx(-1), tidx(-1) { }
        int vidx, nidx, tidx;
      };

      struct Face
      {
        int smggroup;
        int mtlgroup;

        std::vector<Vertex> corner;
      };

      struct Group
      {
        Group (std::string n) : name(n) { }

        std::string      name;
        std::vector<int> fidx_v;
      };

      OBJReader ();

      void read (std::string);
      void buildMesh (Mesh&) const;

      /*!
       * Convenience function to load a file into a mesh.
       **/
      static void load(Mesh& m, const std::string& fn)
      {
        OBJReader objr;
        objr.read(fn);
        objr.buildMesh(m);
      }

    protected:

      void toLowerCase (char*);

      void parseFile (std::string);
      void parseLine (std::vector<char*>&);

      static const unsigned int OBJ_MAXLINELEN = 512;
      static const unsigned int OBJ_MAXARGVLEN = 32;

      static const char* errUnexpectedArgs ()
      {
        return "Unexpected #arguments for directive: ";
      }

      std::vector<gravis::f3Vector>   vertex_v;
      std::vector<gravis::fRGBA>   color_v;
      std::vector<gravis::f3Vector>   normal_v;
      std::vector<Texcoord>   texcrd_v;
      std::vector<Face>       face_v;
      std::vector<Group>      group_v;
      std::vector<Material> mtl_v;

      // whether indices to normals were found
      bool foundNormals;
      // whether indices to texture coordinates were found
      bool foundTexCrds;

      std::string objpath;

      int              active_smggroup;
      int              active_mtlgroup;
      std::vector<int> active_objgroup;
      double my_atof(const char* str);
  };

#include "../../private/OBJReader.hxx"

} // namespace gravis;

#endif
