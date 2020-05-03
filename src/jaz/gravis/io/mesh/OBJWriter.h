/******************************************************************************
**        Title:
**  Description:
******************************************************************************/
#ifndef _OBJWriter_H_
#define _OBJWriter_H_

#include <string>
#include <map>

#include "../../Mesh.h"
#include "../../tMesh.h"
#include <sstream>

namespace gravis
{


  class OBJWriter
  {
      static std::string textureFilename(std::string basename, int i)
      {
        std::stringstream s;
        s << basename << "_" << i << ".png";
        return s.str();
      }

      // TODO: These function belong into libFoundation, but I do not want to have xerces in everything
#ifdef WIN32
#define PATH_SEPARATOR "\\"
#else
#define PATH_SEPARATOR "/"
#endif
      static
      void path_split(std::string& dir, std::string& filename, const std::string& path)
      {
        size_t p = path.rfind(PATH_SEPARATOR);
        if (p == std::string::npos)
        {
          dir = "."PATH_SEPARATOR;
          filename = path;
        }
        else
        {
          dir = path.substr(0, p+1);
          filename = path.substr(p+1);
        }
      }

      static
      std::string path_filename(const std::string& path)
      {
        std::string dir, filename;
        path_split(dir, filename, path);
        return filename;
      };


      static
      void copyfile(const std::string& out, const std::string& in)
      {
        std::ifstream is(in.c_str(), std::ios_base::binary);
        std::ofstream os(out.c_str(), std::ios_base::binary);
        os << is.rdbuf();
      }

      static
      bool file_exists(const std::string& filename)
      {
        bool res = false;
        FILE* f = fopen(filename.c_str(),"r");
        if(f)
        {
          res = true;
          fclose(f);
        }
        return res;
      }

    public:
      OBJWriter();
      void write(std::string filepath, const gravis::Mesh& mesh);
      // It writes the gGenericaMesh data into OBJ format.
      void write( std::ofstream& ofs, const gravis::Mesh* mesh, std::string mtlfn );
      // It writes the material file of an OBJ.
      void writemtl( std::string path, std::string mtlfn, const gravis::Mesh* mesh );

      /*!
       * Convenience function to write a mesh to file without explicitly constructing an OBJWriter
       **/
      static void save( const std::string& fn, const gravis::Mesh& mesh)
      {
        OBJWriter ow;
        ow.write(fn, mesh);
      }

      /*!
       * Convenience function to write a mesh with its textures to file without explicitly constructing an OBJWriter
       **/
      static void save_complete( const std::string& fn, const gravis::Mesh& mesh)
      {
        std::cout << " Writing " << fn << std::endl;
        gravis::Mesh outmesh;
        mesh.clone(outmesh);

        std::string dir, name;
        path_split(dir, name, fn);

        for (size_t i=0; i<mesh.material.size(); ++i)
          if (mesh.material[i].hasTexture)
          {
            if (!(file_exists(mesh.material[i].textureName)))
              GRAVIS_THROW3(Exception,"Did not find texture file",mesh.material[i].textureName);
            outmesh.material[i].textureName = path_filename(mesh.material[i].textureName);
            if (!(file_exists(dir + outmesh.material[i].textureName)))
            {
              copyfile(dir + outmesh.material[i].textureName, mesh.material[i].textureName);
              std::cout << " copied to " << dir + outmesh.material[i].textureName << std::endl;
            }
          }
        OBJWriter ow;
        ow.write(fn, outmesh);
      }
    protected:
      // It writes the gVAMesh data into OBJ format.
      //	void write(std::ofstream& ofs, const gvl::gVAMesh* mesh );
    private: // state variables
  };

#include "../../private/OBJWriter.hxx"

} /* Close Namespace "gvl" */


#endif /* _OBJWriter_H_ */
