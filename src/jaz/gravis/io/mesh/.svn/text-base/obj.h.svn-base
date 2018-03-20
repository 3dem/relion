/******************************************************************************
 **        Title: gravis/io/mesh/obj.h
 **  Description: Implements reader/writer for the .obj mesh file format
 **
 **       Author: Brian Amberg
 **               Computer Science Department, University Basel (CH. **
 ******************************************************************************/
#ifndef __GRAVIS_IO__MESH_OBJ_H__
#define __GRAVIS_IO__MESH_OBJ_H__
#include "OBJReader.h"
#include "OBJWriter.h"

#include <cstdio>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include "../../tImage.h"
#include "../../Exception.h"
#include "absolute_paths.h"

#if defined(HAVE_LIBZ)
#	include <zlib.h>
#endif

namespace gravis
{
  namespace io
  {

    /**
     *
     * saves and load meshes in/from obj format
     * The obj format is a text based represenation for triangle meshes, that
     * can hold a subset of the information in a Mesh object.
     *
     * It is useful to interface different applications, as it is a common format.
     *
     * For internal use the msh format (gravis::io::MeshMSH) might be a better
     * choice, as it is faster, more accurate and covers all of the
     * capabilities of the Mesh datastructure.
     *
     * When implementing loading and saving it is best to use the gravis::io::Mesh
     * and gravis::io::Mesh, which determine the filetype automatically.
     *
     **/
    class MeshOBJ
    {

        static
        inline
        void write(ofstream& ofs, const gravis::fMesh& mesh, string mtlfn)
        {
          ofs << "mtllib " << mtlfn << endl;
          // Write vertices
          for (size_t i = 0; i < mesh.vertex.size(); ++i )
          {
            ofs << "v "
                << mesh.vertex[i].x << " "
                << mesh.vertex[i].y << " "
                << mesh.vertex[i].z << endl;
          }
          ofs << "# " << mesh.vertex.size() << " vertices." << endl << endl;
          // Write texture coordinates
          for (size_t i = 0; i < mesh.texcrd.size(); ++i )
          {
            ofs << "vt "
                << mesh.texcrd[i].x << " "
                << mesh.texcrd[i].y << " "
                << mesh.texcrd[i].z << endl;
          }
          ofs << "# " << mesh.texcrd.size() << " texture coordinates." << endl << endl;

          // Write texture normals
          for (size_t i = 0; i < mesh.normal.size(); ++i )
          {
            ofs << "vn "
                << mesh.normal[i].x << " "
                << mesh.normal[i].y << " "
                << mesh.normal[i].z << " " << endl;
          }
          ofs << "# " << mesh.normal.size() << " normals." << endl << endl;

          // Write faces
          // NOTE: the vertex indices in OBJ format start from 1
          int current_mtl = -1;
          for (size_t i = 0; i < mesh.tvi.size(); ++i )
          {
            if(i<mesh.tmi.size())
            {
              int mtl = mesh.tmi[i];
              if (current_mtl != mtl)
              {
                ofs << "usemtl " << mesh.material[mtl].name << endl;
                current_mtl = mtl;
              }
            }

            ofs << "f ";
            for (size_t c=0 ; c < 3; ++c )
            {
              ofs << (mesh.tvi[i][c])+1;
              if(i<mesh.tti.size())
              {
                ofs 	<< "/" << (mesh.tti[i][c])+1;
              }
              if(i<mesh.tni.size())
              {
                if (i>=mesh.tti.size())
                  ofs << "/";
                ofs 	<< "/" << (mesh.tni[i][c])+1 ;
              }
              ofs 	<< " ";
            }
            ofs << endl;
          }
          ofs << "# " << mesh.tvi.size() << " faces." << endl << endl;
        }

        static
        inline
        void writemtl( ofstream& ofs, const fMesh& mesh, const boost::filesystem::path objdir  )
        {
          for (size_t i=0; i<mesh.material.size(); ++i )
          {
            const fMaterial& material = (mesh.material[i]);
            ofs << "newmtl " << material.name << endl
                << "ns " << material.shininess << endl
                << "ka " << material.ambient.r << " " << material.ambient.g << " " << material.ambient.b << " " << endl
                << "kd " << material.diffuse.r << " " << material.diffuse.g << " " << material.diffuse.b << " " << endl
                << "ks " << material.specular.r << " " << material.specular.g << " " << material.specular.b << " " << endl;
            if ( material.texture.isSet())
            {
              ofs << "map_kd " << mesh_helper::makeRelative(material.texture.getFilenameNative(), objdir).string() << endl;
            }
            if ( material.envMap.isSet())
            {
              ofs << "map_refl " << mesh_helper::makeRelative(material.envMap.getFilenameNative(), objdir).string() << endl;
            }
            if ( material.normalMap.isSet())
            {
              ofs << "map_norm " << mesh_helper::makeRelative(material.normalMap.getFilenameNative(), objdir).string() << endl;
            }
          }
        }


        static
        inline
        void write(ofstream& ofs, const gravis::Mesh& mesh, string mtlfn)
        {
          ofs << "mtllib " << mtlfn << endl;
          // Write vertices
          for (size_t i = 0; i < mesh.vertex.size(); ++i )
          {
            ofs << "v "
                << mesh.vertex[i].x << " "
                << mesh.vertex[i].y << " "
                << mesh.vertex[i].z << endl;
          }
          ofs << "# " << mesh.vertex.size() << " vertices." << endl << endl;
          // Write texture coordinates
          for (size_t i = 0; i < mesh.texcrd.size(); ++i )
          {
            ofs << "vt "
                << mesh.texcrd[i].x << " "
                << mesh.texcrd[i].y << " "
                << mesh.texcrd[i].z << endl;
          }
          ofs << "# " << mesh.texcrd.size() << " texture coordinates." << endl << endl;

          // Write texture normals
          for (size_t i = 0; i < mesh.normal.size(); ++i )
          {
            ofs << "vn "
                << mesh.normal[i].x << " "
                << mesh.normal[i].y << " "
                << mesh.normal[i].z << " " << endl;
          }
          ofs << "# " << mesh.normal.size() << " normals." << endl << endl;

          // Write faces
          // NOTE: the vertex indices in OBJ format start from 1
          int current_mtl = -1;
          for (size_t i = 0; i < mesh.tvi.size(); ++i )
          {
            if(i<mesh.tmi.size())
            {
              int mtl = mesh.tmi[i];
              if (current_mtl != mtl)
              {
                ofs << "usemtl " << mesh.material[mtl].name << endl;
                current_mtl = mtl;
              }
            }

            ofs << "f ";
            for (size_t c=0 ; c < 3; ++c )
            {
              ofs << (mesh.tvi[i][c])+1;
              if(i<mesh.tti.size())
              {
                ofs 	<< "/" << (mesh.tti[i][c])+1;
              }
              if(i<mesh.tni.size())
              {
                if (i>=mesh.tti.size())
                  ofs << "/";
                ofs 	<< "/" << (mesh.tni[i][c])+1 ;
              }
              ofs 	<< " ";
            }
            ofs << endl;
          }
          ofs << "# " << mesh.tvi.size() << " faces." << endl << endl;
        }

        static
        inline
        void writemtl( ofstream& ofs, const Mesh& mesh, const boost::filesystem::path objdir )
        {
          for (size_t i=0; i<mesh.material.size(); ++i )
          {
            const Material& material = (mesh.material[i]);
            ofs << "newmtl " << material.name << endl
                << "ns " << material.shininess << endl
                << "ka " << material.ambient.r << " " << material.ambient.g << " " << material.ambient.b << " " << endl
                << "kd " << material.diffuse.r << " " << material.diffuse.g << " " << material.diffuse.b << " " << endl
                << "ks " << material.specular.r << " " << material.specular.g << " " << material.specular.b << " " << endl;
            if ( material.hasTexture)
            {
              ofs << "map_kd " << mesh_helper::makeRelative(material.textureName, objdir).string() << endl;
            }
            if ( material.hasEnvMap)
            {
              ofs << "map_refl " << mesh_helper::makeRelative(material.envMapName, objdir).string() << endl;
            }
            if ( material.hasNormalMap)
            {
              ofs << "map_norm " << mesh_helper::makeRelative(material.normalMapName, objdir).string() << endl;
            }
          }
        }
        /**
         * Assuming all paths are given in the native format
         **/
        static
        inline
        void write(const string& filename, const gravis::fMesh& mesh)
        {
          const boost::filesystem::path obj_path    (filename);
          const boost::filesystem::path obj_dir     (obj_path.branch_path());
          const std::string             obj_basename( boost::filesystem::basename(obj_path) );
          const boost::filesystem::path mtl_path    (obj_dir / (obj_basename + ".mtl"));

          ofstream obj_fs( obj_path.string().c_str(), ios::out );
          ofstream mtl_fs( mtl_path.string().c_str(), ios::out );
          if (!obj_fs.good()) GRAVIS_THROW3(Exception,"Can not open file for writing!",std::string("File: ")+obj_path.string());
          if (!mtl_fs.good()) GRAVIS_THROW3(Exception,"Can not open file for writing!",std::string("File: ")+mtl_path.string());

          writemtl(mtl_fs, mesh, obj_dir);
          write(obj_fs, mesh, mtl_path.string().c_str());
        }

        /**
         * Assuming all paths are given in the native format
         **/
        static
        inline
        void write(const string& filename, const gravis::Mesh& mesh)
        {
          const boost::filesystem::path obj_path    (filename);
          const boost::filesystem::path obj_dir     (obj_path.branch_path());
          const std::string             obj_basename( boost::filesystem::basename(obj_path) );
          const boost::filesystem::path mtl_path    (obj_dir / (obj_basename + ".mtl"));

          ofstream obj_fs( obj_path.string().c_str(), ios::out );
          ofstream mtl_fs( mtl_path.string().c_str(), ios::out );
          if (!obj_fs.good()) GRAVIS_THROW3(Exception,"Can not open file for writing!",std::string("File: ")+obj_path.string());
          if (!mtl_fs.good()) GRAVIS_THROW3(Exception,"Can not open file for writing!",std::string("File: ")+mtl_path.string());

          writemtl(mtl_fs, mesh, obj_dir);
          write(obj_fs, mesh, mtl_path.string().c_str());
        }

      public:
        static
        inline
        void load(gravis::fMesh& out, const std::string& filename)
        {
          gravis::Mesh m;
          load(m, filename);
          out = m;
          mesh_helper::makePathsAbsolute(out, filename);
        }

        static
        inline
        void save(const std::string& filename, const gravis::fMesh& mesh)
        {
          write(filename, mesh);
        }

        static
        inline
        void save_complete(const std::string& filename, const gravis::fMesh& mesh)
        {
          gravis::fMesh outmesh(mesh);
          boost::filesystem::path msh(filename);
          boost::filesystem::path dir = msh.branch_path();

          for (size_t i=0; i<outmesh.material.size(); ++i)
            if (outmesh.material[i].texture.isSet())
            {
              try
              {
                outmesh.material[i].texture.load();
              }
              catch(const gravis::Exception& e)
              {
                std::cerr << __FILE__ << ":" << __LINE__ << " I was not able to load the texture: " << mesh.material[i].texture.getFilenameNative() << std::endl;
                std::cerr << e << std::endl;
              }
              catch(...)
              {
                std::cerr << __FILE__ << ":" << __LINE__ << " I was not able to load the texture: " << mesh.material[i].texture.getFilenameNative() << std::endl;
              }
              outmesh.material[i].texture.setFilename(boost::filesystem::absolute(dir)/(mesh.material[i].texture.getFilename().leaf()));
              try
              {
                outmesh.material[i].texture.getImage().writePNG(outmesh.material[i].texture.getFilenameNative());
              }
              catch(const gravis::Exception& e)
              {
                std::cerr << __FILE__ << ":" << __LINE__ << " I was not able to save the texture: " << outmesh.material[i].texture.getFilenameNative() << std::endl;
                std::cerr << e << std::endl;
              }
              catch(...)
              {
                std::cerr << __FILE__ << ":" << __LINE__ << " I was not able to save the texture: " << outmesh.material[i].texture.getFilenameNative() << std::endl;
              }
            }
          save(filename, outmesh);
        }

        static
        inline
        void load(gravis::Mesh& out, const std::string& filename)
        {
          OBJReader::load(out, filename);
          mesh_helper::makePathsAbsolute(out, filename);
        }

        static
        inline
        void save(const std::string& filename, const gravis::Mesh& mesh)
        {
          write(filename, mesh);
        }

        static
        inline
        void save_complete(const std::string& filename, const gravis::Mesh& mesh)
        {
          gravis::fMesh outmesh(mesh);
          save_complete(filename, outmesh);
        }

    };
  }
}
#endif
