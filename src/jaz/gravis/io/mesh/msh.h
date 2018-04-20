/******************************************************************************
**        Title: gravis/io/mesh/msh.h
**  Description: Implements reader/writer for the .msh mesh file format
**
**       Author: Brian Amberg
**               Computer Science Department, University Basel (CH)
**
******************************************************************************/
#ifndef __GRAVIS_IO__MESH_MSH_H__
#define __GRAVIS_IO__MESH_MSH_H__

#include <zlib.h>

#include "../../Mesh.h"
#include "../../tMesh.h"
#include <iostream>
#include <sstream>

#ifdef WIN32
#include <malloc.h>
#define alloca _alloca
#else
#include <alloca.h>
#endif

#include "absolute_paths.h"
#include <boost/filesystem/exception.hpp>

#define MEM_ALLOC(type, n, size) ((type)alloca((n)*(size)))

namespace gravis
{
  namespace io
  {

    namespace MSH_CONSTANTS
    {
      static const std::string id_start("GRAVIS::MSH2::BINARY::START\n");
      static const std::string id_end("\nGRAVIS::MSH2::BINARY::END\n");
    }

    /**
     * Saves meshes in msh format
     * The msh format is a direct reflection of our Mesh datastructure, which
     * is extremely fast to save and load while being a lot more accurate and
     * compact than meshes in obj format.
     *
     * When implementing loading and saving it is best to use the gravis::io::Mesh
     * and gravis::io::Mesh, which determine the filetype automatically.
     **/
    class MeshMSH
    {

      private:

        /*** Functions to load and save data ***/
        template <class T>
        static
        inline
        void get(T& v, gzFile& fin)
        {
          if (gzread(fin, (char*)(&v), sizeof(v)) <= 0)
            GRAVIS_THROW2(gravis::Exception, "File ends prematurely");
        }

        template <class T>
        static
        inline
        void put(std::ostream& os, const T& v)
        {
          os.write((char*)(&v), sizeof(v));
        }

        static
        inline
        size_t get_size(gzFile& fin)
        {
          uint32_t size;
          get(size, fin);
          return size;
        }

        static
        inline
        void put_size(std::ostream& os, size_t size)
        {
          uint32_t sz = size;
          os.write((char*)(&sz), sizeof(sz));
        }

        template <class T>
        static
        inline
        void get(std::vector<T> &out, gzFile& fin)
        {
          size_t size = get_size(fin);
          out.resize(size);
          if (size > 0)
            if (gzread(fin, (char*)&(out[0]), sizeof(out[0]) * size) <= 0)
              GRAVIS_THROW2(gravis::Exception, "File ends prematurely");
        }

        template <class T>
        static
        inline
        void get(gravis::tArray<T> &out, gzFile& fin)
        {
          size_t size = get_size(fin);
          out.resize(size);
          if (size > 0)
            if (gzread(fin, (char*)&(out[0]), sizeof(out[0]) * size) <= 0)
              GRAVIS_THROW2(gravis::Exception, "File ends prematurely");
        }

        template <class T>
        static
        inline
        void put(std::ostream& os, const std::vector<T> &in)
        {
          put_size(os, in.size());
          if ( in.size() > 0 )
            os.write((char*)&(in[0]), sizeof(in[0]) * in.size());
        }

        template <class T>
        static
        inline
        void put(std::ostream& os, const gravis::tArray<T> &in)
        {
          put_size(os, in.size());
          if ( in.size() > 0 ) // invalid access of in[0] if omitted!!
            os.write((char*)&(in[0]), sizeof(in[0]) * in.size());
        }

        static
        inline
        void get(std::string& out, gzFile& fin)
        {
          size_t size = get_size(fin);
          out.resize(size);
          if (size == 0)
            return;
          char* temp = MEM_ALLOC(char*,size+1,1);
          temp[size] = 0;
          if (gzread(fin, temp, size) <= 0)
            GRAVIS_THROW2(gravis::Exception, "File ends prematurely");
          out=temp;
        }

        static
        inline
        void put(std::ostream& os, const std::string& in)
        {
          put_size(os, in.size());
          os.write(in.data(), in.size());
        }

      public:

        /**
         * Load a mesh from the binary msh format
         **/
        static
        inline
        void load(gravis::fMesh& mesh, const std::string& filename)
        {

          gzFile fin = gzopen(filename.c_str(), "rb");
          if (0 == fin)
          {
            GRAVIS_THROW3(Exception, "Unable to open file: ", filename.c_str());
            return;
          }

          try
          {
            char _id_start[MSH_CONSTANTS::id_start.size()+1];
            char _id_end[MSH_CONSTANTS::id_end.size()+1];
            _id_end[MSH_CONSTANTS::id_end.size()] = 0;
            _id_start[MSH_CONSTANTS::id_start.size()] = 0;

            gzread(fin, _id_start, MSH_CONSTANTS::id_start.size());
            if (MSH_CONSTANTS::id_start != std::string(_id_start))
              GRAVIS_THROW3(Exception, "File is not in gravis msh format.", filename.c_str());

            mesh.material.resize(get_size(fin));
            for (unsigned int i=0; i<mesh.material.size(); i++)
            {
              get(mesh.material[i].name, fin);
              get(mesh.material[i].ambient, fin);
              get(mesh.material[i].diffuse, fin);
              get(mesh.material[i].specular, fin);
              get(mesh.material[i].shininess, fin);
              bool hasTexture;
              get(hasTexture, fin);
              std::string textureName;
              get(textureName, fin);
              if (hasTexture)
              {
                try
                {
                  mesh.material[i].texture.setFilename(textureName);
                }
                catch(...)
                {
                  std::cout << "Could not open texture " << textureName << std::endl;
                }
              }
            }
            get(mesh.vertex, fin);
            get(mesh.normal, fin);
            get(mesh.texcrd, fin);
            get(mesh.color, fin);
            get(mesh.tvi, fin);
            get(mesh.tni, fin);
            get(mesh.tti, fin);
            get(mesh.tci, fin);
            get(mesh.tmi, fin);
            get(mesh.lvi, fin);
            get(mesh.lti, fin);
            get(mesh.lci, fin);
            get(mesh.pvi, fin);
            get(mesh.pci, fin);
            gzread(fin, _id_end, MSH_CONSTANTS::id_start.size());
            if (MSH_CONSTANTS::id_end != std::string(_id_end))
              GRAVIS_THROW3(Exception, "File is not in gravis msh format.", filename.c_str());

          }
          catch(gravis::Exception& e)
          {
            gzclose(fin);
            throw(e);
          }
          gzclose(fin);
          mesh_helper::makePathsAbsolute(mesh, filename);
        }

        //! Deprecated, you should use gravis::fMesh instead
        static
        inline
        void load(gravis::Mesh& mesh, const std::string& filename)
        {

          gzFile fin = gzopen(filename.c_str(), "rb");
          if (0 == fin)
          {
            GRAVIS_THROW3(Exception, "Unable to open file: ", filename.c_str());
            return;
          }

          try
          {
            char* _id_start = MEM_ALLOC(char*,MSH_CONSTANTS::id_start.size()+1,1);
            char* _id_end   = MEM_ALLOC(char*,MSH_CONSTANTS::id_end.size()+1,1);
            _id_end[MSH_CONSTANTS::id_end.size()] = 0;
            _id_start[MSH_CONSTANTS::id_start.size()] = 0;

            gzread(fin, _id_start, MSH_CONSTANTS::id_start.size());
            if (MSH_CONSTANTS::id_start != std::string(_id_start))
              GRAVIS_THROW3(Exception, "File is not in gravis msh format.", filename.c_str());

            mesh.material.resize(get_size(fin));
            for (unsigned int i=0; i<mesh.material.size(); i++)
            {
              get(mesh.material[i].name, fin);
              get(mesh.material[i].ambient, fin);
              get(mesh.material[i].diffuse, fin);
              get(mesh.material[i].specular, fin);
              get(mesh.material[i].shininess, fin);
              get(mesh.material[i].hasTexture, fin);
              get(mesh.material[i].textureName, fin);
            }
            get(mesh.vertex, fin);
            get(mesh.normal, fin);
            get(mesh.texcrd, fin);
            get(mesh.color, fin);
            get(mesh.tvi, fin);
            get(mesh.tni, fin);
            get(mesh.tti, fin);
            get(mesh.tci, fin);
            get(mesh.tmi, fin);
            get(mesh.lvi, fin);
            get(mesh.lti, fin);
            get(mesh.lci, fin);
            get(mesh.pvi, fin);
            get(mesh.pci, fin);
            gzread(fin, _id_end, MSH_CONSTANTS::id_start.size());
            if (MSH_CONSTANTS::id_end != std::string(_id_end))
              GRAVIS_THROW3(Exception, "File is not in gravis msh format.", filename.c_str());

          }
          catch(gravis::Exception& e)
          {
            gzclose(fin);
            throw(e);
          }
          gzclose(fin);
          mesh_helper::makePathsAbsolute(mesh, filename);
        }

        /**
         * Save a gravis mesh into a fast and complete binary format
         **/
        static
        inline
        void save(const std::string& filename, const gravis::fMesh& mesh)
        {
          const boost::filesystem::path msh_path    (filename);
          const boost::filesystem::path msh_dir     (msh_path.branch_path());
          std::fstream os(filename.c_str(), std::ios::out | std::ios::binary);
          if (!os.good())
          {
            GRAVIS_THROW2(Exception, "Could not save to file");
          }
          os << MSH_CONSTANTS::id_start;
          put_size(os, mesh.material.size());
          for (unsigned int i=0; i<mesh.material.size(); i++)
          {
            put(os, mesh.material[i].name);
            put(os, mesh.material[i].ambient);
            put(os, mesh.material[i].diffuse);
            put(os, mesh.material[i].specular);
            put(os, mesh.material[i].shininess);
            put(os, mesh.material[i].texture.isSet());
            put(os, mesh_helper::makeRelative(mesh.material[i].texture.getFilename(), msh_dir).string());
          }
          put(os, mesh.vertex);
          put(os, mesh.normal);
          put(os, mesh.texcrd);
          put(os, mesh.color);
          put(os, mesh.tvi);
          put(os, mesh.tni);
          put(os, mesh.tti);
          put(os, mesh.tci);
          put(os, mesh.tmi);
          put(os, mesh.lvi);
          put(os, mesh.lti);
          put(os, mesh.lci);
          put(os, mesh.pvi);
          put(os, mesh.pci);
          os << MSH_CONSTANTS::id_end;
        }

        //! Deprecated, you should use gravis::tMesh objects
        static
        inline
        void save(const std::string& filename, const gravis::Mesh& mesh)
        {
          const boost::filesystem::path msh_path    (filename);
          const boost::filesystem::path msh_dir     (msh_path.branch_path());
          std::fstream os(filename.c_str(), std::ios::out | std::ios::binary);
          if (!os.good())
          {
            GRAVIS_THROW2(Exception, "Could not save to file");
          }
          os << MSH_CONSTANTS::id_start;
          put_size(os, mesh.material.size());
          for (unsigned int i=0; i<mesh.material.size(); i++)
          {
            put(os, mesh.material[i].name);
            put(os, mesh.material[i].ambient);
            put(os, mesh.material[i].diffuse);
            put(os, mesh.material[i].specular);
            put(os, mesh.material[i].shininess);
            put(os, mesh.material[i].hasTexture);
            put(os, mesh_helper::makeRelative(boost::filesystem::path(mesh.material[i].textureName), msh_dir).string());
          }
          put(os, mesh.vertex);
          put(os, mesh.normal);
          put(os, mesh.texcrd);
          put(os, mesh.color);
          put(os, mesh.tvi);
          put(os, mesh.tni);
          put(os, mesh.tti);
          put(os, mesh.tci);
          put(os, mesh.tmi);
          put(os, mesh.lvi);
          put(os, mesh.lti);
          put(os, mesh.lci);
          put(os, mesh.pvi);
          put(os, mesh.pci);
          os << MSH_CONSTANTS::id_end;
        }

        static
        inline
        void save_complete(const std::string& filename, const gravis::Mesh& mesh)
        {
          gravis::fMesh outmesh(mesh);
          save_complete(filename, outmesh);
        }

        static
        inline
        void save_complete(const std::string& filename, const gravis::fMesh& mesh)
        {
          gravis::fMesh outmesh(mesh);
          try
          {
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
                  std::cerr << __FILE__ << ":" << __LINE__ << " I was not able to save the texture: " << mesh.material[i].texture.getFilenameNative() << std::endl;
                  std::cerr << e << std::endl;
                }
                catch(...)
                {
                  std::cerr << __FILE__ << ":" << __LINE__ << " I was not able to save the texture: " << mesh.material[i].texture.getFilenameNative() << std::endl;
                }
              }
          }
          catch(const boost::filesystem::filesystem_error& e)
          {
            GRAVIS_THROW3(Exception, e.what(), filename);
          }
          save(filename, outmesh);
        }


    };
  }
}
#endif

