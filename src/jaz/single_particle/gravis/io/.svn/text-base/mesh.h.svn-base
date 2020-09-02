/******************************************************************************
**        Title: gravis/io/mesh/obj.h
**  Description: Implements reader/writer for different mesh file formats.
**
**       Author: Brian Amberg
**               Computer Science Department, University Basel (CH)
**
******************************************************************************/
#ifndef __GRAVIS_IO_MESH__
#define __GRAVIS_IO_MESH__

#include "mesh/obj.h"
#include "mesh/msh.h"

namespace gravis
{
  namespace io
  {

    /**
     * The main functionality of this class is to forward the commands to
     * MeshMSH and MeshOBJ based on the file ending.
     *
     * Using this instead of MeshMSH and MeshOBJ makes your program a lot more
     * flexible.
     **/
    class Mesh
    {

      private:
        static
        inline
        bool has_ending(const std::string& filename, const std::string& ending)
        {
          return (filename.size() >= ending.size()) && (filename.substr(filename.size() - ending.size()) == ending);
        }

        static
        inline
        bool is_obj(const std::string& filename)
        {
          return
            (has_ending(filename, ".obj") ||
             has_ending(filename, ".obj.gz"));
        }

        static
        inline
        bool is_msh(const std::string& filename)
        {
          return
            (has_ending(filename, ".msh") ||
             has_ending(filename, ".msh.gz"));
        }

      public:
        /**
         * Load mesh from a file. The filetype is determined automatically
         **/
        static
        inline
        void load(gravis::Mesh& out, const std::string& filename)
        {
          if (is_msh(filename))
            MeshMSH::load(out, filename);
          else
            MeshOBJ::load(out, filename);
        }

        /**
         * Load mesh from a file. The filetype is determined automatically
         **/
        static
        inline
        void load(gravis::fMesh& out, const std::string& filename)
        {
          if (is_msh(filename))
            MeshMSH::load(out, filename);
          else
          {
            MeshOBJ::load(out, filename);
          }
        }

        /**
         * save mesh to a file. The filetype is determined from the ending
         **/
        static
        inline
        void save(const std::string& filename, const gravis::Mesh& mesh)
        {
          if (is_msh(filename))
            MeshMSH::save(filename, mesh);
          else
            MeshOBJ::save(filename, mesh);
        }

        /**
         * save mesh to a file. The filetype is determined from the ending
         **/
        static
        inline
        void save(const std::string& filename, const gravis::fMesh& mesh)
        {
          if (is_msh(filename))
            MeshMSH::save(filename, mesh);
          else
            MeshOBJ::save(filename, mesh);
        }

        /**
         * save mesh to a file. The filetype is determined from the ending.
         * Any texture files are copied into the same directory as the output
         * file unless they already exist.
         **/
        static
        inline
        void save_complete(const std::string& filename, const gravis::Mesh& mesh)
        {
          if (is_msh(filename))
            MeshMSH::save_complete(filename, mesh);
          else
            MeshOBJ::save_complete(filename, mesh);
        }

        /**
         * save mesh to a file. The filetype is determined from the ending.
         * Any texture files are copied into the same directory as the output
         * file unless they already exist.
         **/
        static
        inline
        void save_complete(const std::string& filename, const gravis::fMesh& mesh)
        {
          if (is_msh(filename))
            MeshMSH::save_complete(filename, mesh);
          else
            MeshOBJ::save_complete(filename, mesh);
        }

    };

  }
}

#endif
