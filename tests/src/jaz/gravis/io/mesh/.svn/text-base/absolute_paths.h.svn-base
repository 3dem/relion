#ifndef __GRAVIS__IO__MESH__ABSOLUTE_PATHS__
#define __GRAVIS__IO__MESH__ABSOLUTE_PATHS__

#define BOOST_FILESYSTEM_VERSION 3
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/operations.hpp>
#include "../../tMesh.h"
#include "../../Mesh.h"

namespace gravis
{
  namespace io
  {
    namespace mesh_helper
    {

      /**
       * Find the relative path to p given base.
       **/
      static
      inline
      boost::filesystem::path makeRelative(const boost::filesystem::path& p, const boost::filesystem::path& base)
      {
        using namespace boost::filesystem;
        if (p.string() == "")
          return path();
        path p2(absolute(base));
        if (base.string() == "")
          p2 = initial_path();
        path p1(absolute(p, p2));
        // Remove initial equal parts
        path::iterator ip1 = p1.begin();
        path::iterator ip2 = p2.begin();
        while((ip1!=p1.end()) && (ip2 != p2.end()) && (*ip1 == *ip2))
        {
          ++ip1;
          ++ip2;
        };
        path relative;
        for (; ip2!=p2.end(); ++ip2) relative /= "..";
        for (; ip1!=p1.end(); ++ip1) relative /= *ip1;
        return relative;
      }


      /**
       * After loading, all paths have to be absolute, otherwise we loose track
       * of the files. To do this we check if we can find the file by any means
       * (from obj directory or from current directory) and expand accordingly
       **/
      static inline
      void makePathsAbsolute(gravis::fMesh& m, const std::string& filename)
      {
        using namespace boost::filesystem;

        path dir(path(filename).branch_path());

        for (size_t i=0; i<m.material.size(); ++i)
        {
          if (m.material[i].texture.isSet())
          {
            {
              path texturefn = m.material[i].texture.getFilename();
              if (dir=="")
              {
                dir = initial_path<path>();
              }
              if (!dir.is_complete())
              {
                dir = absolute(dir);
              }
              path complete_path = absolute( texturefn, dir );
              if (exists( complete_path ))
                m.material[i].texture.setFilename( complete_path );
            }
          }
        }
      }

      /**
       * After loading, all paths have to be absolute, otherwise we loose track
       * of the files. To do this we check if we can find the file by any means
       * (from obj directory or from current directory) and expand accordingly
       **/
      static inline
      void makePathsAbsolute(gravis::Mesh& m, const std::string& filename)
      {
        using namespace boost::filesystem;
        path dir(path(filename).branch_path());

        for (size_t i=0; i<m.material.size(); ++i)
        {
          if (m.material[i].hasTexture)
          {
            {
              path texturefn = m.material[i].textureName;
              if (dir=="")
              {
                dir = initial_path<path>();
              }
              if (!dir.is_complete())
              {
                dir = absolute(dir);
              }
              path complete_path = absolute( texturefn, dir );
              if (exists( complete_path ))
              {
                m.material[i].textureName = complete_path.string();
              }
            }
          }
        }
      }

    }
  }
}

#endif
