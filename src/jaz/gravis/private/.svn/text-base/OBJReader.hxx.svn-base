/******************************************************************************
 **        Title: OBJReader.hxx
 **  Description: Class to import a .obj file info a Mesh.
 **
 **       Author: Jean-Sebastien Pierrard, 2005
 **               Computer Science Department, University Basel (CH)
 **
 ******************************************************************************/
#include <stdio.h>
#include <ctype.h>
#include "../Exception.h"

#define HAVE_LIBZ
#if defined(HAVE_LIBZ)
#	include <zlib.h>
#endif

#ifdef __APPLE__
#include <sstream>
inline double OBJReader::my_atof(const char* str)
{
  std::istringstream s(str);
  double d = 0;
  s >> d;
  return d;
}

#else
inline double OBJReader::my_atof(const char* str)
{
  return atof(str);
}
#endif


inline
OBJReader::OBJReader () :
  foundNormals(false),
  foundTexCrds(false)
{
  /*
  ** Create 'fallback' material */
  Material fb_mtl("_fallback_");
  mtl_v.push_back(fb_mtl);

  /*
  ** Create 'default' object group */
  Group fb_grp("_default_");
  group_v.push_back(fb_grp);

  active_smggroup = 0; // 0 is the default ('no smoothing') group
  active_mtlgroup = 0; // Point to fallback
}


inline
void OBJReader::read (std::string filename)
{

  std::string::size_type n = filename.rfind('/');
  if (n == std::string::npos)
  {
    objpath = "./";
  }
  else
  {
    objpath  = filename.substr(0, n);
    objpath += "/";
  }

  parseFile(filename);
}


inline
void OBJReader::toLowerCase (char* str)
{
  while (*str != '\0')
  {
    if ((*str >= 'A') && (*str <= 'Z'))
    {
      *str = (char)((int)(*str) - (int)'A' + (int)'a');
    }
    ++str;
  }
}


inline
void OBJReader::parseFile (std::string filename)
{


#ifdef HAVE_LIBZ
  gzFile fin = gzopen (filename.c_str(), "rb");
  if (0 == fin)
  {
    GRAVIS_THROW3(Exception, "Unable to open file: ", filename);
    return;
  }
#else
  FILE* fin = fopen(filename.c_str(), "rt");
  if (0 == fin)
  {
    GRAVIS_THROW3(Exception, "Unable to open file: ", filename);
    return;
  }
#endif

  int  linecount = 0;
  char line[OBJ_MAXLINELEN];

  do
  {
#ifdef HAVE_LIBZ
    if (Z_NULL == gzgets(fin, &line[0], OBJ_MAXLINELEN)) break;
#else
    if (0 == fgets(&line[0], OBJ_MAXLINELEN, fin)) break;
#endif

    int linelen = strlen(line);
    ++linecount;

    // Ignore comments and empty lines
    if ((line[0] == '#' ) || (line[0] == '\n')) continue;

    // Check for VERY long input lines,
    // TODO: this should be handled differently
    if (linelen == (OBJ_MAXLINELEN-1))
    {
      std::cerr << "possible buffer overflow" << std::endl;
      continue;
    }

    // Tokenize line into argc,argv style
    int i=0;
    std::vector<char*> argv;

    while ((linelen > 0) && (i < linelen))
    {
      while (isspace(line[i])) ++i;    // skip leading spaces

      if (i < linelen)
      {
        argv.push_back(&line[i]);
        while (!isspace(line[i])) ++i; // read over sequence of non-spaces
        line[i++] = '\0';              // terminate each sequ. of non-spaces
      }
    }

    // Check on parse errors
    try
    {
      parseLine(argv);
    }
    catch (Exception& e)
    {
      std::cerr << "Parse error in '" << filename << "', line " << linecount << ": "
                << e.detail() << e.argument() << std::endl;
    }
  }
  while (true);

#ifdef HAVE_LIBZ
  gzclose(fin);
#else
  fclose(fin);
#endif
}


inline
void OBJReader::parseLine (std::vector<char*>& argv)
{

  int argc = argv.size();
  if (argc <= 0) return; // return on empty lines

  // Transform argv[0] (datatype) to lower-case and
  // derive integer key from result (max. first 4 letters)

  char* argv0=argv[0];
  int   tkey=0, tpos=0;

  toLowerCase(argv0);

  while (*argv0 != '\0')
  {
    if (tpos < 4)
    {
      tkey |= ((int)*argv0) << (24-8*tpos);
      ++tpos;
    }
    ++argv0;
  }


  switch (tkey)
  {

    case 0x76000000 :   // tkey = "v", vertex coordinates
    {
      if (argc < 4) GRAVIS_THROW3(Exception, errUnexpectedArgs(), argv[0]);
      float x = my_atof(argv[1]);
      float y = my_atof(argv[2]);
      float z = my_atof(argv[3]);
      vertex_v.push_back(f3Vector(x, y, z));

      if (argc == 7 || argc == 8)
      {
        float r = my_atof(argv[4])/255.0;
        float g = my_atof(argv[5])/255.0;
        float b = my_atof(argv[6])/255.0;
        float a = 0;
        if(argc == 8) a = my_atof(argv[7]);
        color_v.push_back(fRGBA(r, g, b,a));
      }
      break;
    }

    case 0x766E0000 :   // tkey = "vn", vertex normal
    {
      if (argc < 4) GRAVIS_THROW3(Exception, errUnexpectedArgs(), argv[0]);
      float x = my_atof(argv[1]);
      float y = my_atof(argv[2]);
      float z = my_atof(argv[3]);
      normal_v.push_back(f3Vector(x, y, z));
      break;
    }

    case 0x76740000 :   // tkey = "vt", texture coordinates
    {
      if (argc < 3) GRAVIS_THROW3(Exception, errUnexpectedArgs(), argv[0]);
      Texcoord uvw;
      if(argc < 4)
        uvw = Texcoord(my_atof(argv[1]), my_atof(argv[2]), 0.0);
      else
        uvw = Texcoord(atof(argv[1]), my_atof(argv[2]), my_atof(argv[3]));
      texcrd_v.push_back(uvw);
      break;
    }

    case 0x66000000 :   // tkey = "(f)ace", polygon
    {
      Face face;

      face.smggroup = active_smggroup;
      face.mtlgroup = active_mtlgroup;

      // Number of vertices for this face(polygon) = argc - 1
      // Parse each section separately for vertex/texture/normal
      for (int c=0; c<argc-1; ++c)
      {

        char* cs=argv[c+1];      // pointer to continous string for the c'th vertex
        std::vector<char*> sls;  // pointers to the positions of the slashes in *cs

        Vertex vertex;
        vertex.vidx = atoi(cs);

        // Parse each section for the field separator "/"
        while (*cs != '\0')
        {
          if (*cs == '/')
          {
            sls.push_back(cs+1);
            *cs = '\0';
          }
          ++cs;
        }

        switch (sls.size())
        {
          case 0 :   // no slashes: only vertex index defined
          {
            vertex.tidx = 0;
            vertex.nidx = 0;
            break;
          }

          case 1 :   // one slash: vertex and texcrd index
          {
            foundTexCrds = true;
            vertex.tidx = atoi(sls[0]);
            vertex.nidx = 0;
            break;
          }

          case 2 :   // two slashes: vertex, texcrd and normal index
          {
            foundNormals = true;
            if (sls[0] == (sls[1]-1))   // texcrd index ommited ??
            {

              vertex.tidx = 0;
              vertex.nidx = atoi(sls[1]);
            }
            else
            {
              foundTexCrds = true;
              vertex.tidx = atoi(sls[0]);
              vertex.nidx = atoi(sls[1]);
            }
            break;
          }

          default :
          {
            GRAVIS_THROW3(Exception, "Unsupported face-format", argv[c+1]);
          }
        } // switch(number of slashes)


        // negative references must be remapped relative to current position
        if (vertex.vidx < 0) vertex.vidx += vertex_v.size() + 1;
        if (vertex.nidx < 0) vertex.nidx += normal_v.size() + 1;
        if (vertex.tidx < 0) vertex.tidx += texcrd_v.size() + 1;

        // OBJs first index is 1, C indexing start with 0
        // non-set indices become -1, which is good!
        vertex.vidx -= 1;
        vertex.nidx -= 1;
        vertex.tidx -= 1;

        // Error(range) checking
        if ((vertex.vidx < 0) || (vertex.vidx >= int(vertex_v.size())))
          GRAVIS_THROW3(Exception, "Vertex index out of range", argv[c+1]);

        if ((vertex.tidx <-1) || (vertex.tidx >= int(texcrd_v.size())))
          GRAVIS_THROW3(Exception, "Texture index out of range", sls[0]);

        if ((vertex.nidx <-1) || (vertex.nidx >= int(normal_v.size())))
          vertex.nidx = -1;

        face.corner.push_back(vertex);
      }  //for(each corner)

      face_v.push_back(face);

      // Add face to each active group
      for (int gid=0; gid<int(active_objgroup.size()); ++gid)
      {
        int group = active_objgroup[gid];
        group_v[group].fidx_v.push_back(face_v.size()-1);
      }

      break;
    } //case f(ace)

    /*
    ** Grouping related directives --------------------------------------------
    */

    case 0x6F000000 :   // tkey = "o", object name
    {
      // Ignored
      break;
    }

    case 0x73000000 :   // tkey = "s", smoothing group
    {
      if (argc != 2) GRAVIS_THROW3(Exception, errUnexpectedArgs(), argv[0]);

      toLowerCase(argv[1]);
      if (!strcmp(argv[1], "off"))
        active_smggroup = 0;
      else
        active_smggroup = atoi(argv[1]);
      break;
    }

    case 0x67000000 :   //tkey = "g", object grouping
    {
      // Clear current list of active groups
      active_objgroup.clear();

      // For each string: test if the groupname is known. If yes, just reuse it,
      // otherwise allocate a new group structure. Add each group to the list of
      // active groups. The faces defined later will be referenced by each active
      // group.

      // If groupname is ommited, do nothing.

      for (int c=1; c<argc; ++c)
      {
        bool        isnew = true;
        std::string grpname(argv[c]);

        for (int gid=0; gid < (int)group_v.size(); ++gid)
        {
          if (group_v[gid].name == grpname)
          {
            isnew = false;
            active_objgroup.push_back(gid);
            break;
          }
        }

        if (isnew)
        {
          Group ng(grpname);
          group_v.push_back(ng);
          active_objgroup.push_back(group_v.size()-1);
        }
      }

      break;
    }

    /*
    ** Material related directives --------------------------------------------
    */

    case 0x6D746C6C :   // tkey = "(mtll)ib", material library
    {
      if (argc < 2) GRAVIS_THROW3(Exception, errUnexpectedArgs(), argv[0]);

      // If some of the material libraries in the argument cannot be loaded
      // we still want to continue looking for the remaining ones. Therefore
      // we collect error informations until each argument has be processed.
      // Then, if at least one file failed to load we throw an exception. In
      // this way we load as much information as possible

      bool        failed       = false;
      std::string failed_files = "";

      for (int c=1; c<argc; ++c)
      {

        std::string mtlpath(argv[c]);
        if ( mtlpath.find('/') == std::string::npos &&
             mtlpath.find('\\') == std::string::npos )
        {
          mtlpath = objpath;
          mtlpath.append(argv[c]);
        }

        try
        {
          parseFile(mtlpath);
          continue;
        }
        catch (Exception& e)
        {
          failed = true;
          failed_files.append(argv[c]).append(" ");
        }

      }

      if (failed)
        GRAVIS_THROW3(Exception,
                      "Unable to load following material libraries: ",
                      failed_files
                     );

      break;
    } // case mtllib


    case 0x6E65776D :   // tkey = "(newm)tl", define new material
    {
      if (argc != 2) GRAVIS_THROW3(Exception, errUnexpectedArgs(), argv[0]);

      bool        isnew = true;
      std::string mtlname(argv[1]);

      // Search name in listed materials
      for (int i=0; i < (int)mtl_v.size(); ++i)
      {
        if (mtl_v[i].name == mtlname)
        {
          isnew = false;
          active_mtlgroup = i;
          break;
        }
      }

      if (isnew)
      {
        Material m(mtlname);
        mtl_v.push_back(m);
        active_mtlgroup = mtl_v.size()-1;
      }

      break;
    }

    case 0x7573656D :   // tkey = "(usem)tl, define current material
    {
      if (argc != 2) GRAVIS_THROW3(Exception, errUnexpectedArgs(), argv[0]);

      bool        isnew = true;
      std::string mtlname(argv[1]);

      // Search name in listed materials
      for (int i=0; i < (int)mtl_v.size(); ++i)
      {
        if (mtl_v[i].name == mtlname)
        {
          isnew = false;
          active_mtlgroup = i;
          break;
        }
      }

      if (isnew)             // If name is unknown,
      {
        active_mtlgroup = 0; // revert to fallback material
        GRAVIS_THROW3(Exception, "Material not found, using 'fallback': ", mtlname);
      }

      break;
    }

    case 0x64000000 :   // tkey = "d", material: dissolve(transparency)
    {
      if (argc != 2) GRAVIS_THROW3(Exception, errUnexpectedArgs(), argv[0]);
      // Ignored

      //  			mtl_v[active_mtlgroup].opacity = atof(argv[1]);
      break;
    }

    case 0x6E730000 :   // tkey = "ns", material: specular exponent
    {
      if (argc != 2) GRAVIS_THROW3(Exception, errUnexpectedArgs(), argv[0]);

      mtl_v[active_mtlgroup].shininess = my_atof(argv[1]);
      break;
    }

    case 0x6B610000 :   // tkey = "ka", material: ambient color
    {
      if (argc != 4) GRAVIS_THROW3(Exception, errUnexpectedArgs(), argv[0]);

      mtl_v[active_mtlgroup].ambient = fRGBA(
                                         my_atof(argv[1]), my_atof(argv[2]), my_atof(argv[3]), 1.f
                                       );
      break;
    }

    case 0x6B640000 :   // tkey = "kd", material: diffuse color
    {
      if (argc != 4) GRAVIS_THROW3(Exception, errUnexpectedArgs(), argv[0]);

      mtl_v[active_mtlgroup].diffuse = fRGBA(
                                         my_atof(argv[1]), my_atof(argv[2]), my_atof(argv[3]), 1.f
                                       );
      break;
    }

    case 0x6B730000 :   // tkey = "ks", material: specular color
    {
      if (argc != 4) GRAVIS_THROW3(Exception, errUnexpectedArgs(), argv[0]);

      mtl_v[active_mtlgroup].specular = fRGBA(
                                          my_atof(argv[1]), my_atof(argv[2]), my_atof(argv[3]), 1.f
                                        );
      break;
    }

    case 0x73686172 :   // tkey = "(shar)pness", material: strength of refl.
    {
      if (argc != 2) GRAVIS_THROW3(Exception, errUnexpectedArgs(), argv[0]);
      // Ignored. Need more info.
      //    mtl_v[active_mtlgroup].reflect = atof(argv[1]);
      break;
    }

    case 0x6E690000 :   // tkey = "ni", material: optical density
    {
      // Ignored
      break;
    }

    case 0x74660000 :   // tkey = "tf", material: ??? what is this ???
    {
      // Ignored
      break;
    }

    case 0x696C6C75 :   // tkey = "(illu)m", material: illumination model
    {
      // Ignored; we let the renderer and/or user decide !
      break;
    }

    case 0x6D61705F :   // tkey = "map_", material: texture maps
    {
      if (argc != 2) GRAVIS_THROW3(Exception, errUnexpectedArgs(), argv[0]);

      std::string texpath(argv[1]);

      if (!strcmp(argv[0], "map_Kd") || !strcmp(argv[0], "map_kd"))   //diffuse texture
      {

        mtl_v[active_mtlgroup].textureName = objpath + texpath;
        mtl_v[active_mtlgroup].hasTexture = true;
        break;
      }

      if (!strcmp(argv[0], "map_refl") )   //environment map
      {

        mtl_v[active_mtlgroup].envMapName = objpath + texpath;
        mtl_v[active_mtlgroup].hasEnvMap = true;
        break;
      }

      if (!strcmp(argv[0], "map_norm") )   //normal map
      {

        mtl_v[active_mtlgroup].normalMapName = objpath + texpath;
        mtl_v[active_mtlgroup].hasNormalMap = true;
        break;
      }

      break;
    }

    default :
    {
      GRAVIS_THROW3(Exception, "Unsupported directive" , argv[0]);
    }
  }
}


/*! \brief Builds a mesh from the parsed object file.
 *
 * Only faces are considered (points and lines are dropped).
 *
 * Note, that the size() of tvi and tmi are guaranteed
 * to be equal. If no material was specified or referenced, a fallback
 * matrerial is used.
 *
 * tti.size() and tni.size() may be zero.
 */
inline
void OBJReader::buildMesh (Mesh& mesh) const
{
  /*
  ** Transfer materials ---------------------------------------------------- */
  bool foundColor = color_v.size() > 0 && color_v.size() == vertex_v.size();

  mesh.material.resize(mtl_v.size());
  for (int i=0; i < (int)mtl_v.size(); ++i)
  {
    mesh.material[i] = mtl_v[i];
  }

  /*
  ** Transfer triangle data ------------------------------------------------ */
  // Compute number of triangles after breaking down non-triangular faces.

  int tri_faces = face_v.size();

  for (int i=0; i<(int)face_v.size(); ++i)
  {
    tri_faces += face_v[i].corner.size() - 3;
  }

  mesh.tvi.resize(tri_faces);
  if (foundNormals) mesh.tni.resize(tri_faces);
  if (foundTexCrds) mesh.tti.resize(tri_faces);
  if (foundColor)   mesh.tci.resize(tri_faces);
  mesh.tmi.resize(tri_faces);

  int t = 0;
  for (int i=0; i<(int)face_v.size(); ++i)
  {
    const Face& face = face_v[i];

    for (int c=0; c<=((int)face.corner.size()-3); ++c)
    {
      mesh.tvi[t].c0 = face.corner[0    ].vidx;
      mesh.tvi[t].c1 = face.corner[1 + c].vidx;
      mesh.tvi[t].c2 = face.corner[2 + c].vidx;

      if (foundNormals)
      {
        mesh.tni[t].c0 = face.corner[0    ].nidx;
        mesh.tni[t].c1 = face.corner[1 + c].nidx;
        mesh.tni[t].c2 = face.corner[2 + c].nidx;
      }

      if (foundTexCrds)
      {
        mesh.tti[t].c0 = face.corner[0    ].tidx;
        mesh.tti[t].c1 = face.corner[1 + c].tidx;
        mesh.tti[t].c2 = face.corner[2 + c].tidx;
      }
      if (foundColor)
      {
        mesh.tci[t].c0 = face.corner[0    ].vidx;
        mesh.tci[t].c1 = face.corner[1 + c].vidx;
        mesh.tci[t].c2 = face.corner[2 + c].vidx;
      }
      mesh.tmi[t] = face.mtlgroup;
      ++t;
    }
  }

  /*
  ** Transfer vertex data -------------------------------------------------- */

  mesh.vertex.resize(vertex_v.size());
  mesh.normal.resize(normal_v.size());
  mesh.texcrd.resize(texcrd_v.size());

  memcpy(
    mesh.vertex.data(), &(vertex_v[0]), 3*sizeof(float)*vertex_v.size()
  );
  memcpy(
    mesh.normal.data(), &(normal_v[0]), 3*sizeof(float)*normal_v.size()
  );
  memcpy(
    mesh.texcrd.data(), &(texcrd_v[0]), 3*sizeof(float)*texcrd_v.size()
  );

  if(foundColor)
  {
    mesh.color.resize(color_v.size());
    memcpy(mesh.color.data(), &(color_v[0]), 4*sizeof(float)*color_v.size());
  }
}

