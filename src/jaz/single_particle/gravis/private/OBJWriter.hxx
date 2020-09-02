
#include <cstdio>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include "../tImage.h"
#include "../Exception.h"

#if defined(HAVE_LIBZ)
#	include <zlib.h>
#endif

#ifdef WIN32
#define PATH_SEPARATOR "\\"
#else
#define PATH_SEPARATOR "/"
#endif

using namespace std;

inline
OBJWriter::OBJWriter()
{
}

inline
void OBJWriter::write(string filename, const gravis::Mesh& mesh)
{
  // 	File::addSearchPath(objfile.pathName());
  // 	File::addSearchPath(".");
  int pos = filename.rfind(PATH_SEPARATOR);
  std::string base = "";
  string path      = "."PATH_SEPARATOR;
  if(pos<=0) pos=-1;
  if(pos>0)
    path = filename.substr(0,pos+1);

  base = filename.substr(pos+1,filename.size());

#ifdef DEBUG
  cout << "Writing OBJ File: " << filename << endl;
  cout << "--      pathname: " << path << endl;
  cout << "--      basename: " << base << endl;
#endif

  ofstream ofs( filename.c_str(), ios::out );
  if(!ofs) GRAVIS_THROW3(Exception,"Can not open file for writing!",std::string("File: ")+filename);

  base.replace(base.find("obj"), 3, "mtl");

  write(ofs, &mesh, base);
  //writemtl(path, base, &mesh);
}
/*
void OBJWriter::write( ofstream& ofs, const gVAMesh* mesh ) {
	// Write vertices
	for ( size_t i = 0; i < mesh->getNofVertices(); ++i ) {
		ofs << "v "
				<< mesh->vertex_v[i].x << " "
				<< mesh->vertex_v[i].y << " "
				<< mesh->vertex_v[i].z << endl;
	}
	ofs << "# " << mesh->getNofVertices() << " vertices." << endl << endl;
	// Write faces
	// NOTE: the vertex indices in OBJ format start from 1
	for ( size_t i = 0; i < mesh->getNofTriangles(); ++i ) {
		ofs << "f "
				<< mesh->triang_v[i].index[0]+1 << " "
				<< mesh->triang_v[i].index[1]+1 << " "
				<< mesh->triang_v[i].index[2]+1 << endl;
	}
	ofs << "# " << mesh->getNofTriangles() << " triangles." << endl << endl;
}
*/
inline
void OBJWriter::write(ofstream& ofs, const gravis::Mesh* mesh, string mtlfn)
{
  //ofs << "mtllib " << mtlfn << endl;
  // Write vertices
  for (size_t i = 0; i < mesh->vertex.size(); ++i )
  {
    ofs << "v "
        << mesh->vertex[i].x << " "
        << mesh->vertex[i].y << " "
        << mesh->vertex[i].z << endl;
  }
  ofs << "# " << mesh->vertex.size() << " vertices." << endl << endl;
  // Write texture coordinates
  for (size_t i = 0; i < mesh->texcrd.size(); ++i )
  {
    ofs << "vt "
        << mesh->texcrd[i].x << " "
        << mesh->texcrd[i].y << " "
        << mesh->texcrd[i].z << endl;
  }
  ofs << "# " << mesh->texcrd.size() << " texture coordinates." << endl << endl;

  // Write texture normals
  for (size_t i = 0; i < mesh->normal.size(); ++i )
  {
    ofs << "vn "
        << mesh->normal[i].x << " "
        << mesh->normal[i].y << " "
        << mesh->normal[i].z << " " << endl;
  }
  ofs << "# " << mesh->normal.size() << " normals." << endl << endl;

  // Write faces
  // NOTE: the vertex indices in OBJ format start from 1
  int current_mtl = -1;
  for (size_t i = 0; i < mesh->tvi.size(); ++i )
  {
    if(i<mesh->tmi.size())
    {
      int mtl = mesh->tmi[i];
      if (current_mtl != mtl)
      {
        ofs << "usemtl " << mesh->material[mtl].name << endl;
        current_mtl = mtl;
      }
    }

    ofs << "f ";
    for (size_t c=0 ; c < 3; ++c )
    {
      ofs << (mesh->tvi[i][c])+1;
      if(i<mesh->tti.size())
      {
        ofs 	<< "/" << (mesh->tti[i][c])+1;
      }
      if(i<mesh->tni.size())
      {
        if (i>=mesh->tti.size())
          ofs << "/";
        ofs 	<< "/" << (mesh->tni[i][c])+1 ;
      }
      ofs 	<< " ";
    }
    ofs << endl;
  }
  ofs << "# " << mesh->tvi.size() << " faces." << endl << endl;
}

inline
void OBJWriter::writemtl( string path, string mtlfn, const Mesh* mesh )
{
#ifdef DEBUG
  cout <<"Writing material file: "<< string(path).append(mtlfn).c_str() << endl;
#endif
  ofstream ofs( string(path).append(mtlfn).c_str(), ios::out );
  //int current_mtl = -1;
  vector<int> written;
  for (size_t i = 0; i < mesh->tmi.size(); ++i )
  {
    int mtl = mesh->tmi[i];
    //if (current_mtl != mtl &&  current_mtl != 0) {
    if (find(written.begin(),written.end(),mtl) == written.end() )
    {
      const Material* material = &(mesh->material[mtl]);
      ofs << "newmtl " << material->name << endl
          //						<< "d " << material->opacity << endl
          << "ns " << material->shininess << endl
          << "ka " << material->ambient.r << " " << material->ambient.g << " " << material->ambient.b << " " << endl
          << "kd " << material->diffuse.r << " " << material->diffuse.g << " " << material->diffuse.b << " " << endl
          << "ks " << material->specular.r << " " << material->specular.g << " " << material->specular.b << " " << endl;
      if ( material->hasTexture)
      {
        string dmap_file = material->textureName;
        int pos = dmap_file.rfind(PATH_SEPARATOR);
        if (pos <= 0) pos = -1;
        std::string base = dmap_file.substr(pos+1, dmap_file.size());
        ofs << "map_kd " << base << endl;
      }
      if ( material->hasEnvMap)
      {
        string dmap_file = material->envMapName;
        int pos = dmap_file.rfind(PATH_SEPARATOR);
        if (pos <= 0) pos = -1;
        std::string base = dmap_file.substr(pos+1, dmap_file.size());
        ofs << "map_refl " << base << endl;
      }
      if ( material->hasNormalMap)
      {
        string dmap_file = material->normalMapName;
        int pos = dmap_file.rfind(PATH_SEPARATOR);
        if (pos <= 0) pos = -1;
        std::string base = dmap_file.substr(pos+1, dmap_file.size());
        ofs << "map_norm " << base << endl;
      }
      written.push_back( mtl );
      std::cout << "material " << i << std::endl;
    }
    //current_mtl = mtl;

    //}
  }
}

