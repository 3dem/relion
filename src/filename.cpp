/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres"
 * MRC Laboratory of Molecular Biology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/
/***************************************************************************
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include "src/filename.h"
#include "src/funcs.h"
#include <unistd.h>

// Constructor with root, number and extension .............................
void FileName::compose(const std::string &str, long int no, const std::string &ext, int numberlength)
{
    *this = (FileName) str;
    if (no != -1)
    {

        char aux_str[numberlength+1];
        std::string tmp_fileformat;
        tmp_fileformat = (std::string) "%0" +
                         integerToString(numberlength)+
                         (std::string)"d";
        sprintf(aux_str, tmp_fileformat.c_str(), no);
        *this += aux_str;
    }

    if (ext != "")
        *this += (std::string)"." + ext;
}

// Constructor: prefix number and filename, mainly for selfiles..
void FileName::compose(long int no , const std::string &str, int numberlength)
{
    *this = (FileName) str;
    if (no != -1)
    {

        char aux_str[numberlength+1];
        std::string tmp_fileformat;
        tmp_fileformat = (std::string) "%0" +
                         integerToString(numberlength)+
                         (std::string)"d@";
        sprintf(aux_str, tmp_fileformat.c_str(), no);
        *this = aux_str + str;
    }
    else
        *this = str;


}

// Is in stack ............................................................
bool FileName::isInStack() const
{
    return find("@") != std::string::npos;
}

// Decompose ..............................................................
void FileName::decompose(long int &no, std::string &str) const
{
    size_t idx = find('@');
    if(idx != std::string::npos)
    {
        no = textToInteger(substr(0,idx));
        str = substr(idx+1,length()-idx);
    }
    else{
      no=-1;
      str = *this;
    }
}

// Convert to lower case characters .........................................
FileName FileName::toLowercase() const
{
    FileName result = *this;
    for(unsigned int i=0;i<result.length();i++)
        result[i] = tolower(result[i]);
    return result;
}

// Convert to upper case characters .........................................
FileName FileName::toUppercase() const
{
    FileName result = *this;
    for(unsigned int i=0;i<result.length();i++)
        result[i] = toupper(result[i]);
    return result;
}

// Is substring present?
bool FileName::contains(const std::string& str) const
{
    int point = rfind(str);
    if (point > -1)
        return true;
    else
        return false;
}

// Get substring before first instance of str
FileName FileName::beforeFirstOf(const std::string& str) const
{
    int point = find(str);
    if (point > -1)
        return substr(0, point);
    else
        return *this;
}

// Get substring before last instance of str
FileName FileName::beforeLastOf(const std::string& str) const
{
    int point = rfind(str);
    if (point > -1)
        return substr(0, point);
    else
        return *this;
}

// Get substring after first instance of str
FileName FileName::afterFirstOf(const std::string& str) const
{
    int point = find(str);
    if (point > -1)
        return substr(point + 1);
    else
        return *this;
}

// Get substring after last instance of str
FileName FileName::afterLastOf(const std::string& str) const
{
    int point = rfind(str);
    if (point > -1)
        return substr(point + 1);
    else
        return *this;
}

// Get the base name of a filename .........................................
std::string FileName::getBaseName() const
{
    std::string basename = "";
    std::string myname = *this;
    int myindex = 0;
    for (int p = myname.size() - 1; p >= 0; p--)
    {
        if (myname[p] == '/')
        {
            myindex = p + 1;
            break;
        }
    }
    for (int p = myindex; p < myname.size(); p++)
    {
        if (myname[p] != '.')
            basename += myname[p];
        else
            break;
    }
    return basename;
}

// Get the extension of a filename .........................................
std::string FileName::getExtension() const
{
    int last_point = find_last_of(".");
    if (last_point == -1)
        return "";
    else
        return substr(last_point + 1);
}

// Init random .............................................................
void FileName::initRandom(int length)
{
    randomize_random_generator();
    *this = "";
    for (int i = 0; i < length; i++)
        *this += 'a' + FLOOR(rnd_unif(0, 26));
}

// Add at beginning ........................................................
FileName FileName::addPrefix(const std::string &prefix) const
{
    FileName retval = *this;
    int skip_directories = find_last_of("/") + 1;
    return retval.insert(skip_directories, prefix);
}

// Add at the end ..........................................................
FileName FileName::addExtension(const std::string &ext) const
{
    if (ext == "")
        return *this;
    else
    {
        FileName retval = *this;
        retval = retval.append((std::string)"." + ext);
        return retval;
    }
}

// Remove last extension ...................................................
FileName FileName::withoutExtension() const
{
    FileName retval = *this;
    return retval.substr(0, rfind("."));
}

// Insert before extension .................................................
FileName FileName::insertBeforeExtension(const std::string &str) const
{
    return  withoutExtension() + str + "." + getExtension();
}

// Remove an extension wherever it is ......................................
FileName FileName::removeExtension(const std::string &ext) const
{
    int first = find((std::string)"." + ext);
    if (first == -1)
        return *this;
    else
    {
        FileName retval = *this;
        return retval.erase(first, 1 + ext.length());
    }
}

// Remove all extensions....................................................
FileName FileName::removeAllExtensions() const
{
    int first = rfind("/");
    first = find(".", first + 1);
    if (first == -1)
        return *this;
    else
        return substr(0, first);
}

// Replace all substrings
void FileName::replaceAllSubstrings(std::string from, std::string to)
{
	FileName result;
    size_t start_pos = 0;
    while((start_pos = (*this).find(from, start_pos)) != std::string::npos)
    {
        (*this).replace(start_pos, from.length(), to);
        start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
    }

}

FileName FileName::getFileFormat() const
{
    int first;
    FileName result;
    if (find("#") != std::string::npos)
        return "raw";
    else if ( (first = rfind(":"))!=std::string::npos)
        result = substr(first + 1) ;
    else if ( (first = rfind("."))!=std::string::npos)
        result = substr(first + 1);
    else
        result="spi";
    return result.toLowercase();

}

FileName FileName::removeFileFormat() const
{
    if ( find("#", 0) > -1 )
        REPORT_ERROR("Not implemented for raw data");
    size_t found=rfind(":");
    if (found!=std::string::npos)
        return substr(0, found);
    return *this;
}

bool FileName::isStarFile() const
{
    //file names containing @, : or % are not metadatas
    size_t found=this->find('@');
    if (found!=std::string::npos)
        return false;
    // Allow :star to indicate that file really is a STAR file!
    //found=this->find(':');
    //if (found!=std::string::npos)
    //    return false;
    found=this->find('#');
    if (found!=std::string::npos)
        return false;

    FileName ext = getFileFormat();
    if (ext=="star")
    {
        return true;
    }
    else
    {
    	return false;
    }
}

// Substitute one extension by other .......................................
FileName FileName::substituteExtension(const std::string &ext1,
                                       const std::string &ext2) const
{
    int first = find((std::string)"." + ext1);
    if (first == -1)
        return *this;
    else
    {
        FileName retval = *this;
        return retval.replace(first, 1 + ext1.length(), (std::string)"." + ext2);
    }
}

// Remove a substring ......................................................
FileName FileName::without(const std::string &str) const
{
    if (str.length() == 0)
        return *this;
    int pos = find(str);
    if (pos == -1)
        return *this;
    else
    {
        FileName retval = *this;
        return retval.erase(pos, str.length());
    }
}

// Remove until prefix .....................................................
FileName FileName::removeUntilPrefix(const std::string &str) const
{
    if (str.length() == 0)
        return *this;
    int pos = find(str);
    if (pos == -1)
        return *this;
    else
    {
        FileName retval = *this;
        return retval.erase(0, pos + str.length());
    }
}

// Remove directories ......................................................
FileName FileName::removeDirectories(int keep) const
{
    int last_slash = rfind("/");
    int tokeep = keep;
    while (tokeep > 0)
    {
        last_slash = rfind("/", last_slash - 1);
        tokeep--;
    }
    if (last_slash == -1)
        return *this;
    else
        return substr(last_slash + 1, length() - last_slash);
}

size_t FileName::getFileSize() const
{
	struct stat filestatus;
	stat( this->c_str(), &filestatus );
	return filestatus.st_size;

}

int FileName::globFiles(std::vector<FileName> &files, bool do_clear) const
{
	if (do_clear)
		files.clear();

	glob_t glob_result;
	glob((*this).c_str(), GLOB_TILDE, NULL, &glob_result);
	for(unsigned  int  i = 0; i < glob_result.gl_pathc; ++i)
	{
		files.push_back(std::string(glob_result.gl_pathv[i]));
	}
	globfree(&glob_result);

	return files.size();
}

bool exists(const FileName &fn)
{

    struct stat buffer;
    return (stat (fn.c_str(), &buffer) == 0);
}

void touch(const FileName &fn)
{
	std::ofstream  fh;
    fh.open(fn.c_str(), std::ios::out);
    if (!fh)
        REPORT_ERROR( (std::string)"Filename::touch ERROR: Cannot open file: " + fn);
    fh.close();

}

void copy(const FileName &fn_src, const FileName &fn_dest)
{
    std::ifstream srce( fn_src.c_str(), std::ios::binary ) ;
    std::ofstream dest( fn_dest.c_str(), std::ios::binary ) ;
    dest << srce.rdbuf() ;

}

int mktree(const FileName &fn_dir, mode_t mode)
{
    std::string s = fn_dir;
	size_t pre=0,pos;
    std::string dir;
    int mdret;

    // force trailing / so we can handle everything in loop
    if(s[s.size()-1]!='/')
        s+='/';

    while((pos=s.find_first_of('/',pre))!=std::string::npos)
    {
        dir=s.substr(0,pos++);
        pre=pos;
        // if leading / first time is 0 length
        if (dir.size() == 0)
        	continue;

        if ((mdret = mkdir(dir.c_str(), mode)) && errno != EEXIST)
        {
            return mdret;
        }
    }

    return mdret;
}

bool decomposePipelineFileName(FileName fn_in, FileName &fn_pre, FileName &fn_jobnr, FileName &fn_post)
{

	size_t slashpos = 0;
	int i = 0;
	while (slashpos < fn_in.length())
	{
		i++;
		slashpos = fn_in.find("/", slashpos+1);
		if (fn_in[slashpos+1]=='j' && fn_in[slashpos+2]=='o' && fn_in[slashpos+3]=='b' &&
			std::isdigit(fn_in[slashpos+4]) && std::isdigit(fn_in[slashpos+5]) && std::isdigit(fn_in[slashpos+6]))
		{
			// find the second slash
			size_t slashpos2 = fn_in.find("/", slashpos+6);
			if (slashpos2 == std::string::npos)
				slashpos2 = fn_in.length() - 1;
			fn_pre = fn_in.substr(0, slashpos+1); // this has the first slash
			fn_jobnr = fn_in.substr(slashpos+1, slashpos2-slashpos); // this has the second slash
			fn_post = fn_in.substr(slashpos2+1); // this has the rest
			return true;
		}
	    // TODO: temporary check for - in uniq filename for backward compatibility with early alpha version. Remove in near future!!
		else if (std::isdigit(fn_in[slashpos+1]) && std::isdigit(fn_in[slashpos+2]) && std::isdigit(fn_in[slashpos+3]) &&
		    std::isdigit(fn_in[slashpos+4]) && std::isdigit(fn_in[slashpos+5]) && std::isdigit(fn_in[slashpos+6]) &&
		    (fn_in[slashpos+7] == '.' || fn_in[slashpos+7] == '-') &&
		    std::isdigit(fn_in[slashpos+8]) && std::isdigit(fn_in[slashpos+9]) && std::isdigit(fn_in[slashpos+10]) &&
		    std::isdigit(fn_in[slashpos+11]) && std::isdigit(fn_in[slashpos+12]) && std::isdigit(fn_in[slashpos+13]) )
		{
			fn_pre = fn_in.substr(0, slashpos+1); // this has the first slash
			fn_jobnr = fn_in.substr(slashpos+1,14); // this has the second slash
			fn_post = fn_in.substr(slashpos+15); // this has the rest
			return true;
		}
		if (i>20)
			REPORT_ERROR("decomposePipelineFileName: BUG or found more than 20 directories deep structure for pipeline filename: " + fn_in);
	}
	// This was not a pipeline filename
	fn_pre="";
	fn_jobnr="";
	fn_post=fn_in;
	return false;

}

bool decomposePipelineSymlinkName(FileName fn_in, FileName &fn_pre, FileName &fn_jobnr, FileName &fn_post)
{

	// Symlinks are always in the second directory....
	size_t slashpos = 0;
	int i = 0;
	while (slashpos < fn_in.length())
	{
		i++;
		slashpos = fn_in.find("/", slashpos+1);
		if (i==2)
			break;
	}

	// Check whether this is a symbol link
	char linkname[100];
	ssize_t len = ::readlink(fn_in.substr(0, slashpos).c_str(), linkname, sizeof(linkname)-1);
	if (len != -1)
    {
    	// This is a symbolic link!
    	linkname[len] = '\0';
    	FileName fn_link = std::string(linkname);
    	fn_link = fn_link.afterFirstOf("../") + fn_in.substr(slashpos+1);
    	// So dereference the link, BUT only if the second directory started with "job"!
    	if (decomposePipelineFileName(fn_link, fn_pre, fn_jobnr, fn_post))
    	{
    		return true;
    	}
    	else
    	{
    		fn_pre = fn_jobnr = "";
    		fn_post = fn_in;
    		return false;
    	}
    }

	// If it is not a symlink, just decompose the filename
    return decomposePipelineFileName(fn_in, fn_pre, fn_jobnr, fn_post);

}

FileName getOutputFileWithNewUniqueDate(FileName fn_input, FileName fn_new_outputdir)
{
	FileName fn_pre, fn_jobnr, fn_post;
	decomposePipelineFileName(fn_input, fn_pre, fn_jobnr, fn_post);
	return fn_new_outputdir + fn_post;
}

