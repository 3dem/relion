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
 * Authors: Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *
 *
 * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
 * 02111-1307 USA
 *
 * All comments concerning this program package may be sent to the
 * e-mail address 'xmipp@cnb.csic.es'
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

		char aux_str[numberlength+2];
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
	if (idx != std::string::npos)
	{
		bool ok = true;
		for (int i = 0; i < idx; i++)
		{
			if ((*this)[i] < '0' || (*this)[i] > '9')
			{
				ok = false;
				break;
			}
		}

		if (ok)
		{
			no = textToInteger(substr(0, idx));
			str = substr(idx + 1,length() - idx);
			return;
		}
	}

	no = -1;
	str = *this;
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
		return substr(point + str.length());
	else
		return *this;
}

// Get substring after last instance of str
FileName FileName::afterLastOf(const std::string& str) const
{
	int point = rfind(str);
	if (point > -1)
		return substr(point + str.length());
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
	return	withoutExtension() + str + "." + getExtension();
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
		result="";
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
	//	  return false;
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
	for(unsigned long int i = 0; i < glob_result.gl_pathc; ++i)
	{
		files.push_back(std::string(glob_result.gl_pathv[i]));
	}
	globfree(&glob_result);

	return files.size();
}

bool FileName::getTheOtherHalf(FileName &fn_out) const
{
	FileName ret = this->afterLastOf("/");

	if (ret.contains("half1"))
		ret.replaceAllSubstrings("half1", "half2");
	else if (ret.contains("half2"))
		ret.replaceAllSubstrings("half2", "half1");
	else
		return false;

	if (this->contains("/"))
		ret = this->beforeLastOf("/") + "/" + ret;

	fn_out = ret;
	return true;
}

bool FileName::getHalf(FileName &fn_out, int halfset) const
{
	FileName ret = this->afterLastOf("/");

	if (ret.contains("half1") && halfset == 2)
		ret.replaceAllSubstrings("half1", "half2");
	else if (ret.contains("half2") && halfset == 1)
		ret.replaceAllSubstrings("half2", "half1");
	else
		return false;

	if (this->contains("/"))
		ret = this->beforeLastOf("/") + "/" + ret;

	fn_out = ret;
	return true;
}

bool FileName::validateCharactersStrict(bool do_allow_double_dollar) const
{

	FileName myname = *this;
	if (do_allow_double_dollar)
	{
		myname.replaceAllSubstrings("$$","");
	}
	for (int p = 0, len = myname.size(); p < len; p++)
	{
		char c = myname[p];
		if (!(c >= '0' && c <= '9') && // 0-9
		    !(c >= 'A' && c <= 'Z') && // A-Z
		    !(c >= 'a' && c <= 'z') && // a-z
		      c != '-')
			return false;
	}

	return true;
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

void move(const FileName &fn_src, const FileName &fn_dest)
{
	copy(fn_src, fn_dest);
	remove(fn_src.c_str());
}

int mktree(const FileName &fn_dir, mode_t mode)
{
	if (exists(fn_dir)) return 0;

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

FileName realpath(const FileName &fn, bool allow_nonexisting_path)
{
	char retval[PATH_MAX];

	if (allow_nonexisting_path)
	{
		char *dont_care = realpath(fn.c_str(), retval);
		if (retval == NULL)
		{
			REPORT_ERROR(std::string("Failed to resolve realpath of ") + fn);
		}
	}
	else
	{
		if (realpath(fn.c_str(), retval) == NULL)
		{
			REPORT_ERROR(std::string("Failed to resolve realpath of ") + fn);
		}
	}

	return FileName(retval);
}

void symlink(const FileName &src, const FileName &dst)
{
	// If the output file exists, remove it first
	if (exists(dst))
	{
		if (remove(dst.c_str()) != 0)
		{
			REPORT_ERROR(std::string("Failed to remove existing file ") + dst.c_str() + " to make symlink ");
		}
	}

	if (symlink(src.c_str(), dst.c_str()) != 0)
	{
		REPORT_ERROR(std::string("Failed to make a symlink from ") + src.c_str() + " to " + dst.c_str());
	}
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
	bool dont_expand = false;

	// Symlinks are always in the second directory. (e.g. Refine3D/JOB_ALIAS_AS_LINK/...)
	size_t slashpos = 0;
	int i = 0;
	while (slashpos < fn_in.length())
	{
		i++;
		slashpos = fn_in.find("/", slashpos+1);
		if (i==2)
			break;
	}

	// We ignore links in the first directory (link/XXX)
	if (i != 2)
		dont_expand = true;

	FileName second_dir = fn_in.substr(0, slashpos);
//	std::cout << second_dir << std::endl;
	// We also ignore jobXXX even if it is a symbolic link.
	//  e.g. MotionCorr/job003 ==> /path/to/online_processing/MotionCorr/job003
	// This is safe because we don't allow an alias name to start from 'job'.
	if (second_dir.afterLastOf("/").substr(0, 3) == "job")
		dont_expand = true;

	// Check whether this is a symbol link
	char linkname[4096];
	ssize_t len_max = sizeof(linkname) - 1;
	ssize_t len = ::readlink(second_dir.c_str(), linkname, len_max);
	if (len == len_max)
		REPORT_ERROR("Too long path in decomposePipelineSymlinkName.");
	if (!dont_expand && len != -1)
	{
			// This is a symbolic link!

			if (linkname[len - 1] == '/')
				linkname[len - 1] = '\0'; // remove trailing '/'
			linkname[len] = '\0';

			FileName fn_link = std::string(linkname);
			// TODO: FIXME: This condition is still not perfect. For example,
			// Micrograph/mic001.mrc -> ../../../storage/mic001.mrc breaks the code.
			// Meanwhile one can circumvent this case by using an absolute path in the symlink.
			if (fn_link.substr(0, 3) == "../")
			{
				fn_link = fn_link.substr(3) + fn_in.substr(slashpos);
//			std::cout << "fn_link=" << fn_link << std::endl;
				return decomposePipelineFileName(fn_link, fn_pre, fn_jobnr, fn_post);
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

