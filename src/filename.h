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
 * 02111-1307  USA
 *
 * All comments concerning this program package may be sent to the
 * e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#ifndef FILENAME_H_
#define FILENAME_H_

#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <climits>
#include <algorithm>
#include <vector>
#include <typeinfo>
#include <sys/stat.h>
#include <glob.h>
#include <errno.h>
#include "src/numerical_recipes.h"
#include "src/macros.h"
#include "src/error.h"
#include "src/strings.h"

#define FILENAMENUMBERLENGTH 6

//@{
/** Filenames.
 *
 * This class allows you a lot of usual and common manipulations with filenames.
 * See filename conventions for a detailed explanation of the Filenames dealed
 * here, although most of the functions work with the more general model
 * "name.extension"
 */
class FileName: public std::string
{
public:
	/// @name Filename constructors
	/// @{

	/** Empty constructor
	 *
	 * The empty constructor is inherited from the string class, so an empty
	 * FileName is equal to "".
	 *
	 * @code
	 * FileName fn_blobs;
	 * @endcode
	 */
	FileName(): std::string("") {}

	/* Destructor
	 */
	~FileName() {}

	/** Constructor from string
	 *
	 * The constructor from a string allows building complex expressions based
	 * on the string class. Notice that in the following example the type
	 * casting to string is very important, if not, the operation is just a
	 * pointer movement instead of a string concatenation.
	 *
	 * @code
	 * FileName fn_blobs((std::string) "art00001" + ".blobs");
	 * @endcode
	 */
	FileName(const std::string& str): std::string(str) {}

	/** Constructor from char*
	 */
	FileName(const char* str): std::string(str) {}

	/** Copy constructor
	 */
	FileName(const FileName& fn): std::string(fn) {}

	/** Assignment constructor
	 */
	FileName& operator=(const FileName& op) {
		return (FileName&) std::string::operator=(op);
	}

	/** Constructor from root, number and extension
	 *
	 * The number and extension are optional.
	 *
	 * @code
	 * FileName fn_proj("g1ta000001.xmp"); // fn_proj = "g1ta000001.xmp"
	 * FileName fn_proj("g1ta",1,"xmp"); // fn_proj = "g1ta000001.xmp"
	 * FileName fn_proj("g1ta",1); // fn_proj = "g1ta000001"
	 * @endcode
	 */
	FileName(const char* str, long int no, const std::string& ext = "")
	{
		compose(str, no, ext);
	}

	/** Constructor from root and extension
	 *
	 * None of the parameters is optional
	 *
	 * @code
	 * FileName fn_proj("g1ta00001", "xmp"); // fn_proj = "g1ta00001.xmp"
	 * @endcode
	 */
	FileName(const char* str, const std::string& ext): std::string(str + ext)
	{}
	//@}

	/// @name Composing/Decomposing the filename
	/// @{

	/** Compose from root, number and extension
	 *
	 * @code
	 * fn_proj.compose("g1ta", 1, "xmp");  // fn_proj = "g1ta000001.xmp"
	 * @endcode
	 */
	void compose(const std::string& str, long int no, const std::string& ext, int numberlength = FILENAMENUMBERLENGTH);

	/** Prefix with number @. Mainly for selfiles
	 *
	 * @code
	 * fn_proj.compose(1,"g1ta.xmp");  // fn_proj = "000001@g1ta.xmp"
	 * @endcode
	 */
	void compose(long int no, const std::string& str, int numberlength = FILENAMENUMBERLENGTH);

	/** True if this filename belongs to a stack
	 */
	bool isInStack() const;

	/** Decompose filenames with @. Mainly from selfiles
	 *
	 * @code
	 * fn_proj.decompose(no,filename); // fn_proj = "000001@g1ta000001.xmp"
	 *                                 // no=1
	 *                                 // filename = "g1ta000001.xmp"
	 * @endcode
	 */
	void decompose(long int &no, std::string& str) const;

	/** Get the base name from a filename
	 */
	std::string getBaseName() const;

	/** Get the last extension from filename
	 *
	 * The extension is returned without the dot. If there is no extension "" is
	 * returned.
	 *
	 * @code
	 * std::string ext = fn_proj.get_extension();
	 * @endcode
	 */
	std::string getExtension() const;

	/** Get image format identifier (as in Bsoft)
	 *
	 * @code
	 * fn_proj = "g1ta00001.xmp";
	 * fn_proj = fn_proj.get_file_format(); // fn_proj == "xmp"
	 * fn_proj = "g1ta00001.nor:spi";
	 * fn_proj = fn_proj.get_file_format(); // fn_proj == "spi"
	 * fn_proj = "input.file#d=f#x=120,120,55#h=1024";
	 * fn_proj = fn_proj.get_file_format(); // fn_proj == "raw"
	 * @endcode
	 */
	FileName getFileFormat() const;

	/** Random name
	 *
	 * Generate a random name of the desired length.
	 */
	void initRandom(int length);
	//@}

	///@name Filename utilities
	//@{
	/** Change all characters for lowercases
	 *
	 * @code
	 * FileName fn_proj("g1tA00001");
	 * fn_proj = fn_proj.to_lowercase(); // fn_proj = "g1ta00001"
	 * @endcode
	 */
	FileName toLowercase() const;

	/** Change all characters for uppercases
	 *
	 * @code
	 * FileName fn_proj("g1tA00001");
	 * fn_proj = fn_proj.to_uppercase(); // fn_proj = "G1Ta00001"
	 * @endcode
	 */
	FileName toUppercase() const;

	/** Check whether the filename contains the argument substring
	 *
	 * @code
	 * FileName fn_proj("g1ta00001.raw#d=f");
	 * if (fn_proj.contains("raw) ) // true
	 * @endcode
	 */
	bool contains(const std::string& str) const;

	/** Return substring before first instance of argument (as in Bsoft)
	 *
	  * @code
	 * FileName fn_proj("g1ta00001.raw#d=f");
	 * fn_proj = fn_proj.before_first_of("#"); // fn_proj = "g1ta00001.raw"
	 * @endcode
	 */
	FileName beforeFirstOf(const std::string& str) const;

	/** Return substring before last instance of argument (as in Bsoft)
	 *
	  * @code
	 * FileName fn_proj("g1ta00001.raw#d=f");
	 * fn_proj = fn_proj.before_last_of("#"); // fn_proj = "g1ta00001.raw"
	 * @endcode
	 */
	FileName beforeLastOf(const std::string& str) const;

	/** Return substring after first instance of argument (as in Bsoft)
	 *
	  * @code
	 * FileName fn_proj("g1ta00001.raw#d=f");
	 * fn_proj = fn_proj.after_first_of("#"); // fn_proj = "d=f"
	 * @endcode
	 */
	FileName afterFirstOf(const std::string& str) const;

	/** Return substring after last instance of argument (as in Bsoft)
	 *
	  * @code
	 * FileName fn_proj("g1ta00001.raw#d=f");
	 * fn_proj = fn_proj.after_last_of("#"); // fn_proj = "d=f"
	 * @endcode
	 */
	FileName afterLastOf(const std::string& str) const;

	/** Add string at the beginning
	 *
	 * If there is a path then the prefix is added after the path.
	 *
	 * @code
	 * fn_proj = "imgs/g1ta00001";
	 * fn_proj.add_prefix("h"); // fn_proj == "imgs/hg1ta00001"
	 *
	 * fn_proj = "g1ta00001";
	 * fn_proj.add_prefix("h"); // fn_proj == "hg1ta00001"
	 * @endcode
	 */
	FileName addPrefix(const std::string& prefix) const;

	/** Add extension at the end.
	 *
	 * The "." is added. If teh input extension is "" then the same name is
	 * returned, with nothing added.
	 *
	 * @code
	 * fn_proj = "g1ta00001";
	 * fn_proj.add_extension("xmp"); // fn_proj == "g1ta00001.xmp"
	 * @endcode
	 */
	FileName addExtension(const std::string& ext) const;

	/** Remove last extension, if any
	 *
	 * @code
	 * fn_proj = "g1ta00001.xmp";
	 * fn_proj = fn_proj.without_extension(); // fn_proj == "g1ta00001"
	 *
	 * fn_proj = "g1ta00001";
	 * fn_proj = fn_proj.without_extension(); // fn_proj == "g1ta00001"
	 * @endcode
	 */
	FileName withoutExtension() const;

	/** Insert before first extension
	 *
	 * If there is no extension, the insertion is performed at the end.
	 *
	 * @code
	 * fn_proj = "g1ta00001.xmp";
	 * fn_proj = fn_proj.insert_before_extension("pp");
	 * // fn_proj == "g1ta00001pp.xmp"
	 *
	 * fn_proj = "g1ta00001";
	 * fn_proj = fn_proj.insert_before_extension("pp");
	 * // fn_proj=="g1ta00001pp"
	 * @endcode
	 */
	FileName insertBeforeExtension(const std::string& str) const;

	/** Remove a certain extension
	 *
	 * It doesn't matter if there are several extensions and the one to be
	 * removed is in the middle. If the given extension is not present in the
	 * filename nothing is done.
	 *
	 * @code
	 * fn_proj = "g1ta00001.xmp.bak";
	 * fn_proj = fn_proj.remove_extension("xmp");
	 * // fn_proj == "g1ta00001.bak"
	 * @endcode
	 */
	FileName removeExtension(const std::string& ext) const;

	/** Remove all extensions
	 */
	FileName removeAllExtensions() const;

	/**
	 *	Replace all substrings
	 */
	void replaceAllSubstrings(std::string from, std::string to);

	/** Remove file format
	 * @code
	 * fn_proj = "g1ta00001.xmp";
	 * fn_proj = fn_proj.get_file_format(); // fn_proj == "xmp"
	 * fn_proj = "g1ta00001.nor:spi";
	 * fn_proj = fn_proj.get_file_format(); // fn_proj == "spi"
	 * fn_proj = "input.file#d=f#x=120,120,55#h=1024";
	 * fn_proj = fn_proj.get_file_format(); // fn_proj == "raw"
	 * @endcode
	 */
	FileName removeFileFormat() const;

	/** Is this file a MetaData file?
	 * Returns false if the filename contains "@", ":" or "#"
	 * Returns true if the get_file_format extension == "star"
	 */
	bool isStarFile() const;

	/** Clean image FileName (as in Bsoft)
	 *
	 * @code
	 * fn_proj = "g1ta00001.xmp";
	 * fn_proj = fn_proj.get_file_format(); // fn_proj == "g1ta00001.xmp"
	 * fn_proj = "g1ta00001.nor:spi";
	 * fn_proj = fn_proj.clean_image_name(); // fn_proj == "g1ta00001.nor"
	 * fn_proj = "input.file#d=f#x=120,120,55#h=1024";
	 * fn_proj = fn_proj.clean_image_name(); // fn_proj == "input.file"
	 * @endcode
	 */
	//FileName clean_image_name() const;

	/** Substitute ext1 by ext2
	 *
	 * It doesn't matter if ext1 is in the middle of several extensions. If ext1
	 * is not present in the filename nothing is done.
	 *
	 * @code
	 * fn_proj = "g1ta00001.xmp.bak";
	 * fn_proj = fn_proj.substitute_extension("xmp", "bor");
	 * // fn_proj == "g1ta00001.bor.bak"
	 *
	 * fn_proj = "g1ta00001.xmp.bak";
	 * fn_proj = fn_proj.substitute_extension("tob", "bor");
	 * // fn_proj=="g1ta00001.xmp.bak"
	 * @endcode
	 */
	FileName substituteExtension(const std::string& ext1,
	                             const std::string& ext2) const;

	/** Without a substring
	 *
	 * If the substring is not present the same FileName is returned, if it is
	 * there the substring is removed.
	 */
	FileName without(const std::string& str) const;

	/** Remove until prefix
	 *
	 * Remove the starting string until the given prefix, inclusively. For
	 * instance /usr/local/data/ctf-image00001.fft with ctf- yields
	 * image00001.fft. If the prefix is not found nothing is done.
	 */
	FileName removeUntilPrefix(const std::string& str) const;

	/** Remove all directories
	 *
	 * Or if keep>0, then keep the lowest keep directories
	 */
	FileName removeDirectories(int keep = 0) const;

	/*
	 * Gets the filesize (in bytes)
	 */
	size_t getFileSize() const;

	// Get the other half map by swapping half1 and half2
	bool getTheOtherHalf(FileName &fn_out) const;

	// Get my half map by swapping half# if required
	bool getHalf(FileName &fn_out, int halfset) const;

	bool validateCharactersStrict(bool do_allow_double_dollar = false) const;

	/** From a wild-card containing filename get a vector with all existing filenames,
	 * return number of existing filenames
	 * If do_clear, the output vector will be clear when starting, when false, files will just be added to the vector
	 *
	 */
	int globFiles(std::vector<FileName> &files, bool do_clear = true) const;
	//@}
};

/** This class is used for comparing filenames.
 *
 * Example: "g0ta00001.xmp" is less than "g0ta00002.xmp"
 *
 * This class is needed to define a std::map as
 * map<FileName,FileName,FileNameComparison> myMap;
 *
 * This function is not ported to Python.
 */
class FileNameComparison
{
public:
	inline bool operator ()(const FileName &fn1, const FileName &fn2)
	{
		return fn1<fn2;
	}
};

// The following 2 functions are for the pipelining of RELION-2.0
// Finds the 160101.093245 or job001 UNIQ-ID substring from a larger string, as used in relion-2.0 pipeline
// It returns the part before that id (+trailing slash), the id itself (plus trailing slash) and the part after the id
// The function returns true if the input filename was a pipeline one, and false otherwise
// If the input filename wasn't a pipeline one, then fn_post is set to the input filename and the fn_pre and fn_jobnr are left empty
bool decomposePipelineFileName(FileName fn_in, FileName &fn_pre, FileName &fn_jobnr, FileName &fn_post);
// Also same for possible symlink
bool decomposePipelineSymlinkName(FileName fn_in, FileName &fn_pre, FileName &fn_jobnr, FileName &fn_post);

// Replaces the UNIQDATE substring and its preceding Directory-structure from fn_input, and adds fn_new_outputdir in front of it
FileName getOutputFileWithNewUniqueDate(FileName fn_input, FileName fn_new_outputdir);

/** True if the file exists in the current directory
 *
 * @code
 * if (exists("g1ta00001"))
 *	   std::cout << "The file exists" << std::endl;
 * @endcode
 */
bool exists(const FileName& fn);

/** Touch a file on the file system. */
void touch(const FileName& fn);

/** Copy a file */
void copy(const FileName &fn_src, const FileName &fn_dest);

/** Move a file */
void move(const FileName &fn_src, const FileName &fn_dest);

/** Make a directory tree. mode=0777 is OK, because it is modified by umask.*/
int mktree(const FileName &fn_dir, mode_t mode = 0777);

/* wrapper to realpath in stdlib */
FileName realpath(const FileName &fn, bool allow_nonexisting_path = false);

/* wrapper to symlink in stdlib */
void symlink(const FileName &src, const FileName &dst);

/** True if the path is a directory */
bool isDirectory (const FileName &fn);

/** True if the file exists in the current directory
 *	Remove leading xx@ and tailing :xx
 *
 * @code
 * if (exists("g1ta00001"))
 *	   std::cout << "The file exists" << std::endl;
 * @endcode
 */
int existsTrim(const FileName& fn);

/** Return the list of files within a directory. */
void getdir(const std::string &dir, std::vector<FileName> &files);

/** This function raised an ERROR if the filename if not empty and if
 * the corresponding file does not exist.
 * This may be useful to have a better (killing) control on (mpi-regulated) jobs
 *
 * @code
 *	 exit_if_not_exists("control_file.txt");
 * @endcode
 *
 * This function is not ported to Python.
 */
void exit_if_not_exists(const FileName &fn);

/** Waits until the given filename has a stable size
 *
 * The stable size is defined as having the same size within two samples
 * separated by time_step (microsecs).
 *
 * An exception is throw if the file exists but its size cannot be stated.
 */
void wait_until_stable_size(const FileName& fn,
                            unsigned long time_step = 250000);

/** Write a zero filled file with the desired size.
 *
 * The file is written by blocks to speed up, you can modify the block size.
 * An exception is thrown if any error happens
 */
void create_empty_file(const FileName& fn,
                       unsigned long long size,
                       unsigned long long block_size = 102400);

/** Returns the base directory of the Xmipp installation
 */
FileName xmippBaseDir();
//@}

#endif /* FILENAME_H_ */
