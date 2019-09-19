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
 *
 * Authors:     J.R. Bilbao-Castro (jrbcast@ace.ual.es)
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

#ifndef XMIPP_STRINGS_H
#define XMIPP_STRINGS_H

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <string.h>
#include <stdio.h>
#include "src/macros.h"

/// @defgroup StringUtilities String utilities
/// @ingroup DataLibrary
//@{

//@name String processing
//@{

/** Removes all occurrences of 'character' from the string no matter
where they are */
std::string removeChar(const std::string &str, char character);

/** Removes escaped symbols ESC+n, t, v, b, r, f, and a
 * Note that tabs are replaced by spaces.
 * */
std::string unescape(const std::string &str);

void escapeStringForSTAR(std::string &value);

/** Best precision for a float number.
 *
 * This function returns the best precision to be used in a "printf" format if
 * this number is to fit in a given width. It returns -1 if the exponential
 * format is advised.
 *
 * @code
 * template<typename T>
 * std::ostream& operator<<(std::ostream& out, const T& val)
 * {
 *     int i,j;
 *
 *     if (val.xdim == 0)
 *         out << "NULL matrix" << std::endl;
 *     else
 *         out << std::endl;
 *
 *     T aux = ABSnD(val);
 *     int prec = bestPrecision(aux.max(), 10);
 *
 *     for (i=STARTINGY(val); i<=FINISHINGY(val); i++)
 *     {
 *         for (j=STARTINGX(val); j<=FINISHINGX(val); j++)
 *         {
 *             out << floatToString((float) val(i,j), 10, prec) << ' ';
 *         }
 *         out << std::endl;
 *     }
 *
 *     return out;
 * }
 *
 * @endcode
 */
int bestPrecision(float F, int _width);

bool isNumber(std::string);

/** String (char*) to double conversion.
 *
 * @code
 * RFLOAT key = textToDouble(firstToken(line), 1602, "Error reading key");
 * @endcode
 */
double textToDouble(const char* str,
                    int _errno = 2101,
                    std::string errmsg = "Error in textToDouble");

/** String (STL) to double conversion.
 *
 * @code
 * double key = textToDouble(str, 1602, "Error reading key");
 * @endcode
 */
inline double textToDouble(const std::string& str,
                         int _errno = 2101,
                         std::string errmsg = "Error in textToDouble")
{
    return textToDouble(str.c_str(), _errno, errmsg);
}

/** String (char*) to float conversion.
 *
 * @code
 * float key = textToFloat(firstToken(line), 1602, "Error reading key");
 * @endcode
 */
float textToFloat(const char* str,
                  int _errno = 2101,
                  std::string errmsg = "Error in textToFloat");

/** String (STL) to float conversion.
 *
 * @code
 * float key = textToFloat(str, 1602, "Error reading key");
 * @endcode
 */
inline float textToFloat(const std::string& str,
                         int _errno = 2101,
                         std::string errmsg = "Error in textToFloat")
{
    return textToFloat(str.c_str(), _errno, errmsg);
}

/** String (char*) to integer conversion.
 *
 * @code
 * int param_no = textToInteger(nextToken(), 1602, "Error reading number parameters")
 * @endcode
 */
int textToInteger(const char* str,
                  int _errno = 2102,
                  std::string errmsg = "Error in textToInteger");

/** String (STL) to integer conversion.
 *
 * @code
 * int param_no = textToInteger(str, 1602, "Error reading number parameters")
 * @endcode
 */
inline int textToInteger(const std::string& str,
                         int _errno = 2102,
                         std::string errmsg = "Error in textToInteger")
{
    return textToInteger(str.c_str(), _errno, errmsg);
}

/** String (char*) to boolean conversion.
 *
 * @code
 * bool param_no = textToBool(nextToken(), 1602, "Error reading boolean parameters")
 * @endcode
 */
bool textToBool(const char* str,
                  int _errno = 2102,
                  std::string errmsg = "Error in textToBool");

/** String (STL) to boolean conversion.
 *
 * @code
 * int param_no = textToBool(str, 1602, "Error reading boolean parameters")
 * @endcode
 */
inline int textToBool(const std::string& str,
                         int _errno = 2102,
                         std::string errmsg = "Error in textToBool")
{
    return textToBool(str.c_str(), _errno, errmsg);
}

/** String (char*) to long long integer conversion.
 *
 * @code
 * long long param_no = textToLongLong(nextToken(), 1602, "Error reading number
 *     parameters")
 * @endcode
 */
long long textToLongLong(const char* str,
                         int _errno = 2102,
                         std::string errmsg = "Error in AtoL");

/** String (STL) to long long conversion.
 *
 * @code
 * int param_no = textToLongLong(str, 1602, "Error reading number parameters")
 * @endcode
 */
inline int textToLongLong(const std::string& str,
                         int _errno = 2102,
                         std::string errmsg = "Error in textToInteger")
{
    return textToLongLong(str.c_str(), _errno, errmsg);
}


/** Float to string conversion.
 *
 * If precision==0 the precision is automatically computed in such a way that
 * the number fits the width (the exponential format might be chosen). If
 * precision==-1 then the exponential format is forced. If width==0 then the
 * minimum width is used.
 *
 * @code
 * REPORT_ERROR(1602, "Value not recognised " + floatToString(val));
 * @endcode
 */
std::string floatToString(float F, int _width = 0, int _prec = 0);

/** Integer to string conversion.
 *
 * If width==0 then writes the number with the number of digits needed. The
 * fill_with field indicates which is the filling character for the left
 * positions.
 *
 * @code
 * REPORT_ERROR(1602, "Error reading key " + integerToString(key));
 * @endcode
 */
std::string integerToString(int I, int _width = 0, char fill_with = '0');

/** Character to integer conversion.
 *
 * Takes a character and produces a number according to its ASCII code minus 48.
 * For instance, ASCII=48 produces number 0, ASCII=49 produces 1, ..., ASCII=57
 * produces 9, ASCII=58 produces 10!!, ... This is used when you have codified
 * numbers greater than 9 in a single character.
 *
 * @code
 * int param_no = textToInt(token, 1602, "Error reading number parameters");
 * @endcode
 */
int textToInt(const char* str,
              int _errno = 2103,
              std::string errmsg = "Error in textToInt");

/** String to string with given length conversion.
 *
 * The output string will have the information of the input one with the given
 * width. If the width is smaller than the string length then the string is
 * truncated and if it is greater the string is right padded with spaces. If
 * width==0 then the same string is returned.
 */
std::string stringToString(const std::string& str, int _width = 0);

/** Check angle.
 *
 * If the argument is not "rot", "tilt" nor "psi" an exception is thrown
 */
void checkAngle(const std::string& str);

/** To lower.
 *
 * All characters between A-Z are brought to a-z. Result is rewritten on input
 * string
 */
void toLower(char* _str);

/** To lower, for STL strings.
 */
void toLower(std::string& _str);

/** Removes white spaces from the beginning and the end of the string
as well as escaped characters
and simplifies the rest of groups of white spaces of the string to
a single white space */
std::string simplify( const std::string& str );

/** Remove trailing spaces */
void trim(std::string& str);

/** Remove consecutive spaces.
 *
 * All consecutive spaces are replaced by a single one and starting and
 * finishing spaces are removed
 */
std::string removeSpaces(const std::string& _str);

/** Remove quotes.
 *
 * This function removes the first character if it is a RFLOAT or single quote,
 * as well as the last character. The char pointer might be moved.
 *
 * @code
 * char str[10] = "\"Hello\"";
 * (&str);
 * @endcode
 */
void removeQuotes(char** _str);
//@}

/** @name Tokenization
 *
 * These functions allow to split a string into small pieces separated by blank
 * spaces, giving you a pointer to the different word each time. The different
 * elements from the string are selected using strtok, so after the application
 * of this function to the input string, this is modified and NULL characters
 * are introduced as delimiters of the elements. This is useful in most
 * situations since after reading a list you might go on reading more things,
 * but you must be aware of it. Here goes an example of doing so:
 *
 * @code
 * std::cout << "Whole  line: " << line << std::endl;
 * std::cout << "First  word: " << firstToken(line) << std::endl;
 * std::cout << "Second word: " << nextToken() << std::endl;
 * std::cout << "Third  word: " << nextToken() << std::endl;
 * ...
 * @endcode
 *
 * When there are no more words, both functions return a NULL pointer. Here we
 * make a distinction between tokens (words that might be empty) and words
 * (words that cannot be empty, if they are then an exception or an exit error
 * is thrown).
 *
 * For STL there is another way. You supply a string object and a vector of
 * strings is returned with all the elements
 */
//@{
/** Split a STL string given some delimiter.
 *
 * Returns a the number of tokens found. The tokens are in the variable results.
 */
int splitString(const std::string& input,
                const std::string& delimiter,
                std::vector< std::string >& results,
                bool includeEmpties = false);

/** Returns first token (char*).
 *
 * @code
 * char line[80];
 *
 * std::cout << "First  word: " << firstToken(line) << std::endl;
 * @endcode
 */
inline char* firstToken(const char* str)
{
    return strtok((char*) str, " \t\n");
}

/** Returns first token (STL).
 *
 * @code
 * std::string line;
 *
 * std::cout << "First  word: " << firstToken(line) << std::endl;
 * @endcode
 */
inline char* firstToken(const std::string& str)
{
    return strtok((char*) str.c_str(), " \t\n");
}

/** Returns next token.
 *
 * This functions returns the next word of the line we have given last as
 * parameter to firstToken.
 *
 * @code
 * char line[80];
 * ...
 * firstToken(line);
 * std::cout << "Second  word: " << nextToken(line) << std::endl;
 *
 * stf::string line;
 * ...
 * firstToken(line);
 * std::cout << "Second  word: " << nextToken(line) << std::endl;
 * @endcode
 */
inline char* nextToken()
{
    return strtok((char*) NULL, " \t\n");
}

/** Returns next token.
 *
 * It reads from position i. Returns (in i) the following position to search on.
 * When there are no more tokens. It returns "".
 */
std::string nextToken(const std::string& str, int& i);

bool nextTokenInSTAR(const std::string& str, int& i, std::string &token);

/** Get non empty string (char*).
 *
 * This function returns the first word found in the given line disregarding the
 * leading blanks. If no word is found then an exception or an exit error is
 * produced. After calling this function the first blank character after the
 * word is substituted by a NULL character (as it uses the function firstToken.
 * Further word readings should use the function read_nextWord
 */
char* firstWord(char* str,
                int _errno = 2106,
                const std::string & errmsg = "first word: String not found");

/** Get non empty string (STL).
 *
 * Same as the previous function but for STL strings
 */
inline char* firstWord(std::string& str,
                       int _errno = 2106,
                       std::string errmsg = "first word: String not found")
{
    // FIXME C-style cast
    return firstWord((char*) str.c_str(), _errno, errmsg);
}

/** Get next non empty string.
 *
 * This is the same as the nextToken, but an exception is thrown or an exit
 * error produced if the word is empty
 */
inline char* nextWord(int _errno = 2106,
                      std::string errmsg = "next word: String not found")
{
    return firstWord((char*) NULL, _errno, errmsg);
}

/** Tokenize a string and return a list of tokens
 *
 */
void tokenize(const std::string& str,
              std::vector< std::string >& tokens,
              const std::string& delimiters = " \t");
//@}
//@}
#endif
