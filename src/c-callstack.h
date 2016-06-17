/*
 * c-callstack.h
 *
 *  Created on: Feb 23, 2016
 *      Author: bjornf
 */

#ifndef C_CALLSTACK_H_
#define C_CALLSTACK_H_

/**************************************************************************
//
// Copyright 2013 Kangmo Kim, Nanolat Software.
//
// e-mail : kangmo@nanolat.com
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//        http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// c-callstack.h : Show Java-like callstack in C/C++ projects.
//
***************************************************************************/

#if defined(NDEBUG) /* Release Mode */

#  define NL_RETURN(rc) return (rc)

#else /* Debug Mode */

#  define NL_RETURN(rc)                                       \
   {                                                          \
     if ((rc)) {                                              \
       fprintf( stderr,                                       \
                "Error(code:%d) at : %s (%s:%d)\n",           \
                (rc), __FUNCTION__, __FILE__, __LINE__);      \
     }                                                        \
   }

#endif /* NDEBUG */


#endif /* C_CALLSTACK_H_ */
