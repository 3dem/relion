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
#include "src/memory.h"

char*  askMemory(unsigned long memsize) 
{ 
    char*		ptr = NULL;
    
    if ( memsize == 0 ) {
        REPORT_ERROR("Error in askMemory: Memory allocation size requested is zero!");
        return(NULL);
    }
	
    if ( ( ptr = (char *) calloc(1,memsize*sizeof(char)) ) == NULL ) {
        std::cerr<<"Memory allocation of %ld bytes failed, memsize= "<< memsize<<std::endl;
        REPORT_ERROR("Error in askMemory");
        return(NULL); 
    }
	
    //memset(ptr, 0, memsize); 	 

    return(ptr); 
}

int  freeMemory(void* ptr, unsigned long memsize)
{
    if ( ptr == NULL ) 
        return(0);
	
    if ( memsize < 1 ) 
    {
        return(-1);
    }

    free(ptr);
    ptr = NULL;
    return(0);
}

