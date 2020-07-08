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
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

#include <stdio.h>

#include "src/symmetries.h"

// Read Symmetry file ======================================================
// crystal symmetry matices from http://cci.lbl.gov/asu_gallery/
int SymList::read_sym_file(FileName fn_sym)
{
    int i, j;
    FILE *fpoii;
    char line[80];
    char *auxstr;
    RFLOAT ang_incr, rot_ang;
    int  fold;
    Matrix2D<RFLOAT> L(4, 4), R(4, 4);
    Matrix1D<RFLOAT> axis(3);
    int pgGroup = 0, pgOrder = 0;
    std::vector<std::string> fileContent;

    //check if reserved word

    // Open file ---------------------------------------------------------
    if ((fpoii = fopen(fn_sym.c_str(), "r")) == NULL)
    {
        //check if reserved word and return group and order
        if (isSymmetryGroup(fn_sym, pgGroup, pgOrder))
        {
        	fill_symmetry_class(fn_sym, pgGroup, pgOrder, fileContent);
        }
        else
            REPORT_ERROR((std::string)"SymList::read_sym_file:Can't open file: "
                     + " or do not recognize symmetry group" + fn_sym);
    }
    else
    {
        while (fgets(line, 79, fpoii) != NULL)
        {
            if (line[0] == ';' || line[0] == '#' || line[0] == '\0')
            	continue;
			fileContent.push_back(line);
        }
        fclose(fpoii);
    }

    // Count the number of symmetries ------------------------------------
    true_symNo = 0;
    // count number of axis and mirror planes. It will help to identify
    // the crystallographic symmetry

    int no_axis, no_mirror_planes, no_inversion_points;
    no_axis = no_mirror_planes = no_inversion_points = 0;

    for (int n=0; n<fileContent.size(); n++)
    {
    	strcpy(line,fileContent[n].c_str());
        auxstr = firstToken(line);
        if (auxstr == NULL)
        {
            std::cout << line;
            std::cout << "Wrong line in symmetry file, the line is skipped\n";
            continue;
        }
        if (strcmp(auxstr, "rot_axis") == 0)
        {
            auxstr = nextToken();
            fold = textToInteger(auxstr);
            true_symNo += (fold - 1);
            no_axis++;
        }
        else if (strcmp(auxstr, "mirror_plane") == 0)
        {
            true_symNo++;
            no_mirror_planes++;
        }
        else if (strcmp(auxstr, "inversion") == 0)
        {
            true_symNo += 1;
            no_inversion_points = 1;
        }
    }
    // Ask for memory
    __L.resize(4*true_symNo, 4);
    __R.resize(4*true_symNo, 4);
    __chain_length.resize(true_symNo);
    __chain_length.initConstant(1);

    // Read symmetry parameters
    i = 0;
    for (int n=0; n<fileContent.size(); n++)
    {
        strcpy(line,fileContent[n].c_str());
        auxstr = firstToken(line);
        // Rotational axis ---------------------------------------------------
        if (strcmp(auxstr, "rot_axis") == 0)
        {
            auxstr = nextToken();
            fold = textToInteger(auxstr);
            auxstr = nextToken();
            XX(axis) = textToDouble(auxstr);
            auxstr = nextToken();
            YY(axis) = textToDouble(auxstr);
            auxstr = nextToken();
            ZZ(axis) = textToDouble(auxstr);
            ang_incr = 360. / fold;
            L.initIdentity();
            for (j = 1, rot_ang = ang_incr; j < fold; j++, rot_ang += ang_incr)
            {
                rotation3DMatrix(rot_ang, axis, R);
                R.setSmallValuesToZero();
                set_matrices(i++, L, R.transpose());
            }
            __sym_elements++;
            // inversion ------------------------------------------------------
        }
        else if (strcmp(auxstr, "inversion") == 0)
        {
            L.initIdentity();
            L(2, 2) = -1;
            R.initIdentity();
            R(0, 0) = -1.;
            R(1, 1) = -1.;
            R(2, 2) = -1.;
            set_matrices(i++, L, R);
            __sym_elements++;
            // mirror plane -------------------------------------------------------------
        }
        else if (strcmp(auxstr, "mirror_plane") == 0)
        {
            auxstr = nextToken();
            XX(axis) = textToFloat(auxstr);
            auxstr = nextToken();
            YY(axis) = textToFloat(auxstr);
            auxstr = nextToken();
            ZZ(axis) = textToFloat(auxstr);
            L.initIdentity();
            L(2, 2) = -1;
            Matrix2D<RFLOAT> A;
            alignWithZ(axis,A);
            A = A.transpose();
            R = A * L * A.inv();
            L.initIdentity();
            set_matrices(i++, L, R);
            __sym_elements++;
        }
    }

    compute_subgroup();

    return pgGroup;
}

// Get matrix ==============================================================
void SymList::get_matrices(int i, Matrix2D<RFLOAT> &L, Matrix2D<RFLOAT> &R)
const
{
    int k, l;
    L.initZeros(4, 4);
    R.initZeros(4, 4);
    for (k = 4 * i; k < 4*i + 4; k++)
        for (l = 0; l < 4; l++)
        {
            L(k - 4*i, l) = __L(k, l);
            R(k - 4*i, l) = __R(k, l);
        }
}

// Set matrix ==============================================================
void SymList::set_matrices(int i, const Matrix2D<RFLOAT> &L,
                           const Matrix2D<RFLOAT> &R)
{
    int k, l;
    for (k = 4 * i; k < 4*i + 4; k++)
        for (l = 0; l < 4; l++)
        {
            __L(k, l) = L(k - 4 * i, l);
            __R(k, l) = R(k - 4 * i, l);
        }
}

// Add matrix ==============================================================
void SymList::add_matrices(const Matrix2D<RFLOAT> &L, const Matrix2D<RFLOAT> &R,
                           int chain_length)
{
    if (MAT_XSIZE(L) != 4 || MAT_YSIZE(L) != 4 || MAT_XSIZE(R) != 4 || MAT_YSIZE(R) != 4)
        REPORT_ERROR( "SymList::add_matrix: Transformation matrix is not 4x4");
    if (TrueSymsNo() == SymsNo())
    {
        __L.resize(MAT_YSIZE(__L) + 4, 4);
        __R.resize(MAT_YSIZE(__R) + 4, 4);
        __chain_length.resize(__chain_length.size() + 1);
    }

    set_matrices(true_symNo, L, R);
    __chain_length(__chain_length.size() - 1) = chain_length;
    true_symNo++;
}

// Compute subgroup ========================================================
bool found_not_tried(const Matrix2D<int> &tried, int &i, int &j,
                     int true_symNo)
{
    i = j = 0;
    int n = 0;
    while (n != MAT_YSIZE(tried))
    {
        if (tried(i, j) == 0 && !(i >= true_symNo && j >= true_symNo))
            return true;
        if (i != n)
        {
            // Move downwards
            i++;
        }
        else
        {
            // Move leftwards
            j--;
            if (j == -1)
            {
                n++;
                j = n;
                i = 0;
            }
        }
    }
    return false;
}

//#define DEBUG
void SymList::compute_subgroup()
{
    Matrix2D<RFLOAT> I(4, 4);
    I.initIdentity();
    Matrix2D<RFLOAT> L1(4, 4), R1(4, 4), L2(4, 4), R2(4, 4), newL(4, 4), newR(4, 4);
    Matrix2D<int>    tried(true_symNo, true_symNo);
    int i, j;
    int new_chain_length;
    while (found_not_tried(tried, i, j, true_symNo))
    {
        tried(i, j) = 1;

        get_matrices(i, L1, R1);
        get_matrices(j, L2, R2);
        newL = L1 * L2;
        newR = R1 * R2;
        new_chain_length = __chain_length(i) + __chain_length(j);
        Matrix2D<RFLOAT> newR3 = newR;
        newR3.resize(3,3);
        if (newL.isIdentity() && newR3.isIdentity()) continue;

        // Try to find it in current ones
        bool found;
        found = false;
        for (int l = 0; l < SymsNo(); l++)
        {
        	get_matrices(l, L1, R1);
            if (newL.equal(L1) && newR.equal(R1))
            {
                found = true;
                break;
            }
        }

        if (!found)
        {
//#define DEBUG
#ifdef DEBUG
           std::cout << "Matrix size " << tried.Xdim() << " "
            << "trying " << i << " " << j << " "
            << "chain length=" << new_chain_length << std::endl;
            std::cout << "Result R Sh\n" << newR;
#endif
//#undef DEBUG
            newR.setSmallValuesToZero();
            newL.setSmallValuesToZero();
            add_matrices(newL, newR, new_chain_length);
            tried.resize(MAT_YSIZE(tried) + 1, MAT_XSIZE(tried) + 1);
        }
    }
}

/** translate string fn_sym to symmetry group, return false
    is translation is not possible. See URL
    http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry
    for details  */
bool SymList::isSymmetryGroup(FileName fn_sym, int &pgGroup, int &pgOrder)
{
   char G1,G2,G3,G4;
   char auxChar[3];
   //each case check lenght, check first letter, second, is number
   //Non a point group

   //remove path
   FileName fn_sym_tmp;
   fn_sym_tmp=fn_sym.removeDirectories();
   int mySize=fn_sym_tmp.size();
   bool return_true;
   return_true=false;
   auxChar[2]='\0';
   //size maybe 4 because n maybe a 2 digit number
   if(mySize>4 || mySize<1)
   {
      pgGroup=-1;
      pgOrder=-1;
      return false;
   }
   //get the group character by character
   G1=toupper((fn_sym_tmp.c_str())[0]);
   G2=toupper((fn_sym_tmp.c_str())[1]);
   if (mySize > 2)
   {   G3=toupper((fn_sym_tmp.c_str())[2]);
       if(mySize > 3)
           G4=toupper((fn_sym.c_str())[3]);
   }
   else
       G4='\0';
   //CN
   if (mySize==2 && G1=='C' && isdigit(G2))
   {
       pgGroup=pg_CN;
       pgOrder=int(G2)-48;
       return_true=true;
   }
   if (mySize==3 && G1=='C' && isdigit(G2) && isdigit(G3))
   {
       pgGroup=pg_CN;
       auxChar[0]=G2;
       auxChar[1]=G3;
       pgOrder=atoi(auxChar);
       return_true=true;
   }
   //CI
   else if (mySize==2 && G1=='C' && G2=='I')
   {
       pgGroup=pg_CI;
       pgOrder=-1;
       return_true=true;
   }
   //CS
   else if (mySize==2 && G1=='C' && G2=='S')
   {
       pgGroup=pg_CS;
       pgOrder=-1;
       return_true=true;
   }
   //CNH
   else if (mySize==3 && G1=='C' && isdigit(G2) && G3=='H')
   {
       pgGroup=pg_CNH;
       pgOrder=int(G2)-48;
       return_true=true;
   }
   else if (mySize==4 && G1=='C' && isdigit(G2) && isdigit(G3) && G4=='H')
   {
       pgGroup=pg_CNH;
       auxChar[0]=G2;
       auxChar[1]=G3;
       pgOrder=atoi(auxChar);
       return_true=true;
   }
   //CNV
   else if (mySize==3 && G1=='C' && isdigit(G2) && G3=='V')
   {
       pgGroup=pg_CNV;
       pgOrder=int(G2)-48;
       return_true=true;
   }
   else if (mySize==4 && G1=='C' && isdigit(G2) && isdigit(G3) && G4=='V')
   {
       pgGroup=pg_CNV;
       auxChar[0]=G2;
       auxChar[1]=G3;
       pgOrder=atoi(auxChar);
       return_true=true;
   }
   //SN
   else if (mySize==2 && G1=='S' && isdigit(G2) )
   {
       pgGroup=pg_SN;
       pgOrder=int(G2)-48;
       return_true=true;
   }
   else if (mySize==3 && G1=='S' && isdigit(G2) && isdigit(G3) )
   {
       pgGroup=pg_SN;
       auxChar[0]=G2;
       auxChar[1]=G3;
       pgOrder=atoi(auxChar);
       return_true=true;
   }
   //DN
   else if (mySize==2 && G1=='D' && isdigit(G2) )
   {
       pgGroup=pg_DN;
       pgOrder=int(G2)-48;
       return_true=true;
   }
   if (mySize==3 && G1=='D' && isdigit(G2) && isdigit(G3))
   {
       pgGroup=pg_DN;
       auxChar[0]=G2;
       auxChar[1]=G3;
       pgOrder=atoi(auxChar);
       return_true=true;
   }
   //DNV
   else if (mySize==3 && G1=='D' && isdigit(G2) && G3=='V')
   {
       pgGroup=pg_DNV;
       pgOrder=int(G2)-48;
       return_true=true;
   }
   else if (mySize==4 && G1=='D' && isdigit(G2) && isdigit(G3) && G4=='V')
   {
       pgGroup=pg_DNV;
       auxChar[0]=G2;
       auxChar[1]=G3;
       pgOrder=atoi(auxChar);
       return_true=true;
   }
   //DNH
   else if (mySize==3 && G1=='D' && isdigit(G2) && G3=='H')
   {
       pgGroup=pg_DNH;
       pgOrder=int(G2)-48;
       return_true=true;
   }
   else if (mySize==4 && G1=='D' && isdigit(G2) && isdigit(G3) && G4=='H')
   {
       pgGroup=pg_DNH;
       auxChar[0]=G2;
       auxChar[1]=G3;
       pgOrder=atoi(auxChar);
       return_true=true;
   }
   //T
   else if (mySize==1 && G1=='T')
   {
       pgGroup=pg_T;
       pgOrder=-1;
       return_true=true;
   }
   //TD
   else if (mySize==2 && G1=='T' && G2=='D')
   {
       pgGroup=pg_TD;
       pgOrder=-1;
       return_true=true;
   }
   //TH
   else if (mySize==2 && G1=='T' && G2=='H')
   {
       pgGroup=pg_TH;
       pgOrder=-1;
       return_true=true;
   }
   //O
   else if (mySize==1 && G1=='O')
   {
       pgGroup=pg_O;
       pgOrder=-1;
       return_true=true;
   }
   //OH
   else if (mySize==2 && G1=='O'&& G2=='H')
   {
       pgGroup=pg_OH;
       pgOrder=-1;
       return_true=true;
   }
   //I
   else if (mySize==1 && G1=='I')
   {
       pgGroup=pg_I;
       pgOrder=-1;
       return_true=true;
   }
   //I1
   else if (mySize==2 && G1=='I'&& G2=='1')
   {
       pgGroup=pg_I1;
       pgOrder=-1;
       return_true=true;
   }
   //I2
   else if (mySize==2 && G1=='I'&& G2=='2')
   {
       pgGroup=pg_I2;
       pgOrder=-1;
       return_true=true;
   }
   //I3
   else if (mySize==2 && G1=='I'&& G2=='3')
   {
       pgGroup=pg_I3;
       pgOrder=-1;
       return_true=true;
   }
   //I4
   else if (mySize==2 && G1=='I'&& G2=='4')
   {
       pgGroup=pg_I4;
       pgOrder=-1;
       return_true=true;
   }
   //I5
   else if (mySize==2 && G1=='I'&& G2=='5')
   {
       pgGroup=pg_I5;
       pgOrder=-1;
       return_true=true;
   }
   //IH
   else if (mySize==2 && G1=='I'&& G2=='H')
   {
       pgGroup=pg_IH;
       pgOrder=-1;
       return_true=true;
   }
   //I1H
   else if (mySize==3 && G1=='I'&& G2=='1'&& G3=='H')
   {
       pgGroup=pg_I1H;
       pgOrder=-1;
       return_true=true;
   }
   //I2H
   else if (mySize==3 && G1=='I'&& G2=='2'&& G3=='H')
   {
       pgGroup=pg_I2H;
       pgOrder=-1;
       return_true=true;
   }
   //I3H
   else if (mySize==3 && G1=='I'&& G2=='3'&& G3=='H')
   {
       pgGroup=pg_I3H;
       pgOrder=-1;
       return_true=true;
   }
   //I4H
   else if (mySize==3 && G1=='I'&& G2=='4'&& G3=='H')
   {
       pgGroup=pg_I4H;
       pgOrder=-1;
       return_true=true;
   }
   //I5H
   else if (mySize==3 && G1=='I'&& G2=='5'&& G3=='H')
   {
       pgGroup=pg_I5H;
       pgOrder=-1;
       return_true=true;
   }
//#define DEBUG7
#ifdef DEBUG7
std::cerr << "pgGroup" << pgGroup << " pgOrder " << pgOrder << std::endl;
#endif
#undef DEBUG7

   return return_true;
}
void SymList::fill_symmetry_class(const FileName symmetry, int pgGroup, int pgOrder,
   std::vector<std::string> &fileContent)
{

	fileContent.clear();
	if (pgGroup == pg_CN)
    {
    	fileContent.push_back("rot_axis " + integerToString(pgOrder) + " 0 0 1");
    }
    else if (pgGroup == pg_CI)
    {
    	fileContent.push_back("inversion ");
    }
    else if (pgGroup == pg_CS)
    {
    	fileContent.push_back("mirror_plane 0 0 1");
    }
    else if (pgGroup == pg_CNV)
    {
    	fileContent.push_back("rot_axis " + integerToString(pgOrder) + " 0 0 1");
    	fileContent.push_back("mirror_plane 0 1 0");
    }
    else if (pgGroup == pg_CNH)
    {
    	fileContent.push_back("rot_axis " + integerToString(pgOrder) + " 0 0 1");
    	fileContent.push_back("mirror_plane 0 0 1");
    }
    else if (pgGroup == pg_SN)
    {
        int order = pgOrder / 2;
		if(2*order != pgOrder)
		{
				std::cerr << "ERROR: order for SN group must be even" << std::endl;
				exit(0);
		}
        fileContent.push_back("rot_axis " + integerToString(order) + " 0 0 1");
        fileContent.push_back("inversion ");
    }
    else if (pgGroup == pg_DN)
    {
    	if (pgOrder > 1)
    		fileContent.push_back("rot_axis " + integerToString(pgOrder) + " 0 0 1");
    	fileContent.push_back("rot_axis 2 1 0 0");
    }
    else if (pgGroup == pg_DNV)
    {
    	if (pgOrder > 1)
    		fileContent.push_back("rot_axis " + integerToString(pgOrder) + " 0 0 1");
        fileContent.push_back("rot_axis 2 1 0 0");
        fileContent.push_back("mirror_plane 1 0 0");
    }
    else if (pgGroup == pg_DNH)
    {
    	if (pgOrder > 1)
    		fileContent.push_back("rot_axis " + integerToString(pgOrder) + " 0 0 1");
        fileContent.push_back("rot_axis 2 1 0 0");
        fileContent.push_back("mirror_plane 0 0 1");
    }
    else if (pgGroup == pg_T)
    {
        fileContent.push_back("rot_axis 3  0. 0. 1.");
        fileContent.push_back("rot_axis 2 0. 0.816496 0.577350");
    }
    else if (pgGroup == pg_TD)
    {
        fileContent.push_back("rot_axis 3  0. 0. 1.");
        fileContent.push_back("rot_axis 2 0. 0.816496 0.577350");
        fileContent.push_back("mirror_plane 1.4142136 2.4494897 0.0000000");
    }
    else if (pgGroup == pg_TH)
    {
        fileContent.push_back("rot_axis 3  0. 0. 1.");
        fileContent.push_back("rot_axis 2 0. -0.816496 -0.577350");
        fileContent.push_back("inversion");
    }
    else if (pgGroup == pg_O)
    {
        fileContent.push_back("rot_axis 3  .5773502  .5773502 .5773502");
        fileContent.push_back("rot_axis 4 0 0 1");
    }
    else if (pgGroup == pg_OH)
    {
        fileContent.push_back("rot_axis 3  .5773502  .5773502 .5773502");
        fileContent.push_back("rot_axis 4 0 0 1");
        fileContent.push_back("mirror_plane 0 1 1");
    }
    else if (pgGroup == pg_I || pgGroup == pg_I2)
    {
        fileContent.push_back("rot_axis 2  0 0 1");
        fileContent.push_back("rot_axis 5  0.525731114  0 0.850650807");
        fileContent.push_back("rot_axis 3  0 0.356822076 0.934172364");
    }
    else if (pgGroup == pg_I1)
    {
        fileContent.push_back("rot_axis 2  1  	   0	       0");
        fileContent.push_back("rot_axis 5 0.85065080702670 0 -0.5257311142635");
        fileContent.push_back("rot_axis 3 0.9341723640 0.3568220765 0");
    }
    else if (pgGroup == pg_I3)
    {
        fileContent.push_back("rot_axis 2  -0.5257311143 0 0.8506508070");
        fileContent.push_back("rot_axis 5  0. 0. 1.");
        fileContent.push_back("rot_axis 3  -0.4911234778630044, 0.3568220764705179, 0.7946544753759428");
    }
    else if (pgGroup == pg_I4)
    {
        fileContent.push_back("rot_axis 2  0.5257311143 0 0.8506508070");
        fileContent.push_back("rot_axis 5  0.8944271932547096 0 0.4472135909903704");
        fileContent.push_back("rot_axis 3  0.4911234778630044 0.3568220764705179 0.7946544753759428");
    }
    else if (pgGroup == pg_I5)
    {
        std::cerr << "ERROR: Symmetry pg_I5 not implemented" << std::endl;
        exit(0);
    }
    else if (pgGroup == pg_IH || pgGroup == pg_I2H)
    {
        fileContent.push_back("rot_axis 2  0 0 1");
        fileContent.push_back("rot_axis 5  0.525731114  0 0.850650807");
        fileContent.push_back("rot_axis 3  0 0.356822076 0.934172364");
        fileContent.push_back("mirror_plane 1 0 0");
    }
    else if (pgGroup == pg_I1H)
    {
        fileContent.push_back("rot_axis 2  1  	   0	       0");
        fileContent.push_back("rot_axis 5 0.85065080702670 0 -0.5257311142635");
        fileContent.push_back("rot_axis 3 0.9341723640 0.3568220765 0");
        fileContent.push_back("mirror_plane 0 0 -1");
    }
    else if (pgGroup == pg_I3H)
    {
        fileContent.push_back("rot_axis 2  -0.5257311143 0 0.8506508070");
        fileContent.push_back("rot_axis 5  0. 0. 1.");
        fileContent.push_back("rot_axis 3  -0.4911234778630044, 0.3568220764705179, 0.7946544753759428");
        fileContent.push_back("mirror_plane 0.850650807 0  0.525731114");
    }
    else if (pgGroup == pg_I4H)
    {
        fileContent.push_back("rot_axis 2  0.5257311143 0 0.8506508070");
        fileContent.push_back("rot_axis 5  0.8944271932547096 0 0.4472135909903704");
        fileContent.push_back("rot_axis 3  0.4911234778630044 0.3568220764705179 0.7946544753759428");
        fileContent.push_back("mirror_plane 0.850650807 0 -0.525731114");
    }
    else if (pgGroup == pg_I5H)
    {
        std::cerr << "ERROR: Symmetry pg_I5H not implemented" << std::endl;
        exit(0);
    }
    else
    {
        std::cerr << "ERROR: Symmetry " << symmetry  << "is not known" << std::endl;
        exit(0);
    }

//#define DEBUG5
#ifdef DEBUG5
    for (int n=0; n<fileContent.size(); n++)
    	std::cerr << fileContent[n] << std::endl;
	std::cerr << "fileContent.size()" << fileContent.size() << std::endl;
#endif
#undef DEBUG5
}

void SymList::writeDefinition(std::ostream &outstream, FileName fn_sym)
{
	read_sym_file(fn_sym);
	Matrix2D<RFLOAT> L(3,3), R(3,3);
	outstream << " ++++ Using symmetry group " << fn_sym << ", with the following " << SymsNo()+1 << " transformation matrices:"<< std::endl;
    R.initIdentity();
    outstream << " R(1)= " << R;
    for (int isym = 0; isym < SymsNo(); isym++)
    {
        get_matrices(isym, L, R);
        R.resize(3, 3);
        L.resize(3, 3);
        if (!L.isIdentity())
        	outstream << " L("<< isym+2<<")= "<<L;
        outstream << " R("<< isym+2<<")= "<<R;
        RFLOAT alpha, beta, gamma;
        Euler_matrix2angles(R, alpha, beta, gamma);
        outstream << "     Euler angles: " << alpha << " " << beta << " " << gamma << std::endl;
    }

}

RFLOAT SymList::non_redundant_ewald_sphere(int pgGroup, int pgOrder)
{
    if (pgGroup == pg_CN)
    {
        return 4.*PI/pgOrder;
    }
    else if (pgGroup == pg_CI)
    {
        return 4.*PI/2.;
    }
    else if (pgGroup == pg_CS)
    {
        return 4.*PI/2.;
    }
    else if (pgGroup == pg_CNV)
    {
        return 4.*PI/pgOrder/2;
    }
    else if (pgGroup == pg_CNH)
    {
        return 4.*PI/pgOrder/2;
    }
    else if (pgGroup == pg_SN)
    {
        return 4.*PI/pgOrder;
    }
    else if (pgGroup == pg_DN)
    {
        return 4.*PI/pgOrder/2;
    }
    else if (pgGroup == pg_DNV)
    {
        return 4.*PI/pgOrder/4;
    }
    else if (pgGroup == pg_DNH)
    {
        return 4.*PI/pgOrder/4;
    }
    else if (pgGroup == pg_T)
    {
        return 4.*PI/12;
    }
    else if (pgGroup == pg_TD)
    {
        return 4.*PI/24;
    }
    else if (pgGroup == pg_TH)
    {
        return 4.*PI/24;
    }
    else if (pgGroup == pg_O)
    {
        return 4.*PI/24;
    }
    else if (pgGroup == pg_OH)
    {
        return 4.*PI/48;
    }
    else if (pgGroup == pg_I || pgGroup == pg_I2)
    {
        return 4.*PI/60;
    }
    else if (pgGroup == pg_I1)
    {
        return 4.*PI/60;
    }
    else if (pgGroup == pg_I3)
    {
        return 4.*PI/60;
    }
    else if (pgGroup == pg_I4)
    {
        return 4.*PI/60;
    }
    else if (pgGroup == pg_I5)
    {
        return 4.*PI/60;
    }
    else if (pgGroup == pg_IH || pgGroup == pg_I2H)
    {
        return 4.*PI/120;
    }
    else if (pgGroup == pg_I1H)
    {
        return 4.*PI/120;
    }
    else if (pgGroup == pg_I3H)
    {
        return 4.*PI/120;
    }
    else if (pgGroup == pg_I4H)
    {
        return 4.*PI/120;
    }
    else if (pgGroup == pg_I5H)
    {
        return 4.*PI/120;
    }
    else
    {
        std::cerr << "ERROR: Symmetry group, order=" << pgGroup
                                                     << " "
                                                     <<  pgOrder
                                                     << "is not known"
                                                     << std::endl;
        exit(0);
    }
}

void symmetriseMap(MultidimArray<RFLOAT> &img, FileName &fn_sym, bool do_wrap)
{

	if (img.getDim() != 3)
		REPORT_ERROR("symmetriseMap ERROR: symmetriseMap can only be run on 3D maps!");

	img.setXmippOrigin();

	SymList SL;
	SL.read_sym_file(fn_sym);

	Matrix2D<RFLOAT> L(4, 4), R(4, 4); // A matrix from the list
    MultidimArray<RFLOAT> sum, aux;
    sum = img;
    aux.resize(img);

	for (int isym = 0; isym < SL.SymsNo(); isym++)
    {
        SL.get_matrices(isym, L, R);
        applyGeometry(img, aux, R, IS_INV, do_wrap);
        sum += aux;
    }

	// Overwrite the input
	img = sum / (SL.SymsNo() + 1);

}
