//
//  CPlot2D.cpp
//
//  A simple class for plotting xy values and
//  generating a PostScript output, aimed at
//  easy integration into C/C++ programs.
//
//  You are granted use of the code, but please
//  be a nice guy and acknowledge where you got
//  it from by giving credit!
//
//  If you have any comments, questions, bug
//  reports, you can contact me at
//  attila AT amzsaki.com
//
//  Created by Attila Michael Zsaki on 14-03-17.
//  Copyright (c) 2014 AMZs. All rights reserved.
//

#include <math.h>
#include <iostream>
#include <float.h>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "CPlot2D.h"

void joinMultipleEPSIntoSinglePDF(FileName fn_pdf, std::vector<FileName> fn_eps)
{

    FileName fn_list = fn_pdf + ".lst";
    std::string command = "gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -dDEVICEWIDTHPOINTS=800 -dDEVICEHEIGHTPOINTS=800 -sOutputFile=";
    command += fn_pdf + " @" + fn_list;
    std::ofstream filelist(fn_pdf + ".lst");
    bool have_at_least_one = false;
    for (int i = 0; i < fn_eps.size(); i++)
    {
        // fn_eps[i] could be a Linux wildcard...
    	std::vector<FileName> all_eps_files;
        fn_eps[i].globFiles(all_eps_files);
        for (long int j= 0; j < all_eps_files.size(); j++)
        {
        	if (exists(all_eps_files[j]))
        	{
        		filelist << all_eps_files[j] << "\n";
        		have_at_least_one = true;
        	}
        }
    }
    filelist.close();

    bool have_error_in_gs = false;
    if (have_at_least_one)
    {
		command += " > /dev/null";

		if (system(command.c_str()))
		{
			std::cerr << " ERROR in executing: " << command << "\n";
			have_error_in_gs = true;
		}
    }
    else
    {
    	std::cerr << " Did not find any of the expected EPS files to generate a PDF file" << "\n";
    }

    // std::remove(fn_list.c_str()); // don't know why but Ghostscript fails with this line.
    // system() should wait the termination of the program, so this is very strange...
    if (!have_at_least_one || have_error_in_gs)
    {
    	std::cerr << " + Will make an empty PDF-file in " << fn_pdf << "\n";
    	touch(fn_pdf);
    }

}
bool concatenatePDFfiles(FileName fn_pdf_out, FileName pdf1, FileName pdf2)
{
	std::vector<FileName> fn_pdfs;
	fn_pdfs.push_back(pdf1);
	fn_pdfs.push_back(pdf2);
	return concatenatePDFfiles(fn_pdf_out, fn_pdfs);

}

bool concatenatePDFfiles(FileName fn_pdf_out, std::vector<FileName> fn_pdfs)
{

	FileName fn_comb=fn_pdf_out;
	// check if fn_pdf_out occurs in fn_pdfs
	if (std::find(fn_pdfs.begin(), fn_pdfs.end(), fn_pdf_out) != fn_pdfs.end())
	{
	  // Element in vector.
		fn_comb.withoutExtension() + "_tmp_combine.pdf";
	}

	std::string command="gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=";
	command += fn_comb;
	command += " -dBATCH ";
	for (int i = 0; i < fn_pdfs.size(); i++)
		command += fn_pdfs[i] + " ";

	command += " > /dev/null";

	if (system(command.c_str()))
	{
		std::cerr << " ERROR in executing: " << command << "\n";
		return false;
	}

	if (fn_comb != fn_pdf_out)
	{
		// First remove fn_pdf_out
		if (std::remove(fn_pdf_out.c_str()))
		{
			std::cerr << "ERROR in removing pre-existing " << fn_pdf_out << ".\n";
			return false;
		}
		// Then move fn_comb to fn_pdf_out
		if (std::rename(fn_comb.c_str(), fn_pdf_out.c_str()))
		{
			std::cerr << "ERROR in renaming " << fn_comb << " to " << fn_pdf_out << ".\n";
			return false;
		}
	}

	return true;
}

CPlot2D::CPlot2D(std::string title)
{
    //m_dXAxisSize=809.0; // Golden Ratio
    m_dXAxisSize=600.0; // Golden Ratio
    m_dYAxisSize=500.0;

    m_dBottomFrameSize=90.0;
    m_dRightFrameSize=40.0;
    m_dTopFrameSize=40.0;
    m_dLeftFrameSize=75.0;

    m_dXTotalSize=m_dLeftFrameSize+m_dXAxisSize+m_dRightFrameSize;
    m_dYTotalSize=m_dBottomFrameSize+m_dYAxisSize+m_dTopFrameSize;
    m_dDiagonalSize=sqrt(m_dXAxisSize*m_dXAxisSize+m_dYAxisSize*m_dYAxisSize);

    m_dFrameLineWidth=1.0;
    m_dFrameColor[0]=0.0;
    m_dFrameColor[1]=0.0;
    m_dFrameColor[2]=0.0;

    m_dGridLineWidth=1.0;
    m_dGridColor[0]=0.8;
    m_dGridColor[1]=0.8;
    m_dGridColor[2]=0.8;

    m_dLineDotSpacing=0.0;

    m_bDrawXAxisTickMarks=true;
    m_bDrawYAxisTickMarks=true;
    m_dXAxisNumbersSpacing=1.0;
    m_dYAxisNumbersSpacing=1.0;
    m_iXAxisNumberOfTicks=10;
    m_iYAxisNumberOfTicks=10;
    m_bDrawXAxisGridLines=true;
    m_bDrawYAxisGridLines=true;
    m_bDrawGridLinesDashed=true;

    m_bFlipY = false;
    m_dFlipYOffset = 0;

    // Sjors Scheres 22mar2016: changed all fonts to Times
    m_strXAxisLabelFont="Times";
    m_dXAxisLabelFontSize=12.0;
    m_strYAxisLabelFont="Times";
    m_dYAxisLabelFontSize=12.0;

    m_strXAxisTitleFont="Times";
    m_dXAxisTitleFontSize=16.0;
    m_strYAxisTitleFont="Times";
    m_dYAxisTitleFontSize=16.0;
    m_dXAxisTitleColor[0]=0.0;
    m_dXAxisTitleColor[1]=0.0;
    m_dXAxisTitleColor[2]=0.0;
    m_dYAxisTitleColor[0]=0.0;
    m_dYAxisTitleColor[1]=0.0;
    m_dYAxisTitleColor[2]=0.0;
    m_strXAxisTitle="";
    m_strYAxisTitle="";
    // Sjors Scheres 22mar2016: insert title of plot here
    m_strPlotTitle=title;

    m_dXAxisLabelColor[0]=0.0;
    m_dXAxisLabelColor[1]=0.0;
    m_dXAxisLabelColor[2]=0.0;
    m_dYAxisLabelColor[0]=0.0;
    m_dYAxisLabelColor[1]=0.0;
    m_dYAxisLabelColor[2]=0.0;

    m_bDrawLegend=true;

	m_bSizeSetExternally = false;
}

CPlot2D::~CPlot2D()
{

}

void CPlot2D::OutputPostScriptPlot(std::string fileName)
{

    outputFile.open(fileName.c_str());

    // precompute plot dimensions
    PrecomputeDimensions();

    // header
    outputFile << "%!PS-Adobe-2.0 EPSF-1.2" << "\n";
    outputFile << "%%BoundingBox: 0 0 " << (int)m_dXTotalSize << " " <<  (int)m_dYTotalSize << "\n";
    outputFile << "%%Pages: 1" << "\n";
    outputFile << "%%EndComments" << "\n";

    // grid lines
    if (m_bDrawXAxisGridLines) {
        DrawXAxisGridLinesPostScript();
    }
    if (m_bDrawYAxisGridLines) {
        DrawYAxisGridLinesPostScript();
    }

    // draw plot frame
    DrawFramePostScript();

    // draw axis tick marks
    if (m_bDrawXAxisTickMarks) {
        DrawXAxisTickMarksPostScript();
    }
    if (m_bDrawYAxisTickMarks) {
        DrawYAxisTickMarksPostScript();
    }

    // draw axis labels
    // might be separate from drawing tick marks one day...
    if (m_bDrawXAxisTickMarks) {
        DrawXAxisLabelsPostScript();
    }
    if (m_bDrawYAxisTickMarks) {
        DrawYAxisLabelsPostScript();
    }

    // draw axis titles
    DrawXAxisTitlePostScript();
    DrawYAxisTitlePostScript();


    // draw data
    DrawDataPostScript();

    // draw legend
    if (m_bDrawLegend) {
        DrawLegendPostScript();
    }


    outputFile.close();

}

void CPlot2D::DrawFramePostScript()
{
    outputFile << "newpath" << "\n";
    outputFile << m_dLeftFrameSize << " " << m_dBottomFrameSize << " moveto" << "\n";
    outputFile << m_dXAxisSize << " " << 0 << " rlineto" << "\n";
    outputFile << 0 << " " << m_dYAxisSize << " rlineto" << "\n";
    outputFile << -m_dXAxisSize << " " << 0 << " rlineto" << "\n";
    outputFile << "closepath" << "\n";

    outputFile << m_dFrameLineWidth << " setlinewidth" << "\n";
    outputFile << m_dFrameColor[0] << " " << m_dFrameColor[1] << " " << m_dFrameColor[2] << " setrgbcolor" << "\n";
    outputFile << "stroke" << "\n";

    double labelXCoordinate,labelYCoordinate;

    labelXCoordinate=m_dLeftFrameSize+m_dXAxisSize*0.5;
    labelYCoordinate=m_dBottomFrameSize + m_dYAxisSize + 10;

    outputFile << "/" << m_strXAxisTitleFont << " findfont" << "\n";
    outputFile << m_dYAxisTitleFontSize << " scalefont" << "\n";
    outputFile << "setfont" << "\n";
    outputFile << labelXCoordinate << " " << labelYCoordinate << " moveto" << "\n";
    outputFile << m_dXAxisTitleColor[0] << " " << m_dXAxisTitleColor[1] << " " << m_dXAxisTitleColor[2] << " setrgbcolor" << "\n";

    // let PostScript handle the final adjustment based on the width of the string
    outputFile << "(" << m_strPlotTitle << ")" << " dup stringwidth pop 2 div neg 0 rmoveto show"  << "\n";

}

void CPlot2D::PrecomputeDimensions()
{

    m_dXTotalSize=m_dXAxisSize+m_dRightFrameSize+m_dLeftFrameSize;
    m_dYTotalSize=m_dYAxisSize+m_dTopFrameSize+m_dBottomFrameSize;

    m_dMaxXExtent=0.0;
    m_dMaxYExtent=0.0;
    m_dMinXStartPoint=DBL_MAX;
    m_dMinYStartPoint=DBL_MAX;
    m_dMaxXEndPoint=-DBL_MAX;
    m_dMaxYEndPoint=-DBL_MAX;

	if (!m_bSizeSetExternally) {
		// for all data sets
		for (int i=0;i<m_dataSets.size();++i) {
			if (m_dMinXStartPoint>m_dataSets[i].GetXMinValue()) {
				m_dMinXStartPoint=m_dataSets[i].GetXMinValue();
			}
			if (m_dMinYStartPoint>m_dataSets[i].GetYMinValue()) {
				m_dMinYStartPoint=m_dataSets[i].GetYMinValue();
			}
			if (m_dMaxXEndPoint<m_dataSets[i].GetXMaxValue()) {
				m_dMaxXEndPoint=m_dataSets[i].GetXMaxValue();
			}
			if (m_dMaxYEndPoint<m_dataSets[i].GetYMaxValue()) {
				m_dMaxYEndPoint=m_dataSets[i].GetYMaxValue();
			}
		}
	}
	else {
		m_dMinXStartPoint = m_dMinXStartPointOverride;
		m_dMinYStartPoint = m_dMinYStartPointOverride;
		m_dMaxXEndPoint = m_dMaxXEndPointOverride;
		m_dMaxYEndPoint = m_dMaxYEndPointOverride;
	}

    // Sjors 20Apr2016: prevent zero width of x,y axes
    if (fabs(m_dMinXStartPoint - m_dMaxXEndPoint) < 1e-10)
    	m_dMaxXEndPoint += 1;
    if (fabs(m_dMinYStartPoint - m_dMaxYEndPoint) < 1e-10)
    	m_dMaxYEndPoint += 1;

    // compute a 'nice' set of tickmark locations
    double xGraphMin,xGraphMax;
    double yGraphMin,yGraphMax;

    ComputeLabelTickSpacing(m_dMinXStartPoint,m_dMaxXEndPoint,&xGraphMin,
                            &xGraphMax,&m_dXAxisNumbersSpacing,m_iXAxisNumberOfTicks,"x");
    m_dMinXStartPoint=xGraphMin;
    m_dMaxXEndPoint=xGraphMax;

    ComputeLabelTickSpacing(m_dMinYStartPoint,m_dMaxYEndPoint,&yGraphMin,
                            &yGraphMax,&m_dYAxisNumbersSpacing,m_iYAxisNumberOfTicks,"y");
    m_dMinYStartPoint=yGraphMin;
    m_dMaxYEndPoint=yGraphMax;


    m_dMaxXExtent=m_dMaxXEndPoint-m_dMinXStartPoint;
    m_dMaxYExtent=m_dMaxYEndPoint-m_dMinYStartPoint;
    m_dXScale=m_dXAxisSize/m_dMaxXExtent;
    m_dYScale=m_dYAxisSize/m_dMaxYExtent;

    if (m_bFlipY) {
        m_dMinYStartPoint = m_dMaxYEndPoint;
        m_dYScale *= -1;
        m_dFlipYOffset = m_dYAxisSize;
    }

    // establish a pleasing spacing between dots
    // for a dashed line made up of dots... :)
    m_dLineDotSpacing=m_dDiagonalSize*0.01;

    // tick mark length
    m_dTickMarkLength=m_dDiagonalSize*0.01;

}

void CPlot2D::DrawDataPostScript()
{

    // for all data sets
    for (int i=0;i<m_dataSets.size();++i) {
        if (m_dataSets[i].GetNumberOfDataPointsInSet() == 0) continue;

        CDataPoint point;

        if (m_dataSets[i].GetDrawLine()) {

            if (!m_dataSets[i].GetDashedLine()) {
                outputFile << "newpath" << "\n";

                point=m_dataSets[i].GetDataPoint(0);

                outputFile << m_dLeftFrameSize+(-m_dMinXStartPoint+point.GetX())*m_dXScale
                           << " " << m_dBottomFrameSize+(-m_dMinYStartPoint+point.GetY())*m_dYScale << " moveto" << "\n";
                for (int j=1;j<m_dataSets[i].GetNumberOfDataPointsInSet();++j) {
                    point=m_dataSets[i].GetDataPoint(j);
                    outputFile << m_dLeftFrameSize+(-m_dMinXStartPoint+point.GetX())*m_dXScale
                               << " " << m_dBottomFrameSize+(-m_dMinYStartPoint+point.GetY())*m_dYScale << " lineto" << "\n";
                }
                outputFile << m_dataSets[i].GetLineWidth() << " setlinewidth" << "\n";

                double r,g,b;
                m_dataSets[i].GetDatasetColor(&r,&g,&b);

                outputFile << r << " " << g << " " << b << " setrgbcolor" << "\n";
                outputFile << "stroke" << "\n";
            }
            else {
                if (m_dataSets[i].GetDashedLinePattern()=="dot") {

                    for (int j=0;j<m_dataSets[i].GetNumberOfDataPointsInSet()-1;++j) {

                        point=m_dataSets[i].GetDataPoint(j);
                        CDataPoint point1=m_dataSets[i].GetDataPoint(j+1);

                        double segX=(point1.GetX()-point.GetX())*m_dXScale;
                        double segY=(point1.GetY()-point.GetY())*m_dYScale;

                        double segmentLength=sqrt(segX*segX+segY*segY);

                        int numDots=floor(segmentLength/m_dLineDotSpacing);

                        double delX=(point1.GetX()-point.GetX())/numDots;
                        double delY=(point1.GetY()-point.GetY())/numDots;

                        double r,g,b;
                        m_dataSets[i].GetDatasetColor(&r,&g,&b);

                        for (int k=0;k<numDots+1;++k) {
                            outputFile << m_dLeftFrameSize+(-m_dMinXStartPoint+point.GetX()+k*delX)*m_dXScale << " "
                                       << m_dBottomFrameSize+(-m_dMinYStartPoint+point.GetY()+k*delY)*m_dYScale << " "
                                       << m_dataSets[i].GetLineWidth() << " 0 360 arc closepath" << "\n";

                            outputFile << r << " " << g << " " << b << " setrgbcolor" << "\n";
                            outputFile << "fill" << "\n";
                            outputFile << "stroke" << "\n";
                        }
                    }
                }
                if (m_dataSets[i].GetDashedLinePattern()=="dash") {

                    for (int j=0;j<m_dataSets[i].GetNumberOfDataPointsInSet()-1;++j) {

                        point=m_dataSets[i].GetDataPoint(j);
                        CDataPoint point1=m_dataSets[i].GetDataPoint(j+1);

                        double segX=(point1.GetX()-point.GetX())*m_dXScale;
                        double segY=(point1.GetY()-point.GetY())*m_dYScale;

                        double segmentLength=sqrt(segX*segX+segY*segY);

                        int numDashes=floor(segmentLength/m_dLineDotSpacing);

                        if (numDashes%2==0) {
                            ++numDashes;
                        }

                        double delX=(point1.GetX()-point.GetX())/numDashes;
                        double delY=(point1.GetY()-point.GetY())/numDashes;

                        double r,g,b;
                        m_dataSets[i].GetDatasetColor(&r,&g,&b);

                        for (int k=0;k<numDashes+1;k+=2) {

                            outputFile << m_dLeftFrameSize+(-m_dMinXStartPoint+point.GetX())*m_dXScale+(k)*delX*m_dXScale << " "
                            << m_dBottomFrameSize+(-m_dMinYStartPoint+point.GetY())*m_dYScale+(k)*delY*m_dYScale << " moveto" << "\n";

                            outputFile << delX*m_dXScale << " " << delY*m_dYScale << " rlineto" << "\n";

                            outputFile << r << " " << g << " " << b << " setrgbcolor" << "\n";
                            outputFile << "stroke" << "\n";
                        }
                    }
                }
                if (m_dataSets[i].GetDashedLinePattern()=="dash_dot") {

                    for (int j=0;j<m_dataSets[i].GetNumberOfDataPointsInSet()-1;++j) {

                        point=m_dataSets[i].GetDataPoint(j);
                        CDataPoint point1=m_dataSets[i].GetDataPoint(j+1);

                        double segX=(point1.GetX()-point.GetX())*m_dXScale;
                        double segY=(point1.GetY()-point.GetY())*m_dYScale;

                        double segmentLength=sqrt(segX*segX+segY*segY);

                        int numDashes=floor(segmentLength/(m_dLineDotSpacing*1.5));

                        if (numDashes%2==0) {
                            ++numDashes;
                        }

                        double delX=(point1.GetX()-point.GetX())/numDashes;
                        double delY=(point1.GetY()-point.GetY())/numDashes;

                        double r,g,b;
                        m_dataSets[i].GetDatasetColor(&r,&g,&b);

                        for (int k=0;k<numDashes+1;k+=2) {

                            outputFile << m_dLeftFrameSize+(-m_dMinXStartPoint+point.GetX())*m_dXScale+(k)*delX*m_dXScale << " "
                            << m_dBottomFrameSize+(-m_dMinYStartPoint+point.GetY())*m_dYScale+(k)*delY*m_dYScale << " moveto" << "\n";

                            outputFile << delX*m_dXScale << " " << delY*m_dYScale << " rlineto" << "\n";

                            outputFile << r << " " << g << " " << b << " setrgbcolor" << "\n";
                            outputFile << "stroke" << "\n";

                            outputFile << m_dLeftFrameSize+(-m_dMinXStartPoint+point.GetX())*m_dXScale+(k+1.5)*delX*m_dXScale << " "
                            << m_dBottomFrameSize+(-m_dMinYStartPoint+point.GetY())*m_dYScale+(k+1.5)*delY*m_dYScale << " moveto" << "\n";

                            if (k<numDashes-1) {
                                outputFile << "newpath" << "\n";
                                outputFile << m_dLeftFrameSize+(-m_dMinXStartPoint+point.GetX())*m_dXScale+(k+1.5)*delX*m_dXScale << " "
                                << m_dBottomFrameSize+(-m_dMinYStartPoint+point.GetY())*m_dYScale+(k+1.5)*delY*m_dYScale << " "
                                << m_dataSets[i].GetLineWidth()*0.5 << " 0 360 arc closepath" << "\n";
                                outputFile << "fill" << "\n";
                            }
                        }
                    }
                }
            }
        }

        if (m_dataSets[i].GetDrawMarker()) {
            for (int j=0;j<m_dataSets[i].GetNumberOfDataPointsInSet();++j) {
                point=m_dataSets[i].GetDataPoint(j);
                DrawMarker(m_dataSets[i].GetMarkerSymbol(),m_dataSets[i].GetMarkerSize(),m_dataSets[i].GetDrawMarkerFilled(),
                           m_dLeftFrameSize+(-m_dMinXStartPoint+point.GetX())*m_dXScale,m_dBottomFrameSize+(-m_dMinYStartPoint+point.GetY())*m_dYScale,i);
            }
        }
    }
}

void CPlot2D::AddDataSet(int numPoints, double *xValues, double *yValues)
{
    CDataSet dataSet;
    for (int i=0;i<numPoints;++i) {
        CDataPoint point=CDataPoint(xValues[i],yValues[i]);
        dataSet.AddDataPoint(point);
    }
    m_dataSets.push_back(dataSet);
}

void CPlot2D::AddDataSet(std::vector<RFLOAT> xValues, std::vector<RFLOAT> yValues)
{
    CDataSet dataSet;
    if (m_dataSets.size() == 0)
    	dataSet.SetDatasetColor(1., 0., 0.);
    else if (m_dataSets.size() == 1)
    	dataSet.SetDatasetColor(0., 1., 0.);
    else
    	dataSet.SetDatasetColor(0., 0., 1.);

    dataSet.SetDrawMarker(false);

    if (xValues.size() != yValues.size())
    {
    	REPORT_ERROR("ERROR: xValues and yValues vectors do not have identical sizes.");
    }

    for (long int i = 0; i < yValues.size(); i++)
    {
    	CDataPoint point=CDataPoint(xValues[i],yValues[i]);
        dataSet.AddDataPoint(point);
    }
    m_dataSets.push_back(dataSet);
}

void CPlot2D::AddDataSet(std::vector<RFLOAT> yValues)
{
    CDataSet dataSet;
    if (m_dataSets.size() == 0)
    	dataSet.SetDatasetColor(1., 0., 0.);
    else if (m_dataSets.size() == 1)
    	dataSet.SetDatasetColor(0., 1., 0.);
    else
    	dataSet.SetDatasetColor(0., 0., 1.);

   dataSet.SetDrawMarker(false);

    for (long int i = 0; i < yValues.size(); i++)
    {
    	CDataPoint point=CDataPoint(i+1,yValues[i]);
        dataSet.AddDataPoint(point);
    }
	m_dataSets.push_back(dataSet);
}

void CPlot2D::SetViewArea(double start_x, double start_y, double end_x, double end_y)
{
	m_dMinXStartPointOverride = start_x;
	m_dMaxXEndPointOverride = end_x;
	m_dMinYStartPointOverride = start_y;
	m_dMaxYEndPointOverride = end_y;
	m_bSizeSetExternally = true;
}

void CPlot2D::DrawMarker(std::string symbol, double size, bool filled, double xLocation, double yLocation, int dataSet)
{
    double r,g,b;
    m_dataSets[dataSet].GetDatasetColor(&r,&g,&b);


    if (symbol=="o") {
        outputFile << xLocation << " " << yLocation << " moveto" << "\n";
        outputFile << xLocation << " " << yLocation << " " << size*0.5 << " 0 360 arc closepath" << "\n";
        outputFile << r << " " << g << " " << b << " setrgbcolor" << "\n";
        if (filled) {
            outputFile << "fill" << "\n";
        }
        outputFile << "stroke" << "\n";
    }
    if (symbol=="x" || symbol=="*") {
        double halfSize=0.5*size;
        outputFile << xLocation-halfSize << " " << yLocation-halfSize << " moveto" << "\n";
        outputFile << size << " " << size << " rlineto" << "\n";
        outputFile << -size << " " << 0 << " rmoveto" << "\n";
        outputFile << size << " " << -size << " rlineto" << "\n";
        outputFile << r << " " << g << " " << b << " setrgbcolor" << "\n";
        outputFile << "stroke" << "\n";
    }
    if (symbol=="+" || symbol=="*") {
        double halfSize=0.5*size;
        outputFile << xLocation-halfSize << " " << yLocation << " moveto" << "\n";
        outputFile << size << " " << 0 << " rlineto" << "\n";
        outputFile << -halfSize << " " << halfSize << " rmoveto" << "\n";
        outputFile << 0 << " " << -size << " rlineto" << "\n";
        outputFile << r << " " << g << " " << b << " setrgbcolor" << "\n";
        outputFile << "stroke" << "\n";
    }
    if (symbol=="diamond") {
        double halfSize=0.5*size;
        outputFile << "newpath" << "\n";
        outputFile << xLocation << " " << yLocation-halfSize << " moveto" << "\n";
        outputFile << halfSize << " " << halfSize << " rlineto" << "\n";
        outputFile << -halfSize << " " << halfSize << " rlineto" << "\n";
        outputFile << -halfSize << " " << -halfSize << " rlineto" << "\n";
        outputFile << "closepath" << "\n";
        outputFile << r << " " << g << " " << b << " setrgbcolor" << "\n";
        if (filled) {
            outputFile << "fill" << "\n";
        }
        outputFile << "stroke" << "\n";
    }
    if (symbol=="square") {
        double halfSize=0.5*size;
        outputFile << "newpath" << "\n";
        outputFile << xLocation-halfSize << " " << yLocation-halfSize << " moveto" << "\n";
        outputFile << size << " " << 0 << " rlineto" << "\n";
        outputFile << 0 << " " << size << " rlineto" << "\n";
        outputFile << -size << " " << 0 << " rlineto" << "\n";
        outputFile << "closepath" << "\n";
        outputFile << r << " " << g << " " << b << " setrgbcolor" << "\n";
        if (filled) {
            outputFile << "fill" << "\n";
        }
        outputFile << "stroke" << "\n";
    }
    if (symbol=="triangle") {
        double halfSize=0.5*size;
        double halfEdgeLength=1.5/sqrt(3.0)*halfSize;

        outputFile << "newpath" << "\n";
        outputFile << xLocation-halfEdgeLength << " " << yLocation-halfSize*0.5 << " moveto" << "\n";
        outputFile << halfEdgeLength*2.0 << " " << 0 << " rlineto" << "\n";
        outputFile << -halfEdgeLength << " " << 1.5*halfSize << " rlineto" << "\n";
        outputFile << "closepath" << "\n";
        outputFile << r << " " << g << " " << b << " setrgbcolor" << "\n";
        if (filled) {
            outputFile << "fill" << "\n";
        }
        outputFile << "stroke" << "\n";
    }
}

void CPlot2D::DrawXAxisTickMarksPostScript()
{
    for (int i=0;i<m_iXAxisNumberOfLabels;++i) {
        outputFile << m_dLeftFrameSize+(i*m_dXAxisNumbersSpacing)*m_dXScale << " "
                   << m_dBottomFrameSize << " moveto" << "\n";
        outputFile << 0 << " " << m_dTickMarkLength << " rlineto" << "\n";
        outputFile << m_dFrameLineWidth << " setlinewidth" << "\n";
        outputFile << m_dFrameColor[0] << " " << m_dFrameColor[1] << " " << m_dFrameColor[2] << " setrgbcolor" << "\n";
        outputFile << "stroke" << "\n";
    }
}

void CPlot2D::DrawYAxisTickMarksPostScript()
{
    for (int i=0;i<m_iYAxisNumberOfLabels;++i) {
        outputFile << m_dLeftFrameSize << " "
                   << m_dBottomFrameSize+m_dFlipYOffset+(i*m_dYAxisNumbersSpacing)*m_dYScale << " moveto" << "\n";
        outputFile << m_dTickMarkLength << " " << 0 << " rlineto" << "\n";
        outputFile << m_dFrameLineWidth << " setlinewidth" << "\n";
        outputFile << m_dFrameColor[0] << " " << m_dFrameColor[1] << " " << m_dFrameColor[2] << " setrgbcolor" << "\n";
        outputFile << "stroke" << "\n";
    }
}

void CPlot2D::DrawXAxisGridLinesPostScript()
{
    if (!m_bDrawGridLinesDashed) {
        for (int i=0;i<m_iXAxisNumberOfLabels;++i) {
            outputFile << m_dLeftFrameSize+(i*m_dXAxisNumbersSpacing)*m_dXScale << " "
                       << m_dBottomFrameSize << " moveto" << "\n";
            outputFile << 0 << " " << m_dYAxisSize << " rlineto" << "\n";
            outputFile << m_dGridLineWidth << " setlinewidth" << "\n";
            outputFile << m_dGridColor[0] << " " << m_dGridColor[1] << " " << m_dGridColor[2] << " setrgbcolor" << "\n";
            outputFile << "stroke" << "\n";
        }
    }
    else {

        int numDashes=floor(m_dYAxisSize/m_dLineDotSpacing);

        if (numDashes%2==0) {
            ++numDashes;
        }

        double delY=m_dYAxisSize/numDashes;

        for (int i=0;i<m_iXAxisNumberOfLabels;++i) {

            outputFile << m_dLeftFrameSize+(i*m_dXAxisNumbersSpacing)*m_dXScale << " "
                       << m_dBottomFrameSize << " moveto" << "\n";

            for (int k=0;k<numDashes+1;k+=2) {

                outputFile << m_dLeftFrameSize+(i*m_dXAxisNumbersSpacing)*m_dXScale << " "
                           << m_dBottomFrameSize+k*delY << " moveto" << "\n";
                outputFile << 0 << " " << delY << " rlineto" << "\n";

                outputFile << m_dGridLineWidth << " setlinewidth" << "\n";
                outputFile << m_dGridColor[0] << " " << m_dGridColor[1] << " " << m_dGridColor[2] << " setrgbcolor" << "\n";
                outputFile << "stroke" << "\n";
            }
        }
    }
}

void CPlot2D::DrawYAxisGridLinesPostScript()
{
    if (!m_bDrawGridLinesDashed) {
        for (int i=0;i<m_iYAxisNumberOfLabels;++i) {
            outputFile << m_dLeftFrameSize << " " << m_dBottomFrameSize+m_dFlipYOffset+(i*m_dYAxisNumbersSpacing)*m_dYScale << " moveto" << "\n";
            outputFile << m_dXAxisSize << " " << 0 << " rlineto" << "\n";
            outputFile << m_dGridLineWidth << " setlinewidth" << "\n";
            outputFile << m_dGridColor[0] << " " << m_dGridColor[1] << " " << m_dGridColor[2] << " setrgbcolor" << "\n";
            outputFile << "stroke" << "\n";
        }
    }
    else {

        int numDashes=floor(m_dXAxisSize/m_dLineDotSpacing);

        if (numDashes%2==0) {
            ++numDashes;
        }

        double delX=m_dXAxisSize/numDashes;

        for (int i=0;i<m_iYAxisNumberOfLabels;++i) {

            outputFile << m_dLeftFrameSize << " " << m_dBottomFrameSize+m_dFlipYOffset+(i*m_dYAxisNumbersSpacing)*m_dYScale << " moveto" << "\n";

            for (int k=0;k<numDashes+1;k+=2) {

                outputFile << m_dLeftFrameSize+k*delX  << " " << m_dBottomFrameSize+m_dFlipYOffset+(i*m_dYAxisNumbersSpacing)*m_dYScale << " moveto" << "\n";
                outputFile << delX << " " << 0 << " rlineto" << "\n";

                outputFile << m_dGridLineWidth << " setlinewidth" << "\n";
                outputFile << m_dGridColor[0] << " " << m_dGridColor[1] << " " << m_dGridColor[2] << " setrgbcolor" << "\n";
                outputFile << "stroke" << "\n";
            }
        }
    }
}

//
//  Axis tick mark spacing and labeling algorithm inspired by
//  an article in Graphics Gems I by P. Heckbert,
//  "Nice Numbers for Graph Labels"
//  Code available from: http://www.cs.cmu.edu/~ph/
//

void CPlot2D::ComputeLabelTickSpacing(double dataMin, double dataMax, double *plotMin, double *plotMax, double *tickSpacing, int numTicks, std::string axis)
{

    int nfrac;
    double d;
    double range;


    range=NiceNum(dataMax-dataMin, 0);
    d=NiceNum(range/(numTicks-1),1);
    *plotMin=floor(dataMin/d)*d;
    *plotMax=ceil(dataMax/d)*d;
    *tickSpacing=d;
    nfrac=(int)fmax(-floor(log10(d)),0);


    if (axis=="x") {
        snprintf(m_cXAxisLabelFormat,20,"%%.%df",nfrac);
        char temp[20];
        m_iXAxisNumberOfLabels=0;
        for (double x=*plotMin; x<*plotMax+.5*d; x+=d) {
            snprintf(temp,20,m_cXAxisLabelFormat,x);
            m_strXAxisLabels.push_back(temp);
            m_iXAxisNumberOfLabels++;
        }
    }
    else if (axis=="y") {
        snprintf(m_cYAxisLabelFormat,20,"%%.%df",nfrac);
        char temp[20];
        m_iYAxisNumberOfLabels=0;
        for (double x=*plotMin; x<*plotMax+.5*d; x+=d) {
            snprintf(temp,20,m_cYAxisLabelFormat,x);
            m_strYAxisLabels.push_back(temp);
            m_iYAxisNumberOfLabels++;
        }
    }
}

//
//  Axis tick mark spacing and labeling algorithm inspired by
//  an article in Graphics Gems I by P. Heckbert,
//  "Nice Numbers for Graph Labels"
//  Code available from: http://www.cs.cmu.edu/~ph/
//

double CPlot2D::NiceNum(double x, int round)
{
    int expv;
    double f;
    double nf;

    expv=floor(log10(x));
    f=x/pow(10.0,expv);

    if (round) {
        if (f<1.5) {
            nf=1.0;
        }
        else if (f<3.0) {
            nf=2.0;
        }
        else if (f<7.0) {
            nf=5.0;
        }
        else {
            nf=10.0;
        }
    }
    else {
        if (f<=1.0) {
            nf=1.0;
        }
        else if (f<=2.0) {
            nf=2.0;
        }
        else if (f<=5.0) {
            nf=5.0;
        }
        else {
            nf=10.0;
        }
    }

    return (nf*pow(10.0, expv));
}

void CPlot2D::DrawXAxisLabelsPostScript()
{
    double labelXCoordinate,labelYCoordinate;


    for (int i=0;i<m_iXAxisNumberOfLabels;++i) {

        // adjustment for the label
        labelXCoordinate=m_dLeftFrameSize+(i*m_dXAxisNumbersSpacing)*m_dXScale;
        labelYCoordinate=m_dBottomFrameSize-1.25*m_dXAxisLabelFontSize;

        outputFile << "/" << m_strXAxisLabelFont << " findfont" << "\n";
        outputFile << m_dXAxisLabelFontSize << " scalefont" << "\n";
        outputFile << "setfont" << "\n";
        outputFile << labelXCoordinate << " " << labelYCoordinate << " moveto" << "\n";
        outputFile << m_dXAxisLabelColor[0] << " " << m_dXAxisLabelColor[1] << " " << m_dXAxisLabelColor[2] << " setrgbcolor" << "\n";

        // let PostScript handle the final adjustment based on the width of the string
        outputFile << "(" << m_strXAxisLabels[i] << ")" << " dup stringwidth pop 2 div neg 0 rmoveto show"  << "\n";
    }
}

void CPlot2D::DrawYAxisLabelsPostScript()
{
    double labelXCoordinate,labelYCoordinate;

    for (int i=0;i<m_iYAxisNumberOfLabels;++i) {

        labelXCoordinate=m_dLeftFrameSize;
        labelYCoordinate=m_dBottomFrameSize+m_dFlipYOffset+(i*m_dYAxisNumbersSpacing)*m_dYScale;

        // adjustment for the label
        labelXCoordinate-=m_dYAxisLabelFontSize*0.5;
        labelYCoordinate-=m_dYAxisLabelFontSize*0.25;

        outputFile << "/" << m_strYAxisLabelFont << " findfont" << "\n";
        outputFile << m_dYAxisLabelFontSize << " scalefont" << "\n";
        outputFile << "setfont" << "\n";
        outputFile << labelXCoordinate << " " << labelYCoordinate << " moveto" << "\n";
        outputFile << m_dYAxisLabelColor[0] << " " << m_dYAxisLabelColor[1] << " " << m_dYAxisLabelColor[2] << " setrgbcolor" << "\n";

        // let PostScript handle the final adjustment based on the width of the string
        outputFile << "(" << m_strYAxisLabels[i] << ")" << " dup stringwidth pop neg 0 rmoveto show"  << "\n";
    }
}

void CPlot2D::DrawXAxisTitlePostScript()
{
    double labelXCoordinate,labelYCoordinate;

    labelXCoordinate=m_dLeftFrameSize+m_dXAxisSize*0.5;
    labelYCoordinate=m_dBottomFrameSize-m_dXAxisTitleFontSize*2.0;

    outputFile << "/" << m_strXAxisTitleFont << " findfont" << "\n";
    outputFile << m_dYAxisTitleFontSize << " scalefont" << "\n";
    outputFile << "setfont" << "\n";
    outputFile << labelXCoordinate << " " << labelYCoordinate << " moveto" << "\n";
    outputFile << m_dXAxisTitleColor[0] << " " << m_dXAxisTitleColor[1] << " " << m_dXAxisTitleColor[2] << " setrgbcolor" << "\n";

    // let PostScript handle the final adjustment based on the width of the string
    outputFile << "(" << m_strXAxisTitle << ")" << " dup stringwidth pop 2 div neg 0 rmoveto show"  << "\n";
}

void CPlot2D::DrawYAxisTitlePostScript()
{
    double labelXCoordinate,labelYCoordinate;

    labelXCoordinate=10.0+m_dYAxisTitleFontSize;
    labelYCoordinate=m_dBottomFrameSize+m_dYAxisSize*0.5;

    outputFile << "/" << m_strYAxisTitleFont << " findfont" << "\n";
    outputFile << m_dYAxisTitleFontSize << " scalefont" << "\n";
    outputFile << "setfont" << "\n";
    outputFile << labelXCoordinate << " " << labelYCoordinate << " moveto" << "\n";
    outputFile << "90 rotate" << "\n";
    outputFile << m_dYAxisTitleColor[0] << " " << m_dYAxisTitleColor[1] << " " << m_dYAxisTitleColor[2] << " setrgbcolor" << "\n";

    // let PostScript handle the final adjustment based on the width of the string
    outputFile << "(" << m_strYAxisTitle << ")" << " dup stringwidth pop 2 div neg 0 rmoveto show"  << "\n";
    outputFile << "-90 rotate" << "\n";
}

void CPlot2D::DrawLegendPostScript()
{
    int numDataSets=(int)m_dataSets.size();

    double widthPerLegend=m_dXTotalSize*0.8/numDataSets; // account for margins

    double symbolWidth=30.0;
    double widthPerLegendText=widthPerLegend-symbolWidth;

    int maxNumberOfCharactersInDatasetTitles=0;

    for (int i=0;i<numDataSets;++i) {
        if (m_dataSets[i].GetDatasetTitle().size()>maxNumberOfCharactersInDatasetTitles) {
            maxNumberOfCharactersInDatasetTitles=(int)m_dataSets[i].GetDatasetTitle().size();
        }
    }

    double maxFontHeight=widthPerLegendText/(maxNumberOfCharactersInDatasetTitles*0.75);

    if (maxFontHeight>14.0) {
        maxFontHeight=14.0;
    }
    if (maxFontHeight<6.0) {
        maxFontHeight=6.0;
    }


    double legendXCoordinate,legendYCoordinate;

    for (int i=0;i<numDataSets;++i) {

        double r,g,b;
        m_dataSets[i].GetDatasetColor(&r,&g,&b);

		if (m_dataSets[i].GetDatasetTitle() == "") continue;

        if (m_dataSets[i].GetDrawLine()) {
            // draw the line
            if (!m_dataSets[i].GetDashedLine()) {
                legendXCoordinate=m_dXTotalSize*0.1+i*widthPerLegend+symbolWidth*0.25;
                legendYCoordinate=m_dBottomFrameSize*0.25+maxFontHeight*0.3;
                outputFile << legendXCoordinate << " " << legendYCoordinate << " moveto" << "\n";
                outputFile << symbolWidth*0.5 << " " << 0 << " rlineto" << "\n";
                outputFile << m_dataSets[i].GetLineWidth() << " setlinewidth" << "\n";
                outputFile << r << " " << g << " " << b << " setrgbcolor" << "\n";
                outputFile << "stroke" << "\n";
            }
            else {
                if (m_dataSets[i].GetDashedLinePattern()=="dot") {

                    for (int j=0;j<5;++j) {

                        legendXCoordinate=m_dXTotalSize*0.1+i*widthPerLegend+symbolWidth*0.25+j*(symbolWidth*0.5/4.0);
                        legendYCoordinate=m_dBottomFrameSize*0.25+maxFontHeight*0.3;

                        outputFile << legendXCoordinate << " " << legendYCoordinate << " "
                        << m_dataSets[i].GetLineWidth() << " 0 360 arc closepath" << "\n";

                        outputFile << r << " " << g << " " << b << " setrgbcolor" << "\n";
                        outputFile << "fill" << "\n";
                        outputFile << "stroke" << "\n";
                    }
                }

                if (m_dataSets[i].GetDashedLinePattern()=="dash") {

                    for (int j=0;j<5;j+=2) {

                        legendXCoordinate=m_dXTotalSize*0.1+i*widthPerLegend+symbolWidth*0.25+j*(symbolWidth*0.5/5.0);
                        legendYCoordinate=m_dBottomFrameSize*0.25+maxFontHeight*0.3;

                        outputFile << legendXCoordinate << " " << legendYCoordinate << " moveto" << "\n";

                        outputFile << symbolWidth*0.5/5.0 << " " << 0 << " rlineto" << "\n";
                        outputFile << r << " " << g << " " << b << " setrgbcolor" << "\n";
                        outputFile << "stroke" << "\n";
                    }
                }
                if (m_dataSets[i].GetDashedLinePattern()=="dash_dot") {

                    legendXCoordinate=m_dXTotalSize*0.1+i*widthPerLegend+symbolWidth*0.25;
                    legendYCoordinate=m_dBottomFrameSize*0.25+maxFontHeight*0.3;


                    outputFile << legendXCoordinate << " " << legendYCoordinate << " moveto" << "\n";
                    outputFile << symbolWidth*0.5/2.0 << " " << 0 << " rlineto" << "\n";
                    outputFile << r << " " << g << " " << b << " setrgbcolor" << "\n";
                    outputFile << "stroke" << "\n";

                    legendXCoordinate=m_dXTotalSize*0.1+i*widthPerLegend+symbolWidth*0.75;
                    legendYCoordinate=m_dBottomFrameSize*0.25+maxFontHeight*0.3;

                    outputFile << legendXCoordinate << " " << legendYCoordinate << " "
                    << m_dataSets[i].GetLineWidth()*0.5 << " 0 360 arc closepath" << "\n";
                    outputFile << "fill" << "\n";
                }
            }
        }

        // draw marker
        if (m_dataSets[i].GetDrawMarker()) {
            legendXCoordinate=m_dXTotalSize*0.1+i*widthPerLegend+symbolWidth*0.5;
            legendYCoordinate=m_dBottomFrameSize*0.25+maxFontHeight*0.3;
            DrawMarker(m_dataSets[i].GetMarkerSymbol(),m_dataSets[i].GetMarkerSize(),m_dataSets[i].GetDrawMarkerFilled(),
                       legendXCoordinate,legendYCoordinate,i);
        }

        // write the dataset title
        legendXCoordinate=m_dXTotalSize*0.1+i*widthPerLegend+symbolWidth;
        legendYCoordinate=m_dBottomFrameSize*0.25;

        outputFile << legendXCoordinate << " " << legendYCoordinate << " moveto" << "\n";

        outputFile << "/" << m_dataSets[i].GetDatasetLegendFont() << " findfont" << "\n";
        outputFile << maxFontHeight << " scalefont" << "\n";
        outputFile << "setfont" << "\n";
        outputFile << legendXCoordinate << " " << legendYCoordinate << " moveto" << "\n";

        outputFile << 0.0 << " " << 0.0 << " " << 0.0 << " setrgbcolor" << "\n";

        outputFile << "(" <<  m_dataSets[i].GetDatasetTitle() << ")" << " show"  << "\n";
    }

}
