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


	std::string command = "gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -dDEVICEWIDTHPOINTS=800 -dDEVICEHEIGHTPOINTS=800 -sOutputFile=";
	command += fn_pdf + " ";
	for (int i = 0; i < fn_eps.size(); i++)
		command += fn_eps[i] + " ";

	command += " > /dev/null &";
	if (system(command.c_str()))
	{
		std::cerr << " ERROR in executing: " << command << std::endl;
		std::cerr << " + Will make an empty PDF-file in " << fn_pdf << std::endl;
		std::cerr << " + Solve your issue with the ps-command to get better PDF logfiles." << std::endl;
		touch(fn_pdf);
	}
}

CPlot2D::CPlot2D(std::string title)
{
    m_dXAxisSize=809.0; // Golden Ratio
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
    outputFile << "%!PS-Adobe-2.0 EPSF-1.2" << std::endl;
    outputFile << "%%BoundingBox: 0 0 " << (int)m_dXTotalSize << " " <<  (int)m_dYTotalSize << std::endl;
    outputFile << "%%Pages: 1" << std::endl;
    outputFile << "%%EndComments" << std::endl;

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
    outputFile << "newpath" << std::endl;
    outputFile << m_dLeftFrameSize << " " << m_dBottomFrameSize << " moveto" << std::endl;
    outputFile << m_dXAxisSize << " " << 0 << " rlineto" << std::endl;
    outputFile << 0 << " " << m_dYAxisSize << " rlineto" << std::endl;
    outputFile << -m_dXAxisSize << " " << 0 << " rlineto" << std::endl;
    outputFile << "closepath" << std::endl;

    outputFile << m_dFrameLineWidth << " setlinewidth" << std::endl;
    outputFile << m_dFrameColor[0] << " " << m_dFrameColor[1] << " " << m_dFrameColor[2] << " setrgbcolor" << std::endl;
    outputFile << "stroke" << std::endl;

    double labelXCoordinate,labelYCoordinate;

    labelXCoordinate=m_dLeftFrameSize+m_dXAxisSize*0.5;
    labelYCoordinate=m_dBottomFrameSize + m_dYAxisSize + 10;

    outputFile << "/" << m_strXAxisTitleFont << " findfont" << std::endl;
    outputFile << m_dYAxisTitleFontSize << " scalefont" << std::endl;
    outputFile << "setfont" << std::endl;
    outputFile << labelXCoordinate << " " << labelYCoordinate << " moveto" << std::endl;
    outputFile << m_dXAxisTitleColor[0] << " " << m_dXAxisTitleColor[1] << " " << m_dXAxisTitleColor[2] << " setrgbcolor" << std::endl;

    // let PostScript handle the final adjustment based on the width of the string
    outputFile << "(" << m_strPlotTitle << ")" << " dup stringwidth pop 2 div neg 0 rmoveto show"  << std::endl;

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

        CDataPoint point;

        if (m_dataSets[i].GetDrawLine()) {

            if (!m_dataSets[i].GetDashedLine()) {
                outputFile << "newpath" << std::endl;

                point=m_dataSets[i].GetDataPoint(0);

                outputFile << m_dLeftFrameSize+(-m_dMinXStartPoint+point.GetX())*m_dXScale
                           << " " << m_dBottomFrameSize+(-m_dMinYStartPoint+point.GetY())*m_dYScale << " moveto" << std::endl;
                for (int j=1;j<m_dataSets[i].GetNumberOfDataPointsInSet();++j) {
                    point=m_dataSets[i].GetDataPoint(j);
                    outputFile << m_dLeftFrameSize+(-m_dMinXStartPoint+point.GetX())*m_dXScale
                               << " " << m_dBottomFrameSize+(-m_dMinYStartPoint+point.GetY())*m_dYScale << " lineto" << std::endl;
                }
                outputFile << m_dataSets[i].GetLineWidth() << " setlinewidth" << std::endl;

                double r,g,b;
                m_dataSets[i].GetDatasetColor(&r,&g,&b);

                outputFile << r << " " << g << " " << b << " setrgbcolor" << std::endl;
                outputFile << "stroke" << std::endl;
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
                                       << m_dataSets[i].GetLineWidth() << " 0 360 arc closepath" << std::endl;

                            outputFile << r << " " << g << " " << b << " setrgbcolor" << std::endl;
                            outputFile << "fill" << std::endl;
                            outputFile << "stroke" << std::endl;
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
                            << m_dBottomFrameSize+(-m_dMinYStartPoint+point.GetY())*m_dYScale+(k)*delY*m_dYScale << " moveto" << std::endl;

                            outputFile << delX*m_dXScale << " " << delY*m_dYScale << " rlineto" << std::endl;

                            outputFile << r << " " << g << " " << b << " setrgbcolor" << std::endl;
                            outputFile << "stroke" << std::endl;
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
                            << m_dBottomFrameSize+(-m_dMinYStartPoint+point.GetY())*m_dYScale+(k)*delY*m_dYScale << " moveto" << std::endl;

                            outputFile << delX*m_dXScale << " " << delY*m_dYScale << " rlineto" << std::endl;

                            outputFile << r << " " << g << " " << b << " setrgbcolor" << std::endl;
                            outputFile << "stroke" << std::endl;

                            outputFile << m_dLeftFrameSize+(-m_dMinXStartPoint+point.GetX())*m_dXScale+(k+1.5)*delX*m_dXScale << " "
                            << m_dBottomFrameSize+(-m_dMinYStartPoint+point.GetY())*m_dYScale+(k+1.5)*delY*m_dYScale << " moveto" << std::endl;

                            if (k<numDashes-1) {
                                outputFile << "newpath" << std::endl;
                                outputFile << m_dLeftFrameSize+(-m_dMinXStartPoint+point.GetX())*m_dXScale+(k+1.5)*delX*m_dXScale << " "
                                << m_dBottomFrameSize+(-m_dMinYStartPoint+point.GetY())*m_dYScale+(k+1.5)*delY*m_dYScale << " "
                                << m_dataSets[i].GetLineWidth()*0.5 << " 0 360 arc closepath" << std::endl;
                                outputFile << "fill" << std::endl;
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

void CPlot2D::DrawMarker(std::string symbol, double size, bool filled, double xLocation, double yLocation, int dataSet)
{
    double r,g,b;
    m_dataSets[dataSet].GetDatasetColor(&r,&g,&b);


    if (symbol=="o") {
        outputFile << xLocation << " " << yLocation << " moveto" << std::endl;
        outputFile << xLocation << " " << yLocation << " " << size*0.5 << " 0 360 arc closepath" << std::endl;
        outputFile << r << " " << g << " " << b << " setrgbcolor" << std::endl;
        if (filled) {
            outputFile << "fill" << std::endl;
        }
        outputFile << "stroke" << std::endl;
    }
    if (symbol=="x" || symbol=="*") {
        double halfSize=0.5*size;
        outputFile << xLocation-halfSize << " " << yLocation-halfSize << " moveto" << std::endl;
        outputFile << size << " " << size << " rlineto" << std::endl;
        outputFile << -size << " " << 0 << " rmoveto" << std::endl;
        outputFile << size << " " << -size << " rlineto" << std::endl;
        outputFile << r << " " << g << " " << b << " setrgbcolor" << std::endl;
        outputFile << "stroke" << std::endl;
    }
    if (symbol=="+" || symbol=="*") {
        double halfSize=0.5*size;
        outputFile << xLocation-halfSize << " " << yLocation << " moveto" << std::endl;
        outputFile << size << " " << 0 << " rlineto" << std::endl;
        outputFile << -halfSize << " " << halfSize << " rmoveto" << std::endl;
        outputFile << 0 << " " << -size << " rlineto" << std::endl;
        outputFile << r << " " << g << " " << b << " setrgbcolor" << std::endl;
        outputFile << "stroke" << std::endl;
    }
    if (symbol=="diamond") {
        double halfSize=0.5*size;
        outputFile << "newpath" << std::endl;
        outputFile << xLocation << " " << yLocation-halfSize << " moveto" << std::endl;
        outputFile << halfSize << " " << halfSize << " rlineto" << std::endl;
        outputFile << -halfSize << " " << halfSize << " rlineto" << std::endl;
        outputFile << -halfSize << " " << -halfSize << " rlineto" << std::endl;
        outputFile << "closepath" << std::endl;
        outputFile << r << " " << g << " " << b << " setrgbcolor" << std::endl;
        if (filled) {
            outputFile << "fill" << std::endl;
        }
        outputFile << "stroke" << std::endl;
    }
    if (symbol=="square") {
        double halfSize=0.5*size;
        outputFile << "newpath" << std::endl;
        outputFile << xLocation-halfSize << " " << yLocation-halfSize << " moveto" << std::endl;
        outputFile << size << " " << 0 << " rlineto" << std::endl;
        outputFile << 0 << " " << size << " rlineto" << std::endl;
        outputFile << -size << " " << 0 << " rlineto" << std::endl;
        outputFile << "closepath" << std::endl;
        outputFile << r << " " << g << " " << b << " setrgbcolor" << std::endl;
        if (filled) {
            outputFile << "fill" << std::endl;
        }
        outputFile << "stroke" << std::endl;
    }
    if (symbol=="triangle") {
        double halfSize=0.5*size;
        double halfEdgeLength=1.5/sqrt(3.0)*halfSize;

        outputFile << "newpath" << std::endl;
        outputFile << xLocation-halfEdgeLength << " " << yLocation-halfSize*0.5 << " moveto" << std::endl;
        outputFile << halfEdgeLength*2.0 << " " << 0 << " rlineto" << std::endl;
        outputFile << -halfEdgeLength << " " << 1.5*halfSize << " rlineto" << std::endl;
        outputFile << "closepath" << std::endl;
        outputFile << r << " " << g << " " << b << " setrgbcolor" << std::endl;
        if (filled) {
            outputFile << "fill" << std::endl;
        }
        outputFile << "stroke" << std::endl;
    }
}

void CPlot2D::DrawXAxisTickMarksPostScript()
{
    for (int i=0;i<m_iXAxisNumberOfLabels;++i) {
        outputFile << m_dLeftFrameSize+(i*m_dXAxisNumbersSpacing)*m_dXScale << " "
                   << m_dBottomFrameSize << " moveto" << std::endl;
        outputFile << 0 << " " << m_dTickMarkLength << " rlineto" << std::endl;
        outputFile << m_dFrameLineWidth << " setlinewidth" << std::endl;
        outputFile << m_dFrameColor[0] << " " << m_dFrameColor[1] << " " << m_dFrameColor[2] << " setrgbcolor" << std::endl;
        outputFile << "stroke" << std::endl;
    }
}

void CPlot2D::DrawYAxisTickMarksPostScript()
{
    for (int i=0;i<m_iYAxisNumberOfLabels;++i) {
        outputFile << m_dLeftFrameSize << " "
                   << m_dBottomFrameSize+(i*m_dYAxisNumbersSpacing)*m_dYScale << " moveto" << std::endl;
        outputFile << m_dTickMarkLength << " " << 0 << " rlineto" << std::endl;
        outputFile << m_dFrameLineWidth << " setlinewidth" << std::endl;
        outputFile << m_dFrameColor[0] << " " << m_dFrameColor[1] << " " << m_dFrameColor[2] << " setrgbcolor" << std::endl;
        outputFile << "stroke" << std::endl;
    }
}

void CPlot2D::DrawXAxisGridLinesPostScript()
{
    if (!m_bDrawGridLinesDashed) {
        for (int i=0;i<m_iXAxisNumberOfLabels;++i) {
            outputFile << m_dLeftFrameSize+(i*m_dXAxisNumbersSpacing)*m_dXScale << " "
                       << m_dBottomFrameSize << " moveto" << std::endl;
            outputFile << 0 << " " << m_dYAxisSize << " rlineto" << std::endl;
            outputFile << m_dGridLineWidth << " setlinewidth" << std::endl;
            outputFile << m_dGridColor[0] << " " << m_dGridColor[1] << " " << m_dGridColor[2] << " setrgbcolor" << std::endl;
            outputFile << "stroke" << std::endl;
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
                       << m_dBottomFrameSize << " moveto" << std::endl;

            for (int k=0;k<numDashes+1;k+=2) {

                outputFile << m_dLeftFrameSize+(i*m_dXAxisNumbersSpacing)*m_dXScale << " "
                           << m_dBottomFrameSize+k*delY << " moveto" << std::endl;
                outputFile << 0 << " " << delY << " rlineto" << std::endl;

                outputFile << m_dGridLineWidth << " setlinewidth" << std::endl;
                outputFile << m_dGridColor[0] << " " << m_dGridColor[1] << " " << m_dGridColor[2] << " setrgbcolor" << std::endl;
                outputFile << "stroke" << std::endl;
            }
        }
    }
}

void CPlot2D::DrawYAxisGridLinesPostScript()
{
    if (!m_bDrawGridLinesDashed) {
        for (int i=0;i<m_iYAxisNumberOfLabels;++i) {
            outputFile << m_dLeftFrameSize << " " << m_dBottomFrameSize+(i*m_dYAxisNumbersSpacing)*m_dYScale << " moveto" << std::endl;
            outputFile << m_dXAxisSize << " " << 0 << " rlineto" << std::endl;
            outputFile << m_dGridLineWidth << " setlinewidth" << std::endl;
            outputFile << m_dGridColor[0] << " " << m_dGridColor[1] << " " << m_dGridColor[2] << " setrgbcolor" << std::endl;
            outputFile << "stroke" << std::endl;
        }
    }
    else {

        int numDashes=floor(m_dXAxisSize/m_dLineDotSpacing);

        if (numDashes%2==0) {
            ++numDashes;
        }

        double delX=m_dXAxisSize/numDashes;

        for (int i=0;i<m_iYAxisNumberOfLabels;++i) {

            outputFile << m_dLeftFrameSize << " " << m_dBottomFrameSize+(i*m_dYAxisNumbersSpacing)*m_dYScale << " moveto" << std::endl;

            for (int k=0;k<numDashes+1;k+=2) {

                outputFile << m_dLeftFrameSize+k*delX  << " " << m_dBottomFrameSize+(i*m_dYAxisNumbersSpacing)*m_dYScale << " moveto" << std::endl;
                outputFile << delX << " " << 0 << " rlineto" << std::endl;

                outputFile << m_dGridLineWidth << " setlinewidth" << std::endl;
                outputFile << m_dGridColor[0] << " " << m_dGridColor[1] << " " << m_dGridColor[2] << " setrgbcolor" << std::endl;
                outputFile << "stroke" << std::endl;
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
        sprintf(m_cXAxisLabelFormat,"%%.%df",nfrac);
        char temp[20];
        m_iXAxisNumberOfLabels=0;
        for (double x=*plotMin; x<*plotMax+.5*d; x+=d) {
            sprintf(temp,m_cXAxisLabelFormat,x);
            m_strXAxisLabels.push_back(temp);
            m_iXAxisNumberOfLabels++;
        }
    }
    else if (axis=="y") {
        sprintf(m_cYAxisLabelFormat,"%%.%df",nfrac);
        char temp[20];
        m_iYAxisNumberOfLabels=0;
        for (double x=*plotMin; x<*plotMax+.5*d; x+=d) {
            sprintf(temp,m_cYAxisLabelFormat,x);
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

        outputFile << "/" << m_strXAxisLabelFont << " findfont" << std::endl;
        outputFile << m_dXAxisLabelFontSize << " scalefont" << std::endl;
        outputFile << "setfont" << std::endl;
        outputFile << labelXCoordinate << " " << labelYCoordinate << " moveto" << std::endl;
        outputFile << m_dXAxisLabelColor[0] << " " << m_dXAxisLabelColor[1] << " " << m_dXAxisLabelColor[2] << " setrgbcolor" << std::endl;

        // let PostScript handle the final adjustment based on the width of the string
        outputFile << "(" << m_strXAxisLabels[i] << ")" << " dup stringwidth pop 2 div neg 0 rmoveto show"  << std::endl;
    }
}

void CPlot2D::DrawYAxisLabelsPostScript()
{
    double labelXCoordinate,labelYCoordinate;

    for (int i=0;i<m_iYAxisNumberOfLabels;++i) {

        labelXCoordinate=m_dLeftFrameSize;
        labelYCoordinate=m_dBottomFrameSize+(i*m_dYAxisNumbersSpacing)*m_dYScale;

        // adjustment for the label
        labelXCoordinate-=m_dYAxisLabelFontSize*0.5;
        labelYCoordinate-=m_dYAxisLabelFontSize*0.25;

        outputFile << "/" << m_strYAxisLabelFont << " findfont" << std::endl;
        outputFile << m_dYAxisLabelFontSize << " scalefont" << std::endl;
        outputFile << "setfont" << std::endl;
        outputFile << labelXCoordinate << " " << labelYCoordinate << " moveto" << std::endl;
        outputFile << m_dYAxisLabelColor[0] << " " << m_dYAxisLabelColor[1] << " " << m_dYAxisLabelColor[2] << " setrgbcolor" << std::endl;

        // let PostScript handle the final adjustment based on the width of the string
        outputFile << "(" << m_strYAxisLabels[i] << ")" << " dup stringwidth pop neg 0 rmoveto show"  << std::endl;
    }
}

void CPlot2D::DrawXAxisTitlePostScript()
{
    double labelXCoordinate,labelYCoordinate;

    labelXCoordinate=m_dLeftFrameSize+m_dXAxisSize*0.5;
    labelYCoordinate=m_dBottomFrameSize-m_dXAxisTitleFontSize*2.0;

    outputFile << "/" << m_strXAxisTitleFont << " findfont" << std::endl;
    outputFile << m_dYAxisTitleFontSize << " scalefont" << std::endl;
    outputFile << "setfont" << std::endl;
    outputFile << labelXCoordinate << " " << labelYCoordinate << " moveto" << std::endl;
    outputFile << m_dXAxisTitleColor[0] << " " << m_dXAxisTitleColor[1] << " " << m_dXAxisTitleColor[2] << " setrgbcolor" << std::endl;

    // let PostScript handle the final adjustment based on the width of the string
    outputFile << "(" << m_strXAxisTitle << ")" << " dup stringwidth pop 2 div neg 0 rmoveto show"  << std::endl;
}

void CPlot2D::DrawYAxisTitlePostScript()
{
    double labelXCoordinate,labelYCoordinate;

    labelXCoordinate=10.0+m_dYAxisTitleFontSize;
    labelYCoordinate=m_dBottomFrameSize+m_dYAxisSize*0.5;

    outputFile << "/" << m_strYAxisTitleFont << " findfont" << std::endl;
    outputFile << m_dYAxisTitleFontSize << " scalefont" << std::endl;
    outputFile << "setfont" << std::endl;
    outputFile << labelXCoordinate << " " << labelYCoordinate << " moveto" << std::endl;
    outputFile << "90 rotate" << std::endl;
    outputFile << m_dYAxisTitleColor[0] << " " << m_dYAxisTitleColor[1] << " " << m_dYAxisTitleColor[2] << " setrgbcolor" << std::endl;

    // let PostScript handle the final adjustment based on the width of the string
    outputFile << "(" << m_strYAxisTitle << ")" << " dup stringwidth pop 2 div neg 0 rmoveto show"  << std::endl;
    outputFile << "-90 rotate" << std::endl;
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


        if (m_dataSets[i].GetDrawLine()) {
            // draw the line
            if (!m_dataSets[i].GetDashedLine()) {
                legendXCoordinate=m_dXTotalSize*0.1+i*widthPerLegend+symbolWidth*0.25;
                legendYCoordinate=m_dBottomFrameSize*0.25+maxFontHeight*0.3;
                outputFile << legendXCoordinate << " " << legendYCoordinate << " moveto" << std::endl;
                outputFile << symbolWidth*0.5 << " " << 0 << " rlineto" << std::endl;
                outputFile << m_dataSets[i].GetLineWidth() << " setlinewidth" << std::endl;
                outputFile << r << " " << g << " " << b << " setrgbcolor" << std::endl;
                outputFile << "stroke" << std::endl;
            }
            else {
                if (m_dataSets[i].GetDashedLinePattern()=="dot") {

                    for (int j=0;j<5;++j) {

                        legendXCoordinate=m_dXTotalSize*0.1+i*widthPerLegend+symbolWidth*0.25+j*(symbolWidth*0.5/4.0);
                        legendYCoordinate=m_dBottomFrameSize*0.25+maxFontHeight*0.3;

                        outputFile << legendXCoordinate << " " << legendYCoordinate << " "
                        << m_dataSets[i].GetLineWidth() << " 0 360 arc closepath" << std::endl;

                        outputFile << r << " " << g << " " << b << " setrgbcolor" << std::endl;
                        outputFile << "fill" << std::endl;
                        outputFile << "stroke" << std::endl;
                    }
                }

                if (m_dataSets[i].GetDashedLinePattern()=="dash") {

                    for (int j=0;j<5;j+=2) {

                        legendXCoordinate=m_dXTotalSize*0.1+i*widthPerLegend+symbolWidth*0.25+j*(symbolWidth*0.5/5.0);
                        legendYCoordinate=m_dBottomFrameSize*0.25+maxFontHeight*0.3;

                        outputFile << legendXCoordinate << " " << legendYCoordinate << " moveto" << std::endl;

                        outputFile << symbolWidth*0.5/5.0 << " " << 0 << " rlineto" << std::endl;
                        outputFile << r << " " << g << " " << b << " setrgbcolor" << std::endl;
                        outputFile << "stroke" << std::endl;
                    }
                }
                if (m_dataSets[i].GetDashedLinePattern()=="dash_dot") {

                    legendXCoordinate=m_dXTotalSize*0.1+i*widthPerLegend+symbolWidth*0.25;
                    legendYCoordinate=m_dBottomFrameSize*0.25+maxFontHeight*0.3;


                    outputFile << legendXCoordinate << " " << legendYCoordinate << " moveto" << std::endl;
                    outputFile << symbolWidth*0.5/2.0 << " " << 0 << " rlineto" << std::endl;
                    outputFile << r << " " << g << " " << b << " setrgbcolor" << std::endl;
                    outputFile << "stroke" << std::endl;

                    legendXCoordinate=m_dXTotalSize*0.1+i*widthPerLegend+symbolWidth*0.75;
                    legendYCoordinate=m_dBottomFrameSize*0.25+maxFontHeight*0.3;

                    outputFile << legendXCoordinate << " " << legendYCoordinate << " "
                    << m_dataSets[i].GetLineWidth()*0.5 << " 0 360 arc closepath" << std::endl;
                    outputFile << "fill" << std::endl;
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

        outputFile << legendXCoordinate << " " << legendYCoordinate << " moveto" << std::endl;

        outputFile << "/" << m_dataSets[i].GetDatasetLegendFont() << " findfont" << std::endl;
        outputFile << maxFontHeight << " scalefont" << std::endl;
        outputFile << "setfont" << std::endl;
        outputFile << legendXCoordinate << " " << legendYCoordinate << " moveto" << std::endl;

        outputFile << 0.0 << " " << 0.0 << " " << 0.0 << " setrgbcolor" << std::endl;

        outputFile << "(" <<  m_dataSets[i].GetDatasetTitle() << ")" << " show"  << std::endl;
    }

}
