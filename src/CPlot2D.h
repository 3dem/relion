 //
 //  CPlot2D.h
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

#ifndef __CPlot2D__
#define __CPlot2D__

#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string.h>
#include <vector>
#include <float.h>
#include "src/filename.h"


/* SHWS: join multiple eps files into a single pdf
 *
 */
void joinMultipleEPSIntoSinglePDF(FileName fn_pdf, std::vector<FileName> fn_eps);

/* SHWS: concatenate multiple PDF files into a single one
 *
 */
bool concatenatePDFfiles(FileName fn_pdf_out, FileName pdf1, FileName pdf2);
bool concatenatePDFfiles(FileName fn_pdf_out, std::vector<FileName> fn_pdfs);

/* EL: Including all of the std namespace on the global level both in this and
 *     every single file that includes this header will lead to ambiguous
 *     definitions with Cuda. Fixed by adding std:: prefixes to string, vector,
 *     and ofstream.
 */
// using namespace std;

 /*!
  A simple container class to hold a data point, comprised of an x and y value
  stored in a double. Simple accessors are implemented to set and get the data
  stored.
  */
 class CDataPoint
 {
 public:
     /*!
      Constructor for the class. The x and y values are initialized to zero.
      */
     CDataPoint() {
         m_dDataX=0.0;
         m_dDataY=0.0;
     };
     /*!
      Nothing really to destruct...
      */
     ~CDataPoint() {};

     /*!
      Constructor for the class with initialization of x and y.
      */
     CDataPoint(double x, double y) {
         m_dDataX=x;
         m_dDataY=y;
     };

     /*!
      Another way to set the x and y values.
      */
     void SetValues(double x, double y) {
         m_dDataX=x;
         m_dDataY=y;
     };

     /*!
      Get the x and y values as a pair.
      */
     void GetValues(double *x, double *y) {
         *x=m_dDataX;
         *y=m_dDataY;
     };
     /*!
      Get the x value individually.
      */
     double GetX() {
         return (m_dDataX);
     };
     /*!
      Get the y value individually.
      */
     double GetY() {
         return (m_dDataY);
     };

 protected:

     double m_dDataX; /*!< Storage for the x value of a data point. */
     double m_dDataY; /*!< Storage for the y value of a data point. */
 };


 /*!
  A container class to hold a dataset. This class stores all the data points
  in a vector along with attributes specific to representing a dataset, such
  as the color of the line and markers or the title of the dataset. */
 class CDataSet
 {
 public:

     /*!
      Constructor for the class. The member variables are initialized to
      common values resulting in a, perhaps, pleasing representation of the
      data without the need to set anything up.
      */
     CDataSet() {
         m_dLineWidth=1.0;
         m_dColor[0]=0.0;
         m_dColor[1]=0.0;
         m_dColor[2]=0.0;
         m_strMarker="o";
         m_dMarkerSize=7.0;
         m_bDrawLine=true;
         m_bDrawMarker=true;
         m_bDrawMarkerFilled=true;
         m_bDashedLine=false;
         m_iDashedLinePattern="dash";
         m_strDatasetTitle="";
         m_strDatasetLegendFont="Times";
     };
     /*!
      Nothing really to destruct...
      */
     ~CDataSet() {};

     /*!
      Set a data point at the index location. This function
      neccessitates the presence of the data location pointed
      to by the index variable. Generally, the preferred method
      to set or add data is to use the AddDataPoint() function,
      which appends the new data point to the dataset.
      */
     void SetDataPoint(int index, CDataPoint point) {
         m_dDataPoints[index]=point;
     };

     /*!
      Get the data point at the location indicated by the variable
      index.
      */
     CDataPoint GetDataPoint(int index) {
         return (m_dDataPoints[index]);
     };

     /*!
      Add a new data point by appending it to the end of the data
      vector.
      */
     void AddDataPoint(CDataPoint point) {
         m_dDataPoints.push_back(point);
     };

     /*!
      Returns the number of data points comprisig the data set.
      */
     int GetNumberOfDataPointsInSet() {
         return ((int)m_dDataPoints.size());
     };

     /*!
      Returns the minimum value of the x variable in the dataset
      by iterating over all data points in the set.
      */
     double GetXMinValue() {
         double min=DBL_MAX;
         for (int i=0;i<m_dDataPoints.size();++i) {
             if (m_dDataPoints[i].GetX()<min) {
                 min=m_dDataPoints[i].GetX();
             }
         }
         return (min);
     };

     /*!
      Returns the maximum value of the x variable in the dataset
      by iterating over all data points in the set.
      */
     double GetXMaxValue() {
         double max=-DBL_MAX;
         for (int i=0;i<m_dDataPoints.size();++i) {
             if (m_dDataPoints[i].GetX()>max) {
                 max=m_dDataPoints[i].GetX();
             }
         }
         return (max);

     };

     /*!
      Returns the minimum value of the y variable in the dataset
      by iterating over all data points in the set.
      */
     double GetYMinValue() {
         double min=DBL_MAX;
         for (int i=0;i<m_dDataPoints.size();++i) {
             if (m_dDataPoints[i].GetY()<min) {
                 min=m_dDataPoints[i].GetY();
             }
         }
         return (min);
     };

     /*!
      Returns the maximum value of the y variable in the dataset
      by iterating over all data points in the set.
      */
     double GetYMaxValue() {
         double max=-DBL_MAX;
         for (int i=0;i<m_dDataPoints.size();++i) {
             if (m_dDataPoints[i].GetY()>max) {
                 max=m_dDataPoints[i].GetY();
             }
         }
         return (max);
     };

     /*!
      Returns the extent (maximum minus minimum) of the
      x values in the dataset.
      */
     double GetXExtent() {
         return (GetXMaxValue()-GetXMinValue());
     };

     /*!
      Returns the extent (maximum minus minimum) of the
      y values in the dataset.
      */
     double GetYExtent() {
         return (GetYMaxValue()-GetYMinValue());
     };

     /*!
      Sets the line width that will be used in drawing the
      line representing the dataset.
      */
     void SetLineWidth(double lineWidth)
     {
         m_dLineWidth=lineWidth;
     };

     /*!
      Returns the line width that will be used in drawing the
      line representing the dataset.
      */
     double GetLineWidth()
     {
         return (m_dLineWidth);
     };

     /*!
      Sets the RGB color used for drawing the dataset. The range
      is 0.0-1.0.
      */
     void SetDatasetColor(double r, double g, double b) {
         m_dColor[0]=r;
         m_dColor[1]=g;
         m_dColor[2]=b;
     };

     /*!
      Gets the RGB color used for drawing the dataset.
      */
     void GetDatasetColor(double *r, double *g, double *b) {
         *r=m_dColor[0];
         *g=m_dColor[1];
         *b=m_dColor[2];
     };

     /*!
      Sets the marker symbol used in drawing the dataset.
      See the declaration of variable for the available types.
      */
     void SetMarkerSymbol(std::string symbol)
     {
         m_strMarker=symbol;
     };

     /*!
      Sets the marker symbol size used in drawing the dataset.
      */
     void SetMarkerSize(double size)
     {
         m_dMarkerSize=size;
     };

     /*!
      A flag to enable/disable drawing a line connecting
      the data points.
      */
     void SetDrawLine(bool flag)
     {
         m_bDrawLine=flag;
     };

     /*!
      A flag to enable/disable drawing the marker symbol.
      */
     void SetDrawMarker(bool flag)
     {
         m_bDrawMarker=flag;
     };

     /*!
      A flag to enable/disable filling the interior of
      the marker symbol.
      */
     void SetDrawMarkerFilled(bool flag)
     {
         m_bDrawMarkerFilled=flag;
     };

     /*!
      Returns a string describing the marker symbol.
      */
     std::string GetMarkerSymbol()
     {
         return (m_strMarker);
     };

     /*!
      Returns the size of the marker symbol.
      */
     double GetMarkerSize()
     {
         return (m_dMarkerSize);
     };

     /*!
      Returns a boolean representing if the line spanning
      the data points will be drawn or not.
      */
     bool GetDrawLine()
     {
         return (m_bDrawLine);
     };

     /*!
      Returns a boolean representing if the marker symbol
      will be drawn or not.
      */
     bool GetDrawMarker()
     {
         return (m_bDrawMarker);
     };

     /*!
      Returns a boolean representing if filling the interior of the
      marker symbol is enabled or disabled.
      */
     bool GetDrawMarkerFilled()
     {
         return (m_bDrawMarkerFilled);
     };

     /*!
      Sets the dased line pattern style used in drawing the line
      spanning data points. For the availabel styles see the
      declaration of the variable.
      */
     void SetDashedLinePattern(std::string pattern)
     {
         m_iDashedLinePattern=pattern;
     };

     /*!
      Gets the dased line pattern style used in drawing the line
      spanning data points.
      */
     std::string GetDashedLinePattern()
     {
         return (m_iDashedLinePattern);
     };

     /*!
      Sets the drawing style for the lines spanning data points
      to be dashed.
      */
     void SetDashedLine(bool dashed)
     {
         m_bDashedLine=dashed;
     };

     /*!
      Gets a boolean representing if the dashed drawing style for
      the lines spanning data points is active or not.
      */
     bool GetDashedLine()
     {
         return (m_bDashedLine);
     };

     /*!
      Sets the string used in the plot legend denoting this
      dataset.
      */
     void SetDatasetTitle(std::string title)
     {
         m_strDatasetTitle=title;
     };

     /*!
      Gets the string used in the plot legend denoting this
      dataset.
      */
     std::string GetDatasetTitle()
     {
         return (m_strDatasetTitle);
     };

     /*!
      Sets the font used in the plot legend denoting this
      dataset.
      */
     void SetDatasetLegendFont(std::string font)
     {
         m_strDatasetLegendFont=font;
     };

     /*!
      Gets the font used in the plot legend denoting this
      dataset.
      */
     std::string GetDatasetLegendFont()
     {
         return (m_strDatasetLegendFont);
     }

 protected:

     std::string m_strDatasetTitle; /*!< Storage for the dataset's title. */
     std::string m_strDatasetLegendFont; /*!< Size of the font used to show the dataset's title in the legend. */

     double m_dLineWidth; /*!< The width of the line drawn to connect data points. */
     double m_dColor[3]; /*!< The color of the line drawn to connect data points. */
     std::string m_iDashedLinePattern; /*!< The dashed line pattern drawn to connect data points, the possibilites are (strings): "dot", "dash" or "dash_dot". */

     std::string m_strMarker; /*!< The style for the marker drawn at each data point, possibilities are (strings): "x", "o", "*", "diamond", "square", "triangle". */
     double m_dMarkerSize; /*!< The size of the marker drawn at each data point. */

     // flags
     bool m_bDrawLine; /*!< Boolean flag to enable/disable drawing of the line connecting data points. */
     bool m_bDrawMarker; /*!< Boolean flag to enable/disable drawing of the marker at the location of data points. */
     bool m_bDrawMarkerFilled; /*!< Boolean flag to enable/disable infilling of the marker at the location of data points. */
     bool m_bDashedLine; /*!< Boolean flag to enable/disable drawing of a dashed line connecting data points. */

     // data storage
     std::vector<CDataPoint> m_dDataPoints; /*!< A vector storage for the data points. */
 };

 /*!
  The class responsible for storing and rendering all the data sets it contains. Although the currently implemented
  method for rendering of a plot is into a PostScript file, it is envisioned that other rendering methods, such as OpenGL
  can be implemented as well. The accessors are created for both setting and retrieving plot parameters. */
 class CPlot2D
 {
 public:

     /*!
      The constructor for the class. Member variables are initialized to a set of values,
      which results in a, hopefully, pleasing plot without changing any of them.
      */
     CPlot2D(std::string title = "");

     /*!
      Nothing really to destruct...
      */
     ~CPlot2D();

     // accessors

     /*!
      Sets the total size of the plot in the x dimension.
      */
     void SetTitle(std::string);

    /*!
      Sets the total size of the plot in the x dimension.
      */
     void SetXTotalSize(double value);

     /*!
      Gets the total size of the plot in the x dimension.
      */
     double GetXTotalSize();

     /*!
      Sets the total size of the plot in the y dimension.
      */
     void SetYTotalSize(double value);

     /*!
      Gets the total size of the plot in the y dimension.
      */
     double GetYTotalSize();

     /*!
      Sets the x axis size of the plot.
      */
     void SetXAxisSize(double value);

     /*!
      Gets the x axis size of the plot.
      */
     double GetXAxisSize();

     /*!
      Sets the y axis size of the plot.
      */
     void SetYAxisSize(double value);

     /*!
      Gets the y axis size of the plot.
      */
     double GetYAxisSize();

     /*!
      Sets the bottom frame of the plot (e.g. the distance from the bottom of the image to the
      plot frame.
      */
     void SetBottomFrameSize(double value);

     /*!
      Gets the bottom frame of the plot (e.g. the distance from the bottom of the image to the
      plot frame.
      */
     double GetBottomFrameSize();

     /*!
      Sets the right frame of the plot (e.g. the distance from the right of the image to the
      plot frame.
      */
     void SetRightFrameSize(double value);

     /*!
      Gets the right frame of the plot (e.g. the distance from the right of the image to the
      plot frame.
      */
     double GetRightFrameSize();

     /*!
      Sets the top frame of the plot (e.g. the distance from the top of the image to the
      plot frame.
      */
     void SetTopFrameSize(double value);

     /*!
      Gets the top frame of the plot (e.g. the distance from the top of the image to the
      plot frame.
      */
     double GetTopFrameSize();

     /*!
      Sets the left frame of the plot (e.g. the distance from the left of the image to the
      plot frame.
      */
     void SetLeftFrameSize(double value);

     /*!
      Gets the left frame of the plot (e.g. the distance from the left of the image to the
      plot frame.
      */
     double GetLeftFrameSize();

     /*!
      Sets the thickness of the line that the plot frame is drawn with.
      */
     void SetFrameLineWidth(double value);

     /*!
      Gets the thickness of the line that the plot frame is drawn with.
      */
     double GetFrameLineWidth();

     /*!
      Sets the thickness of the line that the grid is drawn with.
      */
     void SetGridLineWidth(double value);

     /*!
      Gets the thickness of the line that the grid is drawn with.
      */
     double GetGridLineWidth();

     /*!
      Sets the color of the frame as an RGB triplet.
      */
     void SetFrameColor(double r, double g, double b);

     /*!
      Gets the color of the frame as an RGB triplet.
      */
     void GetFrameColor(double *r, double *g, double *b);

     /*!
      Sets the color of the grid as an RGB triplet.
      */
     void SetGridColor(double r, double g, double b);

     /*!
      Gets the color of the grid as an RGB triplet.
      */
     void GetGridColor(double *r, double *g, double *b);

     /*!
     Sets the spacing of tick marks and labels along the x axis.
      */
     void SetXAxisNumbersSpacing(double spacing);

     /*!
      Gets the spacing of tick marks and labels along the x axis.
      */
     double GetXAxisNumbersSpacing();

     /*!
      Sets the spacing of tick marks and labels along the y axis.
      */
     void SetYAxisNumbersSpacing(double spacing);

     /*!
      Gets the spacing of tick marks and labels along the y axis.
      */
     double GetYAxisNumbersSpacing();

     /*!
      Sets a flag to enable/disable drawing the x axis labels.
      */
     void SetDrawXAxisTickMarks(bool flag);

     /*!
      Gets the state of a flag which enables/disables drawing the x axis labels.
      */
     bool GetDrawXAxisTickMarks();

     /*!
      Sets a flag to enable/disable drawing the y axis labels.
      */
     void SetDrawYAxisTickMarks(bool flag);

     /*!
      Gets the state of a flag which enables/disables drawing the y axis labels.
      */
     bool GetDrawYAxisTickMarks();

     /*!
      Sets the number of tick marks and labels along the x axis.
      */
     void SetXAxisNumberOfTicks(int number);

     /*!
      Gets the number of tick marks and labels along the x axis.
      */
     int GetXAxisNumberOfTicks();

     /*!
      Sets the number of tick marks and labels along the y axis.
      */
     void SetYAxisNumberOfTicks(int number);

     /*!
      Gets the number of tick marks and labels along the y axis.
      */
     int GetYAxisNumberOfTicks();

     /*!
      Sets a flag that enables/disables the drawing of grid lines for the x axis.
      */
     void SetDrawXAxisGridLines(bool flag);

     /*!
      Gets the state of a flag that enables/disables the drawing of grid lines for the x axis.
      */
     bool GetDrawXAxisGridLines();

     /*!
      Sets a flag that enables/disables the drawing of grid lines for the y axis.
      */
     void SetDrawYAxisGridLines(bool flag);

     /*!
      Gets the state of a flag that enables/disables the drawing of grid lines for the y axis.
      */
     bool GetDrawYAxisGridLines();

     /*!
      Sets a flag that enables/disables the drawing of dashed grid lines.
      */
     void SetDrawGridLinesDashed(bool flag);

     /*!
      Gets the state of a flag that enables/disables the drawing of dashed grid lines.
      */
     bool GetDrawGridLinesDashed();

     /*!
      Sets the font (as a string) that is used for the labels on the x axis.
      */
     void SetXAxisLabelFont(std::string font);

     /*!
      Gets the font (as a string) that is used for the labels on the x axis.
      */
     std::string GetXAxisLabelFont();

     /*!
      Sets the font size for the x axis labels.
      */
     void SetXAxisLabelFontSize(double value);

     /*!
      Gets the font size for the x axis labels.
      */
     double GetXAxisLabelFontSize();

     /*!
      Sets the font (as a string) that is used for the labels on the y axis.
      */
     void SetYAxisLabelFont(std::string font);

     /*!
      Gets the font (as a string) that is used for the labels on the y axis.
      */
     std::string GetYAxisLabelFont();

     /*!
      Set the font size for the y axis labels.
      */
     void SetYAxisLabelFontSize(double value);

     /*!
      Gets the font size for the y axis labels.
      */
     double GetYAxisLabelFontSize();

     /*!
      Sets the font (as a string) that is used for the title on the legend for the x axis.
      */
     void SetXAxisTitleFont(std::string font);

     /*!
      Gets the font (as a string) that is used for the title on the legend for the x axis.
      */
     std::string GetXAxisTitleFont();

     /*!
      Sets the font size that is used for the title on the legend for the x axis.
      */
     void SetXAxisTitleFontSize(double value);

     /*!
      Gets the font size that is used for the title on the legend for the x axis.
      */
     double GetXAxisTitleFontSize();

     /*!
      Sets the font (as a string) that is used for the title on the legend for the y axis.
      */
     void SetYAxisTitleFont(std::string font);

     /*!
      Gets the font (as a string) that is used for the title on the legend for the x axis.
      */
     std::string GetYAxisTitleFont();

     /*!
      Sets the font size that is used for the title on the legend for the y axis.
      */
     void SetYAxisTitleFontSize(double value);

     /*!
      Gets the font size that is used for the title on the legend for the y axis.
      */
     double GetYAxisTitleFontSize();

     /*!
      Sets the title for the x axis. Usually this is where the units or quantity represented by
      the data should be displayed goes.
      */
     void SetXAxisTitle(std::string title);

     /*!
      Gets the title for the x axis.
      */
     std::string GetXAxisTitle();

     /*!
      Sets the title for the y axis. Usually this is where the units or quantity represented by
      the data should be displayed goes.
      */
     void SetYAxisTitle(std::string title);

     /*!
      Gets the title for the y axis.
      */
     std::string GetYAxisTitle();

     /*!
      Sets the color, as an RGB triplet, for the x axis title.
      */
     void SetXAxisTitleColor(double r, double g, double b);

     /*!
      Gets the color, as an RGB triplet, for the x axis title.
      */
     void GetXAxisTitleColor(double *r, double *g, double *b);

     /*!
      Sets the color, as an RGB triplet, for the y axis title.
      */
     void SetYAxisTitleColor(double r, double g, double b);

     /*!
      Gets the color, as an RGB triplet, for the y axis title.
      */
     void GetYAxisTitleColor(double *r, double *g, double *b);

     /*!
      Sets the color, as an RGB triplet, for the x axis labels.
      */
     void SetXAxisLabelColor(double r, double g, double b);

     /*!
      Gets the color, as an RGB triplet, for the x axis labels.
      */
     void GetXAxisLabelColor(double *r, double *g, double *b);

     /*!
      Sets the color, as an RGB triplet, for the y axis labels.
      */
     void SetYAxisLabelColor(double r, double g, double b);

     /*!
      Gets the color, as an RGB triplet, for the y axis labels.
      */
     void GetYAxisLabelColor(double *r, double *g, double *b);

     /*!
      Set a flag to enable/disable drawing the legend.
      */
     void SetDrawLegend(bool flag);

     /*!
      Get the state of a flag that enables/disables drawing the legend.
      */
     bool GetDrawLegend();

     /*!
      Get and Set a flag that flips the orientation of the Y axis.
     */
     bool GetFlipY();
     void SetFlipY(bool flag);

     // outputs
     /*!
      The function, which is responsible for generating the PostScript output of the plot.
      */
     void OutputPostScriptPlot(std::string fileName);

     // data set functions
     /*!
      Adds a new dataset to the plot.
      */
     void AddDataSet(CDataSet dataSet);

     /*!
      Adds a new dataset to the plot as two arrays comprised of the x and y values. The length of the
      array is supplied as well.
      */
     void AddDataSet(int numPoints, double *xValues, double *yValues);

     /*!
      Adds a new dataset to the plot as two vectors comprised of the x and y values
      */
     void AddDataSet(std::vector<RFLOAT> xValues, std::vector<RFLOAT> yValues);
     void AddDataSet(std::vector<RFLOAT> yValues);

	 void SetViewArea(double origin_x, double origin_y, double width, double height);


 protected:

     /*!
      A function, which precomputes parameters of the plot, such as the overal size of the plot,
      the spacing of tick marks and labels, the length of the dashes in a dashed line.
      */
     void PrecomputeDimensions();

     /*!
      A function to compute label spacing, see the comments at the head of the implementation
      regarding the source of the algorithm.
      */
     void ComputeLabelTickSpacing(double dataMin, double dataMax, double *plotMin,
                                  double *plotMax, double *tickSpacing, int numTicks, std::string axis);
     /*!
      A function to aid the computation of label spacing, see the comments at the head of the implementation
      regarding the source of the algorithm.
      */
     double NiceNum(double x, int round);

     // outputs

     /*!
      A function to draw the frame surrounding the plot in a PostScript format.
      */
     void DrawFramePostScript();

     /*!
      A function to draw the data sets in a PostScript format.
      */
     void DrawDataPostScript();

     /*!
      A function to draw a single marker in a PostScript format.
      */
     void DrawMarker(std::string symbol, double size, bool filled, double xLocation, double yLocation, int dataSet);

     /*!
      A function to draw the x axis tick marks in a PostScript format.
      */
     void DrawXAxisTickMarksPostScript();

     /*!
      A function to draw the y axis tick marks in a PostScript format.
      */
     void DrawYAxisTickMarksPostScript();

     /*!
      A function to draw the x axis grid lines in a PostScript format.
      */
     void DrawXAxisGridLinesPostScript();

     /*!
      A function to draw the y axis grid lines in a PostScript format.
      */
     void DrawYAxisGridLinesPostScript();

     /*!
      A function to draw the x axis labels in a PostScript format.
      */
     void DrawXAxisLabelsPostScript();

     /*!
      A function to draw the y axis labels in a PostScript format.
      */
     void DrawYAxisLabelsPostScript();

     /*!
      A function to draw the x axis title in a PostScript format.
      */
     void DrawXAxisTitlePostScript();

     /*!
      A function to draw the y axis title in a PostScript format.
      */
     void DrawYAxisTitlePostScript();

     /*!
      A function to draw the legend in a PostScript format.
      */
     void DrawLegendPostScript();



 protected:

     // general plot sizes
     double m_dXTotalSize; /*!< Total size of the plot in the x direction. Including the frames around it. */
     double m_dYTotalSize; /*!< Total size of the plot in the y direction. Including the frames around it.*/

     double m_dXAxisSize; /*!< The size of the plot in the x direction. Not including the frames around it. */
     double m_dYAxisSize; /*!< The size of the plot in the y direction. Not including the frames around it. */
     double m_dDiagonalSize; /*!< The size of the plot along its diagonal direction. */

     double m_dFlipYOffset; /*!< Used to shift the Y origin when flipping the Y axis. */

     double m_dBottomFrameSize; /*!< The size of the bottom frame. */
     double m_dRightFrameSize; /*!< The size of the right frame. */
     double m_dTopFrameSize; /*!< The size of the top frame. */
     double m_dLeftFrameSize; /*!< The size of the left frame. */

     int m_iXAxisNumberOfTicks; /*!< Number of tick marks along the x axis. */
     int m_iYAxisNumberOfTicks; /*!< Number of tick marks along the y axis. */
     double m_dXAxisNumbersSpacing; /*!< Spacing between tick marks along the x axis. */
     double m_dYAxisNumbersSpacing; /*!< Spacing between tick marks along the y axis. */
     char m_cXAxisLabelFormat[20]; /*!< Format of labels along the x axis. */
     char m_cYAxisLabelFormat[20]; /*!< Format of labels along the y axis. */
     std::vector<std::string> m_strXAxisLabels; /*!< Labels along the x axis. */
     std::vector<std::string> m_strYAxisLabels; /*!< Labels along the y axis. */
     //Sjors Scheres 22mar2016: insert PlotTitle
     std::string m_strPlotTitle; /*!< Title of plot. */
     std::string m_strXAxisTitle; /*!< Title of x axis. */
     std::string m_strYAxisTitle; /*!< Title of y axis. */
     int m_iXAxisNumberOfLabels; /*!< Number of labels along the x axis. */
     int m_iYAxisNumberOfLabels; /*!< Number of labels along the y axis. */

     double m_dXScale; /*!< Scale along x axis that converts dataset values to the plot space. */
     double m_dYScale; /*!< Scale along yx axis that converts dataset values to the plot spac. */
     double m_dMaxXExtent; /*!< Maximum extent of x axis. */
     double m_dMaxYExtent; /*!< Maximum extent of y axis. */
     double m_dMinXStartPoint; /*!< Minimum starting point along x axis over all data sets. */
     double m_dMinYStartPoint; /*!< Minimum starting point along y axis over all data sets. */
     double m_dMaxXEndPoint; /*!< Maximum end point along x axis over all data sets. */
	 double m_dMaxYEndPoint; /*!< Maximum end point along y axis over all data sets. */

     // line widths
     double m_dFrameLineWidth; /*!< Width (thickness) of line used for the frame. */
     double m_dGridLineWidth; /*!< Width (thickness) of line used for the grid. */

     // colors
     double m_dFrameColor[3]; /*!< Frame color as an RGB triplet. */
     double m_dGridColor[3]; /*!< Grid color as an RGB triplet. */
     double m_dXAxisTitleColor[3]; /*!<  X axis title color as an RGB triplet. */
     double m_dYAxisTitleColor[3]; /*!< Y axis title color as an RGB triplet. */
     double m_dXAxisLabelColor[3]; /*!< X axis label color as an RGB triplet. */
     double m_dYAxisLabelColor[3]; /*!< Y axis label color as an RGB triplet. */

     // pleasing dot spacing, tick mark length
     double m_dLineDotSpacing; /*!< Spacing of lines and dots for a dashed line. */
     double m_dTickMarkLength; /*!< Length of tick marks. */

     // fonts
     std::string m_strXAxisLabelFont; /*!< Font for x axis labels. */
     double m_dXAxisLabelFontSize; /*!< Font size for x axis label. */
     std::string m_strYAxisLabelFont; /*!< Font for yx axis label. */
     double m_dYAxisLabelFontSize; /*!< Font size for y axis label. */
     std::string m_strXAxisTitleFont; /*!< Font for x axis title. */
     double m_dXAxisTitleFontSize; /*!< Font size for x axis title. */
     std::string m_strYAxisTitleFont; /*!< Font for y axis title. */
     double m_dYAxisTitleFontSize; /*!< Font size for y axis title. */

     // flags
     bool m_bDrawXAxisTickMarks; /*!< Flag for enabling/disabling the drawing of x axis tick marks. */
     bool m_bDrawYAxisTickMarks; /*!< Flag for enabling/disabling the drawing of y axis tick marks. */
     bool m_bDrawXAxisGridLines; /*!< Flag for enabling/disabling the drawing of x axis grid lines. */
     bool m_bDrawYAxisGridLines; /*!< Flag for enabling/disabling the drawing of y axis grid lines. */
     bool m_bDrawGridLinesDashed; /*!< Flag for enabling/disabling the drawing dashed grid lines. */
     bool m_bDrawLegend; /*!< Flag for enabling/disabling the drawing of the legend. */
     bool m_bFlipY; /*!< Flag for flipping the Y axis. */

     // output
     std::ofstream outputFile; /*!< The output stream. */

     // data storage
     std::vector<CDataSet> m_dataSets; /*!< Storage for the datasets, implemented as a vector. */

	 // JZ, 6-2020: allow the user to control the view area manually
	 bool m_bSizeSetExternally;
	 double m_dMinXStartPointOverride;
	 double m_dMinYStartPointOverride;
	 double m_dMaxXEndPointOverride;
	 double m_dMaxYEndPointOverride;
 };

 inline void CPlot2D::SetTitle(std::string title)
 {
	 m_strPlotTitle=title;
 }

 inline void CPlot2D::SetXTotalSize(double value)
 {
     m_dXTotalSize=value;
 }

 inline void CPlot2D::SetYTotalSize(double value)
 {
     m_dYTotalSize=value;
 }

 inline void CPlot2D::SetXAxisSize(double value)
 {
     m_dXAxisSize=value;
 }

 inline void CPlot2D::SetYAxisSize(double value)
 {
     m_dYAxisSize=value;
 }

 inline void CPlot2D::SetBottomFrameSize(double value)
 {
     m_dBottomFrameSize=value;
 }

 inline void CPlot2D::SetRightFrameSize(double value)
 {
     m_dRightFrameSize=value;
 }

 inline void CPlot2D::SetTopFrameSize(double value)
 {
     m_dTopFrameSize=value;
 }

 inline void CPlot2D::SetLeftFrameSize(double value)
 {
     m_dLeftFrameSize=value;
 }

 inline double CPlot2D::GetXTotalSize()
 {
     return (m_dXTotalSize);
 }

 inline double CPlot2D::GetYTotalSize()
 {
     return (m_dYTotalSize);
 }

 inline double CPlot2D::GetXAxisSize()
 {
     return (m_dXAxisSize);
 }

 inline double CPlot2D::GetYAxisSize()
 {
     return (m_dYAxisSize);
 }

 inline double CPlot2D::GetBottomFrameSize()
 {
     return (m_dBottomFrameSize);
 }

 inline double CPlot2D::GetRightFrameSize()
 {
     return (m_dRightFrameSize);
 }

 inline double CPlot2D::GetTopFrameSize()
 {
     return (m_dTopFrameSize);
 }

 inline double CPlot2D::GetLeftFrameSize()
 {
     return (m_dLeftFrameSize);
 }

 inline void CPlot2D::SetFrameLineWidth(double value)
 {
     m_dFrameLineWidth=value;
 }

 inline double CPlot2D::GetFrameLineWidth()
 {
     return (m_dFrameLineWidth);
 }

 inline void CPlot2D::SetGridLineWidth(double value)
 {
     m_dGridLineWidth=value;
 }

 inline double CPlot2D::GetGridLineWidth()
 {
     return (m_dGridLineWidth);
 }

 inline void CPlot2D::SetFrameColor(double r, double g, double b)
 {
     m_dFrameColor[0]=r;
     m_dFrameColor[1]=g;
     m_dFrameColor[2]=b;
 }

 inline void CPlot2D::GetFrameColor(double *r, double *g, double *b)
 {
     *r=m_dFrameColor[0];
     *g=m_dFrameColor[1];
     *b=m_dFrameColor[2];
 }

 inline void CPlot2D::SetGridColor(double r, double g, double b)
 {
     m_dGridColor[0]=r;
     m_dGridColor[1]=g;
     m_dGridColor[2]=b;
 }

 inline void CPlot2D::GetGridColor(double *r, double *g, double *b)
 {
     *r=m_dGridColor[0];
     *g=m_dGridColor[1];
     *b=m_dGridColor[2];
 }

 inline void CPlot2D::AddDataSet(CDataSet dataSet)
 {
     m_dataSets.push_back(dataSet);
 }

 inline void CPlot2D::SetXAxisNumbersSpacing(double spacing)
 {
     m_dXAxisNumbersSpacing=spacing;
 }
 inline void CPlot2D::SetYAxisNumbersSpacing(double spacing)
 {
     m_dYAxisNumbersSpacing=spacing;
 }
 inline void CPlot2D::SetDrawXAxisTickMarks(bool flag)
 {
     m_bDrawXAxisTickMarks=flag;
 }
 inline void CPlot2D::SetDrawYAxisTickMarks(bool flag)
 {
     m_bDrawYAxisTickMarks=flag;
 }

 inline void CPlot2D::SetXAxisNumberOfTicks(int number)
 {
     m_iXAxisNumberOfTicks=number;
 }

 inline void CPlot2D::SetYAxisNumberOfTicks(int number)
 {
     m_iYAxisNumberOfTicks=number;
 }

 inline void CPlot2D::SetDrawXAxisGridLines(bool flag)
 {
     m_bDrawXAxisGridLines=flag;
 }

 inline void CPlot2D::SetDrawYAxisGridLines(bool flag)
 {
     m_bDrawYAxisGridLines=flag;
 }

 inline void CPlot2D::SetDrawGridLinesDashed(bool flag)
 {
     m_bDrawGridLinesDashed=flag;
 }

 inline double CPlot2D::GetXAxisNumbersSpacing()
 {
     return (m_dXAxisNumbersSpacing);
 }
 inline double CPlot2D::GetYAxisNumbersSpacing()
 {
     return (m_dYAxisNumbersSpacing);
 }
 inline bool CPlot2D::GetDrawXAxisTickMarks()
 {
     return (m_bDrawXAxisTickMarks);
 }
 inline bool CPlot2D::GetDrawYAxisTickMarks()
 {
     return (m_bDrawYAxisTickMarks);
 }

 inline int CPlot2D::GetXAxisNumberOfTicks()
 {
     return (m_iXAxisNumberOfTicks);
 }

 inline int CPlot2D::GetYAxisNumberOfTicks()
 {
     return (m_iYAxisNumberOfTicks);
 }

 inline bool CPlot2D::GetDrawXAxisGridLines()
 {
     return (m_bDrawXAxisGridLines);
 }

 inline bool CPlot2D::GetDrawYAxisGridLines()
 {
     return (m_bDrawYAxisGridLines);
 }

 inline bool CPlot2D::GetDrawGridLinesDashed()
 {
     return (m_bDrawGridLinesDashed);
 }

 inline void CPlot2D::SetXAxisLabelFont(std::string font)
 {
     m_strXAxisLabelFont=font;
 }

 inline std::string CPlot2D::GetXAxisLabelFont()
 {
     return (m_strXAxisLabelFont);
 }

 inline void CPlot2D::SetXAxisLabelFontSize(double value)
 {
     m_dXAxisLabelFontSize=value;
 }

 inline double CPlot2D::GetXAxisLabelFontSize()
 {
     return (m_dXAxisLabelFontSize);
 }

 inline void CPlot2D::SetYAxisLabelFont(std::string font)
 {
     m_strYAxisLabelFont=font;
 }

 inline std::string CPlot2D::GetYAxisLabelFont()
 {
     return (m_strYAxisLabelFont);
 }

 inline void CPlot2D::SetYAxisLabelFontSize(double value)
 {
     m_dYAxisLabelFontSize=value;
 }

 inline double CPlot2D::GetYAxisLabelFontSize()
 {
     return (m_dYAxisLabelFontSize);
 }

 inline void CPlot2D::SetXAxisTitleFont(std::string font)
 {
     m_strXAxisTitleFont=font;
 }

 inline std::string CPlot2D::GetXAxisTitleFont()
 {
     return (m_strXAxisTitleFont);
 }

 inline void CPlot2D::SetXAxisTitleFontSize(double value)
 {
     m_dXAxisTitleFontSize=value;
 }

 inline double CPlot2D::GetXAxisTitleFontSize()
 {
     return (m_dXAxisTitleFontSize);
 }

 inline void CPlot2D::SetYAxisTitleFont(std::string font)
 {
     m_strYAxisTitleFont=font;
 }

 inline std::string CPlot2D::GetYAxisTitleFont()
 {
     return (m_strYAxisTitleFont);
 }

 inline void CPlot2D::SetYAxisTitleFontSize(double value)
 {
     m_dYAxisTitleFontSize=value;
 }

 inline double CPlot2D::GetYAxisTitleFontSize()
 {
     return (m_dYAxisTitleFontSize);
 }

 inline void CPlot2D::SetXAxisTitle(std::string title)
 {
     m_strXAxisTitle=title;
 }

 inline std::string CPlot2D::GetXAxisTitle()
 {
     return (m_strXAxisTitle);
 }

 inline void CPlot2D::SetYAxisTitle(std::string title)
 {
     m_strYAxisTitle=title;
 }

 inline std::string CPlot2D::GetYAxisTitle()
 {
     return (m_strYAxisTitle);
 }

 inline void CPlot2D::SetXAxisTitleColor(double r, double g, double b)
 {
     m_dXAxisTitleColor[0]=r;
     m_dXAxisTitleColor[1]=g;
     m_dXAxisTitleColor[2]=b;
 }

 inline void CPlot2D::GetXAxisTitleColor(double *r, double *g, double *b)
 {
     *r=m_dXAxisTitleColor[0];
     *g=m_dXAxisTitleColor[1];
     *b=m_dXAxisTitleColor[2];
 }

 inline void CPlot2D::SetYAxisTitleColor(double r, double g, double b)
 {
     m_dYAxisTitleColor[0]=r;
     m_dYAxisTitleColor[1]=g;
     m_dYAxisTitleColor[2]=b;
 }

 inline void CPlot2D::GetYAxisTitleColor(double *r, double *g, double *b)
 {
     *r=m_dYAxisTitleColor[0];
     *g=m_dYAxisTitleColor[1];
     *b=m_dYAxisTitleColor[2];
 }

 inline void CPlot2D::SetXAxisLabelColor(double r, double g, double b)
 {
     m_dXAxisLabelColor[0]=r;
     m_dXAxisLabelColor[1]=g;
     m_dXAxisLabelColor[2]=b;
 }

 inline void CPlot2D::GetXAxisLabelColor(double *r, double *g, double *b)
 {
     *r=m_dXAxisLabelColor[0];
     *g=m_dXAxisLabelColor[1];
     *b=m_dXAxisLabelColor[2];
 }

 inline void CPlot2D::SetYAxisLabelColor(double r, double g, double b)
 {
     m_dYAxisLabelColor[0]=r;
     m_dYAxisLabelColor[1]=g;
     m_dYAxisLabelColor[2]=b;
 }

 inline void CPlot2D::GetYAxisLabelColor(double *r, double *g, double *b)
 {
     *r=m_dYAxisLabelColor[0];
     *g=m_dYAxisLabelColor[1];
     *b=m_dYAxisLabelColor[2];
 }

 inline void CPlot2D::SetDrawLegend(bool flag)
 {
     m_bDrawLegend=flag;
 }

 inline bool CPlot2D::GetDrawLegend()
 {
     return (m_bDrawLegend);
 }

 inline void CPlot2D::SetFlipY(bool flag)
 {
     m_bFlipY=flag;
 }

 inline bool CPlot2D::GetFlipY()
 {
     return (m_bFlipY);
 }

 #endif /* defined(__CPlot2D__) */
