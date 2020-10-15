#include <src/CPlot2D.h>
#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/gravis/t2Matrix.h>
#include <src/jaz/gravis/t3Vector.h>
#include <src/jaz/gravis/tRGB.h>
#include <src/jaz/atomic/pdb_helper.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/math/tensor2x2.h>
#include <iostream>
#include <fstream>
#include <set>

using namespace gravis;

#define LINE_BUFFER_LEN 1024

bool isAtomInteresting(
        std::string atom_name, 
        bool main_chain_only,
        const std::set<std::string>& desired_elements)
{
	const std::string element = atom_name.substr(0,1);
	const std::string position_code = atom_name.substr(1);
	
	if (main_chain_only
	    && element != "S" 
		&& atom_name != "CA" 
		&& position_code.length() > 0) 
	{
		return false;
	}
	
	if (desired_elements.size() != 0 && desired_elements.find(element) == desired_elements.end())
	{
		return false;
	}
	
	return true;
}
        
int main(int argc, char *argv[])
{
	std::string fn_data = "processed_SNR1_sig1_res2_norot/amp_over_phase_sharp_all.dat";
	std::string fn_legend = "processed_SNR1_sig1_res2_norot/amp_over_phase_sharp_all_legend.txt";
	std::string fn_out = "DC_all_CNOS_rot5.eps";
	
	const bool main_chain_only = false;
	const bool rotate = true;
	const double angle = DEG2RAD(-5);
	
	std::set<std::string> desired_elements;
	desired_elements.insert("C");
	desired_elements.insert("O");
	desired_elements.insert("N");
	desired_elements.insert("S");
	
	std::ifstream ifs_data(fn_data);
	std::ifstream ifs_legend(fn_legend);
	
	const int line_buffer_length = 1024;
	char* data_line = new char[line_buffer_length];
	
	
	std::map<std::string, std::vector<d2Vector>> coordinates, averages, inner_ellipses, outer_ellipses;
	std::vector<std::string> atom_names;
	
	while (ifs_legend)
	{
		std::string atom_name, atom_count_string;
		ifs_legend >> atom_name;
		ifs_legend >> atom_count_string;
		
		atom_names.push_back(atom_name);
		
		for (int figure = 0; figure < 4; figure++)
		{
			while (true)
			{
				ifs_data.getline(data_line, line_buffer_length);
				
				std::string data_string(data_line);
				
				if (data_string.length() == 0)
				{
					break;
				}
				else
				{
					std::stringstream sts(data_string);
					
					d2Vector c;
					sts >> c.x;
					sts >> c.y;
					
					c *= 1e8;
					
					if (rotate)
					{
						const d2Matrix A
							( cos(angle), -sin(angle),
							  sin(angle),  cos(angle) );
					
						c = A * c;
					}
					
					switch (figure)
					{
						case 0: 
							coordinates[atom_name].push_back(c);							
							break;
						case 1: 
							averages[atom_name].push_back(c);
							break;
						case 2: 
							inner_ellipses[atom_name].push_back(c);
							break;
						case 3: 
							outer_ellipses[atom_name].push_back(c);
							break;
					}
				}
			}
		}
	}
	
	std::map<std::string, dRGB> element_colours;
	element_colours["C"] = dRGB(0.2,0.2,0.2);
	element_colours["N"] = dRGB(0.0,0.3,1.0);
	element_colours["O"] = dRGB(1.0,0.2,0.0);
	element_colours["S"] = dRGB(0.8,0.8,0.0);
	
	CPlot2D plot2D("");
	
	//plot2D.SetTitle("amplitude over phase");
	plot2D.SetDrawLegend(true);
	
	const int anc = atom_names.size();
	
	std::map<std::string, d2Vector> element_means;
	std::map<std::string, int> element_atom_count;
	
	for (std::set<std::string>::iterator it = desired_elements.begin();
	     it != desired_elements.end(); it++)
	{	
		std::string element = *it;	
		element_means[element] = d2Vector(0,0);	
		element_atom_count[element] = 0;	
	}
	
	for (int an = 0; an < anc; an++)
	{
		const std::string atom_name = atom_names[an];
		
		if (atom_name == "") continue;
		
		const std::string element = atom_name.substr(0,1);
		
		if (!isAtomInteresting(atom_name, main_chain_only, desired_elements))
		{
			continue;
		}
		
		const int ac = coordinates[atom_name].size();
		
		for (int a = 0; a < ac; a++)
		{
			element_means[element] += coordinates[atom_name][a];
		}
		
		element_atom_count[element] += ac;		
	}
	
	std::map<std::string, Tensor2x2<double>> element_covariance_matrices;
	
	for (std::set<std::string>::iterator it = desired_elements.begin();
	     it != desired_elements.end(); it++)
	{	
		std::string element = *it;
		element_means[element] /= element_atom_count[element];
		element_covariance_matrices[element] = Tensor2x2<double>(0,0,0);
		
		std::cout << element << ": mean = " << element_means[element] << std::endl;
	}
	
	for (int an = 0; an < anc; an++)
	{
		const std::string atom_name = atom_names[an];
		
		if (atom_name == "") continue;
		
		const std::string element = atom_name.substr(0,1);		
		
		if (!isAtomInteresting(atom_name, main_chain_only, desired_elements))
		{
			continue;
		}
		
		const int ac = coordinates[atom_name].size();

		for (int a = 0; a < ac; a++)
		{
			const d2Vector d = coordinates[atom_name][a] - element_means[element];
			
			element_covariance_matrices[element].xx += d.x * d.x;
			element_covariance_matrices[element].xy += d.x * d.y;
			element_covariance_matrices[element].yy += d.y * d.y;
		}		
	}

	for (std::set<std::string>::iterator it = desired_elements.begin();
	     it != desired_elements.end(); it++)
	{	
		std::string element = *it;
		
		element_covariance_matrices[element] /= element_atom_count[element];
		
		std::cout << element << ": S = " << element_covariance_matrices[element] << std::endl;
		
		const dRGB colour = element_colours[element];
		
		CDataSet ellipse_points;
		ellipse_points.SetDrawMarker(false);
		ellipse_points.SetDatasetColor(colour.r, colour.g, colour.b);
		ellipse_points.SetDrawLine(true);
		ellipse_points.SetLineWidth(0.5);
		
		const d2Vector mu = element_means[element];
		
		const int ellipse_samples = 100;
		
		Tensor2x2<double> C = element_covariance_matrices[element];
		d2Vector evals, eig0, eig1;
		C.diagonalize(evals, eig0, eig1);
		
		eig0 *= sqrt(evals[0]);
		eig1 *= sqrt(evals[1]);
		
		std::cout << "  e0 = " << eig0 << std::endl;
		std::cout << "  e1 = " << eig1 << std::endl;
		
		for (int i = 0; i < ellipse_samples+1; i++)
		{
			const double phi = 2.0 * PI * i / (double) ellipse_samples;
			const double sp = sin(phi);
			const double cp = cos(phi);
			
			const d2Vector c = mu + 2.0 * (cp * eig0 + sp * eig1);
					
			ellipse_points.AddDataPoint(CDataPoint(c.x, c.y));
		}
		
		plot2D.AddDataSet(ellipse_points);
		
		
		CDataSet average_point;
		average_point.SetDrawMarker(true);
		average_point.SetMarkerSize(16);
		average_point.SetDatasetColor(colour.r, colour.g, colour.b);
		average_point.SetDrawLine(false);
		average_point.SetDrawMarkerFilled(true);
		
		average_point.AddDataPoint(CDataPoint(mu.x, mu.y));
		average_point.SetDatasetTitle(element);
		
		plot2D.AddDataSet(average_point);
		
		
		CDataSet line_to_origin_points;
		line_to_origin_points.SetDrawMarker(false);
		line_to_origin_points.SetDatasetColor(colour.r, colour.g, colour.b);
		line_to_origin_points.SetDrawLine(true);
		
		line_to_origin_points.AddDataPoint(CDataPoint(mu.x, mu.y));
		line_to_origin_points.AddDataPoint(CDataPoint(0,0));
		
		plot2D.AddDataSet(line_to_origin_points);
	}
	
	std::map<std::string, std::string> plot_symbol;
	
	plot_symbol["A"] = "o";
	plot_symbol["B"] = "x";
	plot_symbol["G"] = "+";
	plot_symbol["D"] = "diamond";
	plot_symbol["E"] = "square";
	plot_symbol["Z"] = "triangle";
	plot_symbol["H"] = "triangle";
	        
	for (int an = 0; an < anc; an++)
	{
		const std::string atom_name = atom_names[an];
		
		if (atom_name == "") continue;
		
		const std::string element = atom_name.substr(0,1);
		
		if (!isAtomInteresting(atom_name, main_chain_only, desired_elements))
		{
			continue;
		}
		
		std::cout << atom_name << ": " << element << std::endl;
		
		const dRGB colour = element_colours[element];
		const bool is_hollow = main_chain_only && (atom_name == "CA" || atom_name == "SD");
		
		const std::vector<d2Vector>& coords = coordinates[atom_name];
		const std::vector<d2Vector>& average = averages[atom_name];
		const std::vector<d2Vector>& inner_ellipse = inner_ellipses[atom_name];
		const std::vector<d2Vector>& outer_ellipse = outer_ellipses[atom_name];
		
		const int ac = coords.size();
		const int ec = inner_ellipse.size();
		
		CDataSet cloud_points;
		cloud_points.SetDrawMarker(true);
		cloud_points.SetMarkerSize(1);
		cloud_points.SetDatasetColor(colour.r, colour.g, colour.b);
		cloud_points.SetDrawLine(false);
		cloud_points.SetDrawMarkerFilled(!is_hollow);
		
		for (int a = 0; a < ac; a++)
		{
			const d2Vector c = coords[a];
			cloud_points.AddDataPoint(CDataPoint(c.x, c.y));
		}
		
		
		CDataSet average_point;
		average_point.SetDrawMarker(true);
		average_point.SetMarkerSize(8);
		
		if (atom_name.length() > 1)
		{
			average_point.SetMarkerSymbol(plot_symbol[atom_name.substr(1,1)]);
		}
		
		average_point.SetDatasetColor(colour.r, colour.g, colour.b);
		average_point.SetDrawLine(false);
		
		if (main_chain_only)
		{
			average_point.SetDatasetTitle(atom_name);
		}
		
		average_point.SetDrawMarkerFilled(!is_hollow);
		
		{
			const d2Vector c = average[0];
			average_point.AddDataPoint(CDataPoint(c.x, c.y));
		}
		
		
		CDataSet line_to_origin_points;
		line_to_origin_points.SetDrawMarker(false);
		line_to_origin_points.SetDatasetColor(colour.r, colour.g, colour.b);
		line_to_origin_points.SetDrawLine(true);
		
		if (main_chain_only)
		{
			const d2Vector c = average[0];
			line_to_origin_points.AddDataPoint(CDataPoint(c.x, c.y));
			line_to_origin_points.AddDataPoint(CDataPoint(0,0));
		}
		
		CDataSet inner_ellipse_points;
		inner_ellipse_points.SetDrawMarker(false);
		inner_ellipse_points.SetDatasetColor(colour.r, colour.g, colour.b);
		inner_ellipse_points.SetDrawLine(true);
		inner_ellipse_points.SetLineWidth(0.5);
		
		CDataSet outer_ellipse_points;
		outer_ellipse_points.SetDrawMarker(false);
		outer_ellipse_points.SetDatasetColor(colour.r, colour.g, colour.b);
		outer_ellipse_points.SetDrawLine(true);
		outer_ellipse_points.SetLineWidth(0.5);
		
		
		for (int e = 0; e < ec; e++)
		{
			const d2Vector ci = inner_ellipse[e];
			const d2Vector co = outer_ellipse[e];
			
			inner_ellipse_points.AddDataPoint(CDataPoint(ci.x, ci.y));
			outer_ellipse_points.AddDataPoint(CDataPoint(co.x, co.y));
		}	
		
		if (main_chain_only) 
		{	
			//plot2D.AddDataSet(outer_ellipse_points);
			plot2D.AddDataSet(line_to_origin_points);
		}
		
		plot2D.AddDataSet(cloud_points);		
		plot2D.AddDataSet(average_point);
	}
	
	
	if (rotate)
	{
		//plot2D.SetXAxisTitle("amplitude + "+ZIO::itoa(AC)+" * phase");
		//plot2D.SetYAxisTitle("amplitude - "+ZIO::itoa(AC)+" * phase");
	}
	else
	{
		plot2D.SetXAxisTitle("phase");
		plot2D.SetYAxisTitle("amplitude");
	}
	
	plot2D.OutputPostScriptPlot(fn_out);
		
	return 0;
}
