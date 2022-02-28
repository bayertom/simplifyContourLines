#define _CRT_SECURE_NO_WARNINGS
#define _SILENCE_ALL_CXX17_DEPRECATION_WARNINGS

#include <vector>
#include <set>
#include <map>
#include <fstream>
#include <sstream>
#include <iostream>
#include <memory>
#include <format>
#include <algorithm>
#include <string>

#include <string_view>

#include "Exception.h"
#include "TVector.h"
#include "TVector2D.h"
#include "File.h"
#include "ContourLinesSimplify.h"
#include "DXFExport.h"
//#include "Point3D.h"



int main(int argc, char* argv[])
{

	//Initial parameters of the contour lines and the simplification
	bool var_buffer = false, var_delta = false, generalize = true;
	int method = 4, smooth = 0, max_iter = 3000, nc = 20;
	double z_min = 110.0, z_max = 650.0, z_step = 1.0;
	double dist_threshold = 5.0, dz_threshold = 0.10, ang_threshold = 10.0, dens = 0.0;
	double alpha = 0.00001, beta = 0.00001, gamma = 15, delta = 0.1, kappa = 1.0;

	//Path to the folder
	std::string path = "E:\\Tomas\\CPP\\SimplifyContourLines\\SimplifyContourLines\\tests\\ALS\\ZM5\\514655_DIP\\DMR5G_vrstevnice01m\\MIMO\\MIMO62\\";

	//Countour line file mask
	//std::string contours_file_mask = "*v279_smooth*.txt";
	std::string contours_file_mask = "*vrstevnice*final*CL*27*.csv";

	//Buffer file mask
	//std::string buff1_file_mask = "*_B1_0_1*.csv";
	//std::string buff2_file_mask = "*_B2_0_1*.csv";
	std::string buff1_file_mask = "*_B1_0_1*27*.csv";
	std::string buff2_file_mask = "*_B2_0_1*27*.csv";
	std::string output_file_name = "MIMO62_5g.xyz";

	//Process command-line argument:
	while (--argc > 0)
	{
		//Get - (A parameter follows)
		if (*argv[argc] == '-')
		{
			//Process parameter after -
			for (char* parameter = argv[argc]; ; )
			{
				switch (*(++parameter))
				{
					//Set generalization phse to true
				case 'g':
				{
					generalize = true;
					break;
				}

				//Set variable vertical buffer to true
				case 'v':
				{
					var_delta = true;
					break;
				}

				//Terminate character \0 of the argument
				case '\0':
					break;

					//Throw exception
				default:
					throw Exception("Exception: Invalid parameter in command line!");
				}

				//Stop processing of the argument
				break;
			}
		}

		//Set values
		else if (*argv[argc] == '+')
		{
			//Get new command line parameter
			char* attribute = const_cast <char*> (argv[argc] + 1), * value = NULL;
			char* command = attribute;

			//Find splitter: =
			for (; (*command != '=') && (*command != '\0'); command++);

			//We found splitter, trim command and copy to value
			if ((*command == '=') && (*command != '\0'))
			{
				*command = '\0';
				value = command + 1;
			}

			//Throw exception
			if (attribute == NULL || value == NULL)
				throw Exception("Exception: Invalid value in command line!");

			//Set dz threshold
			if (!strcmp("dz", attribute))
			{
				dz_threshold = std::max(std::min(atof(value), 1.0), 0.0);
			}

			//Set distance threshold
			else if (!strcmp("dist", attribute))
			{
				dist_threshold = std::max(std::min(atof(value), 100.0), 0.0);
			}

			//Set z_step
			else if (!strcmp("z_step", attribute))
			{
				z_step = std::max(std::min(atof(value), 10.0), 0.1);
			}


			//Set alpha
			else if (!strcmp("alpha", attribute))
			{
				alpha = std::max(std::min(atof(value), 1.0), 0.0);
			}

			//Set beta
			else if (!strcmp("beta", attribute))
			{
				beta = std::max(std::min(atof(value), 1.0), 0.0);
			}

			//Set gamma
			else if (!strcmp("gamma", attribute))
			{
				gamma = std::max(std::min(atof(value), 1000.0), 10.0);
			}

			//Set delta
			else if (!strcmp("delta", attribute))
			{
				delta = std::max(std::min(atof(value), 1.0), 0.0);
			}

			//Set kappa
			else if (!strcmp("kappa", attribute))
			{
				kappa = std::max(std::min(atof(value), 1.5), 0.5);
			}

			//Set max iterations
			else if (!strcmp("iter", attribute))
			{
				max_iter = std::max(std::min(atoi(value), 5000), 100);
			}

			//Set min contour points
			else if (!strcmp("nc", attribute))
			{
				nc = std::max(std::min(atoi(value), 100), 10);
			}

			//Set densification step
			else if (!strcmp("dens", attribute))
			{
				dens = std::max(std::min(atof(value), 10.0), 0.0);
			}

			//Set buffer 1 file
			else if (!strcmp("buff1", attribute))
			{
				buff1_file_mask = value;
			}

			//Set buffer 2 file
			else if (!strcmp("buff2", attribute))
			{
				buff2_file_mask = value;
			}

			//Set contour file
			else if (!strcmp("cont", attribute))
			{
				contours_file_mask = value;
			}

			//Set path
			else if (!strcmp("file", attribute))
			{
				output_file_name = value;
			}

			//Set path
			else if (!strcmp("path", attribute))
			{
				path = value;
			}

			//Bad argument
			else
			{
				std::cout << attribute << '\n';
				throw Exception("Exception: Invalid attribute in command line!");
			}
		}

		//Process file
		else
		{
			//contours_file_name = argv[argc];
		}
	}

	//Write parameters
	std::cout << "*** SIMPLIFY CONTOUR LINES, MINIMUM-ENERGY SPLINES ***\n\n";
	std::cout << "Input parameters: \n" <<
		"  Buffer height = " << dz_threshold << '\n' <<
		"  Generalization = " << generalize << '\n' <<
		"  Densification step = " << dens << '\n' <<
		"  Alpha = " << alpha << '\n' <<
		"  Beta = " << beta << '\n' <<
		"  Gamma = " << gamma << '\n' <<
		"  Delta = " << delta << '\n' <<
		"  Variable delta = " << var_delta << '\n' <<
		"  Kappa = " << kappa << '\n' <<
		"  Max iterations = " << max_iter << '\n' <<
		"  Min contour points = " << nc << '\n' <<
		"  Contour mask =" << contours_file_mask << '\n' <<
		"  Buffer 1 mask = " << buff1_file_mask << '\n' <<
		"  Buffer 2 mask = " << buff2_file_mask << '\n' <<
		"  Output file = " << output_file_name << '\n' <<
		"  Path = " << path << '\n' << "\n";

	//Find contour line files
	std::cout << ">>> Read input files: ";
	TVector <std::string> cont_files;
	File::findFilesInDirByMask(path, contours_file_mask, 1, cont_files);

	//Load contours one by one
	TVector2D <std::shared_ptr <Point3D > > contours_polylines;

	File::loadContours(cont_files, contours_polylines);

	//Find buffer 1 files
	TVector <std::string> buff1_files;
	File::findFilesInDirByMask(path, buff1_file_mask, 1, buff1_files);

	//Load first buffer one by one
	std::map <double, TVector < std::shared_ptr < Point3D > > > contour_buffers1;
	File::loadBuffers(buff1_files, contour_buffers1);

	//Find buffer 1 files
	TVector <std::string> buff2_files;
	File::findFilesInDirByMask(path, buff2_file_mask, 1, buff2_files);

	//Load second buffer one by one
	std::map <double, TVector < std::shared_ptr < Point3D > > > contour_buffers2;
	File::loadBuffers(buff2_files, contour_buffers2);
	std::cout << "OK \n";

	//Simplify contour lines given by polylines
	std::cout << ">>> Simplify contour lines: \n";

	//Simplify contour lines
	try
	{
		//Create initial approximation	
		TVector2D <std::shared_ptr <Point3D > > contours_polylines_simpl;
		
		if (generalize)
			ContourLinesSimplify::simplifyContourLinesPotentialEBezier(contours_polylines, contour_buffers1, contour_buffers2, z_step, dz_threshold, ang_threshold, 2, 2, contours_polylines_simpl);
		

		else
			contours_polylines_simpl = contours_polylines;

		//Densify contour line
		TVector2D <std::shared_ptr <Point3D > > contours_polylines_dens;
		if (dens > 0)
			ContourLinesSimplify::densifyContourLines(contours_polylines_simpl, dens, contours_polylines_dens);
		else
			contours_polylines_dens = contours_polylines_simpl;

		//Minimum energy spline
		TVector2D <std::shared_ptr <Point3D > > contours_polylines_spline;
		ContourLinesSimplify::simplifyContourLinesMinimumEnergy(contours_polylines_dens, contour_buffers1, contour_buffers2, z_step, dz_threshold, alpha, beta, gamma, delta, var_delta, kappa, nc, max_iter, contours_polylines_spline);
		
		
		//Export contours
		std::string file_name_simp = "contours_" + output_file_name + "_simp_dz_" + std::format("{:.2f}", dz_threshold) + "_delta_"
			+ std::format("{:.2f}", delta) + "_var_delta_" + std::to_string(var_delta) + "_gamma_" + std::format("{:.2f}", gamma)
			+ "_iter_" + std::to_string(max_iter) + "_gener_" + std::to_string(generalize) + "_dens_" + std::format("{:.2f}", dens) + ".dxf";
		DXFExport::exportContourLinesToDXF(file_name_simp, contours_polylines_spline, 10.0);
		
	}

	//Throw exception
	catch (Exception& e)
	{
		e.printException();
	}



	std::cout << "Hello CMake." << std::endl;
	return 0;
}