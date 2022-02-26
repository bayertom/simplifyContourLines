// Description: Simplification of contour lines

// Copyright (c) 2017 - 2018
// Tomas Bayer
// Charles University in Prague, Faculty of Science
// bayertom@natur.cuni.cz

// This library is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#ifndef ContourLinesSimplify_H
#define ContourLinesSimplify_H

//#include <ctime>
//#include <cstdlib>
//#include <iostream>
//#include <queue>

#include <memory>
#include <map>

#include "TVector2D.h"
#include "Point3D.h"


//Contour line simplification
class ContourLinesSimplify
{
	public:

		static void simplifyContourLinesPotentialEBezier(const TVector2D <std::shared_ptr <Point3D > >& contours, std::map <double, TVector < std::shared_ptr < Point3D > > >& contour_points_buffers_dz1, std::map <double, TVector < std::shared_ptr < Point3D > > >& contour_points_buffers_dz2, const double z_step, const double dz_threshold, const double ang_threshold, const unsigned int k, const int iter, TVector2D <std::shared_ptr <Point3D > >& contours_polylines_simplified);
		static std::tuple<double, TVector <std::shared_ptr <Point3D > > > computePotentialEBezier(const int i1, const int i3, const TVector < std::shared_ptr<Point3D > >& contour, const TVector < std::shared_ptr<Point3D > >& contour_points_buffer_dz1, const TVector < std::shared_ptr<Point3D > >& contour_points_buffer_dz2);
		static void densifyContourLines(const TVector2D <std::shared_ptr <Point3D > >& contours, const double dmin, TVector2D <std::shared_ptr <Point3D > >& contours_dens);
		static void simplifyContourLinesMinimumEnergy(const TVector2D <std::shared_ptr <Point3D > >& contours, std::map <double, TVector < std::shared_ptr < Point3D > > >& contour_points_buffers_dz1, std::map <double, TVector < std::shared_ptr < Point3D > > >& contour_points_buffers_dz2, const double z_step, const double dz_threshold, const double alpha, const double beta, const double gamma, const double delta, const bool var_delta, const double kappa, const int nc, const int max_iter, TVector2D <std::shared_ptr <Point3D > >& contours_polylines_simplified);
		static void findAllNN(const TVector <std::shared_ptr <Point3D > >& points, const TVector <std::shared_ptr <Point3D > >& qpoints, TVector <size_t>& knn_ids, TVector <float>& knn_dists);
		static void findAllNNS(const TVector <std::shared_ptr <Point3D> >& points, const TVector <std::shared_ptr <Point3D> >& qpoints, TVector <size_t>& knn_ids, TVector <float>& knn_dists, TVector <std::shared_ptr <Point3D > >& knn_points);
		static double getOE3(const TVector <std::shared_ptr <Point3D > >& cp, const TVector < std::shared_ptr<Point3D > >& contour_points_buffer_dz1, const TVector < std::shared_ptr<Point3D > >& contour_points_buffer_dz2, const TVector <size_t>& knn_id1, const TVector <size_t>& knn_id2);
		static double getOE3Point(const double xc, const double yc, const double xb1, const double yb1, const double xb2, const double yb2);

}; 

#include "ContourLinesSimplify.hpp"

#endif