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


#ifndef ContourLinesSimplify_HPP
#define ContourLinesSimplify_HPP

/*
#include <set>
#include <functional>
#include <chrono>
#include <algorithm>  
#include <cstdlib>
#include <map>
*/

#include <memory>
#include <map>

#include "TVector.h"
#include "TVector2D.h"
#include "Point3D.h"
#include "MinimumEnergySpline.h"
#include "isEqualPointByPlanarCoordinates.h"
#include "Round.h"
#include "Matrix.h"
#include "DeCasteljau.h"
#include "EuclDistance.h"
#include "PointLineDistance.h"
#include "sortPointsByDist.h"

	
void ContourLinesSimplify::simplifyContourLinesPotentialEBezier(const TVector2D <std::shared_ptr <Point3D > > & contours, std::map <double, TVector < std::shared_ptr < Point3D > > >& contour_points_buffers_dz1, std::map <double, TVector < std::shared_ptr < Point3D > > >& contour_points_buffers_dz2, const double z_step, const double dz_threshold, const double ang_threshold, const unsigned int k, const int iter, TVector2D <std::shared_ptr <Point3D > >& contours_polylines_simplified)
{
	//Simplify contour lines using potential analogous to outer energy
	//Increment the step h
	//Vertices shifted using deDasteljau algorithm
	const clock_t begin_time = clock();
	int index = 0;

	std::cout << "\n >>> Generalize contour lines: \n";

	for (auto c : contours)
	{
		//Get height of the contour
		const double z1 = c[0]->getZ() - dz_threshold;
		const double z1r = Round::roundNumber(z1, 2);
		const double z2 = c[0]->getZ() + dz_threshold;
		const double z2r = Round::roundNumber(z2, 2);

		std::cout << ">>> Z = " << c[0]->getZ() << "m, n = " << c.size() << '\n';

		//Find corresponding buffer z - dz
		TVector < std::shared_ptr<Point3D > > contour_points_buffer_dz1;
		typename std::map <double, TVector <std::shared_ptr<Point3D > > >::iterator  it_contour_points_buffers_dz1 = contour_points_buffers_dz1.find(z1r);

		//Buffer z - dz already created
		if (it_contour_points_buffers_dz1 != contour_points_buffers_dz1.end())
			contour_points_buffer_dz1 = it_contour_points_buffers_dz1->second;

		//No buffer found
		else
			continue;

		//Find corresponding buffer z + dz
		TVector < std::shared_ptr<Point3D > > contour_points_buffer_dz2;
		typename std::map <double, TVector <std::shared_ptr<Point3D > > >::iterator it_contour_points_buffers_dz2 = contour_points_buffers_dz2.find(z2r);

		//Buffer z + dz already created
		if (it_contour_points_buffers_dz2 != contour_points_buffers_dz2.end())
			contour_points_buffer_dz2 = it_contour_points_buffers_dz2->second;

		//No buffer found
		else
			continue;

		//Are there enough points?
		int nb1 = contour_points_buffer_dz1.size();
		int nb2 = contour_points_buffer_dz2.size();

		if ((c.size() > k) && (nb1 > k) && (nb2 > k))
		{
			
			//Create list of predecessors and successors
			for (int h = 10; h <= 10; h++)
			{
				std::cout << "h = " << h << ' ';
				int indexii = 0;

				//Repeat until any swap is available
				int i_first = 0, i_last = c.size();
				TVector <double> pots(c.size(), 1.0e16);
				TVector2D <std::shared_ptr <Point3D > > bps(c.size());
				for (;;)
				{
					//Get two consequent edges with the highest replacement ratio
					int i1 = 0;
					int i3 = (i1 + 2*h < c.size() ? i1 + 2*h : -1);

					//Actulize mninimum
					double pot_min = 1.e16;
					int i_min = -1;

					//Compute potential
					while (i3 != -1)
					{
						//std::cout << i1 << " " << i2 << " " << i3  << '\n';

						//Compute potential for all points inside the interval
						if (i1 >= i_first && i1 <= i_last)
						{
							//Compute potential
							auto [pot, bp] = computePotentialEBezier(i1, i3, c, contour_points_buffer_dz1, contour_points_buffer_dz2);
							
							//Store potential and Bezier vertices
							pots[i1] = pot;
							bps[i1] = bp;
						}

						//Actualize maximum replacement potential: c(i, i3) as much as possible in the middle of the buffer
						if (pots[i1] < pot_min)
						{
							i_min = i1;
							pot_min = pots[i1];
						}

						//Increment vertices
						i1 = (i1 + h < c.size() ? i1 + h : -1);
						i3 = (i1 + 2 * h < c.size() ? i1 + 2 * h : -1);
					}

					//A suitable vertex with the maximum replacement potential has not been found
					if (pot_min >= -0.1)
						break;
					//if (pot_min >= -0.01)
					//if (pot_min >= 1.1)

					std::cout << i_min << " " << pot_min << '\n';

					//Amount of shifted points
					int np = pow(2, h + 1) - 1;

					//Get remaining indices
					int i3_min = (i_min + 2 * h < c.size() ? i_min + 2 * h : -1);

					//Add first Bezier control point
					TVector <double> xp, yp;
		
					//Update contour line vertices by Bezier points
					int j = 0;
					for (int i = i_min; i <= i3_min; i++, j++)
					{
						if (j >= bps[i_min].size())
						{
							std::cout << "error";
							std::cout << i_min << " " << i3_min << bps[i_min].size();
						}

						c[i]->setX(bps[i_min][j]->getX());
						c[i]->setY(bps[i_min][j]->getY());
					}

					//Actualize start index to prev-prev element
					i_first = (i_min - 3 * h < 0 ? 0 : i_min - 3 * h);
	
					//Actualize end index to next-next element
					i_last = (i_min + 3 * h < c.size() ? i_min + 3 * h : c.size() - 1); 
				}
			}

			//Add simplified contour to the list of contours
			contours_polylines_simplified.push_back(c);

			index++;
		}
	}

	std::cout << "OK";
	std::cout << float(clock() - begin_time) / CLOCKS_PER_SEC;
}


std::tuple<double, TVector <std::shared_ptr <Point3D > > > ContourLinesSimplify::computePotentialEBezier(const int i1, const int i3, const TVector < std::shared_ptr<Point3D> >& contour, const TVector < std::shared_ptr<Point3D > >& contour_points_buffer_dz1, const TVector < std::shared_ptr<Point3D > >& contour_points_buffer_dz2)
{
	//Compute simplification potential of the created segment using Energy function
	//Approximate buffer, smooth points using deCasteljau algorithm
	TVector <size_t> knn_id1, knn_id2, knn_id3, knn_id4;
	TVector <float> knn_dist1, knn_dist2, knn_dist3, knn_dist4;
	TVector < std::shared_ptr<Point3D > > bps, knn_points1, knn_points2, knn_points3, knn_points4;

	//Compute length of the segment
	const double length = EuclDistance::getEuclDistance2D(contour[i1]->getX(), contour[i1]->getY(), contour[i3]->getX(), contour[i3]->getY());

	//Identical points
	if (length < 0.001)
		return { 1.0e16, bps };

	//Get control points, vertices of the contour line in the given interval
	TVector <std::shared_ptr <Point3D > > cp(contour.cbegin() + i1, contour.cbegin() + i3 + 1);

	//Find NN to newly created vertices
	findAllNNS(contour_points_buffer_dz1, cp, knn_id1, knn_dist1, knn_points1);
	findAllNNS(contour_points_buffer_dz2, cp, knn_id2, knn_dist2, knn_points2);

	//Compute potential of newly created vertices
	double pot_sum_old = getOE3(cp, contour_points_buffer_dz1, contour_points_buffer_dz2, knn_id1, knn_id2);

	//Compute points on the Bezier curve
	for (int i = 0; i < cp.size(); i++)
	{
		//Get u parameter
		double u = (double)(i) / (cp.size() - 1);

		//Compute Bezier point
		std::shared_ptr <Point3D> bp = DeCasteljau::computeBezierPoint(u, cp);
		
		//Add to the list
		bps.push_back(bp);
	}

	//Find NN to newly created vertices
	findAllNNS(contour_points_buffer_dz1, bps, knn_id3, knn_dist3, knn_points3);
	findAllNNS(contour_points_buffer_dz2, bps, knn_id4, knn_dist4, knn_points4);

	//Compute potential of newly created vertices
	double pot_sum_new = getOE3(bps, contour_points_buffer_dz1, contour_points_buffer_dz2, knn_id3, knn_id4);

	//return (pot_sum_new - pot_sum_old)/ pot_sum_old;
	return { pot_sum_new - pot_sum_old, bps };
}


void ContourLinesSimplify::densifyContourLines(const TVector2D <std::shared_ptr <Point3D > >& contours, const double dmin, TVector2D <std::shared_ptr <Point3D > >& contours_dens)
{
	//Densify contour line vertices with a given step
	for (auto c : contours)
	{
		TVector < std::shared_ptr <Point3D > > c_dens;
		for (int i = 0; i < c.size() - 1; i++)
		{
			//Direction vector
			const double ux = c[i + 1]->getX() - c[i]->getX();
			const double uy = c[i + 1]->getY() - c[i]->getY();
			
			//Length of the segment
			const double nu = sqrt(ux * ux + uy * uy);
			const double nx = ux / nu;
			const double ny = uy / nu;

			//Densify contour points
			for (int j = 0; j <= std::max(0.0, nu - 0.75 * dmin); j += dmin)
			{
				std::shared_ptr <Point3D > pc = std::make_shared<Point3D >(c[i]->getX() + j * nx, c[i]->getY() + j * ny, c[0]->getZ());
				c_dens.push_back(pc);
			}
		}

		//Last point of the last segment
		c_dens.push_back(c.back());

		//Add densified contour line to the list
		contours_dens.push_back(c_dens);
	}
}


void ContourLinesSimplify::simplifyContourLinesMinimumEnergy(const TVector2D <std::shared_ptr <Point3D > >& contours, std::map <double, TVector < std::shared_ptr < Point3D > > >& contour_points_buffers_dz1, std::map <double, TVector < std::shared_ptr < Point3D > > >& contour_points_buffers_dz2, const double z_step, const double dz_threshold, const double alpha, const double beta, const double gamma, const double delta, const bool var_delta, const double kappa, const int nc, const int max_iter, TVector2D <std::shared_ptr <Point3D > >& contours_polylines_simplified)
{
	const clock_t begin_time = clock();
	int index = 0;

	std::cout << "\n >>> Minimum energy splines: \n";

	for (auto c : contours)
	{
		//Get height of the contour
		const double z1 = c[0]->getZ() - dz_threshold;
		const double z1r = Round::roundNumber(z1, 2);
		const double z2 = c[0]->getZ() + dz_threshold;
		const double z2r = Round::roundNumber(z2, 2);

		std::cout << ">>> Z = " << c[0]->getZ() << "m, n = " << c.size() << '\n';

		//Find corresponding buffer z - dz
		TVector < std::shared_ptr<Point3D > > contour_points_buffer_dz1;
		typename std::map <double, TVector <std::shared_ptr<Point3D > > >::iterator  it_contour_points_buffers_dz1 = contour_points_buffers_dz1.find(z1r);

		//Buffer z - dz already created
		if (it_contour_points_buffers_dz1 != contour_points_buffers_dz1.end())
			contour_points_buffer_dz1 = it_contour_points_buffers_dz1->second;

		//No buffer found
		else
			continue;

		//Find corresponding buffer z + dz
		TVector < std::shared_ptr<Point3D > > contour_points_buffer_dz2;
		typename std::map <double, TVector <std::shared_ptr<Point3D > > >::iterator it_contour_points_buffers_dz2 = contour_points_buffers_dz2.find(z2r);

		//Buffer z + dz already created
		if (it_contour_points_buffers_dz2 != contour_points_buffers_dz2.end())
			contour_points_buffer_dz2 = it_contour_points_buffers_dz2->second;

		//No buffer found
		else
			continue;

		//Remove duplicate points
		contour_points_buffer_dz1.erase(std::unique(contour_points_buffer_dz1.begin(), contour_points_buffer_dz1.end(), isEqualPointByPlanarCoordinates<std::shared_ptr <Point3D > >()), contour_points_buffer_dz1.end());
		contour_points_buffer_dz2.erase(std::unique(contour_points_buffer_dz2.begin(), contour_points_buffer_dz2.end(), isEqualPointByPlanarCoordinates<std::shared_ptr <Point3D > >()), contour_points_buffer_dz2.end());

		//Are there enough contour line points?
		if ((c.size() > nc) && (contour_points_buffer_dz1.size() > nc) && (contour_points_buffer_dz2.size() > nc))
		{
			//Find nearest neighbors
			TVector <size_t> knn_id1, knn_id2;
			TVector <float> knn_dist1, knn_dist2;

			//Create approximate vertical buffer: find nearest points on both buffers
			findAllNN(contour_points_buffer_dz1, c, knn_id1, knn_dist1);
			findAllNN(contour_points_buffer_dz2, c, knn_id2, knn_dist2);

			//Create buffer 1 matrix
			Matrix<double> xb1(contour_points_buffer_dz1.size(), 1), yb1(contour_points_buffer_dz1.size(), 1);
			for (int j = 0; j < contour_points_buffer_dz1.size(); j++)
			{
				xb1(j, 0) = contour_points_buffer_dz1[j]->getX();
				yb1(j, 0) = contour_points_buffer_dz1[j]->getY();
			}

			//Create buffer 2 matrix
			Matrix<double> xb2(contour_points_buffer_dz2.size(), 1), yb2(contour_points_buffer_dz2.size(), 1);
			for (int j = 0; j < contour_points_buffer_dz2.size(); j++)
			{
				xb2(j, 0) = contour_points_buffer_dz2[j]->getX();
				yb2(j, 0) = contour_points_buffer_dz2[j]->getY();
			}

			//Process each fragment
			const int n_max = 150;
			int i_last = n_max;
			std::shared_ptr<Point3D > p_last;
			for (int i = 0; i < c.size() - 1; i = i_last)
			{
				//Get last index, avoid short segments
				const int i_last_min = std::min((int)c.size() - 1, i + n_max);
				i_last = (c.size() - i_last_min > 10 ? i_last_min : c.size() - 1);
				std::cout << i << ":" << i_last << "  ";

				//Create submatrix of knn
				TVector <size_t> knn_idf1, knn_idf2;
				std::copy(knn_id1.begin() + i, knn_id1.begin() + i_last, std::back_inserter(knn_idf1));
				std::copy(knn_id2.begin() + i, knn_id2.begin() + i_last, std::back_inserter(knn_idf2));

				//Create contour matrix
				Matrix<double> xc(i_last - i + 1, 1), yc(i_last - i + 1, 1);
				for (int j = 0; j < i_last - i + 1; j++)
				{
					xc(j, 0) = c[i + j]->getX();
					yc(j, 0) = c[i + j]->getY();
				}

				//Create minimum energy spline
				auto [xcn, ycn] = MinimumEnergySpline::createSpline(xc, yc, xb1, yb1, xb2, yb2, knn_idf1, knn_idf2, alpha, beta, gamma, delta, var_delta, kappa, max_iter);

				//Convert spline to contour line
				TVector <std::shared_ptr <Point3D > > contour_polyline_simplified;
				for (int j = 0; j < xc.rows(); j++)
				{
					//Create new point
					std::shared_ptr<Point3D > p = std::make_shared<Point3D >(xcn(j, 0), ycn(j, 0), c[0]->getZ());

					//Add point to the contour line
					contour_polyline_simplified.push_back(p);
				}

				//Modify the junction points
				if (i > 0)
				{
					//Last point of the previous segment and first point of the actual segment
					std::shared_ptr<Point3D >& p_last = contours_polylines_simplified.back().back();
					std::shared_ptr<Point3D >& p_first = contour_polyline_simplified.front();

					//Average point
					const double xa = 0.5 * (p_last->getX() + p_first->getX());
					const double ya = 0.5 * (p_last->getY() + p_first->getY());

					//Modify the junction point
					p_last->setX(xa); p_last->setY(ya);
					p_first->setX(xa); p_first->setY(ya);

				}

				//Add simplified contour to the list of contours
				contours_polylines_simplified.push_back(contour_polyline_simplified);
			}
		}

		index++;

	}

	std::cout << "OK";
	std::cout << float(clock() - begin_time) / CLOCKS_PER_SEC;
}


void ContourLinesSimplify::findAllNN(const TVector <std::shared_ptr <Point3D > >& points, const TVector <std::shared_ptr <Point3D > >& qpoints, TVector <size_t>& knn_ids, TVector <float>& knn_dists)
{
	//Find nearest neighbor to any contour line vertex
	const int n = points.size(), nq = qpoints.size();

	knn_ids.clear(); knn_dists.clear();

	for (int i = 0; i < qpoints.size(); i++)
	{
		//Find element closest to q
		sortPointsByDist sd(qpoints[i]);
		auto it_min = min_element(points.begin(), points.end(), sd);

		//Get its distance
		const float d_min = EuclDistance::getEuclDistance2D((*it_min)->getX(), (*it_min)->getY(), qpoints[i]->getX(), qpoints[i]->getY());

		//Get its index
		const size_t i_min = std::distance(points.begin(), it_min);

		//Add to the list
		knn_ids.push_back(i_min);
		knn_dists.push_back(d_min);
	}
}


void ContourLinesSimplify::findAllNNS(const TVector <std::shared_ptr <Point3D> >& points, const TVector <std::shared_ptr <Point3D> >& qpoints, TVector <size_t>& knn_ids, TVector <float>& knn_dists, TVector <std::shared_ptr <Point3D > >& knn_points)
{
	//Find nearest neighbor to any contour line segment
	const double l_min = 50.0;
	for (int i = 0; i < qpoints.size(); i++)
	{
		//Check all line segments
		int i_min = -1;
		double d_min = 1.0e16, xi_min = 0, yi_min = 0;

		for (int j = 0; j < points.size() - 1; j++)
		{
			//Length of the segment
			const double dx = points[j]->getX() - points[j+1]->getX();
			const double dy = points[j]->getY() - points[j + 1]->getY();
			double l = sqrt(dx * dx + dy * dy);

			//Avoid too long segments
			if (l < l_min)
			{
				double xi, yi;
				double d = PointLineDistance::getPointLineSegmentDistance2D(qpoints[i]->getX(), qpoints[i]->getY(), points[j]->getX(), points[j]->getY(), points[j + 1]->getX(), points[j + 1]->getY(), xi, yi);

				//Remember actual minimum
				if (d < d_min)
				{
					d_min = d; i_min = j;
					xi_min = xi; yi_min = yi;
				}
			}
		}

		//Create nearest point
		std::shared_ptr <Point3D> p_min = std::make_shared<Point3D>(xi_min, yi_min);


		//Add to the list
		knn_ids.push_back(i_min);
		knn_dists.push_back(d_min);
		knn_points.push_back(p_min);
	}
}


double ContourLinesSimplify::getOE3(const TVector <std::shared_ptr <Point3D > >& cp, const TVector < std::shared_ptr<Point3D > >& contour_points_buffer_dz1, const TVector < std::shared_ptr<Point3D > >& contour_points_buffer_dz2, const TVector <size_t>& knn_id1, const TVector <size_t>& knn_id2)
{
	//Compute outer energy (potential) of shifted vertices, var 3
	double pot_sum = 0;
	for (int i = 0; i < cp.size(); i++)
	{
		//Compute potential, corresponding to outer enetgy 3
		double pot = getOE3Point(cp[i]->getX(), cp[i]->getY(), contour_points_buffer_dz1[knn_id1[i]]->getX(), contour_points_buffer_dz1[knn_id1[i]]->getY(),
			contour_points_buffer_dz2[knn_id2[i]]->getX(), contour_points_buffer_dz2[knn_id2[i]]->getY());

		//Summarize potential
		pot_sum += pot;
	}

	return pot_sum;
}


double ContourLinesSimplify::getOE3Point(const double xc, const double yc, const double xb1, const double yb1, const double xb2, const double yb2)
{
	//Compute outer energy of a point
	const double dxb1 = xc - xb1;
	const double dyb1 = yc - yb1;
	const double dxb2 = xc - xb2;
	const double dyb2 = yc - yb2;
	const double dxb12 = xb2 - xb1;
	const double dyb12 = yb2 - yb1;
	const double db1 = sqrt(dxb1 * dxb1 + dyb1 * dyb1);
	const double db2 = sqrt(dxb2 * dxb2 + dyb2 * dyb2);
	const double db12 = sqrt(dxb12 * dxb12 + dyb12 * dyb12);

	return fabs(db1 * db1 - db2 * db2) / db12;
}

#endif
