// Description: Contour line simplification using potential and minimum energy splines

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

	
void ContourLinesSimplify::smoothContourLinesByPotential(const TVector2D <std::shared_ptr <Point3D > > & contours, std::multimap <double, TVector < std::shared_ptr < Point3D > > >& contour_points_buffers_dz1, std::multimap <double, TVector < std::shared_ptr < Point3D > > >& contour_points_buffers_dz2, const double z_step, const double dz_threshold, const unsigned int min_points, const int max_iter, TVector2D <std::shared_ptr <Point3D > >& contours_polylines_simplified)
{
	//Simplify contour lines using potential analogous to outer energy
	//Increment the step h
	//Vertices shifted using deDasteljau algorithm
	const clock_t begin_time = clock();
	int index = 0;

	std::cout << "\n>>> PHASE1: Smoothing contour lines \n\n";

	for (auto c : contours)
	{
		//Get height of the contour
		const double z1 = c[0]->getZ() - dz_threshold;
		const double z1r = Round::roundNumber(z1, 2);
		const double z2 = c[0]->getZ() + dz_threshold;
		const double z2r = Round::roundNumber(z2, 2);

		std::cout << ">>> Z = " << c[0]->getZ() << "m, n = " << c.size() << ":\n";

		//Find corresponding buffer z - dz
		TVector2D < std::shared_ptr<Point3D > > contour_points_buffer_dz1;

		//Find corresponding buffer z - dz
		auto res1 = contour_points_buffers_dz1.equal_range(z1r);

		///No buffer found
		if (res1.first == res1.second)
			continue;

		// Buffer z - dz already created
		else
		{
			//Copy buffer segments
			for (auto it = res1.first; it != res1.second; ++it)
			{
				if (it->second.size() > min_points)
					contour_points_buffer_dz1.push_back(it->second);
			}
		}

		//Find corresponding buffer z + dz
		TVector2D < std::shared_ptr<Point3D > > contour_points_buffer_dz2;

		//Find corresponding buffer z - dz
		auto res2 = contour_points_buffers_dz2.equal_range(z2r);

		///No buffer found
		if (res2.first == res2.second)
			continue;

		// Buffer z + dz already created
		else
		{
			//Copy buffer segments
			for (auto it = res2.first; it != res2.second; ++it)
			{
				if (it->second.size() > min_points)
					contour_points_buffer_dz2.push_back(it->second);
			}
		}

		//Are there enough points?
		if ((c.size() > min_points) && (contour_points_buffers_dz1.size() > 0) && (contour_points_buffers_dz2.size() > 0))
		{
			//if (c.size() != 11)
			//	continue;

			int i1 = 0;
			int n = c.size();

			//Find NN to contour line vertices
			TVector<int> buffer_ids1(n, -1), buffer_ids2(n, -1), buffer_ids3(n, -1), buffer_ids4(n, -1);
			auto [nn_dist1, nn_points1] = findNearestNeighbors(c, contour_points_buffer_dz1, i1, buffer_ids1);
			auto [nn_dist2, nn_points2] = findNearestNeighbors(c, contour_points_buffer_dz2, i1, buffer_ids2);
			
			//Create list of predecessors and successors
			for (int h = 1; h <= 10; h++)
			{
				std::cout << "h = " << h << ' ';

				//Repeat until any swap is available
				int i_first = 0, i_last = n;
				TVector <double> pots_diff(n, 1.0e16);
				TVector2D <std::shared_ptr <Point3D > > bps(n), nn_points3(n), nn_points4(n);

				for (;;)
				{
					//Get two consequent edges with the highest replacement ratio
					int i1 = 0;
					int i3 = (i1 + 2*h < c.size() ? i1 + 2*h : -1);

					//Minimum values
					double pot_diff_min = 1.e16;
					int i_min = -1;
					TVector <std::shared_ptr <Point3D > > nn_points3_min, nn_points4_min;

					//Compute potential differeces
					while (i3 != -1)
					{
						//std::cout << i1 << " " << i2 << " " << i3  << '\n';

						//Compute potentiaå difference over the contour interval c(i, i3) 
						if (i1 >= i_first && i1 <= i_last)
						{
							TVector <float> nn_dist3, nn_dist4;
							std::tie(pots_diff[i1], bps[i1], nn_dist3, nn_points3[i1], nn_dist4, nn_points4[i1]) = computePotentialDifferenceBezier(i1, i3, c, contour_points_buffer_dz1, contour_points_buffer_dz2, nn_points1, nn_points2, buffer_ids3, buffer_ids4);
						}

						//Actualize minimum potential difference
						if (pots_diff[i1] < pot_diff_min)
						{
							i_min = i1;
							pot_diff_min = pots_diff[i1];
							nn_points3_min = nn_points3[i1];
							nn_points4_min = nn_points4[i1];
						}

						//Increment interval
						i1 = (i1 + h < c.size() ? i1 + h : -1);
						i3 = (i1 + 2 * h < c.size() ? i1 + 2 * h : -1);
					}

					//A suitable vertex interval (i1, i3) with the maximum potential difference has not been found
					if (pot_diff_min >= -0.1)
						break;

					//Update contour line vertices between (i1, i3) by Bezier points
					std::copy_n(bps[i_min].cbegin(), bps[i_min].size(), &c[i_min]);
	
					//Update nearest neighbors of the actualized contour part between (i1, i3)
					std::copy_n(nn_points3_min.cbegin(), nn_points3_min.size(), &nn_points1[i_min]);
					std::copy_n(nn_points4_min.cbegin(), nn_points3_min.size(), &nn_points2[i_min]);

					//std::cout << i_min << " " << pot_diff_min << '\n';
					//std::cout << ".";

					//Amount of shifted points
					int np = pow(2, h + 1) - 1;

					//Get last point of the interval
					int i3_min = (i_min + 2 * h < c.size() ? i_min + 2 * h : -1);

					//Actualize start index for potential recomputation to prev-prev element
					i_first = (i_min - 3 * h < 0 ? 0 : i_min - 3 * h);
	
					//Actualize end index for potential recomputation to next-next element
					i_last = (i_min + 3 * h < c.size() ? i_min + 3 * h : c.size() - 1); 
				}
			}

			//Add simplified contour to the list of contours
			contours_polylines_simplified.push_back(c);

			index++;

			std::cout << '\n';
		}
	}

	std::cout << "OK";
	std::cout << float(clock() - begin_time) / CLOCKS_PER_SEC;
}


std::tuple<double, TVector <std::shared_ptr <Point3D > >, TVector <float>, TVector <std::shared_ptr <Point3D > >, TVector <float>, TVector <std::shared_ptr <Point3D > > > ContourLinesSimplify::computePotentialDifferenceBezier(const int i1, const int i3, const TVector < std::shared_ptr<Point3D> >& contour, const TVector2D < std::shared_ptr<Point3D > >& contour_points_buffer_dz1, const TVector2D < std::shared_ptr<Point3D > >& contour_points_buffer_dz2, const TVector < std::shared_ptr<Point3D> >& nn_points1, const TVector < std::shared_ptr<Point3D> >& nn_points2, TVector<int>& buffer_ids3, TVector<int>& buffer_ids4)
{
	//Compute simplification potential difference of the segment (i1, i3)
	//Approximate buffer, smooth points using deCasteljau algorithm
	
	//Get part of the contour line and correpsonding nearest points
	const TVector <std::shared_ptr <Point3D > > cp(contour.cbegin() + i1, contour.cbegin() + i3 + 1);
	const TVector <std::shared_ptr <Point3D > > nn_points1p(nn_points1.cbegin() + i1, nn_points1.cbegin() + i3 + 1);
	const TVector <std::shared_ptr <Point3D > > nn_points2p(nn_points2.cbegin() + i1, nn_points2.cbegin() + i3 + 1);

	//Compute potential of newly created vertices
	double pot_sum_old = getOE3(cp, nn_points1p, nn_points2p);

	//Compute contour points on the Bezier curve
	TVector < std::shared_ptr<Point3D > > bps;
	for (int i = 0; i < cp.size(); i++)
	{
		//Get u parameter
		const double u = (double)(i) / (cp.size() - 1);

		//Compute Bezier point
		std::shared_ptr <Point3D> bp = DeCasteljau::computeBezierPoint(u, cp);

		//Add to the list
		bps.push_back(bp);
	}

	//Find NN to the contour points formed by Bezier vertices
	auto [nn_dist3p, nn_points3p] = findNearestNeighbors(bps, contour_points_buffer_dz1, i1, buffer_ids3);
	auto [nn_dist4p, nn_points4p] = findNearestNeighbors(bps, contour_points_buffer_dz2, i1, buffer_ids4);

	//Compute potential of newly created vertices
	const double pot_sum_new = getOE3(bps, nn_points3p, nn_points4p);

	//Compute potential improvement
	return { pot_sum_new - pot_sum_old, bps, nn_dist3p, nn_points3p,  nn_dist4p, nn_points4p};
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


void ContourLinesSimplify::simplifyContourLinesMinimumEnergy(const TVector2D <std::shared_ptr <Point3D > >& contours, std::multimap <double, TVector < std::shared_ptr < Point3D > > >& contour_points_buffers_dz1, std::multimap <double, TVector < std::shared_ptr < Point3D > > >& contour_points_buffers_dz2, const double z_step, const double dz_threshold, const double alpha, const double beta, const double gamma, const double delta, const double kappa, const int min_points, const int max_iter, TVector2D <std::shared_ptr <Point3D > >& contours_polylines_simplified)
{
	//Simplify contour lines using the minimum energy spline
	const clock_t begin_time = clock();
	int index = 0;

	std::cout << "\n\n>>> PHASE2: Minimum energy splines \n\n";

	for (auto c : contours)
	{
		//Get height of the contour
		const double z1 = c[0]->getZ() - dz_threshold;
		const double z1r = Round::roundNumber(z1, 2);
		const double z2 = c[0]->getZ() + dz_threshold;
		const double z2r = Round::roundNumber(z2, 2);

		std::cout << ">>> Z = " << c[0]->getZ() << "m, n = " << c.size() << ":\n";

		//Find corresponding buffer z - dz
		TVector2D < std::shared_ptr<Point3D > > contour_points_buffer_dz1;

		//Find corresponding buffer z - dz
		auto res1 = contour_points_buffers_dz1.equal_range(z1r);

		///No buffer found
		if (res1.first == res1.second)
			continue;

		// Buffer z - dz already created
		else
		{
			//Copy buffer segments
			for (auto it = res1.first; it != res1.second; ++it)
			{
				if (it->second.size() > min_points)
					contour_points_buffer_dz1.push_back(it->second);
			}
		}

		//Find corresponding buffer z + dz
		TVector2D < std::shared_ptr<Point3D > > contour_points_buffer_dz2;

		//Find corresponding buffer z - dz
		auto res2 = contour_points_buffers_dz2.equal_range(z2r);

		///No buffer found
		if (res2.first == res2.second)
			continue;

		// Buffer z + dz already created
		else
		{
			//Copy buffer segments
			for (auto it = res2.first; it != res2.second; ++it)
			{
				if (it->second.size() > min_points)
					contour_points_buffer_dz2.push_back(it->second);
			}
		}

		//Are there enough contour line points?
		if ((c.size() > min_points) && (contour_points_buffers_dz1.size() > 0) && (contour_points_buffers_dz2.size() > 0))
		{
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

				//Create minimum energy spline
				int m = i_last - i + 1;
				TVector <std::shared_ptr <Point3D > > contour_polyline_simplified = MinimumEnergySpline::createSpline(c, contour_points_buffer_dz1, contour_points_buffer_dz2, i, m, alpha, beta, gamma, delta, kappa, max_iter);
				
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


std::tuple<TVector <float>, TVector <std::shared_ptr <Point3D > > > ContourLinesSimplify::findNearestNeighbors(const TVector <std::shared_ptr <Point3D> >& qpoints, const TVector2D <std::shared_ptr <Point3D> >& buffers, const int i1, TVector <int> &buffer_ids)
{
	//Find nearest neighbor to any contour line vertex
	const int n = qpoints.size();

	TVector <float> nn_dists(n, MAX_FLOAT);
	TVector <std::shared_ptr <Point3D > > nn_points (n);

	//Process all query points
	for (int i = 0; i < n; i++)
	{
		//Check specific buffer fragments
		const int j_start_buff = (buffer_ids[i + i1] == -1 ? 0 : buffer_ids[i + i1]);
		const int j_end_buff = (buffer_ids[i + i1] == -1 ? buffers.size() : j_start_buff + 1);

		for (int j = j_start_buff;j < j_end_buff; j++)
		{
			//Find nearest line segment point
			const auto [d_min, xi_min, yi_min] = getNearestLineSegmentPoint(qpoints[i]->getX(), qpoints[i]->getY(), buffers[j]);

			//We found a closer point
			if (d_min < nn_dists[i])
			{
				//Create nearest point
				std::shared_ptr <Point3D> p_min = std::make_shared<Point3D>(xi_min, yi_min);

				//Actualize lists of neighbors
				nn_dists[i] = d_min;
				nn_points[i] = p_min;

				//Actualize nearest buffer fragment
				buffer_ids[i + i1] = j;
			}
		}
	}

	return { nn_dists, nn_points };
}


std::tuple<double, double, double> ContourLinesSimplify::getNearestLineSegmentPoint(const double xq, const double yq, const TVector <std::shared_ptr <Point3D > > &points)
{
	//Check all line segments
	double d_min = 1.0e16, xi_min = 0, yi_min = 0;

	for (int i = 0; i < points.size() - 1; i++)
	{
		double xi, yi;
		double d = PointLineDistance::getPointLineSegmentDistance2D(xq, yq, points[i]->getX(), points[i]->getY(), points[i + 1]->getX(), points[i + 1]->getY(), xi, yi);
		
		//Remember actual minimum
		if (d < d_min)
		{
			d_min = d;
			xi_min = xi; yi_min = yi;
		}
	}

	return { d_min, xi_min, yi_min };
}


double ContourLinesSimplify::getOE3(const TVector <std::shared_ptr <Point3D > >& cp, const TVector <std::shared_ptr <Point3D > >& contour_points_buffer_dz1, const TVector <std::shared_ptr <Point3D > >& contour_points_buffer_dz2)
{
	//Compute outer energy (potential) of shifted vertices, var 3
	double pot_sum = 0;
	for (int i = 0; i < cp.size(); i++)
	{
		//Compute potential, corresponding to outer enetgy 3
		double pot = getOE3Point(cp[i]->getX(), cp[i]->getY(), contour_points_buffer_dz1[i]->getX(), contour_points_buffer_dz1[i]->getY(),
			contour_points_buffer_dz2[i]->getX(), contour_points_buffer_dz2[i]->getY());

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
