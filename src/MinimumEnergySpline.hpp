// Description: Minimum energy spline (snake), discrete approximation

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


#ifndef MinimumEnergySpline_HPP
#define MinimumEnergySpline_HPP

#include <filesystem>

#include "MatrixOperations.h"
#include "PointLineDistance.h"
#include "ContourLinesSimplify.h"

//Set namespaces
using namespace MatrixOperations;


TVector <std::shared_ptr <Point3D > > MinimumEnergySpline::createSpline(const TVector <std::shared_ptr <Point3D > >& contour, const TVector2D < std::shared_ptr<Point3D > > & buffers1, const TVector2D < std::shared_ptr<Point3D > > & buffers2, int i, int n, const double alpha, const double beta, const double gamma, const double delta, const double kappa, const unsigned int max_iter)
{
	//Create minimum energy spline
	TVector <std::shared_ptr <Point3D > > contour_simplified;

	//Create contour coordinate matrix
	Matrix <double> xc0(n, 1), yc0(n, 1);
	for (int j = 0; j < n; j++)
	{
		xc0(j, 0) = contour[i + j]->getX();
		yc0(j, 0) = contour[i + j]->getY();
	}

	//Assign initial solution
	Matrix <double> xc = xc0, yc = yc0;
	Matrix <double> xck = xc0, yck = yc0;

	//Compute differences
	Matrix <double> dx = xc - xc0;
	Matrix <double> dy = yc - yc0;
	
	//Compute average distance
	const Matrix <double> diffx = diff(xc);
	const Matrix <double> diffy = diff(yc);
	const Matrix <double> dist = sqrtm(diffx % diffx + diffy % diffy);
	const double d_aver = mean(dist);
	

	//Create A matrix
	Matrix <double> A = createA(alpha, beta, gamma, d_aver, n);
	//Matrix <double> A = createAFull(xc, yc, alpha, beta, gamma, n);

	//Inverse using LU decomposition
	short sign = 1;
	Matrix <double> L(n, n), U(n, n), P(n, n, 0, 1), I(n, n, 0, 1);
	lu(A + delta * I, L, U, P, sign);
	Matrix <double> Ainv = inv(L) * inv(U);
	
	//NN search: find nearest points on both buffers 
	TVector<int> buffer_ids1(n, -1), buffer_ids2(n, -1);
	auto [nn_dist1, nn_b1] = findNearestNeighbors(xc, yc, buffers1, buffer_ids1);
	auto [nn_dist2, nn_b2] = findNearestNeighbors(xc, yc, buffers2, buffer_ids2);

	//Perform iterations 
	int iter = 0;
	const double res_threshold = 0.2;
	double res_min = 0, res_max = 0, res_mean = 0;
	Matrix <double> Ex(n, 1), Ey(n, 1);
	while (iter < max_iter)
	{
		// Process all affected points of the input spline
		for (int j = 0; j < n; j++)
		{
			//Partial derivatives of the outer energy
			auto [eox, eoy] = getEOxy3(xc(j, 0), yc(j, 0), nn_b1[j]->getX(), nn_b1[j]->getY(), nn_b2[j]->getX(), nn_b2[j]->getY());
		
			//Partial derivatives of the regularization energy
			auto [erx, ery] = getErxy3(xc(std::max(j - 1, 0), 0), yc(std::max(j - 1, 0), 0), xc(j, 0), yc(j, 0), xc(std::min(j + 1, n - 1), 0), yc(std::min(j + 1, n - 1), 0));
			
			//Compute energy impact of the nearest point of the buffer
			//Minimize distances between points of the simplified contour line
			Ex(j, 0) = eox + erx;
			Ey(j, 0) = eoy + ery;

			//std::cout << eox << ' ' << erx << ' ' << eoy << ' ' << ery << '\n';
			
		}
		//std::cout << "\n";

		//Compute right side of the equation
		Matrix <double> ssx = delta * dx - kappa * Ex;
		Matrix <double> ssy = delta * dy - kappa * Ey;

		//Compute new shifs
		dx = Ainv * ssx;
		dy = Ainv * ssy;

		//New position of the vertices
		xc = xc0 + dx;
		yc = yc0 + dy;
		
		//Actualize matrices
		if (iter % 25 == 0)
		{
			//Coordinate differences
			Matrix <double> deltax = xc - xck;
			Matrix <double> deltay = yc - yck;

			//Compute residuals
			Matrix <double> R = sqrtm(deltax % deltax + deltay % deltay);
			res_min = min(R);
			res_mean = mean(R);
			res_max = max(R);

			//std::cout << "[ #" << iter << ": res_min = " << res_min << ", res_mean = " << res_mean << ", res_max = " << res_max << " ]\n";

			//Assign solution
			xck = xc;
			yck = yc;
			
			//Actualize A
			const Matrix <double> diffx = diff(xc);
			const Matrix <double> diffy = diff(yc);
			const Matrix <double> dist = sqrtm(diffx % diffx + diffy % diffy);
			double d_aver = mean(dist);
		
			A = createA(alpha, beta, gamma, d_aver, n);
			//A = createAFull(xc, yc, alpha, beta, gamma, n);

			//Actualize inverze matrix
			Matrix <double> L(n, n), U(n, n), P(n, n, 0, 1), I(n, n, 0, 1);
			lu(A + delta * I, L, U, P, sign);
			Ainv = inv(L) * inv(U);
			
			//Chceck terminal condition
			if ((res_max < res_threshold) && (i > 0) || res_max / res_mean > 100)
				break;

			//NN search: update nearest points on both buffers 
			std::tie(nn_dist1, nn_b1) = findNearestNeighbors(xc, yc, buffers1, buffer_ids1);
			std::tie(nn_dist2, nn_b2) = findNearestNeighbors(xc, yc, buffers2, buffer_ids2);
		}
		
		//Increment iterations
		iter++;
	}

	std::cout << "[ #" << iter << ": res_min = " << res_min << ", res_mean = " << res_mean << ", res_max = " << res_max << " ]\n";

	//Convert spline to contour line
	for (int i = 0; i < xc.rows(); i++)
	{
		//Create new point
		std::shared_ptr<Point3D > p = std::make_shared<Point3D >(xc(i, 0), yc(i, 0), contour[0]->getZ());

		//Add point to the contour line
		contour_simplified.push_back(p);
	}

	return contour_simplified;
}


Matrix<double> MinimumEnergySpline::createA(const double alpha, const double beta, const double gamma, const double h, const unsigned int n)
{
	 //Create pentadiagonal matrix A of the non-closed spline
	 //Use average h step
	 Matrix <double> A(n, n);
	 
	 //Coefficients of the trifiagonal matrix
	 const double a = alpha + 2.0 * beta /(h * h) + 6.0 * gamma / (h * h * h * h);
	 const double b = -beta / (h * h) - 4.0 * gamma / (h * h * h * h);
	 const double c = gamma / (h * h * h * h);

	 //Main diagonal
	 for (int i = 0; i < n; i++)
		 A(i, i) = a;

	 //Subdiagonals 2,3
	 for (int i = 0; i < n - 1; i++)
	 {
		 A(i, i + 1) = b;
		 A(i + 1, i) = b;
	 }

	 //Subdiagonals 4,5
	 for (int i = 0; i < n - 2; i++)
	 {
		 A(i, i + 2) = c;
		 A(i + 2, i) = c;
	 }

	 return A;
}


 Matrix<double> MinimumEnergySpline::createAFull(const Matrix <double>& xc, const Matrix <double>& yc, const double alpha, const double beta, const double gamma, const unsigned int n)
 {
	 //Create pentadiagonal matrix A of the non-closed spline
	 Matrix <double> A(n, n);

	 //Compute H2, H4
	 Matrix <double> H2 = createH(xc, yc, 1);
	 Matrix <double> H4 = createH(xc, yc, 2);

	 //Main diagonal
	 for (int i = 0; i < n; i++)
	 {
		 const double a = alpha + 2.0 * beta / (H2(i, 0) * H2(i, 0)) + 6.0 * gamma / (H4(i, 0) * H4(i, 0) * H4(i, 0) * H4(i, 0));
		 A(i, i) = a;
	 }

	 //Subdiagonals 2,3
	 for (int i = 0; i < n - 1; i++)
	 {
		 const double b = -beta / (H2(i, 0) * H2(i, 0)) - 4.0 * gamma / (H4(i, 0) * H4(i, 0) * H4(i, 0) * H4(i, 0));
		 A(i, i + 1) = b;
		 A(i + 1, i) = b;
	 }

	 //Subdiagonals 4,5
	 for (int i = 0; i < n - 2; i++)
	 {
		 const double c = gamma / (H4(i, 0) * H4(i, 0) * H4(i, 0) * H4(i, 0));
		 A(i, i + 2) = c;
		 A(i + 2, i) = c;
	 }

	 return A;
 }


 Matrix<double> MinimumEnergySpline::createH(const Matrix <double>& xc, const Matrix <double>& yc, const int k)
 {
	 //Create vector h of the length step
	 int n = xc.rows();
	 Matrix <double> H(n, 1);

	 for (int i = 0; i < n; i++)
	 {
		 //Interval of indices
		 const int jmin = std::max(i - k, 0);
		 const int jmax = std::min(i + k, n - 1);

		 //Compute average h
		 double hi = 0;
		 for (int j = jmin + 1; j <= jmax; j++)
		 {
			 //Compute length
			 const double dx = xc(j, 0) - xc(j - 1, 0);
			 const double dy = yc(j, 0) - yc(j - 1, 0);
			 const double d = sqrt(dx * dx + dy * dy);

			 hi = hi + d;
		 }

		 H(i, 0) = hi / (jmax - jmin);
	 }

	 return H;
 }


std::pair<double, double>  MinimumEnergySpline::getEOxy(const double xc, const double yc, const double xb1, const double yb1, const double xb2, const double yb2)
{
	//Get partial derivatives of Eo according to x, y: abs(d1-d2)
	const double dxb1 = xc - xb1;
	const double dyb1 = yc - yb1;
	const double dxb2 = xc - xb2;
	const double dyb2 = yc - yb2;
	const double db1 = sqrt(dxb1 * dxb1 + dyb1 * dyb1);
	const double db2 = sqrt(dxb2 * dxb2 + dyb2 * dyb2);

	//Partial derivatives
	const double eox = (dxb1 / db1 - dxb2 / db2) * sgn(db1 - db2);
	const double eoy = (dyb1 / db1 - dyb2 / db2) * sgn(db1 - db2);

	return { eox, eoy };
}


std::pair<double, double> MinimumEnergySpline::getEOxy2(const double xc, const double yc, const double xb1, const double yb1, const double xb2, const double yb2)
{
	//Get partial derivatives of Eo according to x, y: (d1 - d2)^2
	const double dxb1 = xc - xb1;
	const double dyb1 = yc - yb1;
	const double dxb2 = xc - xb2;
	const double dyb2 = yc - yb2;
	const double db1 = sqrt(dxb1 * dxb1 + dyb1 * dyb1);
	const double db2 = sqrt(dxb2 * dxb2 + dyb2 * dyb2);

	//Partial derivatives
	const double eox = 2*(dxb1 / db1 - dxb2 / db2) * (db1 - db2);
	const double eoy = 2*(dyb1 / db1 - dyb2 / db2) * (db1 - db2);

	return { eox, eoy };
}


std::pair<double, double> MinimumEnergySpline::getEOxy3(const double xc, const double yc, const double xb1, const double yb1, const double xb2, const double yb2)
{
	//Get partial derivatives of Eo according to x, y: (d1^2 - d2^2)/d12
	const double dxb1 = xc - xb1;
	const double dyb1 = yc - yb1;
	const double dxb2 = xc - xb2;
	const double dyb2 = yc - yb2;
	const double dxb12 = xb2 - xb1;
	const double dyb12 = yb2 - yb1;
	const double db1 = sqrt(dxb1 * dxb1 + dyb1 * dyb1);
	const double db2 = sqrt(dxb2 * dxb2 + dyb2 * dyb2);
	const double db12 = sqrt(dxb12 * dxb12 + dyb12 * dyb12);

	//Partial derivatives
	const double eox = 2 * (dxb1 - dxb2) / db12 * sgn(db1 * db1 - db2 * db2);
	const double eoy = 2 * (dyb1 - dyb2) / db12 * sgn(db1 * db1 - db2 * db2);

	return { eox, eoy };
}


std::pair<double, double>  MinimumEnergySpline::getErxy(const double xc0, const double yc0, const double xc1, const double yc1, const double xc2, const double yc2)
{
	//Get partial derivatives of Er according to x, y: abs(dc1 - dc2)
	const double eps = 1.0e-3;

	const double dx1 = xc1 - xc0;
	const double dy1 = yc1 - yc0;
	const double dx2 = xc1 - xc2;
	const double dy2 = yc1 - yc2;
	const double dc1 = sqrt(dx1 * dx1 + dy1 * dy1);
	const double dc2 = sqrt(dx2 * dx2 + dy2 * dy2);

	//Partial derivatives
	const double erx = (dc1 > eps) && (dc2 > eps) ? (dx1 / dc1 - dx2 / dc2) * sgn(dc1 - dc2) : 0;
	const double ery = (dc1 > eps) && (dc2 > eps) ? (dy1 / dc1 - dy2 / dc2) * sgn(dc1 - dc2) : 0;

	return { erx, ery };
}


std::pair<double, double>  MinimumEnergySpline::getErxy2(const double xc0, const double yc0, const double xc1, const double yc1, const double xc2, const double yc2)
{
	//Get partial derivatives of Er according to x, y: (dc1-dc2)^2
	const double eps = 1.0e-3;

	const double dx1 = xc1 - xc0;
	const double dy1 = yc1 - yc0;
	const double dx2 = xc1 - xc2;
	const double dy2 = yc1 - yc2;
	const double dc1 = sqrt(dx1 * dx1 + dy1 * dy1);
	const double dc2 = sqrt(dx2 * dx2 + dy2 * dy2);

	//Partial derivatives
	const double erx = (dc1 > eps) && (dc2 > eps) ? 2.0*(dx1 / dc1 - dx2 / dc2) * (dc1 - dc2) : 0;
	const double ery = (dc1 > eps) && (dc2 > eps) ? 2.0*(dy1 / dc1 - dy2 / dc2) * (dc1 - dc2) : 0;

	return { erx, ery };
}


std::pair<double, double>  MinimumEnergySpline::getErxy3(const double xc0, const double yc0, const double xc1, const double yc1, const double xc2, const double yc2)
{
	//Get partial derivatives of Er according to x, y: (dc1-dc2)^2
	const double eps = 1.0e-3;

	const double dx1 = xc1 - xc0;
	const double dy1 = yc1 - yc0;
	const double dx2 = xc1 - xc2;
	const double dy2 = yc1 - yc2;
	const double dx12 = xc2 - xc0;
	const double dy12 = yc2 - yc0;
	const double dc1 = sqrt(dx1 * dx1 + dy1 * dy1);
	const double dc2 = sqrt(dx2 * dx2 + dy2 * dy2);
	const double dc12 = sqrt(dx12 * dx12 + dy12 * dy12);

	//Partial derivatives
	const double erx = (dc1 > eps) && (dc2 > eps) ? 2.0 * (dx1 - dx2) / dc12 * sgn(dc1 * dc1 - dc2 * dc2) : 0;
	const double ery = (dc1 > eps) && (dc2 > eps) ? 2.0 * (dy1 - dy2) / dc12 * sgn(dc1 * dc1 - dc2 * dc2) : 0;

	return { erx, ery };
}


std::pair<TVector <float>, TVector <std::shared_ptr <Point3D > > > MinimumEnergySpline::findNearestNeighbors(const Matrix <double>& xq, const Matrix <double>& yq, const TVector2D <std::shared_ptr <Point3D> >& buffers, TVector <int>& buffer_ids)
{
	//Find nearest neighbor to any contour line vertex
	const int n = xq.rows();

	TVector <float> nn_dists(n, MAX_FLOAT);
	TVector <std::shared_ptr <Point3D > > nn_points(n);

	//Process all query points
	for (int i = 0; i < n; i++)
	{
		//Check specific buffer fragments
		const int j_start = (buffer_ids[i] == -1 ? 0 : buffer_ids[i]);
		const int j_end = (buffer_ids[i] == -1 ? buffers.size() : j_start + 1);

		for (int j = j_start; j < j_end; j++)
		{
			//Find nearest line segment point
			const auto [d_min, xi_min, yi_min] = ContourLinesSimplify::getNearestLineSegmentPoint(xq(i, 0), yq(i, 0), buffers[j]);

			//We found a closer point
			if (d_min < nn_dists[i])
			{
				//Create nearest point
				std::shared_ptr <Point3D> p_min = std::make_shared<Point3D>(xi_min, yi_min);

				//Actualize lists
				nn_dists[i] = d_min;
				nn_points[i] = p_min;
				buffer_ids[i] = j;
			}
		}
	}

	return { nn_dists, nn_points };
}



#endif