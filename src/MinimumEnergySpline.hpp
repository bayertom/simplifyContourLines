// Description: Minimum energy spline (snake), discrete approximation

// Copyright (c) 2015 - 2021
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


#ifndef MinimumEnergySpline_HPP
#define MinimumEnergySpline_HPP

#include <filesystem>

#include "MatrixOperations.h"
#include "PointLineDistance.h"


//Set namespaces
using namespace MatrixOperations;


std::tuple<Matrix<double>, Matrix<double> > MinimumEnergySpline::createSpline(const Matrix <double> & xc0, const Matrix <double> & yc0, const Matrix <double> & xb1, const Matrix <double> & yb1, const Matrix <double> & xb2, const Matrix <double> & yb2, const TVector <size_t>& knn_id1, const TVector <size_t>& knn_id2, const double alpha, const double beta, const double gamma, const double delta, const bool var_delta, const double kappa, const unsigned int max_iter)
{
	//Create minimum energy spline
	const int n = xc0.rows();

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
	
	//Create matrix
	Matrix <double> D = createDelta(xc, yc, xb1, yb1, xb2, yb2, delta, d_aver, var_delta);

	//Create A matrix
	Matrix <double> A = createA(alpha, beta, D, d_aver, n);
	//Matrix <T> A = createAFull(xc, yc, alpha, beta, D, n);

	//Inverse using LU decomposition
	short sign = 1;
	Matrix <double> L(n, n), U(n, n), P(n, n, 0, 1), I(n, n, 0, 1);
	lu(A + gamma * I, L, U, P, sign);
	Matrix <double> Ainv = inv(L) * inv(U);
	
	//NN search
	TVector <double> nn_dist1(n), nn_dist2(n), xn1(n), yn1(n), xn2(n), yn2(n);
	getNearestLineSegmentPoints(xc, yc, xb1, yb1, xn1, yn1, nn_dist1);
	getNearestLineSegmentPoints(xc, yc, xb2, yb2, xn2, yn2, nn_dist2);

	//Perform iterations 
	int  i = 0;
	const double res_threshold = 0.2;
	double res_min = 0, res_max = 0, res_mean = 0;
	Matrix <double> Ex(n, 1), Ey(n, 1);
	while (i < max_iter)
	{
		// Process all affected points of the input spline
		for (int j = 0; j < n; j++)
		{
			//Find nearest point of the simplified contour line on the first buffer
			unsigned int i1min; //= nn_id1[j];
			double d1min, xb1min = xn1[j], yb1min = yn1[j];
			//getNearestLineSegmentPoint(xc(j, 0), yc(j, 0), xb1, yb1, d1min, i1min, xb1min, yb1min);
			
			//Find nearest point of the simplified contour line on the second buffer
			unsigned int i2min; //= nn_id2[j];
			double d2min, xb2min = xn2[j], yb2min = yn2[j];
			//getNearestLineSegmentPoint(xc(j, 0), yc(j, 0), xb2, yb2, d2min, i2min, xb2min, yb2min);

			//Partial derivatives of the outer energy
			auto [eox, eoy] = getEOxy3(xc(j, 0), yc(j, 0), xb1min, yb1min, xb2min, yb2min);
		
			//Partial derivatives of the regularization energy
			auto [erx, ery] = getErxy3(xc(std::max(j - 1, 0), 0), yc(std::max(j - 1, 0), 0), xc(j, 0), yc(j, 0), xc(std::min(j + 1, n - 1), 0), yc(std::min(j + 1, n - 1), 0));
			
			//Compute energy impact of the nearest point of the buffer
			//Minimize distances between points of the simplified contour line
			Ex(j, 0) = eox + erx;
			Ey(j, 0) = eoy + ery;

			//std::cout << eox << ' ' << erx << ' ' << eoy << ' ' << ery << '\n';
		}

		//Compute right side of the equation
		Matrix <double> ssx = gamma * dx - kappa * Ex;
		Matrix <double> ssy = gamma * dy - kappa * Ey;

		//Compute new shifs
		dx = Ainv * ssx;
		dy = Ainv * ssy;

		//New position of the vertices
		xc = xc0 + dx;
		yc = yc0 + dy;

		//Actualize matrices
		if (i % 25 == 0)
		{
			//Coordinate differences
			Matrix <double> deltax = xc - xck;
			Matrix <double> deltay = yc - yck;

			//Compute residuals
			Matrix <double> R = sqrtm(deltax % deltax + deltay % deltay);
			res_min = min(R);
			res_mean = mean(R);
			res_max = max(R);

			//std::cout << "[ #" << i << ": res_min = " << res_min << ", res_mean = " << res_mean << ", res_max = " << res_max << " ]\n";

			//Assign solution
			xck = xc;
			yck = yc;
			
			//Actualize A
			const Matrix <double> diffx = diff(xc);
			const Matrix <double> diffy = diff(yc);
			const Matrix <double> dist = sqrtm(diffx % diffx + diffy % diffy);
			double d_aver = mean(dist);
		
			A = createA(alpha, beta, D, d_aver, n);
			//A = createAFull(xc, yc, alpha, beta, D, n);

			//Actualize inverze matrix
			Matrix <double> L(n, n), U(n, n), P(n, n, 0, 1), I(n, n, 0, 1);
			lu(A + gamma * I, L, U, P, sign);
			Ainv = inv(L) * inv(U);
			
			//Chceck terminal condition
			if ((res_max < res_threshold) && (i > 0) || res_max / res_mean > 100)
				break;

			//NN search
			getNearestLineSegmentPoints(xc, yc, xb1, yb1, xn1, yn1, nn_dist1);
			getNearestLineSegmentPoints(xc, yc, xb2, yb2, xn2, yn2, nn_dist2);
		}

		//Increment i
		i++;
	}

	std::cout << "[ #" << i << ": res_min = " << res_min << ", res_mean = " << res_mean << ", res_max = " << res_max << " ]\n";

	return { xc, yc };
}


Matrix<double> MinimumEnergySpline::createDelta(const Matrix <double>& xc, const Matrix <double>& yc, const Matrix <double>& xb1, const Matrix <double>& yb1, const Matrix <double>& xb2, const Matrix <double>& yb2, const double delta, const double daver, const bool variable_delta)
{
	//Create delta matrix
	const int n = xc.rows();
	Matrix <double> D(n, n, 0, delta);

	if (variable_delta)
	{
		for (int i = 0; i < n; i++)
		{
			//Nearest point on the first buffer
			auto [d1min, i1min] = getNearestPoint(xc(i, 0), yc(i, 0), xb1, yb1);

			//Nearest point on the second buffer
			auto [d2min, i2min] = getNearestPoint(xc(i, 0), yc(i, 0), xb2, yb2);

			//Compute buffer width
			const double dx = xb1(i1min, 0) - xb2(i2min, 0);
			const double dy = yb1(i1min, 0) - yb2(i2min, 0);
			const double wb = sqrt(dx * dx + dy * dy);

			//Compute distances to the buffer
			const double dx1 = xc(i, 0) - xb1(i1min, 0);
			const double dy1 = yc(i, 0) - yb1(i1min, 0);
			const double dx2 = xc(i, 0) - xb2(i2min, 0);
			const double dy2 = yc(i, 0) - yb2(i2min, 0);
			const double d1 = sqrt(dx1 * dx1 + dy1 * dy1);
			const double d2 = sqrt(dx2 * dx2 + dy2 * dy2);

			//Compute delta dynamically
			const double deltai = log10(std::max(1.0, wb / daver / 3));
			D(i, i) = std::min(std::max(0.1, deltai), 2.0);
		}
	}

	return D;
}


Matrix<double> MinimumEnergySpline::createA(const double alpha, const double beta, const Matrix<double> &D, const double h, const unsigned int n)
{
	 //Create pentadiagonal matrix A of the non-closed spline
	 //Use average h step
	 Matrix <double> A(n, n);
	 
	 //Coefficients of the trifiagonal matrix
	 const double a = 2.0 * alpha /(h * h)+ 6.0 * beta / (h * h * h * h);
	 const double b = -alpha / (h * h) - 4.0 * beta / (h * h * h * h);
	 const double c = beta / (h * h * h * h);

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

	 //Subtract D
	 A = A + D;

	 return A;
}


 Matrix<double> MinimumEnergySpline::createAFull(const Matrix <double>& xc, const Matrix <double>& yc, const double alpha, const double beta, const Matrix <double>& D, const unsigned int n)
 {
	 //Create pentadiagonal matrix A of the non-closed spline
	 Matrix <double> A(n, n);

	 //Compute H2, H4
	 Matrix <double> H2 = createH(xc, yc, 1);
	 Matrix <double> H4 = createH(xc, yc, 2);

	 //Main diagonal
	 for (int i = 0; i < n; i++)
	 {
		 const double a = 2.0 * alpha / (H2(i, 0) * H2(i, 0)) + 6.0 * beta / (H4(i, 0) * H4(i, 0) * H4(i, 0) * H4(i, 0));
		 A(i, i) = a;
	 }

	 //Subdiagonals 2,3
	 for (int i = 0; i < n - 1; i++)
	 {
		 const double b = -alpha / (H2(i, 0) * H2(i, 0)) - 4.0 * beta / (H4(i, 0) * H4(i, 0) * H4(i, 0) * H4(i, 0));
		 A(i, i + 1) = b;
		 A(i + 1, i) = b;
	 }

	 //Subdiagonals 4,5
	 for (int i = 0; i < n - 2; i++)
	 {
		 const double c = beta / (H4(i, 0) * H4(i, 0) * H4(i, 0) * H4(i, 0));
		 A(i, i + 2) = c;
		 A(i + 2, i) = c;
	 }

	 //Subtract D
	 A = A + D;

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


 std::tuple <double, unsigned int> MinimumEnergySpline::getNearestPoint(const double x, const double y, const Matrix <double>& xp, const Matrix <double>& yp)
 {
	 //Get point on the line segment nearest to P[x, y]
	 unsigned int n = xp.rows(), i_min = 0;
	 double d_min = 1.0e16;

	 for (unsigned int i = 0; i < n - 1; i++)
	 {
		 //Compute distance
		 const double dx = xp(i, 0) - x;
		 const double dy = yp(i, 0) - y;
		 const double d = sqrt(dx * dx + dy * dy);

		 //Actualize minimum
		 if (d < d_min)
		 {
			 d_min = d; i_min = i;
		 }
	 }

	 return { d_min, i_min };
 }


 void MinimumEnergySpline::getNearestLineSegmentPoints(const Matrix <double>& xq, const Matrix <double>& yq, const Matrix <double>& xp, const Matrix <double>& yp, TVector <double>& xn, TVector <double>& yn, TVector <double>& nn_dist)
 {
	 //Find nerest neighbors
	 for (int i = 0; i < xq.rows(); i++)
	 {
		 //Find nearest neighbor on the line segment
		 auto [d_min, i_min, x_near, y_near] = getNearestLineSegmentPoint(xq(i, 0), yq(i, 0), xp, yp);

		 //Store nearest neighbor and its properties
		 xn[i] = x_near; yn[i] = y_near; nn_dist[i] = d_min;
	 }
 }


 std::tuple<double, unsigned int, double, double> MinimumEnergySpline::getNearestLineSegmentPoint(const double x, const double y, const Matrix <double>& xp, const Matrix <double>& yp)
 {
	 //Get point on the line segment nearest to P[x, y]
	 unsigned int n = xp.rows(), imin = 0;
	 double dmin = 1.0e16, xmin = 0., ymin = 0., dthreshold = 0.01;

	 //Initialize
	 for (unsigned int i = 0; i < n - 1; i++)
	 {
		 double xi = 0, yi = 0;
		 const double d = std::fabs(PointLineDistance::getPointLineSegmentDistance2D(x, y, xp(i, 0), yp(i, 0), xp(i + 1, 0), yp(i + 1, 0), xi, yi));

		 //Avoid identical points
		 if ((d < dmin) && (d > dthreshold))
		 {
			 dmin = d; imin = i; xmin = xi; ymin = yi;
		 }
	 }

	 return { dmin, imin, xmin, ymin };
 }


std::tuple<double, double>  MinimumEnergySpline::getEOxy(const double xc, const double yc, const double xb1, const double yb1, const double xb2, const double yb2)
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


std::tuple<double, double> MinimumEnergySpline::getEOxy2(const double xc, const double yc, const double xb1, const double yb1, const double xb2, const double yb2)
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


std::tuple<double, double> MinimumEnergySpline::getEOxy3(const double xc, const double yc, const double xb1, const double yb1, const double xb2, const double yb2)
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


std::tuple<double, double>  MinimumEnergySpline::getErxy (const double xc0, const double yc0, const double xc1, const double yc1, const double xc2, const double yc2)
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


std::tuple<double, double>  MinimumEnergySpline::getErxy2(const double xc0, const double yc0, const double xc1, const double yc1, const double xc2, const double yc2)
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


std::tuple<double, double>  MinimumEnergySpline::getErxy3(const double xc0, const double yc0, const double xc1, const double yc1, const double xc2, const double yc2)
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




#endif