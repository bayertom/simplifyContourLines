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


#ifndef MinimumEnergySpline_H
#define MinimumEnergySpline_H

#include <cmath>
#include <tuple>

#include "Matrix.h"
#include "MatrixOperations.h"

class MinimumEnergySpline
{
        public:
                template <typename T> 
                static int sgn(T val) {
                        return (T(0) < val) - (val < T(0));
                }

                template <typename T>
                static std::tuple<Matrix<T>, Matrix<T> > createSpline(const Matrix <T>& xc, const Matrix <T>& yc, const Matrix <T>& xb1, const Matrix <T>& yb1, const Matrix <T>& xb2, const Matrix <T>& yb2, const TVector <size_t>& knn_id1, const TVector <size_t>& knn_id2, const T alpha = 0.00001, const T beta = 0.00001, const T gamma = 15, const T delta = 0, const bool var_delta = true, const T kappa = 1.0, const unsigned int max_iter = 600);

        private:

                template <typename T>
                static Matrix<T> createDelta(const Matrix <T>& xc, const Matrix <T>& yc, const Matrix <T>& xb1, const Matrix <T>& yb1, const Matrix <T>& xb2, const Matrix <T>& yb2, const T delta, const T daver, const bool variable_delta);

                template <typename T>
                static Matrix<T> createA(const T alpha, const T beta, const Matrix <T> &D, const T h, const unsigned int n);

                template <typename T>
                static Matrix<T> createAFull(const Matrix <T>& xc, const Matrix <T>& yc, const T alpha, const T beta, const Matrix <T>& D, const unsigned int n);

                static std::tuple <double, unsigned int> getNearestPoint(const double x, const double y, const Matrix <double> &xp, const Matrix <double> &yp);

                template <typename T>
                static void getNearestLineSegmentPoints(const Matrix <T>& xq, const Matrix <T>& yq, const Matrix <T>& xp, const Matrix <T>& yp, TVector <T> &xn, TVector <T> &yn, TVector <T>& nn_dist);

                static std::tuple<double, unsigned int, double, double> getNearestLineSegmentPoint(const double xc, const double yc, const Matrix <double> & xp, const Matrix <double>& yp);

                template <typename T>
                static void getEOxy(const T xc, const T yc, const T xb1, const T yb1, const T xb2, const T yb2, T &ex, T &ey);

                template <typename T>
                static void getEOxy2(const T xc, const T yc, const T xb1, const T yb1, const T xb2, const T yb2, T& ex, T& ey);

                static std::tuple<double, double> getEOxy3(const double xc, const double yc, const double xb1, const double yb1, const double xb2, const double yb2);

                template <typename T>
                static  void getErxy(const T xc0, const T yc0, const T xc1, const T yc1, const T xc2, const T yc2, T& ex, T& ey);

                template <typename T>
                static void getErxy2(const T xc0, const T yc0, const T xc1, const T yc1, const T xc2, const T yc2, T& ex, T& ey);

                
                static std::tuple<double, double>  getErxy3(const double xc0, const double yc0, const double xc1, const double yc1, const double xc2, const double yc2);

                template <typename T>
                static void getAverageForceVector(const T xc, const T yc, const Matrix <T>& xb, const Matrix <T>& yb, T &uxa, T &uya);

                template <typename T>
                static Matrix <T> createH(const Matrix <T>& xc, const Matrix <T>& yc, const int k);


};

#include "MinimumEnergySpline.hpp"

#endif