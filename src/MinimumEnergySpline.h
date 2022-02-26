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

                static std::tuple<Matrix<double>, Matrix<double> > createSpline(const Matrix <double> & xc, const Matrix <double> & yc, const Matrix <double> & xb1, const Matrix <double> & yb1, const Matrix <double> & xb2, const Matrix <double> & yb2, const TVector <size_t>& knn_id1, const TVector <size_t>& knn_id2, const double alpha = 0.00001, const double beta = 0.00001, const double gamma = 15, const double delta = 0, const bool var_delta = true, const double kappa = 1.0, const unsigned int max_iter = 600);

        private:

                static Matrix<double> createDelta(const Matrix <double>& xc, const Matrix <double>& yc, const Matrix <double>& xb1, const Matrix <double>& yb1, const Matrix <double>& xb2, const Matrix <double>& yb2, const double delta, const double daver, const bool variable_delta);
                static Matrix<double> createA(const double alpha, const double beta, const Matrix <double> &D, const double h, const unsigned int n);
                static Matrix<double> createAFull(const Matrix <double>& xc, const Matrix <double>& yc, const double alpha, const double beta, const Matrix <double>& D, const unsigned int n);
                static Matrix <double> createH(const Matrix <double>& xc, const Matrix <double>& yc, const int k);

                static  std::tuple <double, unsigned int> getNearestPoint(const double x, const double y, const Matrix <double> & xp, const Matrix <double> & yp);
                static void getNearestLineSegmentPoints(const Matrix <double>& xq, const Matrix <double>& yq, const Matrix <double>& xp, const Matrix <double>& yp, TVector <double> &xn, TVector <double> &yn, TVector <double>& nn_dist);
                static std::tuple<double, unsigned int, double, double> getNearestLineSegmentPoint(const double xc, const double yc, const Matrix <double> & xp, const Matrix <double> & yp);
                
        private:   
                //Outer and regularization energy, partial derivatives
                static std::tuple<double, double> getEOxy(const double xc, const double yc, const double xb1, const double yb1, const double xb2, const double yb2);
                static std::tuple<double, double> getEOxy2(const double xc, const double yc, const double xb1, const double yb1, const double xb2, const double yb2);
                static std::tuple<double, double> getEOxy3(const double xc, const double yc, const double xb1, const double yb1, const double xb2, const double yb2);
                
                static std::tuple<double, double> getErxy(const double xc0, const double yc0, const double xc1, const double yc1, const double xc2, const double yc2);
                static std::tuple<double, double> getErxy2(const double xc0, const double yc0, const double xc1, const double yc1, const double xc2, const double yc2);
                static std::tuple<double, double> getErxy3(const double xc0, const double yc0, const double xc1, const double yc1, const double xc2, const double yc2);

             

};

#include "MinimumEnergySpline.hpp"

#endif