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

//Minimum energy spline (snake), discrete approximation
class MinimumEnergySpline
{
        public:
                template <typename T> 
                static int sgn(T val) {
                        return (T(0) < val) - (val < T(0));
                }

                static TVector <std::shared_ptr <Point3D > > createSpline(const TVector <std::shared_ptr <Point3D > >& contour, const TVector2D < std::shared_ptr<Point3D > >& buffers1, const TVector2D < std::shared_ptr<Point3D > >& buffers2, int i, int n, const double alpha = 0, const double beta = 0.00001, const double gamma = 0.00001, const double delta = 15, const double kappa = 1.0, const unsigned int max_iter = 600);

        private:
                //Matrices A nad H
                static Matrix <double> createA(const double alpha, const double beta, const double gamma, const double h, const unsigned int n);
                static Matrix <double> createAFull(const Matrix <double>& xc, const Matrix <double>& yc, const double alpha, const double beta, const double gamma, const unsigned int n);
                static Matrix <double> createH(const Matrix <double>& xc, const Matrix <double>& yc, const int k);

        private:   
                //Outer and regularization energy, partial derivatives
                static std::pair<double, double> getEOxy(const double xc, const double yc, const double xb1, const double yb1, const double xb2, const double yb2);
                static std::pair<double, double> getEOxy2(const double xc, const double yc, const double xb1, const double yb1, const double xb2, const double yb2);
                static std::pair<double, double> getEOxy3(const double xc, const double yc, const double xb1, const double yb1, const double xb2, const double yb2);
                
                static std::pair<double, double> getErxy(const double xc0, const double yc0, const double xc1, const double yc1, const double xc2, const double yc2);
                static std::pair<double, double> getErxy2(const double xc0, const double yc0, const double xc1, const double yc1, const double xc2, const double yc2);
                static std::pair<double, double> getErxy3(const double xc0, const double yc0, const double xc1, const double yc1, const double xc2, const double yc2);
        
        private:
                static std::pair<TVector <float>, TVector <std::shared_ptr <Point3D > > > findNearestNeighbors(const Matrix <double>& xq, const Matrix <double>& yq, const TVector2D <std::shared_ptr <Point3D> >& buffers, TVector <int>& buffer_ids);

             
};

#include "MinimumEnergySpline.hpp"

#endif