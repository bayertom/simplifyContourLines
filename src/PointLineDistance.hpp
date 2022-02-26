// Description: Distance point-loine

// Copyright (c) 2010 - 2013
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


#ifndef PointLineDistance_HPP
#define PointLineDistance_HPP

#include <cmath>


#include "EuclDistance.h"

#include "MathZeroDevisionException.h"


template <typename Point>
typename Point::Type PointLineDistance::getPointLineDistance2D ( const Point &p, const Point &p1, const Point &p2 )
{
        //Compute distance point and line
        return getPointLineDistance2D ( p->getX(), p->getY(), p1->getX(), p1->getY(), p2->getX(), p2->getY() );
}


template <typename T>
T PointLineDistance::getPointLineDistance2D(const T xa, const T ya, const T x1, const T y1, const T x2, const  T y2)
{
        //Compute unsigned distance point and line
        T xi, yi;
        return fabs(getPointLineDistance2DSigned(xa, ya, x1, y1, x2, y2, xi, yi));
}


template <typename T>
T PointLineDistance::getPointLineDistance2D ( const T xa, const T ya, const T x1, const T y1, const T x2, const  T y2, T &xi, T &yi )
{
        //Compute unsigned distance point and line
        return fabs ( getPointLineDistance2DSigned ( xa, ya, x1, y1, x2, y2, xi, yi ) );
}

template <typename T>
T PointLineDistance::getPointLineDistance2DSigned(const T xa, const T ya, const T x1, const T y1, const T x2, const  T y2)
{
        //Compute signed distance point and line
        T xi, yi;
        return getPointLineDistance2DSigned(xa, ya, x1, y1, x2, y2, xi, yi);
}


template <typename T>
T PointLineDistance::getPointLineDistance2DSigned ( const T xa, const T ya, const T x1, const T y1, const T x2, const  T y2, T& xi, T & yi)
{
        //Compute signed distance point and line
        //         If distance:
        //		> 0 = point lies in left half plane,
        //		< 0 = point lies in the right half plane,
        //		= 0 = point lies on the line,

        T d12 = sqrt ( ( x2 - x1 ) * ( x2 - x1 ) + ( y2 - y1 ) * ( y2 - y1 ) );

        //Throw exception
        if ( fabs (d12) < MIN_FLOAT )
        {
                //std::cout << "error";
                throw MathZeroDevisionException <T> ( "MathZeroDevisionException: can not compute distance point - line, ", "end points of the line are identical.", d12 );
        }

        const T d =  ( xa * ( y1 - y2 ) + x1 * ( y2 - ya ) + x2 * ( ya - y1 ) ) / d12;

        //Intersection point
        const T da1 = sqrt((xa - x1) * (xa - x1) + (ya - y1) * (ya - y1));
        const T k = sqrt(std::max(da1 * da1 - d * d, 0.0));
        xi = x1 + k * (x2 - x1) / d12;
        yi = y1 + k * (y2 - y1) / d12;

        return d;
}


template <typename Point>
typename Point::Type PointLineDistance::getPointLineSegmentDistance2D ( const Point &p, const Point &p1, const Point &p2 )
{
        //Compute distance point and line segment
        return getPointLineSegmentDistance2D( p->getX(), p->getY(), p1->getX(), p1->getY(), p2->getX(), p2->getY());
}

template <typename T>
T PointLineDistance::getPointLineSegmentDistance2D(const T x, const T y, const T x1, const T y1, const T x2, const T y2)
{
        //Compute distance point and line segment
        T xi, yi;
        return getPointLineSegmentDistance2D(x, y, x1, y1, x2, y2, xi, yi);
}

template <typename T>
T PointLineDistance::getPointLineSegmentDistance2D(const T x, const T y, const T x1, const T y1, const T x2, const T y2, T &xi, T &yi)
{
        //Compute distance point and line segment
        const T nx = y1 - y2, ny = x2 - x1;

        //Point p3 is on normal n1: given by p1 and perpendicular to (p1, p2)
        const T x3 = x1 + nx, y3 =  y1 + ny;

        //Point p4 is on normal n2: given by p2 and perpendicular to (p1, p2)
        const T x4 = x2 + nx, y4 = y2 + ny;

        //Position of the point according to the both normals using signed distance
        T xi3, yi3, xi4, yi4;
        const T dist1 = getPointLineDistance2DSigned(x, y, x1, y1, x3, y3, xi3, yi3);
        const T dist2 = getPointLineDistance2DSigned(x, y, x2, y2, x4, y4, xi4, yi4);

        //Point is between both normals n1 and n2
        if ((dist1 < 0) && (dist2 > 0) || (dist1 > 0) && (dist2 < 0))
        {
                return getPointLineDistance2D(x, y, x1, y1, x2, y2, xi, yi);
        }

        //Point is left to the n1
        if (dist1 >= 0)
        {
                xi = x1; yi = y1;
                return EuclDistance::getEuclDistance(x, y, 0.0, x1, y1, 0.0);
        }

        //Point is right to n2
        xi = x2; yi = y2;
        return EuclDistance::getEuclDistance(x, y, 0.0, x2, y2, 0.0);
}


#endif
