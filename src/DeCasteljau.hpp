//Description: Bezier cubic construction using deCasteljau algorithm

// Copyright (c) 2010 - 2021
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

#ifndef DeCasteljau_HPP
#define DeCasteljau_HPP

template <typename Point>
Point DeCasteljau::computeBezierPoint(const double u, const TVector <Point> &points)
{
	//Compute point on Bezier curve using de Casteljau algorithm
	TVector <Point> pointst;
	for (auto p : points)
		pointst.push_back(std::make_shared<Point3D>(*p));

	///Perform iterations
	for (int i = 1; i < points.size(); i++)
	{
		for (int j = 0; j < points.size() - i; j++)
		{
			pointst[j]->setX((1 - u) * pointst[j]->getX() + u * pointst[j+1]->getX());
			pointst[j]->setY((1 - u) * pointst[j]->getY() + u * pointst[j+1]->getY());
		}
	}

	//Return Bezier point
	return pointst[0];
}


std::tuple <double, double> DeCasteljau::computeBezierPoint(const double u, const TVector <double> &xp, const TVector <double> &yp)
{
	//Compute point on Bezier curve using de Casteljau algorithm
	TVector <double> xt = xp, yt = yp;

	//Perform iterations
	for (int i = 1; i < xp.size(); i++)
	{
		for (int j = 0; j < xp.size() - i; j++)
		{
			xt[j] = (1 - u) * xt[j] + u * xt[j + 1];
			yt[j] = (1 - u) * yt[j] + u * yt[j + 1];
		}
	}
	
	//Return Bezier point
	return { xt[0] , yt[0] };
}


#endif