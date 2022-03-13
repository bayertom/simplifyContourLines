// Description: Sort points stored in list according to their distance from a point

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


#ifndef sortPointsByDist_H
#define sortPointsByDist_H


//Sort points stored in list according to their distance from a poi
class sortPointsByDist
{
        private:
                const std::shared_ptr <Point3D> q;

        public:
                sortPointsByDist( const std::shared_ptr <Point3D> &q_ ) : q ( q_ ) {}

                bool operator() (const std::shared_ptr <Point3D>&p1, const std::shared_ptr <Point3D>&p2) const
                {
                        const double dx1 = p1->getX() - q->getX();
                        const double dy1 = p1->getY() - q->getY();
                        const double dx2 = p2->getX() - q->getX();
                        const double dy2 = p2->getY() - q->getY();
                        
                        return (dx1 * dx1 + dy1 * dy1 < dx2 * dx2 + dy2 * dy2);
                }
};



#endif