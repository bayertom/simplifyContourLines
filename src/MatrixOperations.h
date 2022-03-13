// Description: Basic matrix operations

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

#ifndef MatrixOperations_H
#define MatrixOperations_H

#include <cmath>

#include "Matrix.h"

//Define namespace: basic matrix operations
namespace MatrixOperations
{      
        template <typename T>
	T min(const Matrix <T> &A, unsigned int & row, unsigned int & col);
        
        template <typename T>
	T min(const Matrix <T> &A);
       
        template <typename T>
	T max(const Matrix <T> &A, unsigned int & row, unsigned int & col);
        
        template <typename T>
	T max(const Matrix <T> &A);
       
	template <typename T>
	T mean(const Matrix <T>& A);

	template <typename T>
	Matrix <T> diff(const Matrix <T>& A);
               
        template <typename T>
	Matrix <T> sqrtm(const Matrix <T> &A);
        
        template <typename T>
	void lu(const Matrix <T> &A, Matrix <T> &L, Matrix <T> &U, Matrix <T> &P, short &sign);
       	
        template <typename T>
	Matrix <T> inv(const Matrix <T> &A);
};

#include "MatrixOperations.hpp"

#endif
