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


#ifndef MatrixOperations_HPP
#define MatrixOperations_HPP

#include <cmath>

#include "TVector.h"
#include "TVector2D.h"
#include "Const.h"

#include "MathOverflowException.h"
#include "MathMatrixDifferentSizeException.h"
#include "MathMatrixSingularException.h"
#include "MathMatrixNotSquareException.h"

namespace MatrixOperations
{
	//Find min value, row and col position
	template <typename T>
	T min(const Matrix <T> &A, unsigned int & row, unsigned int & col)
	{
		const unsigned int m = A.rows(), n = A.cols();

		//Initialize minimum
		T min = A(0, 0);

		//Process all items
		for (unsigned int i = 0; i < m; i++)
		{
			for (unsigned int j = 0; j < n; j++)
			{
				if (A(i, j) < min)
				{
					min = A(i, j);
					row = i; col = j;
				}
			}
		}

		return min;
	}


	//Find min value
	template <typename T>
	T min(const Matrix <T> &A)
	{
		unsigned int row_index, column_index;
		return min(A, row_index, column_index);
	}


	//Find max value, row and col index
	template <typename T>
	T max(const Matrix <T> &A, unsigned int & row, unsigned int & col)
	{
		const unsigned int m = A.rows(), n = A.cols();

		//Initialize maximum
		T max = A(0, 0);

		//Process all items
		for (unsigned int i = 0; i < m; i++)
		{
			for (unsigned int j = 0; j < n; j++)
			{
				if (A(i, j) > max)
				{
					max = A(i, j);
					row = i; col = j;
				}
			}
		}

		return max;
	}


	//Find max value
	template <typename T>
	T max(const Matrix <T> &A)
	{
		unsigned int row_index, column_index;
		return max(A, row_index, column_index);
	}


	//Find max value, row and col index
	template <typename T>
	T mean(const Matrix <T>& A)
	{
		const unsigned int m = A.rows(), n = A.cols();

		//Initialize sum
		T sum = 0.0;

		//Process all items
		for (unsigned int i = 0; i < m; i++)
		{
			for (unsigned int j = 0; j < n; j++)
			{
				sum += A(i, j);
			}
		}

		return sum / (m * n);
	}


	//Abolute value of a matrix
	template <typename T>
	Matrix <T> diff(const Matrix <T>& A)
	{
		//Create trans matrix
		const unsigned int m = A.rows(), n = A.cols();
		Matrix <T> A_d(m - 1, n);

		//Copy m to m_trans
		for (unsigned int i = 0; i < m - 1; i++)
		{
			for (unsigned int j = 0; j < n; j++)
			{
				A_d(i, j) = A(i + 1, j) - A(i, j);
			}
		}

		return A_d;
	}


	//Sqrt of the matrix
	template <typename T>
	Matrix <T> sqrtm(const Matrix <T> &A)
	{
		const unsigned int m = A.rows(), n = A.cols();
		Matrix <T> B(m, n);

		for (unsigned int i = 0; i < m; i++)
		{
			for (unsigned int j = 0; j < n; j++)
			{
				B(i, j) = sqrt(A(i, j));
			}
		}

		return B;
	}


	//LU decomposition of the matrix
	template <typename T>
	void lu(const Matrix <T> &A, Matrix <T> &L, Matrix <T> &U, Matrix <T> &P, short &sign)
	{
		//LU decomposition of A, L = lower triangular matrix, U = upper triangular matrix, P = permutation matrix
		const unsigned int m = A.rows(), n = A.cols();

		//Set the determinant det(U) sign to 1
		sign = 1;

		//Is A rectangular matrix ?
		if (m != n)
		{
			throw MathMatrixNotSquareException <Matrix <T> >("MathMatrixNotSquareException: ", " invalid dimension of the matrix (rectangle matrix), can not perform LU decomposition; (rows_count, columns_count):  ", A);
		}

		//Create row permutation vector
		Matrix <unsigned int> PR(1, n);

		//Create scale vector
		Matrix <T> S(1, n);

		//Set diagonal items of L to 1, otherwise to 0
		//Set items of the row permutation matrix to <0; n-1>
		for (unsigned int i = 0; i < m; i++)
		{
			L(i, i) = 1.0;
			PR(0, i) = i;

			for (unsigned int j = 0; j < m; j++)
			{
				if (j != i)  L(i, j) = 0;

				P(i, j) = 0;
			}
		}

		//Initialize U = A
		U = A;

		//Find max item in each row to compute the scale vector
		for (unsigned int i = 0; i < n; i++)
		{
			T max_val = 0.0;

			for (unsigned int j = 0; j < n; j++)
			{
				if (fabs(U(i, j)) > max_val)
					max_val = fabs(U(i, j));
			}

			//Actualize scale vector
			if (max_val > MIN_FLOAT)
				S(0, i) = 1.0 / max_val;
		}

		//Start LU decomposition
		for (unsigned int j = 0; j < n; j++)
		{
			for (unsigned int i = 0; i < j; i++)
			{
				T sum = U(i, j);

				//Compute new U ( i, j ) item: multiply ith row and j-th column
				for (unsigned int k = 0; k < i; k++) sum -= U(i, k) * U(k, j);

				U(i, j) = sum;
			}

			//Initialize max_val and pivot index
			T max_val = 0.0;
			unsigned int i_pivot = n;

			//Find row that will be swapped and actualize row index
			for (unsigned int i = j; i < n; i++)
			{
				T sum = U(i, j);

				//Compute new U ( i, j ) item: multiply ith row and j-th column
				for (unsigned int k = 0; k < j; k++) sum -= U(i, k) * U(k, j);

				//Compute new U (i, j)
				U(i, j) = sum;

				//Compute index of the pivot
				const T val = S(0, i) * fabs(sum);

				if (val >= max_val)
				{
					max_val = val;
					i_pivot = i;
				}
			}

			//Perform row swaps in U,PR: j <-> i_pivot
			if ((j != i_pivot) && (i_pivot < n))
			{
				//Perform swap in U matrix
				const Matrix <T> U_temp = U(i_pivot, i_pivot, 0, n - 1);
				U(U(j, j, 0, n - 1), i_pivot, 0);
				U(U_temp, j, 0);

				//Perform swap in the row permutation matrix
				const unsigned int perm_temp = PR(0, i_pivot);
				PR(0, i_pivot) = PR(0, j);
				PR(0, j) = perm_temp;

				//Actualize also the scale vector
				S(0, i_pivot) = S(0, j);

				//Actualize sign of the determinant det(U)
				sign *= -1;
			}

			//Change diagonal item U ( j, j ) = 0 to "small" value before the devision
			if (U(j, j) == 0.0)
				U(j, j) = MIN_FLOAT;

			//Actualize U (i, j) from diagonal items
			if (j != n - 1)
			{
				const T val = 1.0 / U(j, j);

				for (unsigned int i = j + 1; i < n; i++)
					U(i, j) *= val;
			}
		}

		//Process L matrix together with U matrix
		for (unsigned int i = 0; i < n; i++)
		{
			for (unsigned int j = 0; j < i; j++)
			{
				L(i, j) = U(i, j);
				U(i, j) = 0.0;
			}

			//Actualize permutation matrix from the row permutation matrix
			P(i, PR(0, i)) = 1.0;
		}
	}


	//Inverse matrix calculation using LU decomposition
	template <typename T>
	Matrix <T> inv(const Matrix <T> &A)
	{
		const unsigned int m = A.rows(), n = A.cols();

		//Rectangular matrix
		if (m != n)
		{
			throw MathMatrixNotSquareException <Matrix <T> >("MathMatrixNotSquareException: ", " invalid dimension of the matrix (rectangle matrix), can not compute inverse matrix; (rows_count, columns_count):  ", A);
		}

		//Find maximum
		const T max_val = max(A);

		if (max_val > MAX_DOUBLE)
			throw MathOverflowException <T>("MathOverflowException: bad scaled matrix, can not compute inverse matrix. ", "Max item > MAX_FLOAT.", max_val);

		//Create LU decomposition
		Matrix <T> L(m, m);
		Matrix <T> U(m, m);
		Matrix <T> P(m, m, 0, 1);
		short sign = 1;
		lu(A, L, U, P, sign);

		//Compute X = L^-1 (lower triangular matrix)
		Matrix <T> X(m, m);

		for (unsigned int j = 0; j < m; j++)
		{
			X(j, j) = 1.0;

			for (unsigned int i = j + 1; i < m; i++)
			{
				T sum = 0;

				for (unsigned int k = j; k <= i - 1; k++)
				{
					sum -= L(i, k) * X(k, j);
				}

				X(i, j) = sum;
			}
		}

		//Compute Y = U^-1 (upper triangular matrix)
		Matrix <T> Y(m, m);

		for (unsigned int j = 0; j < m; j++)
		{
			Y(j, j) = 1 / U(j, j);

			for (int i = j - 1; i >= 0; i--)
			{
				T sum = 0.0;

				for (unsigned int k = i + 1; k <= j; k++)
				{
					sum -= U(i, k) * Y(k, j) / U(i, i);
				}

				Y(i, j) = sum;
			}
		}

		//Compute inverse matrix A^-1 = U^-1 * L^-1 = X * Y * P
		return Y * X * P;
	}
}


#endif
