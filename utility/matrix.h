/******************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*******************************************************************************
* Author: Changhe Li
* Email: changhe.lw@gmail.com
* Language: C++
*-------------------------------------------------------------------------------
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*
*  see https://github.com/Changhe160/OFEC for more information
*
*-------------------------------------------------------------------------------
* class matrix generates orthogonal rotation matrixs
*
*********************************************************************************/

#ifndef OFEC_MATRIX_H
#define OFEC_MATRIX_H

#include "vector.h"

namespace OFEC {

	class matrix final {								// *****************Orthogonal rotation matrix***********************
	public:
		matrix(int dim = 0);
		matrix(const int c, const int r);
		~matrix();
		matrix & operator*=(const matrix & m);
		matrix & operator*=(double x);
		matrix & operator=(const matrix & m);
		bool identify();
		bool is_identity();
		void set_rotation_axes(int r, int c, double angle);
		void generate_rotation(normal *rand, double CondiNum);
		void randomize(uniform * rand, real min=-1, real max=1);
		void randomize(normal *rand);
		void orthonormalize();
		void zero();
		void set_row(const double *d, int c, int r = 1);
		void set_col(const double *d, const int r, const int c = 1);
		void set(const double * const * d);
		void diagonalize(normal * rand, double CondiNum);
		void transpose();
		void inverse();
		bool is_inverse();
		Vector & operator[](int idx);
		const Vector & operator[](int idx)const;
		void resize(int row, int col);
		void diagonalize(double alpha);
		bool is_diagonal();
		bool is_orthogonal();
		double determinant_(std::vector<Vector>& temp, int num);
		double determinant();
		std::vector<Vector>& data();
		/* Eigendecomposition of a matrix */
		void eigendecomposition();
		void Householder2(int n, double **V, double *d, double *e);
		void QLalgo2(int n, double *d, double *e, double **V);
		double myhypot(double a, double b);
		double* get_eigen_value();
		double** get_eigen_vector();
		/* * * * * * * * * * * * * * * * * */
		void load(std::ifstream &in);
		void print(std::ofstream & out);
			
	private:
		void clear();
	private:
		int m_col_size, m_row_size;									// matrix size
		std::vector<Vector> m_row;
		bool m_is_det_flag = false;
		double m_det;
		double* m_eigen_value;
		double** m_eigen_vector;
	};

}

#endif