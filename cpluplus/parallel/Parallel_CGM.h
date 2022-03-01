/*
 * Serial_CGM.h
 *
 *  Created on: Apr 15, 2021
 *      Author: ubuntu
 */

#ifndef PARALLEL_CGM_H_
#define PARALLEL_CGM_H_
#include "Serial_CGM.h"
#include "mpi.h"
#include <omp.h>
#define TOL 0.000001
#define _NUM_THREADS 4
typedef float float_t;
typedef int int_t;

float_t Parallel_INNER_PRODUCT(std::vector<float_t> &a,std::vector<float_t>&b);
std::vector<float_t>  Parallel_Residual(std::vector<std::vector<float_t>>&A,std::vector<float_t> &x,std::vector<float_t>& b);
std::vector<float_t> Parallel_CG(std::vector<std::vector<float_t>>& A,std::vector<float_t>& b);
std::vector<float_t>  Parallel_MATRIX_VECTOR_PRODUCT(std::vector<std::vector<float_t>>&A, std::vector<float_t> &x);
//std::vector<std::vector<float_t>> Parallel_MATRIX_MATRIX_MULTIPLY(std::vector<std::vector<float_t>> &a, std::vector<std::vector<float_t>> &b);
std::vector<float_t> Parallel_saxpy(float_t a, std::vector<float_t> &x, std::vector<float_t> &y);

#endif /* SERIAL_SERIAL_CGM_H_ */
