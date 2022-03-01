/*
 * Serial_CGM.h
 *
 *  Created on: Apr 15, 2021
 *      Author: ubuntu
 */

#ifndef SERIAL_CGM_H_
#define SERIAL_CGM_H_
#include<vector>
#include<iostream>
#include<iomanip>
#include<cmath>
#include<stdio.h>
#define TOL 0.000001
typedef float float_t;
typedef int int_t;

float_t INNER_PRODUCT(std::vector<float_t> &a,std::vector<float_t>&b);
std::vector<float_t>  Residual(std::vector<std::vector<float_t>>&A,std::vector<float_t> &x,std::vector<float_t>& b);
std::vector<float_t> Serial_CG(std::vector<std::vector<float_t>>& A,std::vector<float_t>& b);
std::vector<float_t>  MATRIX_VECTOR_PRODUCT(std::vector<std::vector<float_t>>&A, std::vector<float_t> &x);
std::vector<std::vector<float_t>> MATRIX_MATRIX_MULTIPLY(std::vector<std::vector<float_t>> &a, std::vector<std::vector<float_t>> &b);
std::vector<float_t> saxpy(float_t a, std::vector<float_t> &x, std::vector<float_t> &y);
void displayMat(std::vector<std::vector<float_t>>&A);
void displayVec(std::vector<float_t> &b);
#endif /* SERIAL_SERIAL_CGM_H_ */
