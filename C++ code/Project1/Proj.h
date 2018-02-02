#pragma once
#include <math.h>
#include <iostream>
#include <Eigen/Dense>
#include <stdio.h>
#include <fstream>

using namespace std;
using namespace Eigen;

class Proj
{
	double W;
	double L;
	double B;
	double dx;
	double dy;
	int M;
	int N;
	double dt;
	double mu;
	double rho;
	double beta;
	double D;

public:
	Proj(double W, double L, double B, double dx, double dy, double dt, double mu, double rho, double D) {
		Proj::W = W;
		Proj::L = L;
		Proj::B = B;
		Proj::dx = dx;
		Proj::dy = dy;
		M = round(L / dx);
		N = round(B / dy);
		Proj::dt = dt;
		Proj::mu = mu;
		Proj::rho = rho;
		beta = pow(dx / dy, 2);
		Proj::D = D;
	}

	int get_M() {
		return M;
	}
	int get_N() {
		return N;
	}

	double uf_calc(double, double, double, double, double, double, double, double, double, double, double);
	double vf_calc(double, double, double, double, double, double, double, double, double, double, double);
	//tuple<MatrixXd, MatrixXd, MatrixXd, MatrixXd, MatrixXd> solve();
	void solve();
	//double **mymesh(int r, int c, double dist, double dist2);
	//tuple<MatrixXd, MatrixXd> mymesh(int r, int c, double dist, double dist2, double, double);
	double bij(double uij, double uim1j, double vij, double vijm1);  // 
	double c_calc(double uij, double uim1j, double vij, double vijm1,
		double cij, double cip1j, double cim1j, double cijp1, double cijm1);
};