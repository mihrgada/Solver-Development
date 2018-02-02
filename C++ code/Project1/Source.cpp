#include "stdafx.h"

int main()
{
	double W = 0.1e-3;  // m
	double L = 36 * W;
	double B = 14 * W;
	double dx = W / 5;
	double dy = W / 5;
	double dt = 1e-4;  // s
	double mu = 8.9e-4;  // Pa.s
	double rho = 998;  // kg/m3
	double D = 5e-10;  // m2/s

	Proj one(W, L, B, dx, dy, dt, mu, rho, D);
	int M = one.get_M();
	int N = one.get_N();
	/*
	auto z = one.mymesh(M, N, dx / 2, dx, dy / 2, dy);
	MatrixXd Xmat = get<0>(z);
	MatrixXd Ymat = get<1>(z);
	*/
	one.solve();
	//auto z2 = one.solve();
	/*
	MatrixXd u = get<0>(z2);
	MatrixXd v = get<1>(z2);
	MatrixXd P = get<2>(z2);
	MatrixXd u_avg = get<3>(z2);
	MatrixXd v_avg = get<4>(z2);
	*/
	//cout << P << "\n";

	/*
	for (int j = 0; j < N; j++)
	{
	for (int i = 0; i < M; i++)
	{
	printf(" %lf", Xmat(j, i));
	}
	printf("\n");
	}


	for (int j = 0; j < N; j++)
	{
	for (int i = 0; i < M; i++)
	{
	printf(" %lf",Ymat(j,i));
	}
	printf("\n");
	}
	*/

	/*
	ofstream file_P("P.csv");
	ofstream file_X("X.csv");
	ofstream file_Y("Y.csv");
	ofstream file_uav("u_avg.csv");
	ofstream file_vav("v_avg.csv");

	for (int i = 0; i < M; i++) {
	for (int j = 0; j < N; j++) {
	file_P << P(i, j) << ",";
	file_X << Xmat(i, j) << ",";
	file_Y << Ymat(i, j) << ",";
	file_uav << u_avg(i, j) << ",";
	file_vav << v_avg(i, j) << ",";
	}
	file_P << "\n";
	file_X << "\n";
	file_Y << "\n";
	file_uav << "\n";
	file_vav << "\n";
	}

	file_P.close();
	file_X.close();
	file_Y.close();
	file_uav.close();
	file_vav.close();

	ofstream file_u("u.csv");

	for (int i = 0; i < M + 1; i++) {
	for (int j = 0; j < N; j++) {
	file_u << u(i, j) << ",";
	}
	file_u << "\n";
	}
	file_u.close();
	*/
	/*
	ofstream file_P("P_g.csv");
	ofstream file_X("u_g.csv");
	ofstream file_Y("v_g.csv");

	for (int i = 0; i < M+2; i++) {
	for (int j = 0; j < N+1; j++) {
	//file_P << P(i, j) << ",";
	//file_X << Xmat(i, j) << ",";
	file_Y << v(i, j) << ",";
	}
	//file_P << "\n";
	//file_X << "\n";
	file_Y << "\n";
	}

	file_P.close();
	file_X.close();
	file_Y.close();
	*/
	return 0;

}
