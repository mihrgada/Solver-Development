#include "stdafx.h"

//tuple<MatrixXd, MatrixXd, MatrixXd, MatrixXd, MatrixXd> Proj::solve() {
void Proj::solve(){
	
	// convenience variable
	int n1 = round(W / dx);

	// Initialization (including ghost nodes)
	MatrixXd u_g = MatrixXd::Constant(M + 1, N + 2, 0.0);
	MatrixXd v_g = MatrixXd::Constant(M + 2, N + 1, 0.0);
	MatrixXd P_g = MatrixXd::Constant(M + 2, N + 2, 0.0);
	MatrixXd c_g = MatrixXd::Constant(M + 2, N + 2, 0.0);

	// Dirichlet BCs (BCs for u already initialized)
	for (int i = 1; i < n1 + 1; i++) {
		v_g(i, 0) = 0.001; // m/s
		v_g(i, N) = -0.001; // m/s
	}

	for (int i = 1; i < n1 + 1; i++) {
		c_g(i, 0) = 2.0;
	}

	MatrixXd P_up = P_g;  // at n+1
	MatrixXd P_t = P_g;  // iteration variable (temporary variable)
	MatrixXd c_up = c_g;  // at n+1

	// FSM METHOD
	// defining Vf (same size as u_g)
	MatrixXd uf = u_g;
	MatrixXd vf = v_g;
	/*
	MatrixXd u = u_g.block(0, 1, M + 1, N);
	MatrixXd v = v_g.block(1, 0, M, N + 1);
	MatrixXd P = P_g.block(1, 1, M, N);
	*/
	MatrixXd u = MatrixXd::Constant(M + 1, N, NAN);
	MatrixXd v = MatrixXd::Constant(M, N + 1, NAN);
	MatrixXd P = MatrixXd::Constant(M, N, NAN);
	MatrixXd c = MatrixXd::Constant(M, N, NAN);

	MatrixXd u_avg = MatrixXd::Constant(M, N, NAN);
	MatrixXd v_avg = MatrixXd::Constant(M, N, NAN);

	// CALCULATIONS

	for (int n = 0; n < 3e5 + 1; n++) {
	// calulating uf
	{
		// block 1
		for (int j = 1; j < N + 1; j++) {
			for (int i = 1; i < n1; i++) {				
				uf(i, j) = uf_calc(u_g(i, j), u_g(i + 1, j), u_g(i - 1, j), u_g(i, j + 1), u_g(i, j - 1),
					v_g(i, j), v_g(i + 1, j), v_g(i, j - 1), v_g(i + 1, j - 1), P_g(i, j), P_g(i + 1, j));		
			}
		}

		// block A1
		for (int j = 6 * n1 + 1; j < 8 * n1 + 1; j++) {
			for (int i = n1; i < 3 * n1 + 1; i++) {
				uf(i, j) = uf_calc(u_g(i, j), u_g(i + 1, j), u_g(i - 1, j), u_g(i, j + 1), u_g(i, j - 1),
					v_g(i, j), v_g(i + 1, j), v_g(i, j - 1), v_g(i + 1, j - 1), P_g(i, j), P_g(i + 1, j));
			}
		}

		// block C1
		for (int j = 6 * n1 + 1; j < N + 1; j++) {
			for (int i = 3 * n1 + 1; i < 5 * n1; i++) {
				uf(i, j) = uf_calc(u_g(i, j), u_g(i + 1, j), u_g(i - 1, j), u_g(i, j + 1), u_g(i, j - 1),
					v_g(i, j), v_g(i + 1, j), v_g(i, j - 1), v_g(i + 1, j - 1), P_g(i, j), P_g(i + 1, j));
			}
		}

		// block B1
		for (int j = 12 * n1 + 1; j < N + 1; j++) {
			for (int i = 5 * n1; i < 7 * n1 + 1; i++) {
				uf(i, j) = uf_calc(u_g(i, j), u_g(i + 1, j), u_g(i - 1, j), u_g(i, j + 1), u_g(i, j - 1),
					v_g(i, j), v_g(i + 1, j), v_g(i, j - 1), v_g(i + 1, j - 1), P_g(i, j), P_g(i + 1, j));
			}
		}

		// block D1
		for (int j = 6 * n1 + 1; j < N + 1; j++) {
			for (int i = 7 * n1 + 1; i < 9 * n1; i++) {
				uf(i, j) = uf_calc(u_g(i, j), u_g(i + 1, j), u_g(i - 1, j), u_g(i, j + 1), u_g(i, j - 1),
					v_g(i, j), v_g(i + 1, j), v_g(i, j - 1), v_g(i + 1, j - 1), P_g(i, j), P_g(i + 1, j));
			}
		}
		//////////////////////////////////////////////////////////////////////////////////////////////////
		// block A2
		for (int j = 6 * n1 + 1; j < 8 * n1 + 1; j++) {
			for (int i = 9 * n1; i < 11 * n1 + 1; i++) {
				uf(i, j) = uf_calc(u_g(i, j), u_g(i + 1, j), u_g(i - 1, j), u_g(i, j + 1), u_g(i, j - 1),
					v_g(i, j), v_g(i + 1, j), v_g(i, j - 1), v_g(i + 1, j - 1), P_g(i, j), P_g(i + 1, j));
			}
		}

		// block C2
		for (int j = 6 * n1 + 1; j < N + 1; j++) {
			for (int i = 11 * n1 + 1; i < 13 * n1; i++) {
				uf(i, j) = uf_calc(u_g(i, j), u_g(i + 1, j), u_g(i - 1, j), u_g(i, j + 1), u_g(i, j - 1),
					v_g(i, j), v_g(i + 1, j), v_g(i, j - 1), v_g(i + 1, j - 1), P_g(i, j), P_g(i + 1, j));
			}
		}

		// block B2
		for (int j = 12 * n1 + 1; j < N + 1; j++) {
			for (int i = 13 * n1; i < 15 * n1 + 1; i++) {
				uf(i, j) = uf_calc(u_g(i, j), u_g(i + 1, j), u_g(i - 1, j), u_g(i, j + 1), u_g(i, j - 1),
					v_g(i, j), v_g(i + 1, j), v_g(i, j - 1), v_g(i + 1, j - 1), P_g(i, j), P_g(i + 1, j));
			}
		}

// block D2
for (int j = 6 * n1 + 1; j < N + 1; j++) {
	for (int i = 15 * n1 + 1; i < 17 * n1; i++) {
		uf(i, j) = uf_calc(u_g(i, j), u_g(i + 1, j), u_g(i - 1, j), u_g(i, j + 1), u_g(i, j - 1),
			v_g(i, j), v_g(i + 1, j), v_g(i, j - 1), v_g(i + 1, j - 1), P_g(i, j), P_g(i + 1, j));
	}
}
//////////////////////////////////////////////////////////////////////////////////////////////////
// block A3
for (int j = 6 * n1 + 1; j < 8 * n1 + 1; j++) {
	for (int i = 17 * n1; i < 19 * n1 + 1; i++) {
		uf(i, j) = uf_calc(u_g(i, j), u_g(i + 1, j), u_g(i - 1, j), u_g(i, j + 1), u_g(i, j - 1),
			v_g(i, j), v_g(i + 1, j), v_g(i, j - 1), v_g(i + 1, j - 1), P_g(i, j), P_g(i + 1, j));
	}
}

// block C3
for (int j = 6 * n1 + 1; j < N + 1; j++) {
	for (int i = 19 * n1 + 1; i < 21 * n1; i++) {
		uf(i, j) = uf_calc(u_g(i, j), u_g(i + 1, j), u_g(i - 1, j), u_g(i, j + 1), u_g(i, j - 1),
			v_g(i, j), v_g(i + 1, j), v_g(i, j - 1), v_g(i + 1, j - 1), P_g(i, j), P_g(i + 1, j));
	}
}

// block B3
for (int j = 12 * n1 + 1; j < N + 1; j++) {
	for (int i = 21 * n1; i < 23 * n1 + 1; i++) {
		uf(i, j) = uf_calc(u_g(i, j), u_g(i + 1, j), u_g(i - 1, j), u_g(i, j + 1), u_g(i, j - 1),
			v_g(i, j), v_g(i + 1, j), v_g(i, j - 1), v_g(i + 1, j - 1), P_g(i, j), P_g(i + 1, j));
	}
}

// block D3
for (int j = 6 * n1 + 1; j < N + 1; j++) {
	for (int i = 23 * n1 + 1; i < 25 * n1; i++) {
		uf(i, j) = uf_calc(u_g(i, j), u_g(i + 1, j), u_g(i - 1, j), u_g(i, j + 1), u_g(i, j - 1),
			v_g(i, j), v_g(i + 1, j), v_g(i, j - 1), v_g(i + 1, j - 1), P_g(i, j), P_g(i + 1, j));
	}
}
//////////////////////////////////////////////////////////////////////////////////////////////////
// block A4
for (int j = 6 * n1 + 1; j < 8 * n1 + 1; j++) {
	for (int i = 25 * n1; i < 27 * n1 + 1; i++) {
		uf(i, j) = uf_calc(u_g(i, j), u_g(i + 1, j), u_g(i - 1, j), u_g(i, j + 1), u_g(i, j - 1),
			v_g(i, j), v_g(i + 1, j), v_g(i, j - 1), v_g(i + 1, j - 1), P_g(i, j), P_g(i + 1, j));
	}
}

// block C4
for (int j = 6 * n1 + 1; j < N + 1; j++) {
	for (int i = 27 * n1 + 1; i < 29 * n1; i++) {
		uf(i, j) = uf_calc(u_g(i, j), u_g(i + 1, j), u_g(i - 1, j), u_g(i, j + 1), u_g(i, j - 1),
			v_g(i, j), v_g(i + 1, j), v_g(i, j - 1), v_g(i + 1, j - 1), P_g(i, j), P_g(i + 1, j));
	}
}

// block B4
for (int j = 12 * n1 + 1; j < N + 1; j++) {
	for (int i = 29 * n1; i < 31 * n1 + 1; i++) {
		uf(i, j) = uf_calc(u_g(i, j), u_g(i + 1, j), u_g(i - 1, j), u_g(i, j + 1), u_g(i, j - 1),
			v_g(i, j), v_g(i + 1, j), v_g(i, j - 1), v_g(i + 1, j - 1), P_g(i, j), P_g(i + 1, j));
	}
}

// block D4
for (int j = 6 * n1 + 1; j < N + 1; j++) {
	for (int i = 31 * n1 + 1; i < 33 * n1; i++) {
		uf(i, j) = uf_calc(u_g(i, j), u_g(i + 1, j), u_g(i - 1, j), u_g(i, j + 1), u_g(i, j - 1),
			v_g(i, j), v_g(i + 1, j), v_g(i, j - 1), v_g(i + 1, j - 1), P_g(i, j), P_g(i + 1, j));
	}
}
//////////////////////////////////////////////////////////////////////////////////////////////////
// block A5
for (int j = 6 * n1 + 1; j < 8 * n1 + 1; j++) {
	for (int i = 33 * n1; i < 35 * n1 + 1; i++) {
		uf(i, j) = uf_calc(u_g(i, j), u_g(i + 1, j), u_g(i - 1, j), u_g(i, j + 1), u_g(i, j - 1),
			v_g(i, j), v_g(i + 1, j), v_g(i, j - 1), v_g(i + 1, j - 1), P_g(i, j), P_g(i + 1, j));
	}
}
//////////////////////////////////////////////////////////////////////////////////////////////////
/*
// block 2
for (int j = 6 * n1 + 1; j < 8 * n1 + 1; j++) {
	for (int i = n1; i < M + 1 - n1; i++) {
		uf(i, j) = uf_calc(u_g(i, j), u_g(i + 1, j), u_g(i - 1, j), u_g(i, j + 1), u_g(i, j - 1),
			v_g(i, j), v_g(i + 1, j), v_g(i, j - 1), v_g(i + 1, j - 1), P_g(i, j), P_g(i + 1, j));
	}
}
*/

// block 3
for (int j = 1; j < N + 1; j++) {
	for (int i = M - n1 + 1; i < M; i++) {
		uf(i, j) = uf_calc(u_g(i, j), u_g(i + 1, j), u_g(i - 1, j), u_g(i, j + 1), u_g(i, j - 1),
			v_g(i, j), v_g(i + 1, j), v_g(i, j - 1), v_g(i + 1, j - 1), P_g(i, j), P_g(i + 1, j));
	}
}
	}
	// calulating vf
	{

}
	{
		// block 1
		for (int j = 1; j < N; j++) {
			for (int i = 1; i < n1 + 1; i++) {				
				vf(i, j) = vf_calc(v_g(i, j), v_g(i, j + 1), v_g(i, j - 1), v_g(i + 1, j), v_g(i - 1, j)
					, u_g(i, j), u_g(i, j + 1), u_g(i - 1, j), u_g(i - 1, j + 1), P_g(i, j), P_g(i, j + 1));				
			}
		}
		//////////////////////////////////////////////////////////////////////////////////////////////////
		// block A1
		for (int j = 6 * n1 + 1; j < 8 * n1; j++) {
			for (int i = n1 + 1; i < 3 * n1 + 1; i++) {
				vf(i, j) = vf_calc(v_g(i, j), v_g(i, j + 1), v_g(i, j - 1), v_g(i + 1, j), v_g(i - 1, j)
					, u_g(i, j), u_g(i, j + 1), u_g(i - 1, j), u_g(i - 1, j + 1), P_g(i, j), P_g(i, j + 1));
			}
		}

		// block C1
		for (int j = 6 * n1 + 1; j < N; j++) {
			for (int i = 3 * n1 + 1; i < 5 * n1 + 1; i++) {
				vf(i, j) = vf_calc(v_g(i, j), v_g(i, j + 1), v_g(i, j - 1), v_g(i + 1, j), v_g(i - 1, j)
					, u_g(i, j), u_g(i, j + 1), u_g(i - 1, j), u_g(i - 1, j + 1), P_g(i, j), P_g(i, j + 1));
			}
		}

		// block B1
		for (int j = 12 * n1 + 1; j < N; j++) {
			for (int i = 5 * n1 + 1; i < 7 * n1 + 1; i++) {
				vf(i, j) = vf_calc(v_g(i, j), v_g(i, j + 1), v_g(i, j - 1), v_g(i + 1, j), v_g(i - 1, j)
					, u_g(i, j), u_g(i, j + 1), u_g(i - 1, j), u_g(i - 1, j + 1), P_g(i, j), P_g(i, j + 1));
			}
		}

		// block D1
		for (int j = 6 * n1 + 1; j < N; j++) {
			for (int i = 7 * n1 + 1; i < 9 * n1 + 1; i++) {
				vf(i, j) = vf_calc(v_g(i, j), v_g(i, j + 1), v_g(i, j - 1), v_g(i + 1, j), v_g(i - 1, j)
					, u_g(i, j), u_g(i, j + 1), u_g(i - 1, j), u_g(i - 1, j + 1), P_g(i, j), P_g(i, j + 1));
			}
		}
		/////////////////////////////////////////////////////////////////////////////////////////////////
		// block A2
		for (int j = 6 * n1 + 1; j < 8 * n1; j++) {
			for (int i = 9 * n1 + 1; i < 11 * n1 + 1; i++) {
				vf(i, j) = vf_calc(v_g(i, j), v_g(i, j + 1), v_g(i, j - 1), v_g(i + 1, j), v_g(i - 1, j)
					, u_g(i, j), u_g(i, j + 1), u_g(i - 1, j), u_g(i - 1, j + 1), P_g(i, j), P_g(i, j + 1));
			}
		}

		// block C2
		for (int j = 6 * n1 + 1; j < N; j++) {
			for (int i = 11 * n1 + 1; i < 13 * n1 + 1; i++) {
				vf(i, j) = vf_calc(v_g(i, j), v_g(i, j + 1), v_g(i, j - 1), v_g(i + 1, j), v_g(i - 1, j)
					, u_g(i, j), u_g(i, j + 1), u_g(i - 1, j), u_g(i - 1, j + 1), P_g(i, j), P_g(i, j + 1));
			}
		}

		// block B2
		for (int j = 12 * n1 + 1; j < N; j++) {
			for (int i = 13 * n1 + 1; i < 15 * n1 + 1; i++) {
				vf(i, j) = vf_calc(v_g(i, j), v_g(i, j + 1), v_g(i, j - 1), v_g(i + 1, j), v_g(i - 1, j)
					, u_g(i, j), u_g(i, j + 1), u_g(i - 1, j), u_g(i - 1, j + 1), P_g(i, j), P_g(i, j + 1));
			}
		}

		// block D2
		for (int j = 6 * n1 + 1; j < N; j++) {
			for (int i = 15 * n1 + 1; i < 17 * n1 + 1; i++) {
				vf(i, j) = vf_calc(v_g(i, j), v_g(i, j + 1), v_g(i, j - 1), v_g(i + 1, j), v_g(i - 1, j)
						, u_g(i, j), u_g(i, j + 1), u_g(i - 1, j), u_g(i - 1, j + 1), P_g(i, j), P_g(i, j + 1));
			}
		}
		/////////////////////////////////////////////////////////////////////////////////////////////////
		// block A3
		for (int j = 6 * n1 + 1; j < 8 * n1; j++) {
			for (int i = 17 * n1 + 1; i < 19 * n1 + 1; i++) {
				vf(i, j) = vf_calc(v_g(i, j), v_g(i, j + 1), v_g(i, j - 1), v_g(i + 1, j), v_g(i - 1, j)
					, u_g(i, j), u_g(i, j + 1), u_g(i - 1, j), u_g(i - 1, j + 1), P_g(i, j), P_g(i, j + 1));
			}
		}

		// block C3
		for (int j = 6 * n1 + 1; j < N; j++) {
			for (int i = 19 * n1 + 1; i < 21 * n1 + 1; i++) {
				vf(i, j) = vf_calc(v_g(i, j), v_g(i, j + 1), v_g(i, j - 1), v_g(i + 1, j), v_g(i - 1, j)
					, u_g(i, j), u_g(i, j + 1), u_g(i - 1, j), u_g(i - 1, j + 1), P_g(i, j), P_g(i, j + 1));
			}
		}

		// block B3
		for (int j = 12 * n1 + 1; j < N; j++) {
			for (int i = 21 * n1 + 1; i < 23 * n1 + 1; i++) {
				vf(i, j) = vf_calc(v_g(i, j), v_g(i, j + 1), v_g(i, j - 1), v_g(i + 1, j), v_g(i - 1, j)
					, u_g(i, j), u_g(i, j + 1), u_g(i - 1, j), u_g(i - 1, j + 1), P_g(i, j), P_g(i, j + 1));
			}
		}

		// block D3
		for (int j = 6 * n1 + 1; j < N; j++) {
			for (int i = 23 * n1 + 1; i < 25 * n1 + 1; i++) {
				vf(i, j) = vf_calc(v_g(i, j), v_g(i, j + 1), v_g(i, j - 1), v_g(i + 1, j), v_g(i - 1, j)
					, u_g(i, j), u_g(i, j + 1), u_g(i - 1, j), u_g(i - 1, j + 1), P_g(i, j), P_g(i, j + 1));
			}
		}
		/////////////////////////////////////////////////////////////////////////////////////////////////
		// block A4
		for (int j = 6 * n1 + 1; j < 8 * n1; j++) {
			for (int i = 25 * n1 + 1; i < 27 * n1 + 1; i++) {
				vf(i, j) = vf_calc(v_g(i, j), v_g(i, j + 1), v_g(i, j - 1), v_g(i + 1, j), v_g(i - 1, j)
					, u_g(i, j), u_g(i, j + 1), u_g(i - 1, j), u_g(i - 1, j + 1), P_g(i, j), P_g(i, j + 1));
			}
		}

		// block C4
		for (int j = 6 * n1 + 1; j < N; j++) {
			for (int i = 27 * n1 + 1; i < 29 * n1 + 1; i++) {
				vf(i, j) = vf_calc(v_g(i, j), v_g(i, j + 1), v_g(i, j - 1), v_g(i + 1, j), v_g(i - 1, j)
					, u_g(i, j), u_g(i, j + 1), u_g(i - 1, j), u_g(i - 1, j + 1), P_g(i, j), P_g(i, j + 1));
			}
		}

		// block B4
		for (int j = 12 * n1 + 1; j < N; j++) {
			for (int i = 29 * n1 + 1; i < 31 * n1 + 1; i++) {
				vf(i, j) = vf_calc(v_g(i, j), v_g(i, j + 1), v_g(i, j - 1), v_g(i + 1, j), v_g(i - 1, j)
					, u_g(i, j), u_g(i, j + 1), u_g(i - 1, j), u_g(i - 1, j + 1), P_g(i, j), P_g(i, j + 1));
			}
		}

		// block D4
		for (int j = 6 * n1 + 1; j < N; j++) {
			for (int i = 31 * n1 + 1; i < 33 * n1 + 1; i++) {
				vf(i, j) = vf_calc(v_g(i, j), v_g(i, j + 1), v_g(i, j - 1), v_g(i + 1, j), v_g(i - 1, j)
					, u_g(i, j), u_g(i, j + 1), u_g(i - 1, j), u_g(i - 1, j + 1), P_g(i, j), P_g(i, j + 1));
			}
		}
		/////////////////////////////////////////////////////////////////////////////////////////////////
		// block A5
		for (int j = 6 * n1 + 1; j < 8 * n1; j++) {
			for (int i = 33 * n1 + 1; i < 35 * n1 + 1; i++) {
				vf(i, j) = vf_calc(v_g(i, j), v_g(i, j + 1), v_g(i, j - 1), v_g(i + 1, j), v_g(i - 1, j)
					, u_g(i, j), u_g(i, j + 1), u_g(i - 1, j), u_g(i - 1, j + 1), P_g(i, j), P_g(i, j + 1));
			}
		}
		/*
		// block 2 (doubtful)
		for (int j = 6 * n1 + 1; j < 8 * n1; j++) {
			for (int i = n1 + 1; i < M + 1 - n1; i++) {
				vf(i, j) = vf_calc(v_g(i, j), v_g(i, j + 1), v_g(i, j - 1), v_g(i + 1, j), v_g(i - 1, j)
					, u_g(i, j), u_g(i, j + 1), u_g(i - 1, j), u_g(i - 1, j + 1), P_g(i, j), P_g(i, j + 1));
			}
		}
		*/
		// block 3
		for (int j = 1; j < N; j++) {
			for (int i = M + 1 - n1; i < M + 1; i++) {
				vf(i, j) = vf_calc(v_g(i, j), v_g(i, j + 1), v_g(i, j - 1), v_g(i + 1, j), v_g(i - 1, j)
					, u_g(i, j), u_g(i, j + 1), u_g(i - 1, j), u_g(i - 1, j + 1), P_g(i, j), P_g(i, j + 1));
			}
		}
		// uf and vf at boundaries
		// vf only at outlet (rest not required bcoz does not appear in P calculation)
		for (int i = M + 1 - n1; i < M + 1; i++) {
			vf(i, 0) = vf(i, 1);
			vf(i, N) = vf(i, N - 1);
		}
	}

	// P CALCULATION
	{

}
	{
	double tolerance = 2.0;
	int k = 0;
	while (tolerance > 1e-5) {	
		// POINT JACOBI		
		// block 1		
		for (int j = 1; j < N + 1; j++) {
			for (int i = 1; i < n1 + 1; i++) {
				P_up(i, j) = (1.0 / (2.0 * (1 + beta))) * (P_t(i - 1, j) + P_t(i + 1, j) + beta*(P_t(i, j - 1) + P_t(i, j + 1))
						- (pow(dx, 2) * bij(uf(i, j), uf(i - 1, j), vf(i, j), vf(i, j - 1))));
			}
		}
		////////////////////////////////////////////////////////////////////////////////////////
		// block A1
		for (int j = 6 * n1 + 1; j < 8 * n1 + 1; j++) {
			for (int i = n1 + 1; i < 3 * n1 + 1; i++) {
				P_up(i, j) = (1.0 / (2.0 * (1 + beta))) * (P_t(i - 1, j) + P_t(i + 1, j) + beta*(P_t(i, j - 1) + P_t(i, j + 1))
					- (pow(dx, 2) * bij(uf(i, j), uf(i - 1, j), vf(i, j), vf(i, j - 1))));
			}
		}
		// block C1
		for (int j = 6 * n1 + 1; j < N + 1; j++) {
			for (int i = 3 * n1 + 1; i < 5 * n1 + 1; i++) {
				P_up(i, j) = (1.0 / (2.0 * (1 + beta))) * (P_t(i - 1, j) + P_t(i + 1, j) + beta*(P_t(i, j - 1) + P_t(i, j + 1))
					- (pow(dx, 2) * bij(uf(i, j), uf(i - 1, j), vf(i, j), vf(i, j - 1))));
			}
		}
		// block B1
		for (int j = j = 12 * n1 + 1; j < N + 1; j++) {
			for (int i = 5 * n1 + 1; i < 7 * n1 + 1; i++) {
				P_up(i, j) = (1.0 / (2.0 * (1 + beta))) * (P_t(i - 1, j) + P_t(i + 1, j) + beta*(P_t(i, j - 1) + P_t(i, j + 1))
					- (pow(dx, 2) * bij(uf(i, j), uf(i - 1, j), vf(i, j), vf(i, j - 1))));
			}
		}
		// block D1
		for (int j = 6 * n1 + 1; j < N + 1; j++) {
			for (int i = 7 * n1 + 1; i < 9 * n1 + 1; i++) {
				P_up(i, j) = (1.0 / (2.0 * (1 + beta))) * (P_t(i - 1, j) + P_t(i + 1, j) + beta*(P_t(i, j - 1) + P_t(i, j + 1))
					- (pow(dx, 2) * bij(uf(i, j), uf(i - 1, j), vf(i, j), vf(i, j - 1))));
			}
		}
		////////////////////////////////////////////////////////////////////////////////////////
		// block A2
		for (int j = 6 * n1 + 1; j < 8 * n1 + 1; j++) {
			for (int i = 9 * n1 + 1; i < 11 * n1 + 1; i++) {
				P_up(i, j) = (1.0 / (2.0 * (1 + beta))) * (P_t(i - 1, j) + P_t(i + 1, j) + beta*(P_t(i, j - 1) + P_t(i, j + 1))
					- (pow(dx, 2) * bij(uf(i, j), uf(i - 1, j), vf(i, j), vf(i, j - 1))));
			}
		}
		// block C2
		for (int j = 6 * n1 + 1; j < N + 1; j++) {
			for (int i = 11 * n1 + 1; i < 13 * n1 + 1; i++) {
				P_up(i, j) = (1.0 / (2.0 * (1 + beta))) * (P_t(i - 1, j) + P_t(i + 1, j) + beta*(P_t(i, j - 1) + P_t(i, j + 1))
					- (pow(dx, 2) * bij(uf(i, j), uf(i - 1, j), vf(i, j), vf(i, j - 1))));
			}
		}
		// block B2
		for (int j = j = 12 * n1 + 1; j < N + 1; j++) {
			for (int i = 13 * n1 + 1; i < 15 * n1 + 1; i++) {
				P_up(i, j) = (1.0 / (2.0 * (1 + beta))) * (P_t(i - 1, j) + P_t(i + 1, j) + beta*(P_t(i, j - 1) + P_t(i, j + 1))
					- (pow(dx, 2) * bij(uf(i, j), uf(i - 1, j), vf(i, j), vf(i, j - 1))));
			}
		}
		// block D2
		for (int j = 6 * n1 + 1; j < N + 1; j++) {
			for (int i = 15 * n1 + 1; i < 17 * n1 + 1; i++) {
				P_up(i, j) = (1.0 / (2.0 * (1 + beta))) * (P_t(i - 1, j) + P_t(i + 1, j) + beta*(P_t(i, j - 1) + P_t(i, j + 1))
					- (pow(dx, 2) * bij(uf(i, j), uf(i - 1, j), vf(i, j), vf(i, j - 1))));
			}
		}
		////////////////////////////////////////////////////////////////////////////////////////
		// block A3
		for (int j = 6 * n1 + 1; j < 8 * n1 + 1; j++) {
			for (int i = 17 * n1 + 1; i < 19 * n1 + 1; i++) {
				P_up(i, j) = (1.0 / (2.0 * (1 + beta))) * (P_t(i - 1, j) + P_t(i + 1, j) + beta*(P_t(i, j - 1) + P_t(i, j + 1))
					- (pow(dx, 2) * bij(uf(i, j), uf(i - 1, j), vf(i, j), vf(i, j - 1))));
			}
		}
		// block C3
		for (int j = 6 * n1 + 1; j < N + 1; j++) {
			for (int i = 19 * n1 + 1; i < 21 * n1 + 1; i++) {
				P_up(i, j) = (1.0 / (2.0 * (1 + beta))) * (P_t(i - 1, j) + P_t(i + 1, j) + beta*(P_t(i, j - 1) + P_t(i, j + 1))
					- (pow(dx, 2) * bij(uf(i, j), uf(i - 1, j), vf(i, j), vf(i, j - 1))));
			}
		}
		// block B3
		for (int j = j = 12 * n1 + 1; j < N + 1; j++) {
			for (int i = 21 * n1 + 1; i < 23 * n1 + 1; i++) {
				P_up(i, j) = (1.0 / (2.0 * (1 + beta))) * (P_t(i - 1, j) + P_t(i + 1, j) + beta*(P_t(i, j - 1) + P_t(i, j + 1))
					- (pow(dx, 2) * bij(uf(i, j), uf(i - 1, j), vf(i, j), vf(i, j - 1))));
			}
		}
		// block D3
		for (int j = 6 * n1 + 1; j < N + 1; j++) {
			for (int i = 23 * n1 + 1; i < 25 * n1 + 1; i++) {
				P_up(i, j) = (1.0 / (2.0 * (1 + beta))) * (P_t(i - 1, j) + P_t(i + 1, j) + beta*(P_t(i, j - 1) + P_t(i, j + 1))
					- (pow(dx, 2) * bij(uf(i, j), uf(i - 1, j), vf(i, j), vf(i, j - 1))));
			}
		}
		////////////////////////////////////////////////////////////////////////////////////////
		// block A4
		for (int j = 6 * n1 + 1; j < 8 * n1 + 1; j++) {
			for (int i = 25 * n1 + 1; i < 27 * n1 + 1; i++) {
				P_up(i, j) = (1.0 / (2.0 * (1 + beta))) * (P_t(i - 1, j) + P_t(i + 1, j) + beta*(P_t(i, j - 1) + P_t(i, j + 1))
					- (pow(dx, 2) * bij(uf(i, j), uf(i - 1, j), vf(i, j), vf(i, j - 1))));
			}
		}
		// block C4
		for (int j = 6 * n1 + 1; j < N + 1; j++) {
			for (int i = 27 * n1 + 1; i < 29 * n1 + 1; i++) {
				P_up(i, j) = (1.0 / (2.0 * (1 + beta))) * (P_t(i - 1, j) + P_t(i + 1, j) + beta*(P_t(i, j - 1) + P_t(i, j + 1))
					- (pow(dx, 2) * bij(uf(i, j), uf(i - 1, j), vf(i, j), vf(i, j - 1))));
			}
		}
		// block B4
		for (int j = j = 12 * n1 + 1; j < N + 1; j++) {
			for (int i = 29 * n1 + 1; i < 31 * n1 + 1; i++) {
				P_up(i, j) = (1.0 / (2.0 * (1 + beta))) * (P_t(i - 1, j) + P_t(i + 1, j) + beta*(P_t(i, j - 1) + P_t(i, j + 1))
					- (pow(dx, 2) * bij(uf(i, j), uf(i - 1, j), vf(i, j), vf(i, j - 1))));
			}
		}
		// block D4
		for (int j = 6 * n1 + 1; j < N + 1; j++) {
			for (int i = 31 * n1 + 1; i < 33 * n1 + 1; i++) {
				P_up(i, j) = (1.0 / (2.0 * (1 + beta))) * (P_t(i - 1, j) + P_t(i + 1, j) + beta*(P_t(i, j - 1) + P_t(i, j + 1))
					- (pow(dx, 2) * bij(uf(i, j), uf(i - 1, j), vf(i, j), vf(i, j - 1))));
			}
		}
		// block A5
		for (int j = 6 * n1 + 1; j < 8 * n1 + 1; j++) {
			for (int i = 33 * n1 + 1; i < 35 * n1 + 1; i++) {
				P_up(i, j) = (1.0 / (2.0 * (1 + beta))) * (P_t(i - 1, j) + P_t(i + 1, j) + beta*(P_t(i, j - 1) + P_t(i, j + 1))
					- (pow(dx, 2) * bij(uf(i, j), uf(i - 1, j), vf(i, j), vf(i, j - 1))));
			}
		}
		////////////////////////////////////////////////////////////////////////////////////////////
		/*
		// block 2
		for (int j = 6 * n1 + 1; j < 8 * n1 + 1; j++) {
			for (int i = n1 + 1; i < M + 1 - n1; i++) {
				P_up(i, j) = (1.0 / (2.0 * (1 + beta))) * (P_t(i - 1, j) + P_t(i + 1, j) + beta*(P_t(i, j - 1) + P_t(i, j + 1))
					- (pow(dx, 2) * bij(uf(i, j), uf(i - 1, j), vf(i, j), vf(i, j - 1))));
			}
		}
		*/
// block 3
for (int j = 1; j < N + 1; j++) {
	for (int i = M + 1 - n1; i < M + 1; i++) {
		P_up(i, j) = (1.0 / (2.0 * (1 + beta))) * (P_t(i - 1, j) + P_t(i + 1, j) + beta * (P_t(i, j - 1) + P_t(i, j + 1))
			- (pow(dx, 2) * bij(uf(i, j), uf(i - 1, j), vf(i, j), vf(i, j - 1))));
	}
}
// Neumann bcs for P_up
// For Left and Right Boundary
for (int j = 1; j < N + 1; j++) {
	P_up(0, j) = P_up(1, j);
	P_up(M + 1, j) = P_up(M, j);
}
// For Top and Bottom except the outlet
for (int i = 1; i < M - n1 + 1; i++) {
	P_up(i, 0) = P_up(i, 1);
	P_up(i, N + 1) = P_up(i, N);
}
// For Block 1(Bottom Right Walls) and Block 2(Bottom left walls)
for (int j = 1; j < 6 * n1 + 1; j++) {
	P_up(n1 + 1, j) = P_up(n1, j);
	P_up(M - n1, j) = P_up(M + 1 - n1, j);
}
// For Block 1(Top Right Walls) and Block 2(Top Left walls)
for (int j = 8 * n1 + 1; j < N + 1; j++) {
	P_up(n1 + 1, j) = P_up(n1, j);
	P_up(M - n1, j) = P_up(M + 1 - n1, j);
}
// For all the walls at j=60(Ghost Nodes) in Serpentine Geometry
for (int i = n1 + 1; i < M + 1 - n1; i++) {
	P_up(i, 6 * n1) = P_up(i, 6 * n1 + 1);
	//P_up(i, 8 * n1 + 1) = P_up(i, 8 * n1); //Commented because of Serpentine Geometry
}
/////////////////////////////////////////////////////////////////////////////////////
// For all the walls at j=81(Ghost Nodes) in Serpentine Geometry
for (int i = n1 + 1; i < 3 * n1 + 1; i++) {
	P_up(i, 8 * n1 + 1) = P_up(i, 8 * n1);
} // A1 top boundary
for (int i = 9 * n1 + 1; i < 11 * n1 + 1; i++) {
	P_up(i, 8 * n1 + 1) = P_up(i, 8 * n1);
} // A2 top boundary
for (int i = 17 * n1 + 1; i < 19 * n1 + 1; i++) {
	P_up(i, 8 * n1 + 1) = P_up(i, 8 * n1);
} // A3 top boundary
for (int i = 25 * n1 + 1; i < 27 * n1 + 1; i++) {
	P_up(i, 8 * n1 + 1) = P_up(i, 8 * n1);
} // A4 top boundary
for (int i = 33 * n1 + 1; i < 35 * n1 + 1; i++) {
	P_up(i, 8 * n1 + 1) = P_up(i, 8 * n1);
} // A5 top boundary
/////////////////////////////////////////////////////////////////////////////////////
// For all the walls lying between j = 8*n1+1 to N
for (int j = 8 * n1 + 1; j < N + 1; j++) {
	P_up(3 * n1, j) = P_up(3 * n1 + 1, j); //Block C1 left Boundary
	P_up(9 * n1 + 1, j) = P_up(9 * n1, j); //Block D1 right Boundary
	P_up(11 * n1, j) = P_up(11 * n1 + 1, j); //Block C2 left Boundary
	P_up(17 * n1 + 1, j) = P_up(17 * n1, j); //Block D2 right Boundary
	P_up(19 * n1, j) = P_up(19 * n1 + 1, j); //Block C3 left Boundary
	P_up(25 * n1 + 1, j) = P_up(25 * n1, j); //Block D3 right Boundary
	P_up(27 * n1, j) = P_up(27 * n1 + 1, j); //Block C4 left Boundary
	P_up(33 * n1 + 1, j) = P_up(33 * n1, j); //Block D4 right Boundary
}
/////////////////////////////////////////////////////////////////////////////////////
// For all the walls lying between j = 6*n1+1 to 12*n1
for (int j = 6 * n1 + 1; j < 12 * n1 + 1; j++) {
	P_up(5 * n1 + 1, j) = P_up(5 * n1, j); //Block C1 right Boundary
	P_up(7 * n1, j) = P_up(7 * n1 + 1, j); //Block D1 left Boundary
	P_up(13 * n1 + 1, j) = P_up(13 * n1, j); //Block C2 right Boundary
	P_up(15 * n1, j) = P_up(15 * n1 + 1, j); //Block D2 left Boundary
	P_up(21 * n1 + 1, j) = P_up(21 * n1, j); //Block C3 right Boundary
	P_up(23 * n1, j) = P_up(23 * n1 + 1, j); //Block D3 left Boundary
	P_up(29 * n1 + 1, j) = P_up(29 * n1, j); //Block C4 right Boundary
	P_up(31 * n1, j) = P_up(31 * n1 + 1, j); //Block D4 left Boundary
}
/////////////////////////////////////////////////////////////////////////////////////
// For all the walls lying at j = 12*n1
for (int i = 5 * n1 + 1; i < 7 * n1 + 1; i++) {
	P_up(i, 12 * n1) = P_up(i, 12 * n1 + 1);
} // For B1 bottom boundary
for (int i = 13 * n1 + 1; i < 15 * n1 + 1; i++) {
	P_up(i, 12 * n1) = P_up(i, 12 * n1 + 1);
} // For B2 bottom boundary
for (int i = 21 * n1 + 1; i < 23 * n1 + 1; i++) {
	P_up(i, 12 * n1) = P_up(i, 12 * n1 + 1);
} // For B3 bottom boundary
for (int i = 29 * n1 + 1; i < 31 * n1 + 1; i++) {
	P_up(i, 12 * n1) = P_up(i, 12 * n1 + 1);
} // For B4 bottom boundary

///////////////////////// End of Pressure Boundary Conditions implementation
tolerance = (P_up - P_t).cwiseAbs().maxCoeff();

P_t = P_up;
k++;

cout << "iteration count: " << k << "\n";
cout << "tolerance: " << tolerance << "\n";
	}
}

// Vn+1 calculation
	{

}
	{
		// un+1 calc
		// block 1
		for (int j = 1; j < N + 1; j++) {
			for (int i = 1; i < n1; i++) {
				u_g(i, j) = uf(i, j) - dt*((P_up(i + 1, j) - P_up(i, j)) / dx);
			}
		}

		// block A1
		for (int j = 6 * n1 + 1; j < 8 * n1 + 1; j++) {
			for (int i = n1; i < 3 * n1 + 1; i++) {
				u_g(i, j) = uf(i, j) - dt*((P_up(i + 1, j) - P_up(i, j)) / dx);
			}
		}

		// block C1
		for (int j = 6 * n1 + 1; j < N + 1; j++) {
			for (int i = 3 * n1 + 1; i < 5 * n1; i++) {
				u_g(i, j) = uf(i, j) - dt*((P_up(i + 1, j) - P_up(i, j)) / dx);
			}
		}

		// block B1
		for (int j = 12 * n1 + 1; j < N + 1; j++) {
			for (int i = 5 * n1; i < 7 * n1 + 1; i++) {
				u_g(i, j) = uf(i, j) - dt*((P_up(i + 1, j) - P_up(i, j)) / dx);
			}
		}

		// block D1
		for (int j = 6 * n1 + 1; j < N + 1; j++) {
			for (int i = 7 * n1 + 1; i < 9 * n1; i++) {
				u_g(i, j) = uf(i, j) - dt*((P_up(i + 1, j) - P_up(i, j)) / dx);
			}
		}
		//////////////////////////////////////////////////////////////////////////////////////////////////
		// block A2
		for (int j = 6 * n1 + 1; j < 8 * n1 + 1; j++) {
			for (int i = 9 * n1; i < 11 * n1 + 1; i++) {
				u_g(i, j) = uf(i, j) - dt*((P_up(i + 1, j) - P_up(i, j)) / dx);
			}
		}

		// block C2
		for (int j = 6 * n1 + 1; j < N + 1; j++) {
			for (int i = 11 * n1 + 1; i < 13 * n1; i++) {
				u_g(i, j) = uf(i, j) - dt*((P_up(i + 1, j) - P_up(i, j)) / dx);
			}
		}

		// block B2
		for (int j = 12 * n1 + 1; j < N + 1; j++) {
			for (int i = 13 * n1; i < 15 * n1 + 1; i++) {
				u_g(i, j) = uf(i, j) - dt*((P_up(i + 1, j) - P_up(i, j)) / dx);
			}
		}

		// block D2
		for (int j = 6 * n1 + 1; j < N + 1; j++) {
			for (int i = 15 * n1 + 1; i < 17 * n1; i++) {
				u_g(i, j) = uf(i, j) - dt*((P_up(i + 1, j) - P_up(i, j)) / dx);
			}
		}
		//////////////////////////////////////////////////////////////////////////////////////////////////
		// block A3
		for (int j = 6 * n1 + 1; j < 8 * n1 + 1; j++) {
			for (int i = 17 * n1; i < 19 * n1 + 1; i++) {
				u_g(i, j) = uf(i, j) - dt*((P_up(i + 1, j) - P_up(i, j)) / dx);
			}
		}

		// block C3
		for (int j = 6 * n1 + 1; j < N + 1; j++) {
			for (int i = 19 * n1 + 1; i < 21 * n1; i++) {
				u_g(i, j) = uf(i, j) - dt*((P_up(i + 1, j) - P_up(i, j)) / dx);
			}
		}

		// block B3
		for (int j = 12 * n1 + 1; j < N + 1; j++) {
			for (int i = 21 * n1; i < 23 * n1 + 1; i++) {
				u_g(i, j) = uf(i, j) - dt*((P_up(i + 1, j) - P_up(i, j)) / dx);
			}
		}

		// block D3
		for (int j = 6 * n1 + 1; j < N + 1; j++) {
			for (int i = 23 * n1 + 1; i < 25 * n1; i++) {
				u_g(i, j) = uf(i, j) - dt*((P_up(i + 1, j) - P_up(i, j)) / dx);
			}
		}
		//////////////////////////////////////////////////////////////////////////////////////////////////
		// block A4
		for (int j = 6 * n1 + 1; j < 8 * n1 + 1; j++) {
			for (int i = 25 * n1; i < 27 * n1 + 1; i++) {
				u_g(i, j) = uf(i, j) - dt*((P_up(i + 1, j) - P_up(i, j)) / dx);
			}
		}

		// block C4
		for (int j = 6 * n1 + 1; j < N + 1; j++) {
			for (int i = 27 * n1 + 1; i < 29 * n1; i++) {
				u_g(i, j) = uf(i, j) - dt*((P_up(i + 1, j) - P_up(i, j)) / dx);
			}
		}

		// block B4
		for (int j = 12 * n1 + 1; j < N + 1; j++) {
			for (int i = 29 * n1; i < 31 * n1 + 1; i++) {
				u_g(i, j) = uf(i, j) - dt*((P_up(i + 1, j) - P_up(i, j)) / dx);
			}
		}

		// block D4
		for (int j = 6 * n1 + 1; j < N + 1; j++) {
			for (int i = 31 * n1 + 1; i < 33 * n1; i++) {
				u_g(i, j) = uf(i, j) - dt*((P_up(i + 1, j) - P_up(i, j)) / dx);
			}
		}
		//////////////////////////////////////////////////////////////////////////////////////////////////
		// block A5
		for (int j = 6 * n1 + 1; j < 8 * n1 + 1; j++) {
			for (int i = 33 * n1; i < 35 * n1 + 1; i++) {
				u_g(i, j) = uf(i, j) - dt*((P_up(i + 1, j) - P_up(i, j)) / dx);
			}
		}
		//////////////////////////////////////////////////////////////////////////////////////////////////
		/*
		// block 2
		for (int j = 6 * n1 + 1; j < 8 * n1 + 1; j++) {
		for (int i = n1; i < M + 1 - n1; i++) {
		u_g(i, j) = uf(i, j) - dt*((P_up(i + 1, j) - P_up(i, j)) / dx);
		}
		}
		*/

		// block 3
		for (int j = 1; j < N + 1; j++) {
			for (int i = M - n1 + 1; i < M; i++) {
				u_g(i, j) = uf(i, j) - dt*((P_up(i + 1, j) - P_up(i, j)) / dx);
			}
		}
		
		//////////////////////////////////////////////////////////////////////////////////////////////////
		// vn+1 calc
		// Block 1
		for (int j = 1; j < N; j++) {
			for (int i = 1; i < n1 + 1; i++) {
				v_g(i, j) = vf(i, j) - dt*((P_up(i, j + 1) - P_up(i, j)) / dy);
			}
		}
		//////////////////////////////////////////////////////////////////////////////////////////////////
		// block A1
		for (int j = 6 * n1 + 1; j < 8 * n1; j++) {
			for (int i = n1 + 1; i < 3 * n1 + 1; i++) {
				v_g(i, j) = vf(i, j) - dt*((P_up(i, j + 1) - P_up(i, j)) / dy);
			}
		}

		// block C1
		for (int j = 6 * n1 + 1; j < N; j++) {
			for (int i = 3 * n1 + 1; i < 5 * n1 + 1; i++) {
				v_g(i, j) = vf(i, j) - dt*((P_up(i, j + 1) - P_up(i, j)) / dy);
			}
		}

		// block B1
		for (int j = 12 * n1 + 1; j < N; j++) {
			for (int i = 5 * n1 + 1; i < 7 * n1 + 1; i++) {
				v_g(i, j) = vf(i, j) - dt*((P_up(i, j + 1) - P_up(i, j)) / dy);
			}
		}

		// block D1
		for (int j = 6 * n1 + 1; j < N; j++) {
			for (int i = 7 * n1 + 1; i < 9 * n1 + 1; i++) {
				v_g(i, j) = vf(i, j) - dt*((P_up(i, j + 1) - P_up(i, j)) / dy);
			}
		}
		/////////////////////////////////////////////////////////////////////////////////////////////////
		// block A2
		for (int j = 6 * n1 + 1; j < 8 * n1; j++) {
			for (int i = 9 * n1 + 1; i < 11 * n1 + 1; i++) {
				v_g(i, j) = vf(i, j) - dt*((P_up(i, j + 1) - P_up(i, j)) / dy);
			}
		}

		// block C2
		for (int j = 6 * n1 + 1; j < N; j++) {
			for (int i = 11 * n1 + 1; i < 13 * n1 + 1; i++) {
				v_g(i, j) = vf(i, j) - dt*((P_up(i, j + 1) - P_up(i, j)) / dy);
			}
		}

		// block B2
		for (int j = 12 * n1 + 1; j < N; j++) {
			for (int i = 13 * n1 + 1; i < 15 * n1 + 1; i++) {
				v_g(i, j) = vf(i, j) - dt*((P_up(i, j + 1) - P_up(i, j)) / dy);
			}
		}

		// block D2
		for (int j = 6 * n1 + 1; j < N; j++) {
			for (int i = 15 * n1 + 1; i < 17 * n1 + 1; i++) {
				v_g(i, j) = vf(i, j) - dt*((P_up(i, j + 1) - P_up(i, j)) / dy);
			}
		}
		/////////////////////////////////////////////////////////////////////////////////////////////////
		// block A3
		for (int j = 6 * n1 + 1; j < 8 * n1; j++) {
			for (int i = 17 * n1 + 1; i < 19 * n1 + 1; i++) {
				v_g(i, j) = vf(i, j) - dt*((P_up(i, j + 1) - P_up(i, j)) / dy);
			}
		}

		// block C3
		for (int j = 6 * n1 + 1; j < N; j++) {
			for (int i = 19 * n1 + 1; i < 21 * n1 + 1; i++) {
				v_g(i, j) = vf(i, j) - dt*((P_up(i, j + 1) - P_up(i, j)) / dy);
			}
		}

		// block B3
		for (int j = 12 * n1 + 1; j < N; j++) {
			for (int i = 21 * n1 + 1; i < 23 * n1 + 1; i++) {
				v_g(i, j) = vf(i, j) - dt*((P_up(i, j + 1) - P_up(i, j)) / dy);
			}
		}

		// block D3
		for (int j = 6 * n1 + 1; j < N; j++) {
			for (int i = 23 * n1 + 1; i < 25 * n1 + 1; i++) {
				v_g(i, j) = vf(i, j) - dt*((P_up(i, j + 1) - P_up(i, j)) / dy);
			}
		}
		/////////////////////////////////////////////////////////////////////////////////////////////////
		// block A4
		for (int j = 6 * n1 + 1; j < 8 * n1; j++) {
			for (int i = 25 * n1 + 1; i < 27 * n1 + 1; i++) {
				v_g(i, j) = vf(i, j) - dt*((P_up(i, j + 1) - P_up(i, j)) / dy);
			}
		}

		// block C4
		for (int j = 6 * n1 + 1; j < N; j++) {
			for (int i = 27 * n1 + 1; i < 29 * n1 + 1; i++) {
				v_g(i, j) = vf(i, j) - dt*((P_up(i, j + 1) - P_up(i, j)) / dy);
			}
		}

		// block B4
		for (int j = 12 * n1 + 1; j < N; j++) {
			for (int i = 29 * n1 + 1; i < 31 * n1 + 1; i++) {
				v_g(i, j) = vf(i, j) - dt*((P_up(i, j + 1) - P_up(i, j)) / dy);
			}
		}

		// block D4
		for (int j = 6 * n1 + 1; j < N; j++) {
			for (int i = 31 * n1 + 1; i < 33 * n1 + 1; i++) {
				v_g(i, j) = vf(i, j) - dt*((P_up(i, j + 1) - P_up(i, j)) / dy);
			}
		}
		/////////////////////////////////////////////////////////////////////////////////////////////////
		// block A5
		for (int j = 6 * n1 + 1; j < 8 * n1; j++) {
			for (int i = 33 * n1 + 1; i < 35 * n1 + 1; i++) {
				v_g(i, j) = vf(i, j) - dt*((P_up(i, j + 1) - P_up(i, j)) / dy);
			}
		}
		/*
		for (int j = 6 * n1 + 1; j < 8 * n1; j++) {
			for (int i = n1 + 1; i < M + 1 - n1; i++) {
				v_g(i, j) = vf(i, j) - dt*((P_up(i, j + 1) - P_up(i, j)) / dy);
			}
		}
		*/
		// Block 3
		for (int j = 1; j < N; j++) {
			for (int i = M + 1 - n1; i < M + 1; i++) {
				v_g(i, j) = vf(i, j) - dt*((P_up(i, j + 1) - P_up(i, j)) / dy);
			}
		}		
		//////////////////////////////////// Calculation of Vn+1 and Un+1 finished here.

		//////////////////////////////////// Boundary Condition for V and U starts here.
		{
		// neumann at outlet for vn+1
		for (int i = M + 1 - n1; i < M + 1; i++) {
			v_g(i, 0) = v_g(i, 1);
			v_g(i, N) = v_g(i, N - 1);
		}
		// un+1 at ghost
		// block 1
		for (int i = 1; i < n1; i++) {
			u_g(i, 0) = -u_g(i, 1);
			u_g(i, N + 1) = -u_g(i, N);
		}
		/////////////////////////////////////////////////////////////
		for (int i = n1 + 1; i < 5 * n1; i++) {
			u_g(i, 6 * n1) = -u_g(i, 6 * n1 + 1);
		} // For A1 bottom boundary
		for (int i = 7 * n1 + 1; i < 13 * n1; i++) {
			u_g(i, 6 * n1) = -u_g(i, 6 * n1 + 1);
		} // For A2 bottom boundary
		for (int i = 15 * n1 + 1; i < 21 * n1; i++) {
			u_g(i, 6 * n1) = -u_g(i, 6 * n1 + 1);
		} // For A3 bottom boundary
		for (int i = 23 * n1 + 1; i < 29 * n1; i++) {
			u_g(i, 6 * n1) = -u_g(i, 6 * n1 + 1);
		} // For A4 bottom boundary
		for (int i = 31 * n1 + 1; i < 35 * n1; i++) {
			u_g(i, 6 * n1) = -u_g(i, 6 * n1 + 1);
		} // For A5 bottom boundary
		/////////////////////////////////////////////////////////////
		for (int i = n1 + 1; i < 3 * n1; i++) {
			u_g(i, 8 * n1 + 1) = -u_g(i, 8 * n1);
		} // For A1 top boundary
		for (int i = 9 * n1 + 1; i < 11 * n1; i++) {
			u_g(i, 8 * n1 + 1) = -u_g(i, 8 * n1);
		} // For A2 top boundary
		for (int i = 17 * n1 + 1; i < 19 * n1; i++) {
			u_g(i, 8 * n1 + 1) = -u_g(i, 8 * n1);
		} // For A3 top boundary
		for (int i = 25 * n1 + 1; i < 27 * n1; i++) {
			u_g(i, 8 * n1 + 1) = -u_g(i, 8 * n1);
		} // For A4 top boundary
		for (int i = 33 * n1 + 1; i < 35 * n1; i++) {
			u_g(i, 8 * n1 + 1) = -u_g(i, 8 * n1);
		} // For A5 top boundary
		/////////////////////////////////////////////////////////////
		for (int i = 3 * n1 + 1; i < 33 * n1; i++) {
			u_g(i, N + 1) = -u_g(i, N);
		} // For all the top boundaries of B1, B2, B3, B4.
		/////////////////////////////////////////////////////////////
		for (int i = 5 * n1 + 1; i < 7 * n1; i++) {
			u_g(i, 12 * n1) = -u_g(i, 12 * n1 + 1);
		} // For B1 bottom boundary
		for (int i = 13 * n1 + 1; i < 15 * n1; i++) {
			u_g(i, 12 * n1) = -u_g(i, 12 * n1 + 1);
		} // For B2 bottom boundary
		for (int i = 21 * n1 + 1; i < 23 * n1; i++) {
			u_g(i, 12 * n1) = -u_g(i, 12 * n1 + 1);
		} // For B3 bottom boundary
		for (int i = 29 * n1 + 1; i < 31 * n1; i++) {
			u_g(i, 12 * n1) = -u_g(i, 12 * n1 + 1);
		} // For B4 bottom boundary
		/////////////////////////////////////////////////////////////
		/*
		// block 2
		for (int i = n1 + 1; i < M - n1; i++) {
			u_g(i, 6 * n1) = -u_g(i, 6 * n1 + 1);
			u_g(i, 8 * n1 + 1) = -u_g(i, 8 * n1);
		}
		*/
		// block 3
		for (int i = M - n1 + 1; i < M; i++) {
			u_g(i, 0) = -u_g(i, 1);
			u_g(i, N + 1) = -u_g(i, N);
		}
		/////////////////////////////////////////////////////////////
		// vn+1 at ghost
		for (int j = 1; j < N; j++) {
			v_g(0, j) = -v_g(1, j);
			v_g(M + 1, j) = -v_g(M, j);
		} // For left and right boundary

		for (int j = 1; j < 6 * n1; j++) {
			v_g(n1 + 1, j) = -v_g(n1, j);
			v_g(M - n1, j) = -v_g(M - n1 + 1, j);
		} // For Block 1 right boundary and Block 3 left Boundary
		/////////////////////////////////////////////////////////////
		for (int j = 8 * n1 + 1; j < N; j++) {
			v_g(n1 + 1, j) = -v_g(n1, j); // Block 1 Top right boundary
			v_g(3 * n1, j) = -v_g(3 * n1 + 1, j); // C1 left Boundary
			v_g(9 * n1 + 1, j) = -v_g(9 * n1, j); // D1 right boundary
			v_g(11 * n1, j) = -v_g(11 * n1 + 1, j); // C2 left Boundary
			v_g(17 * n1 + 1, j) = -v_g(17 * n1, j); // D2 right boundary
			v_g(19 * n1, j) = -v_g(19 * n1 + 1, j); // C3 left Boundary
			v_g(25 * n1 + 1, j) = -v_g(25 * n1, j); // D3 right boundary
			v_g(27 * n1, j) = -v_g(27 * n1 + 1, j); // C4 left Boundary
			v_g(33 * n1 + 1, j) = -v_g(33 * n1, j); // D4 right boundary
			v_g(M - n1, j) = -v_g(M - n1 + 1, j); // Block 3 Top left Boundary
		}
		/////////////////////////////////////////////////////////////
		for (int j = 6 * n1 + 1; j < 12 * n1; j++) {
			v_g(5 * n1 + 1, j) = -v_g(5 * n1, j); // C1 right boundary
			v_g(7 * n1, j) = -v_g(7 * n1 + 1, j); // D1 left Boundary
			v_g(13 * n1 + 1, j) = -v_g(13 * n1, j); // C2 right boundary
			v_g(15 * n1, j) = -v_g(15 * n1 + 1, j); // D2 left Boundary
			v_g(21 * n1 + 1, j) = -v_g(21 * n1, j); // C3 right boundary
			v_g(23 * n1, j) = -v_g(23 * n1 + 1, j); // D3 left Boundary
			v_g(29 * n1 + 1, j) = -v_g(29 * n1, j); // C4 right boundary
			v_g(31 * n1, j) = -v_g(31 * n1 + 1, j); // D4 left Boundary
		}
		/////////////////////////////////////////////////////////////
		}
	}

	// c calc
	{

}
	{
	// block 1		
	for (int j = 1; j < N + 1; j++) {
		for (int i = 1; i < n1 + 1; i++) {
			c_up(i, j) = c_calc(u_g(i, j), u_g(i - 1, j), v_g(i, j), v_g(i, j - 1),
				c_g(i, j), c_g(i + 1, j), c_g(i - 1, j), c_g(i, j + 1), c_g(i, j - 1));
		}
	}
	////////////////////////////////////////////////////////////////////////////////////////
	// block A1
	for (int j = 6 * n1 + 1; j < 8 * n1 + 1; j++) {
		for (int i = n1 + 1; i < 3 * n1 + 1; i++) {
			c_up(i, j) = c_calc(u_g(i, j), u_g(i - 1, j), v_g(i, j), v_g(i, j - 1),
				c_g(i, j), c_g(i + 1, j), c_g(i - 1, j), c_g(i, j + 1), c_g(i, j - 1));
		}
	}
	// block C1
	for (int j = 6 * n1 + 1; j < N + 1; j++) {
		for (int i = 3 * n1 + 1; i < 5 * n1 + 1; i++) {
			c_up(i, j) = c_calc(u_g(i, j), u_g(i - 1, j), v_g(i, j), v_g(i, j - 1),
				c_g(i, j), c_g(i + 1, j), c_g(i - 1, j), c_g(i, j + 1), c_g(i, j - 1));
		}
	}
	// block B1
	for (int j = j = 12 * n1 + 1; j < N + 1; j++) {
		for (int i = 5 * n1 + 1; i < 7 * n1 + 1; i++) {
			c_up(i, j) = c_calc(u_g(i, j), u_g(i - 1, j), v_g(i, j), v_g(i, j - 1),
				c_g(i, j), c_g(i + 1, j), c_g(i - 1, j), c_g(i, j + 1), c_g(i, j - 1));
		}
	}
	// block D1
	for (int j = 6 * n1 + 1; j < N + 1; j++) {
		for (int i = 7 * n1 + 1; i < 9 * n1 + 1; i++) {
			c_up(i, j) = c_calc(u_g(i, j), u_g(i - 1, j), v_g(i, j), v_g(i, j - 1),
				c_g(i, j), c_g(i + 1, j), c_g(i - 1, j), c_g(i, j + 1), c_g(i, j - 1));
		}
	}
	////////////////////////////////////////////////////////////////////////////////////////
	// block A2
	for (int j = 6 * n1 + 1; j < 8 * n1 + 1; j++) {
		for (int i = 9 * n1 + 1; i < 11 * n1 + 1; i++) {
			c_up(i, j) = c_calc(u_g(i, j), u_g(i - 1, j), v_g(i, j), v_g(i, j - 1),
				c_g(i, j), c_g(i + 1, j), c_g(i - 1, j), c_g(i, j + 1), c_g(i, j - 1));
		}
	}
	// block C2
	for (int j = 6 * n1 + 1; j < N + 1; j++) {
		for (int i = 11 * n1 + 1; i < 13 * n1 + 1; i++) {
			c_up(i, j) = c_calc(u_g(i, j), u_g(i - 1, j), v_g(i, j), v_g(i, j - 1),
				c_g(i, j), c_g(i + 1, j), c_g(i - 1, j), c_g(i, j + 1), c_g(i, j - 1));
		}
	}
	// block B2
	for (int j = j = 12 * n1 + 1; j < N + 1; j++) {
		for (int i = 13 * n1 + 1; i < 15 * n1 + 1; i++) {
			c_up(i, j) = c_calc(u_g(i, j), u_g(i - 1, j), v_g(i, j), v_g(i, j - 1),
				c_g(i, j), c_g(i + 1, j), c_g(i - 1, j), c_g(i, j + 1), c_g(i, j - 1));
		}
	}
	// block D2
	for (int j = 6 * n1 + 1; j < N + 1; j++) {
		for (int i = 15 * n1 + 1; i < 17 * n1 + 1; i++) {
			c_up(i, j) = c_calc(u_g(i, j), u_g(i - 1, j), v_g(i, j), v_g(i, j - 1),
				c_g(i, j), c_g(i + 1, j), c_g(i - 1, j), c_g(i, j + 1), c_g(i, j - 1));
		}
	}
	////////////////////////////////////////////////////////////////////////////////////////
	// block A3
	for (int j = 6 * n1 + 1; j < 8 * n1 + 1; j++) {
		for (int i = 17 * n1 + 1; i < 19 * n1 + 1; i++) {
			c_up(i, j) = c_calc(u_g(i, j), u_g(i - 1, j), v_g(i, j), v_g(i, j - 1),
				c_g(i, j), c_g(i + 1, j), c_g(i - 1, j), c_g(i, j + 1), c_g(i, j - 1));
		}
	}
	// block C3
	for (int j = 6 * n1 + 1; j < N + 1; j++) {
		for (int i = 19 * n1 + 1; i < 21 * n1 + 1; i++) {
			c_up(i, j) = c_calc(u_g(i, j), u_g(i - 1, j), v_g(i, j), v_g(i, j - 1),
				c_g(i, j), c_g(i + 1, j), c_g(i - 1, j), c_g(i, j + 1), c_g(i, j - 1));
		}
	}
	// block B3
	for (int j = j = 12 * n1 + 1; j < N + 1; j++) {
		for (int i = 21 * n1 + 1; i < 23 * n1 + 1; i++) {
			c_up(i, j) = c_calc(u_g(i, j), u_g(i - 1, j), v_g(i, j), v_g(i, j - 1),
				c_g(i, j), c_g(i + 1, j), c_g(i - 1, j), c_g(i, j + 1), c_g(i, j - 1));
		}
	}
	// block D3
	for (int j = 6 * n1 + 1; j < N + 1; j++) {
		for (int i = 23 * n1 + 1; i < 25 * n1 + 1; i++) {
			c_up(i, j) = c_calc(u_g(i, j), u_g(i - 1, j), v_g(i, j), v_g(i, j - 1),
				c_g(i, j), c_g(i + 1, j), c_g(i - 1, j), c_g(i, j + 1), c_g(i, j - 1));
		}
	}
	////////////////////////////////////////////////////////////////////////////////////////
	// block A4
	for (int j = 6 * n1 + 1; j < 8 * n1 + 1; j++) {
		for (int i = 25 * n1 + 1; i < 27 * n1 + 1; i++) {
			c_up(i, j) = c_calc(u_g(i, j), u_g(i - 1, j), v_g(i, j), v_g(i, j - 1),
				c_g(i, j), c_g(i + 1, j), c_g(i - 1, j), c_g(i, j + 1), c_g(i, j - 1));
		}
	}
	// block C4
	for (int j = 6 * n1 + 1; j < N + 1; j++) {
		for (int i = 27 * n1 + 1; i < 29 * n1 + 1; i++) {
			c_up(i, j) = c_calc(u_g(i, j), u_g(i - 1, j), v_g(i, j), v_g(i, j - 1),
				c_g(i, j), c_g(i + 1, j), c_g(i - 1, j), c_g(i, j + 1), c_g(i, j - 1));
		}
	}
	// block B4
	for (int j = j = 12 * n1 + 1; j < N + 1; j++) {
		for (int i = 29 * n1 + 1; i < 31 * n1 + 1; i++) {
			c_up(i, j) = c_calc(u_g(i, j), u_g(i - 1, j), v_g(i, j), v_g(i, j - 1),
				c_g(i, j), c_g(i + 1, j), c_g(i - 1, j), c_g(i, j + 1), c_g(i, j - 1));
		}
	}
	// block D4
	for (int j = 6 * n1 + 1; j < N + 1; j++) {
		for (int i = 31 * n1 + 1; i < 33 * n1 + 1; i++) {
			c_up(i, j) = c_calc(u_g(i, j), u_g(i - 1, j), v_g(i, j), v_g(i, j - 1),
				c_g(i, j), c_g(i + 1, j), c_g(i - 1, j), c_g(i, j + 1), c_g(i, j - 1));
		}
	}
	// block A5
	for (int j = 6 * n1 + 1; j < 8 * n1 + 1; j++) {
		for (int i = 33 * n1 + 1; i < 35 * n1 + 1; i++) {
			c_up(i, j) = c_calc(u_g(i, j), u_g(i - 1, j), v_g(i, j), v_g(i, j - 1),
				c_g(i, j), c_g(i + 1, j), c_g(i - 1, j), c_g(i, j + 1), c_g(i, j - 1));
		}
	}
	////////////////////////////////////////////////////////////////////////////////////////////
	/*
	// block 2
	for (int j = 6 * n1 + 1; j < 8 * n1 + 1; j++) {
	for (int i = n1 + 1; i < M + 1 - n1; i++) {
	P_up(i, j) = (1.0 / (2.0 * (1 + beta))) * (P_t(i - 1, j) + P_t(i + 1, j) + beta*(P_t(i, j - 1) + P_t(i, j + 1))
	- (pow(dx, 2) * bij(uf(i, j), uf(i - 1, j), vf(i, j), vf(i, j - 1))));
	}
	}
	*/
	// block 3
	for (int j = 1; j < N + 1; j++) {
		for (int i = M + 1 - n1; i < M + 1; i++) {
			c_up(i, j) = c_calc(u_g(i, j), u_g(i - 1, j), v_g(i, j), v_g(i, j - 1),
				c_g(i, j), c_g(i + 1, j), c_g(i - 1, j), c_g(i, j + 1), c_g(i, j - 1));
		}
	}
}

// c at ghost nodes
{

}
	{
	// For Left and Right Boundary
	for (int j = 1; j < N + 1; j++) {
		c_up(0, j) = c_up(1, j);
		c_up(M + 1, j) = c_up(M, j);
	}
	// For top and bottom inlet
	for (int i = 1; i < n1 + 1; i++) {
		c_up(i, 0) = 2.0 - c_up(i, 1);
		c_up(i, N + 1) = -c_up(i, N);
	}
	// For Top and Bottom except the inlet
	for (int i = 3 * n1 + 1; i < M + 1; i++) {
		c_up(i, 0) = c_up(i, 1);
		c_up(i, N + 1) = c_up(i, N);
	}
	// For Block 1(Bottom Right Walls) and Block 2(Bottom left walls)
	for (int j = 1; j < 6 * n1 + 1; j++) {
		c_up(n1 + 1, j) = c_up(n1, j);
		c_up(M - n1, j) = c_up(M + 1 - n1, j);
	}
	// For Block 1(Top Right Walls) and Block 2(Top Left walls)
	for (int j = 8 * n1 + 1; j < N + 1; j++) {
		c_up(n1 + 1, j) = c_up(n1, j);
		c_up(M - n1, j) = c_up(M + 1 - n1, j);
	}
	// For all the walls at j=60(Ghost Nodes) in Serpentine Geometry
	for (int i = n1 + 1; i < M + 1 - n1; i++) {
		c_up(i, 6 * n1) = c_up(i, 6 * n1 + 1);
		//c_up(i, 8 * n1 + 1) = c_up(i, 8 * n1); //Commented because of Serpentine Geometry
	}
	/////////////////////////////////////////////////////////////////////////////////////
	// For all the walls at j=81(Ghost Nodes) in Serpentine Geometry
	for (int i = n1 + 1; i < 3 * n1 + 1; i++) {
		c_up(i, 8 * n1 + 1) = c_up(i, 8 * n1);
	} // A1 top boundary
	for (int i = 9 * n1 + 1; i < 11 * n1 + 1; i++) {
		c_up(i, 8 * n1 + 1) = c_up(i, 8 * n1);
	} // A2 top boundary
	for (int i = 17 * n1 + 1; i < 19 * n1 + 1; i++) {
		c_up(i, 8 * n1 + 1) = c_up(i, 8 * n1);
	} // A3 top boundary
	for (int i = 25 * n1 + 1; i < 27 * n1 + 1; i++) {
		c_up(i, 8 * n1 + 1) = c_up(i, 8 * n1);
	} // A4 top boundary
	for (int i = 33 * n1 + 1; i < 35 * n1 + 1; i++) {
		c_up(i, 8 * n1 + 1) = c_up(i, 8 * n1);
	} // A5 top boundary
	  /////////////////////////////////////////////////////////////////////////////////////
	  // For all the walls lying between j = 8*n1+1 to N
	for (int j = 8 * n1 + 1; j < N + 1; j++) {
		c_up(3 * n1, j) = c_up(3 * n1 + 1, j); //Block C1 left Boundary
		c_up(9 * n1 + 1, j) = c_up(9 * n1, j); //Block D1 right Boundary
		c_up(11 * n1, j) = c_up(11 * n1 + 1, j); //Block C2 left Boundary
		c_up(17 * n1 + 1, j) = c_up(17 * n1, j); //Block D2 right Boundary
		c_up(19 * n1, j) = c_up(19 * n1 + 1, j); //Block C3 left Boundary
		c_up(25 * n1 + 1, j) = c_up(25 * n1, j); //Block D3 right Boundary
		c_up(27 * n1, j) = c_up(27 * n1 + 1, j); //Block C4 left Boundary
		c_up(33 * n1 + 1, j) = c_up(33 * n1, j); //Block D4 right Boundary
	}
	/////////////////////////////////////////////////////////////////////////////////////
	// For all the walls lying between j = 6*n1+1 to 12*n1
	for (int j = 6 * n1 + 1; j < 12 * n1 + 1; j++) {
		c_up(5 * n1 + 1, j) = c_up(5 * n1, j); //Block C1 right Boundary
		c_up(7 * n1, j) = c_up(7 * n1 + 1, j); //Block D1 left Boundary
		c_up(13 * n1 + 1, j) = c_up(13 * n1, j); //Block C2 right Boundary
		c_up(15 * n1, j) = c_up(15 * n1 + 1, j); //Block D2 left Boundary
		c_up(21 * n1 + 1, j) = c_up(21 * n1, j); //Block C3 right Boundary
		c_up(23 * n1, j) = c_up(23 * n1 + 1, j); //Block D3 left Boundary
		c_up(29 * n1 + 1, j) = c_up(29 * n1, j); //Block C4 right Boundary
		c_up(31 * n1, j) = c_up(31 * n1 + 1, j); //Block D4 left Boundary
	}
	/////////////////////////////////////////////////////////////////////////////////////
	// For all the walls lying at j = 12*n1
	for (int i = 5 * n1 + 1; i < 7 * n1 + 1; i++) {
		c_up(i, 12 * n1) = c_up(i, 12 * n1 + 1);
	} // For B1 bottom boundary
	for (int i = 13 * n1 + 1; i < 15 * n1 + 1; i++) {
		c_up(i, 12 * n1) = c_up(i, 12 * n1 + 1);
	} // For B2 bottom boundary
	for (int i = 21 * n1 + 1; i < 23 * n1 + 1; i++) {
		c_up(i, 12 * n1) = c_up(i, 12 * n1 + 1);
	} // For B3 bottom boundary
	for (int i = 29 * n1 + 1; i < 31 * n1 + 1; i++) {
		c_up(i, 12 * n1) = c_up(i, 12 * n1 + 1);
	} // For B4 bottom boundary
}
	
	P_g = P_up;
	c_g = c_up;
	cout << "time step: " << n + 1 << "\n";

	u = u_g.block(0, 1, M + 1, N);
	v = v_g.block(1, 0, M, N + 1);
	P = P_g.block(1, 1, M, N);
	c = c_g.block(1, 1, M, N);

	// u_avg
	{
	// block 1
	for (int i = 0; i < n1; i++) {
		for (int j = 0; j < N; j++) {
			u_avg(i, j) = (u(i, j) + u(i + 1, j)) / 2;
			v_avg(i, j) = (v(i, j) + v(i, j + 1)) / 2;
		}
	}

	// block A1
	for (int j = 6 * n1; j < 8 * n1; j++) {
		for (int i = n1; i < 3 * n1; i++) {
			u_avg(i, j) = (u(i, j) + u(i + 1, j)) / 2;
			v_avg(i, j) = (v(i, j) + v(i, j + 1)) / 2;
		}
	}

	// block C1
	for (int j = 6 * n1; j < N; j++) {
		for (int i = 3 * n1; i < 5 * n1; i++) {
			u_avg(i, j) = (u(i, j) + u(i + 1, j)) / 2;
			v_avg(i, j) = (v(i, j) + v(i, j + 1)) / 2;
		}
	}

	// block B1
	for (int j = 12 * n1; j < N; j++) {
		for (int i = 5 * n1; i < 7 * n1; i++) {
			u_avg(i, j) = (u(i, j) + u(i + 1, j)) / 2;
			v_avg(i, j) = (v(i, j) + v(i, j + 1)) / 2;
		}
	}

	// block D1
	for (int j = 6 * n1; j < N; j++) {
		for (int i = 7 * n1; i < 9 * n1; i++) {
			u_avg(i, j) = (u(i, j) + u(i + 1, j)) / 2;
			v_avg(i, j) = (v(i, j) + v(i, j + 1)) / 2;
		}
	}
	//////////////////////////////////////////////////////////////////////////////////////////////////
	// block A2
	for (int j = 6 * n1; j < 8 * n1; j++) {
		for (int i = 9 * n1; i < 11 * n1; i++) {
			u_avg(i, j) = (u(i, j) + u(i + 1, j)) / 2;
			v_avg(i, j) = (v(i, j) + v(i, j + 1)) / 2;
		}
	}

	// block C2
	for (int j = 6 * n1; j < N; j++) {
		for (int i = 11 * n1; i < 13 * n1; i++) {
			u_avg(i, j) = (u(i, j) + u(i + 1, j)) / 2;
			v_avg(i, j) = (v(i, j) + v(i, j + 1)) / 2;
		}
	}

	// block B2
	for (int j = 12 * n1; j < N; j++) {
		for (int i = 13 * n1; i < 15 * n1; i++) {
			u_avg(i, j) = (u(i, j) + u(i + 1, j)) / 2;
			v_avg(i, j) = (v(i, j) + v(i, j + 1)) / 2;
		}
	}

	// block D2
	for (int j = 6 * n1; j < N; j++) {
		for (int i = 15 * n1; i < 17 * n1; i++) {
			u_avg(i, j) = (u(i, j) + u(i + 1, j)) / 2;
			v_avg(i, j) = (v(i, j) + v(i, j + 1)) / 2;
		}
	}
	//////////////////////////////////////////////////////////////////////////////////////////////////
	// block A3
	for (int j = 6 * n1; j < 8 * n1; j++) {
		for (int i = 17 * n1; i < 19 * n1; i++) {
			u_avg(i, j) = (u(i, j) + u(i + 1, j)) / 2;
			v_avg(i, j) = (v(i, j) + v(i, j + 1)) / 2;;
		}
	}

	// block C3
	for (int j = 6 * n1; j < N; j++) {
		for (int i = 19 * n1; i < 21 * n1; i++) {
			u_avg(i, j) = (u(i, j) + u(i + 1, j)) / 2;
			v_avg(i, j) = (v(i, j) + v(i, j + 1)) / 2;
		}
	}

	// block B3
	for (int j = 12 * n1; j < N; j++) {
		for (int i = 21 * n1; i < 23 * n1; i++) {
			u_avg(i, j) = (u(i, j) + u(i + 1, j)) / 2;
			v_avg(i, j) = (v(i, j) + v(i, j + 1)) / 2;
		}
	}

	// block D3
	for (int j = 6 * n1; j < N; j++) {
		for (int i = 23 * n1; i < 25 * n1; i++) {
			u_avg(i, j) = (u(i, j) + u(i + 1, j)) / 2;
			v_avg(i, j) = (v(i, j) + v(i, j + 1)) / 2;
		}
	}
	//////////////////////////////////////////////////////////////////////////////////////////////////
	// block A4
	for (int j = 6 * n1; j < 8 * n1; j++) {
		for (int i = 25 * n1; i < 27 * n1; i++) {
			u_avg(i, j) = (u(i, j) + u(i + 1, j)) / 2;
			v_avg(i, j) = (v(i, j) + v(i, j + 1)) / 2;
		}
	}

	// block C4
	for (int j = 6 * n1; j < N; j++) {
		for (int i = 27 * n1; i < 29 * n1; i++) {
			u_avg(i, j) = (u(i, j) + u(i + 1, j)) / 2;
			v_avg(i, j) = (v(i, j) + v(i, j + 1)) / 2;
		}
	}

	// block B4
	for (int j = 12 * n1; j < N; j++) {
		for (int i = 29 * n1; i < 31 * n1; i++) {
			u_avg(i, j) = (u(i, j) + u(i + 1, j)) / 2;
			v_avg(i, j) = (v(i, j) + v(i, j + 1)) / 2;
		}
	}

	// block D4
	for (int j = 6 * n1; j < N; j++) {
		for (int i = 31 * n1; i < 33 * n1; i++) {
			u_avg(i, j) = (u(i, j) + u(i + 1, j)) / 2;
			v_avg(i, j) = (v(i, j) + v(i, j + 1)) / 2;
		}
	}
	//////////////////////////////////////////////////////////////////////////////////////////////////
	// block A5
	for (int j = 6 * n1; j < 8 * n1; j++) {
		for (int i = 33 * n1; i < 35 * n1; i++) {
			u_avg(i, j) = (u(i, j) + u(i + 1, j)) / 2;
			v_avg(i, j) = (v(i, j) + v(i, j + 1)) / 2;
		}
	}
	//////////////////////////////////////////////////////////////////////////////////////////////////
	/*
	// block 2
	for (int i = n1; i < M - n1; i++) {
		//for (int j = 6 * n1; j < 8 * n1; j++) {
		for (int j = 6 * n1; j < 14 * n1; j++) {
			u_avg(i, j) = (u(i, j) + u(i + 1, j)) / 2;
			v_avg(i, j) = (v(i, j) + v(i, j + 1)) / 2;
		}
	}
	*/
	// block 3 
	for (int i = M - n1; i < M; i++) {
		for (int j = 0; j < N; j++) {
			u_avg(i, j) = (u(i, j) + u(i + 1, j)) / 2;
			v_avg(i, j) = (v(i, j) + v(i, j + 1)) / 2;
		}
	}
	}

	if (n % 200 == 0) {
		ofstream file_uavn("u_avg_" + to_string(n) + ".csv");
		ofstream file_vavn("v_avg_" + to_string(n) + ".csv");
		ofstream file_Pn("P_n_" + to_string(n) + ".csv");
		ofstream file_cn("c_n_" + to_string(n) + ".csv");

		for (int i = 0; i < M; i++) {
			for (int j = 0; j < N; j++) {
				file_uavn << u_avg(i, j) << ",";
				file_vavn << v_avg(i, j) << ",";
				file_Pn << P(i, j) << ",";
				file_cn << c(i, j) << ",";
			}
			file_uavn << "\n";
			file_vavn << "\n";
			file_Pn << "\n";
			file_cn << "\n";
		}
		file_uavn.close();
		file_vavn.close();
		file_Pn.close();
		file_cn.close();
	}
	

	}


	//-----------------------------------------------------------
	/*
	ofstream file_P("P_g.csv");
	ofstream file_X("u_g.csv");
	ofstream file_Y("v_g.csv");
	ofstream file_uf("u_f.csv");
	ofstream file_vf("v_f.csv");
	ofstream file_uav("u_avg.csv");
	ofstream file_vav("v_avg.csv");

	for (int i = 0; i < M + 2; i++) {
		for (int j = 0; j < N + 1; j++) {
			file_Y << v_g(i, j) << ",";
			file_vf << vf(i, j) << ",";
		}
		file_Y << "\n";
		file_vf << "\n";
	}
	for (int i = 0; i < M + 1; i++) {
		for (int j = 0; j < N + 2; j++) {
			file_X << u_g(i, j) << ",";
			file_uf << uf(i, j) << ",";
		}
		file_X << "\n";
		file_uf << "\n";
	}

	for (int i = 0; i < M + 2; i++) {
		for (int j = 0; j < N + 2; j++) {
			file_P << P_g(i, j) << ",";
		}
		file_P << "\n";
	}

	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) {
			file_uav << u_avg(i, j) << ",";
			file_vav << v_avg(i, j) << ",";
		}
		file_uav << "\n";
		file_vav << "\n";
	}

	file_P.close();
	file_X.close();
	file_Y.close();
	file_uf.close();
	file_vf.close();
	file_uav.close();
	file_vav.close();
	*/

	//return { u,v,P,u_avg,v_avg };
	//return { u_g,v_g,P_g,u_avg,v_avg };
}

double Proj::uf_calc(double uij, double uip1j, double uim1j, double uijp1, double uijm1, 
					double vij, double vip1j, double vijm1, double vip1jm1, double Pij, double Pip1j) {
	double Auu = -(1.0 / dx) * ( pow((uij + uip1j ) / 2.0, 2) - pow( (uij + uim1j) / 2.0, 2) );
	double Auv = -(1.0 / dy) * ((vij + vip1j)*(uij + uijp1) / 4.0 - (vijm1 + vip1jm1)*(uijm1 + uij) / 4.0);
	//double Auu = -(uij / dx)*((uij + uip1j) / 2 - (uij + uim1j) / 2);
	//double Auv = -(1.0 / (2.0*dy)) * ((vij + vip1j) / 2 + (vijm1 + vip1jm1) / 2)*((uij + uijp1) / 2 - (uijm1 + uij) / 2);
	double A = Auu + Auv;
	double B = (mu / rho)  *  ( (uip1j - 2 * uij + uim1j) / pow(dx, 2) + (uijp1 - 2 * uij + uijm1) / pow(dy, 2) );
	//double P = (Pip1j - Pij) / dx;
	//double uf = uij + dt*(A + B - P);
	double uf = uij + dt*(A + B);
	return uf;
}

double Proj::vf_calc(double vij, double vijp1, double vijm1, double vip1j, double vim1j, 
					double uij, double uijp1, double uim1j, double uim1jp1, double Pij, double Pijp1) {
	double Avv = -(1.0 / dy) * ( pow((vij + vijp1) / 2.0, 2) - pow((vij + vijm1) / 2.0, 2) );
	double Avu = -(1.0 / dx) * ( (uijp1 + uij)*(vij + vip1j) / 4 - (uim1j + uim1jp1)*(vim1j + vij) / 4 );
	//double Avv = -(vij / dy)*((vij + vijp1) / 2 - (vij + vijm1) / 2);
	//double Avu = -(1.0 / (2.0*dx))*((uij + uijp1) / 2 + (uim1j + uim1jp1) / 2)*((vij + vip1j) / 2 - (vij + vim1j) / 2);
	double A = Avv + Avu;
	double B = (mu / rho)  *  ( (vijp1 - 2 * vij + vijm1) / pow(dy, 2) + (vip1j - 2 * vij + vim1j) / pow(dx, 2) );
	//double P = (Pijp1 - Pij) / dy;
	//double vf = vij + dt*(A + B - P);
	double vf = vij + dt*(A + B);
	return vf;
}

// Calculation of bij term for Pressure Poisson eqn (Point Jacobi Method)
double Proj::bij(double uij, double uim1j, double vij, double vijm1) { 
	return (1 / dt) * (((uij - uim1j) / dx) + ((vij - vijm1) / dy));
}

// Conc function
double Proj::c_calc(double uij, double uim1j, double vij, double vijm1,
	double cij, double cip1j, double cim1j, double cijp1, double cijm1) {
	double Auc = -(1.0 / dx)*(uij*(cip1j + cij) / 2 - uim1j * (cij + cim1j) / 2);
	double Avc = -(1.0 / dy)*(vij*(cijp1 + cij) / 2 - vijm1 * (cij + cijm1) / 2);
	double A = Auc + Avc;
	double Bx = D * (1.0 / pow(dx, 2))*(cip1j - 2 * cij + cim1j);
	double By = D * (1.0 / pow(dy, 2))*(cijp1 - 2 * cij + cijm1);
	double B = Bx + By;
	double cnp1 = cij + dt * (A + B);
	return cnp1;
}

/*
//similar to matlab meshgrid
tuple<MatrixXd, MatrixXd> Proj::mymesh(int r, int c, double distx, double dist2x, double disty, double dist2y)
{
	double *x, *y;
	int element, inputx = 0, inputy = 0, count = 0;
	element = r*c;

	x = new double[element];
	y = new double[element];
	for (int i = 0; i < element; i++)
	{
		++count;
		if (i%c == 0)
		{
			++inputx;
		}
		inputy = count;
		if (count%c == 0)
		{
			count = 0;
		}
		x[i] = dist2x*(inputx - 1) + distx;
		y[i] = dist2y*(inputy - 1) + disty;
	}

	MatrixXd Xmat = MatrixXd::Constant(M, N, 1.233);
	MatrixXd Ymat = MatrixXd::Constant(M, N, -231.233);

	for (int j = 0; j < N; j++)
	{
		for (int i = 0; i < M; i++)
		{
			Xmat(i, j) = x[i*N + j];
			Ymat(i, j) = y[i*N + j];
		}
	}
	return { Xmat, Ymat };
}
*/