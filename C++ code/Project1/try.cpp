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