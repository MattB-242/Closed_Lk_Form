# Closed_Lk_Form
Uses AT formula to calculate the linking number of two straight line segments in R^3 given their endpoints as inputs

Line pairs are instantiated as SEGATAN objects of the form [[[a1,a2,a3][b1,b2,b3]],[[c1,c2,c3][d1,d2,d3]]]
where a_i, b_i are the endpoint co-ordinates of the first segment and c_i, d_i the endpoint co-ordinates of the final segment

Lk calculation is done by the method self.seglink()

Also includes two test calculations for the Hopf link and the periodic linking number between a single line segment and a 
1-periodic lattice
