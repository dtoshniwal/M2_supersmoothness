R = QQ[x,y,z];

M = (a,b,s) -> (ideal(a,b))^(s+1);
I0 = (a,r) -> ideal(a^(r+1));
I1 = (a,b,r,s) -> intersect(I0(b,r), M(a,b,s));
I2 = (a,b,c,r,s0,s1) -> intersect(I0(b,r), intersect(M(a,b,s0),M(b,c,s1)));
I3 = (a,b,c,r,s0,s1,k) -> (1-k)*I2(a,b,a,r,s0,s0) + k*I2(c,b,c,r,s1,s1);

--Triangulation
V = {{47/50,41/100}, {11/20,67/100}, {73/100,93/100}, {29/50,81/100}, {3/100,12/25}, {9/20,19/25}, {13/20,21/50}, {13/25,97/100}, {37/100,99/100}, {47/50,43/50}, {83/100,39/100}, {17/20,9/20}, {37/100,1/4}, {59/100,39/50}, {87/100,22/25}, {93/100,91/100}, {67/100,14/25}, {21/100,3/5}, {13/20,3/20}, {7/100,9/10}};
E = {{0,10}, {0,11}, {1,5}, {1,6}, {1,12}, {1,13}, {1,16}, {1,17}, {2,3}, {2,7}, {2,13}, {2,14}, {2,15}, {3,5}, {3,7}, {3,13}, {4,17}, {5,7}, {5,8}, {5,13}, {5,17}, {5,19}, {6,10}, {6,11}, {6,12}, {6,16}, {6,18}, {9,11}, {9,14}, {9,16}, {10,11}, {10,18}, {11,16}, {12,17}, {13,14}, {13,16}, {14,15}, {14,16}, {17,19}};
V = apply(V, v->append(v,1));
varlist = vars R;
varCol := transpose varlist;

--Linear forms
M := (transpose(matrix(R,V)));
mM := numrows M;
minorList := apply(E, e-> gens gb minors(mM,matrix(M_e)|varCol));
L = apply(#minorList, i -> minorList_i_(0,0));
--Perpendicular linear forms
Lp0 = apply(E, e -> (Me := mutableMatrix(M_e);
	a := Me_(1,1)-Me_(1,0)+Me_(0,0);
	b := -Me_(0,1)+Me_(0,0)+Me_(1,0);
	Me_(0,1) = a;
	Me_(1,1) = b;
	gens gb minors(mM,matrix(Me)|varCol)));
Lp1 = apply(E, e -> (Me := mutableMatrix(M_e);
	a := Me_(1,0)-Me_(1,1)+Me_(0,1);
	b := -Me_(0,0)+Me_(0,1)+Me_(1,1);
	Me_(0,0) = a;
	Me_(1,0) = b;
	gens gb minors(mM,matrix(Me)|varCol)));
Lp = apply(#Lp0, i -> {Lp0_i_(0,0), Lp1_i_(0,0)});
Lp = flatten(Lp);

--Edge connectivity
EV = {0, 10, 0, 11, 1, 5, 1, 6, 1, 12, 1, 13, 1, 16, 1, 17, 2, 3, 2, 7, 2, 13, 2, 14, 2, 15, 3, 5, 3, 7, 3, 13, 4, 17, 5, 7, 5, 8, 5, 13, 5, 17, 5, 19, 6, 10, 6, 11, 6, 12, 6, 16, 6, 18, 9, 11, 9, 14, 9, 16, 10, 11, 10, 18, 11, 16, 12, 17, 13, 14, 13, 16, 14, 15, 14, 16, 17, 19};
--Edge smoothness
rE = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
--Vertex smoothness
rV = {3, 2, 4, 4, 4, 2, 4, 3, 2, 3, 2, 4, 2, 4, 4, 3, 3, 4, 4, 3};

--Compute dimension
--global polynomials
dpoly = d -> binomial(d+2,2);
--edge contribution
dedge = d -> (
	 n = 0;
	 for i from 0 to 38 do (
		 r = rE_i; s0 = rV_(EV_(2*i)); s1 = rV_(EV_(2*i+1));
		 Ii := I2(Lp_(2*i),L_(i),Lp_(2*i+1),r,s0,s1);
		 n = n + hilbertFunction(d,module Ii);
	 );
	 return n;);
--vertex contribution
V = {1, 2, 3, 5, 6, 10, 11, 13, 14, 16, 17};
E = {{2,3,4,5,6,7}, {8,9,10,11,12}, {8,13,14,15}, {2,13,17,18,19,20,21}, {3,22,23,24,25,26}, {0,22,30,31}, {1,23,27,30,32}, {5,10,15,19,34,35}, {11,28,34,36,37}, {6,25,29,32,35,37}, {7,16,20,33,38}};
dvertex = d -> (
	 n = 0;
	 for i from 0 to 10 do (
		 Ii := 0;
		 for j in E_i do (
			 r = rE_j; s0 = rV_(EV_(2*j)); s1 = rV_(EV_(2*j+1));
			 Ii = Ii + I2(Lp_(2*j),L_(j),Lp_(2*j+1),r,s0,s1);
		);
		 n = n + hilbertFunction(d,module Ii);
	 );
	 return n;);
--approximate vertex contribution
Ep = {{0,0,0,0,0,0}, {0,0,0,0,0}, {1,0,0,0}, {1,1,0,0,0,0,0}, {1,0,0,0,0,0}, {1,1,0,0}, {1,1,1,1,0}, {1,1,1,1,0,0}, {1,1,1,0,0}, {1,1,1,1,1,1}, {1,1,1,1,0}};
dvertexa = d -> (
	 n = 0;
	 for i from 0 to 10 do (
		 Ii := 0;
		 for j from 0 to (length(E_i)-1) do (
			 k = E_i_j;
			 r = rE_k; s0 = rV_(EV_(2*k)); s1 = rV_(EV_(2*k+1));
			 Ii = Ii + I3(Lp_(2*k),L_(k),Lp_(2*k+1),r,s0,s1,Ep_i_j);
		);
		 n = n + hilbertFunction(d,module Ii);
	 );
	 return n;);

--Degree range
dmin = 4; dmax=7;

--Dimension
apply(dmin..dmax, d -> max(dpoly(d),dpoly(d)+dedge(d)-dvertex(d)))
print "---------------------------";
apply(dmin..dmax, d -> max(dpoly(d),dpoly(d)+dedge(d)-dvertexa(d)))
