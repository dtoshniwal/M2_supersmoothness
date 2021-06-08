uninstallPackage("SimplicialMDSplines")
installPackage("SimplicialMDSplines",FileName=>"./package/SimplicialMDSplines.m2")
--Degree range
a = 0;b = 20;
--Triangulation
V = {{0/1,0/1}, {18/1,0/1}, {9/1,15/1}, {9/1,3/1}, {11/1,6/1}, {7/1,6/1}, {9/1,3/2}, {113/10,19/5}, {123/10,34/5}, {9/1,38/5}, {57/10,34/5}, {67/10,19/5}, {9/1,49/10}, {9/1,0/1}, {1533/340,511/68}, {63/8,21/8}, {1673/270,239/45}, {4587/340,511/68}, {81/8,21/8}, {3187/270,239/45}, {15087/1405,10086/1405}, {10203/1405,10086/1405}, {4532/455,4041/910}, {3658/455,4041/910}, {9/1,6/1}};
F = {{0,6,13}, {1,6,13}, {1,6,18}, {3,6,18}, {3,6,15}, {0,6,15}, {1,7,19}, {4,7,19}, {4,7,22}, {3,7,22}, {3,7,18}, {1,7,18}, {1,8,17}, {2,8,17}, {2,8,20}, {4,8,20}, {4,8,19}, {1,8,19}, {2,9,21}, {5,9,21}, {5,9,24}, {4,9,24}, {4,9,20}, {2,9,20}, {2,10,14}, {0,10,14}, {0,10,16}, {5,10,16}, {5,10,21}, {2,10,21}, {0,11,15}, {3,11,15}, {3,11,23}, {5,11,23}, {5,11,16}, {0,11,16}, {3,12,22}, {4,12,22}, {4,12,24}, {5,12,24}, {5,12,23}, {3,12,23}};
E = {{0,6}, {0,10}, {0,11}, {0,15}, {0,16}, {1,6}, {1,7}, {1,8}, {1,18}, {1,19}, {2,8}, {2,9}, {2,10}, {2,20}, {2,21}, {3,6}, {3,7}, {3,11}, {3,12}, {3,15}, {3,18}, {3,22}, {3,23}, {4,7}, {4,8}, {4,9}, {4,12}, {4,19}, {4,20}, {4,22}, {4,24}, {5,9}, {5,10}, {5,11}, {5,12}, {5,16}, {5,21}, {5,23}, {5,24}, {6,13}, {6,15}, {6,18}, {7,18}, {7,19}, {7,22}, {8,17}, {8,19}, {8,20}, {9,20}, {9,21}, {9,24}, {10,14}, {10,16}, {10,21}, {11,15}, {11,16}, {11,23}, {12,22}, {12,23}, {12,24}};

--Spline space configuration
DF = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
re = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};
rv = {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};

--Homology dimension table
homologyDimensionTable(a,b,{V,F,{E,re,rv},DF})