uninstallPackage("SimplicialMDSplines")
installPackage("SimplicialMDSplines",FileName=>"./package/SimplicialMDSplines.m2")
--Degree range
a = 0;b = 20;
--Triangulation
V = {{0/1,0/1}, {18/1,0/1}, {9/1,15/1}, {9/1,3/1}, {11/1,6/1}, {7/1,6/1}};
F = {{0,1,3}, {1,3,4}, {1,2,4}, {2,4,5}, {0,2,5}, {0,3,5}, {3,4,5}};
E = {{0,3}, {0,5}, {1,3}, {1,4}, {2,4}, {2,5}, {3,4}, {3,5}, {4,5}};

--Spline space configuration
DF = {0, 0, 0, 0, 0, 0, 0};
re = {1, 1, 1, 1, 1, 1, 1, 1, 1};
rv = {2, 2, 2, 2, 2, 2};

--Homology dimension table
homologyDimensionTable(a,b,{V,F,{E,re,rv},DF})