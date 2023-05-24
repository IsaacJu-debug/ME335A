close all
clear
X = [-1, -1, 0, 0, 1, 1;...
    1, -1, 1, -1, 1, -1];
LV = [1, 3, 4, 3;...
    2, 2, 6, 6;...
    3, 4, 3, 5];

[row1,col1] = find(LV==2) ;[row2,col2] = find(LV==3); 
C = intersect(col1,col2);
d = setdiff(C,1);