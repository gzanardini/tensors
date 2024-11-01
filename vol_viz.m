%% Slicing volume visual
function vol_viz(G_test, title_str)
x=0:1:9;
y=0:1:9;
z=0:1:9;

[x,y,z]=meshgrid(x,y,z);

xslice = [1,2,3,4,5,6,7,8,9];    %[3,3,3]   % location of y-z planes
yslice = [];                     % location of x-z plane
zslice = [];    %[2,0]     % location of x-y planes

figure();
slice(x,y,z,G_test,xslice,yslice,zslice)
xlabel('x')
ylabel('y')
zlabel('z')
title(title_str)