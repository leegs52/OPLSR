function [r1 r2 r3] = gen3Dsample3(n1, n2, n3)


drawfig = false;
msize = 10;
if nargin < 1
    n1 = 15;
    n2 = 15;
    n3 = 15;    
    drawfig = true;
end

mu1 = [1 2 3];
U = [   1   -1  1;
        1    1  1;
        1    1  -1];
Lam = diag([10 1 1]);
SIGMA = U*Lam*U';
r1 = mvnrnd(mu1,SIGMA,n1);

mu2 = [5 5 1];
%rotangle = pi/9;
rotangle = 0;
U1 = U*rotationmat3D(rotangle,[0 0 1]);
Lam1 = diag([10 1 1]);
SIGMA1 = U1*Lam1*U1';
r2 = mvnrnd(mu2,SIGMA1,n2);

mu3 = mu2 + (mu2 - mu1);
rotangle = 0;
U1 = U*rotationmat3D(rotangle,[0 0 1]);
Lam1 = diag([10 1 1]);
SIGMA1 = U1*Lam1*U1';
r3 = mvnrnd(mu3,SIGMA1,n3);


if (drawfig) 
    figure(1); clf;
    plot3(r1(:,1),r1(:,2), r1(:,3), '+'); hold on;
    plot3(r2(:,1),r2(:,2), r2(:,3), 'ro');
    plot3(r3(:,1),r3(:,2), r3(:,3), 'gs');
    xlabel('X'); ylabel('Y'); zlabel('Z');
end
