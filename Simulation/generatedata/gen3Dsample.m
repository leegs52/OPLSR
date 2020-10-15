function [r r1] = gen3Dsample(n, n1)

drawfig = false;

% n = 200;
% n1 = 300;
% n = 25;
% n1 = 15;
% drawfig = true;

mu = [1 2 3];
U = [   1   -1  1;
        1    1  1;
        1    1  -1];
Lam = diag([10 1 1]);
SIGMA = U*Lam*U';
r = mvnrnd(mu,SIGMA,n);

mu1 = [5 5 1];
%rotangle = pi/9;
rotangle = 0;
U1 = U*rotationmat3D(rotangle,[0 0 1]);
Lam1 = diag([10 1 1]);
SIGMA1 = U1*Lam1*U1';
r1 = mvnrnd(mu1,SIGMA1,n1);

if (drawfig) 
    figure(1); clf;
    plot3(r(:,1),r(:,2), r(:,3), '+'); hold on;
    plot3(r1(:,1),r1(:,2), r1(:,3), 'ro');
    xlabel('X'); ylabel('Y'); zlabel('Z');
end
