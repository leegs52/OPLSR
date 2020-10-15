function [r1 r2 r3] = gen2Dsample3(n1, n2, n3)

drawfig = false;
msize = 10;
if nargin < 1
    n1 = 15;
    n2 = 15;
    n3 = 15;    
    drawfig = true;
end

 
A = randn(2, n1);
%B = [1 0.6; 0.6 0.7];
%A = randn(2, 10);
%B = [1 0.4; 0.4 0.4];
B = [1 0.4; 0.4 0.4];
r1 = B * A;


A = randn(2, n2);
%B = [1 0.6; 0.6 0.7];
%A = randn(2, 10);
%B = [1 0.4; 0.4 0.4];
r2 = B * A + repmat([0 0.5]', 1, n2);

A = randn(2, n3);
r3 = B * A + repmat([-0.2 0.9]', 1, n3);

if (drawfig)
    figure(11); clf; hold on;
    plot(r1(1,:), r1(2,:), 'b.', 'MarkerSize', msize);
    plot(r2(1,:), r2(2,:), 'rx', 'MarkerSize', msize);
    plot(r3(1,:), r3(2,:), 'gs', 'MarkerSize', msize);
    axis equal
end

r1 = r1';
r2 = r2';
r3 = r3';