function [r1 r2] = gen2Dsample(n1, n2)

drawfig = false;
% n1 = 25;
% n2 = 15;
% drawfig = true;
% msize = 10;
% 
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

if (drawfig)
    figure(11); clf; hold on;
    plot(r1(1,:), r1(2,:), 'b.', 'MarkerSize', msize);
    plot(r2(1,:), r2(2,:), 'rx', 'MarkerSize', msize);
    axis equal
end

r1 = r1';
r2 = r2';