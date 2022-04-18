clear all, close all
clc
tx = cos(30*pi/180)
ty = sin(30*pi/180)

plot([0 6*tx], [0 6*ty], '--k')
hold on, grid on

u1 = 10*cos(40*pi/180)
v1 = 10*sin(40*pi/180)
plot([0 u1], [0 v1], 'r')

u0 = u1*(tx^2 -ty^2) + 2*v1*tx*ty
v0 = 2*u1*tx*ty + v1 * (ty^2 -tx^2)
% u0 = -u1
% v0 = -v1
plot([0 u0], [0 v0], 'b')
axis equal

V1 = sqrt(u1^2+v1^2)
V0 = sqrt(u0^2+v0^2)
