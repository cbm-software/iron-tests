clc;
clear all;
close all;
format longg

% conductivity tensor in material coordinates
tensor_2D = [2.0 0.0; ...
             0.0 3.0];
tensor_3D = [2.0 0.0 0; ...
             0.0 3.0 0.0; ...
             0.0 0.0 7.0];

% rotation angles (degree, radian) about x-, y-, z-axes
ra  = [30.0, 40.0, 10.0];
rar = ra / 180.0 * pi;

% rotation matrices
% like OpenCMISS, rotate xy-plane counterclockwise about first angle,
% see also: https://en.wikipedia.org/wiki/Rotation_matrix
rm_2D   = [cos(rar(1)), -sin(rar(1)); ...
           sin(rar(1)),  cos(rar(1))];
% like OpenCMISS, we rotate about z-axis, about rotated y-axis, about
% rotated x-axis, see also:
% http://danceswithcode.net/engineeringnotes/rotations_in_3d/rotations_in_3d_part1.html
rm_3D   = [cos(rar(2))*cos(rar(1)), ...
           sin(rar(3))*sin(rar(2))*cos(rar(1))-cos(rar(3))*sin(rar(1)), ...
           sin(rar(3))*sin(rar(1))+cos(rar(3))*sin(rar(2))*cos(rar(1)); ...
           cos(rar(2))*sin(rar(1)), ...
           cos(rar(3))*cos(rar(1))+sin(rar(3))*sin(rar(2))*sin(rar(1)), ...
           cos(rar(3))*sin(rar(2))*sin(rar(1))-sin(rar(3))*cos(rar(1)); ...
           -sin(rar(2)), ...
           sin(rar(3))*cos(rar(2)), ...
           cos(rar(3))*cos(rar(2))]

% rotated conductivity tensor for use in CHeart for deriving reference data
rot_2D = rm_2D' * tensor_2D * rm_2D
rot_3D = rm_3D' * tensor_3D * rm_3D

% print rotated tensors row-major order
rot_2D(:)'
rot_3D(:)'

% check that rotation matrices are orthonormal
fprintf(sprintf('det(rm_2D) = %.g\n', det(rm_2D)));
fprintf(sprintf('det(rm_3D) = %.g\n', det(rm_3D)));
% rm_2D' * rm_2D
% rm_2D * rm_2D'
% rm_3D' * rm_3D
% rm_3D * rm_3D'