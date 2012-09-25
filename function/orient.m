function [crd_new] = orient(crd)
%% orient 
% orient molecule using the principal axes of inertia
%
% function [crd_new] = orient(crd,index)
%
% input: crd (1 x natom*3)
%        index 
%        mass (1 x natom)
%
% output: crd_new (1 x natom*3)
%
% example:
% [crd,b] = readpdb('trap11.pdb',{'CA'});
% crd = orient(crd);
% writepdb('orient_ca.pdb',crd,b);
% 

%% preparation
natom3 = size(crd,2);
natom = natom3/3;
nstep = size(crd,1);
mass = ones(1,natom);
index = 1:natom;

crd_new = zeros(1,natom3);

%% subtract by the geometric center
crd(1,1:3:end) = crd(1,1:3:end) - mean(crd(1,1:3:end));
crd(1,2:3:end) = crd(1,2:3:end) - mean(crd(1,2:3:end));
crd(1,3:3:end) = crd(1,3:3:end) - mean(crd(1,3:3:end));

%% calculate the principal axis of inertia
x = crd(1,1:3:end);
y = crd(1,2:3:end);
z = crd(1,3:3:end);

I = zeros(3,3);

I(1,1) = sum(mass.*(y.^2 + z.^2));
I(2,2) = sum(mass.*(x.^2 + z.^2));
I(3,3) = sum(mass.*(x.^2 + y.^2));

I(1,2) = - sum(mass.*(x.*y));
I(2,1) = I(1,2);

I(1,3) = - sum(mass.*(x.*z));
I(3,1) = I(1,3);

I(2,3) = - sum(mass.*(y.*z));
I(3,2) = I(2,3);

I
[p,q,a] = svd(I); %q is already sorted by descending order
p_axis = a(:,end:-1:1); %z-axis has the largest inertia

% check reflection
if det(p_axis) < 0
  p_axis(:,1) = - p_axis(:,1);
end

%% project onto the principal axis of inertia
proj = [x' y' z']*p_axis;
crd_new(:,1:3:end) = proj(:,1)';
crd_new(:,2:3:end) = proj(:,2)';
crd_new(:,3:3:end) = proj(:,3)';

