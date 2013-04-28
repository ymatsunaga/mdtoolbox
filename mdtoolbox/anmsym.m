function [pcs,f,pcs_cell,f_cell] = anmsym(crd,rc,nsym)
% function [pcs,f,B_xyz,B_sum,B_exp,ind] = anmsym(crd,rc,nsym)
% descriptions about the func
% 
% example:
% [crd,b] = readpdb('orient_ca.pdb');
% rc = 12.0;
% nsym = 11;
% [pcs,f,B_xyz,B_sum,B_exp,ind] = anmsym(crd,rc,nsym)
% 

x = crd(1,1:3:end);
y = crd(1,2:3:end);
z = crd(1,3:3:end);
natom = length(x);
natom3 = natom*3;
nmono  = natom/nsym;
nmono3 = nmono*3;

% make a distance matrix S
S = zeros(natom,natom);
for i = 1:natom
  x_diff = x - x(i);
  y_diff = y - y(i);
  z_diff = z - z(i);
  S(i,:) = sqrt(x_diff.^2 + y_diff.^2 + z_diff.^2);
  S(i,i) = realmax;
end

% make a kirchhoff matrix G
G = zeros(natom,natom);
for i = 1:natom
  G(i,:) = - (S(i,:) < rc);
  G(i,i) = - sum(G(i,:));
end

% make a Hessian matrix H from G
H = zeros(nmono3,natom3);
for i = 1:nmono
  x_diff = x - x(i);
  y_diff = y - y(i);
  z_diff = z - z(i);

  % off-diagonal block
  H(3*(i-1)+1,1:3:end) = G(i,:) .* x_diff .* x_diff ./ (S(i,:).^2); % dV/(dx_i dx_j)
  H(3*(i-1)+1,2:3:end) = G(i,:) .* x_diff .* y_diff ./ (S(i,:).^2); % dV/(dx_i dy_j)
  H(3*(i-1)+1,3:3:end) = G(i,:) .* x_diff .* z_diff ./ (S(i,:).^2); % dV/(dx_i dz_j)
  H(3*(i-1)+2,2:3:end) = G(i,:) .* y_diff .* y_diff ./ (S(i,:).^2); % dV/(dy_i dy_j)
  H(3*(i-1)+2,3:3:end) = G(i,:) .* y_diff .* z_diff ./ (S(i,:).^2); % dV/(dy_i dz_j)
  H(3*(i-1)+3,3:3:end) = G(i,:) .* z_diff .* z_diff ./ (S(i,:).^2); % dV/(dz_i dz_j)

  H(3*(i-1)+2,1:3:end) = H(3*(i-1)+1,2:3:end); % dV/(dy_i dx_j)
  H(3*(i-1)+3,1:3:end) = H(3*(i-1)+1,3:3:end); % dV/(dz_i dx_j)
  H(3*(i-1)+3,2:3:end) = H(3*(i-1)+2,3:3:end); % dV/(dz_i dy_j)

  % diagonal block
  H(3*(i-1)+1,3*(i-1)+1) = 0;
  H(3*(i-1)+1,3*(i-1)+2) = 0;
  H(3*(i-1)+1,3*(i-1)+3) = 0;
  H(3*(i-1)+1,3*(i-1)+1) = - sum(H(3*(i-1)+1,1:3:end));
  H(3*(i-1)+1,3*(i-1)+2) = - sum(H(3*(i-1)+1,2:3:end));
  H(3*(i-1)+1,3*(i-1)+3) = - sum(H(3*(i-1)+1,3:3:end));

  H(3*(i-1)+2,3*(i-1)+1) = 0;
  H(3*(i-1)+2,3*(i-1)+2) = 0;
  H(3*(i-1)+2,3*(i-1)+3) = 0;
  H(3*(i-1)+2,3*(i-1)+1) = - sum(H(3*(i-1)+2,1:3:end));
  H(3*(i-1)+2,3*(i-1)+2) = - sum(H(3*(i-1)+2,2:3:end));
  H(3*(i-1)+2,3*(i-1)+3) = - sum(H(3*(i-1)+2,3:3:end));

  H(3*(i-1)+3,3*(i-1)+1) = 0;
  H(3*(i-1)+3,3*(i-1)+2) = 0;
  H(3*(i-1)+3,3*(i-1)+3) = 0;
  H(3*(i-1)+3,3*(i-1)+1) = - sum(H(3*(i-1)+3,1:3:end));
  H(3*(i-1)+3,3*(i-1)+2) = - sum(H(3*(i-1)+3,2:3:end));
  H(3*(i-1)+3,3*(i-1)+3) = - sum(H(3*(i-1)+3,3:3:end));
end

% make a character table
T    = zeros(nsym,nsym);
w    = exp(2*pi*j./nsym);
for k = 1:nsym
  ww     = w^(k-1);
  T(k,:) = ww.^(0:(nsym-1));
end

% make a rotation matrix
t  = 2*pi/nsym;
r  = [cos(t) sin(t) 0; -sin(t) cos(t) 0; 0 0 1];
R  = zeros(nmono3,nmono3);
for i=1:nmono
  R((3*(i-1)+1):(3*(i-1)+3),(3*(i-1)+1):(3*(i-1)+3)) = r;
end

% 対称座標系からみたHessiansを計算する
V  = zeros(nmono3,nmono3,nsym);
for p=1:nsym % loop for irreps
  for l=1:nsym % loop for symmetry operations
    V(:,:,p) = V(:,:,p) + T(p,l)*H(1:nmono3,(nmono3*(l-1)+1):(nmono3*(l-1)+nmono3)) * R^(l-1);
  end
end

nfold = ceil(nsym/2);
for p = 1:nfold % loop for real irreps
  p
  [u,s,v] = svd(squeeze(V(:,:,p)));
  pcs{p} = zeros(natom3,nmono3);
  for l=1:nsym
    pcs{p}((nmono3*(l-1)+1):(nmono3*(l-1)+nmono3),:) = T(p,l)*R^(l-1)*u;
  end
  pcs{p} = pcs{p}./sqrt(nsym);
  f{p} = sqrt(diag(s));

  %順番を低振動からに並び替える
  pcs{p} = pcs{p}(:,end:-1:1);
  f{p} = f{p}(end:-1:1);

  if p == 1
    %z方向の並進と回転のfrequencyを取り除く
    pcs{p}(:,1:2) = [];
    f{p}(1:2) = [];
  else
    if p == 2
      %xy方向の並進と回転のfrequencyを取り除く
      pcs{p}(:,1:2) = [];
      f{p}(1:2)     = [];
    end
    pcs_tmp = pcs{p};
    pcs{p}    = 2*real(pcs_tmp)./sqrt(2);
    pcs{nsym-p+2} = 2*imag(pcs_tmp)./sqrt(2);
    f{nsym-p+2} = f{p};
  end
end

if nsym == 2
  p = 2;
  p
  [u,s,v] = svd(squeeze(V(:,:,p)));
  pcs{p} = zeros(natom3,nmono3);
  for l=1:nsym
    pcs{p}((nmono3*(l-1)+1):(nmono3*(l-1)+nmono3),:) = T(p,l)*R^(l-1)*u;
  end
  pcs{p} = real(pcs{p}./sqrt(nsym));
  f{p} = sqrt(diag(s));

  %順番を低振動からに並び替える
  pcs{p} = pcs{p}(:,end:-1:1);
  f{p} = f{p}(end:-1:1);

  %xy方向の並進と回転のfrequencyを取り除く
  pcs{p}(:,1:4) = [];
  f{p}(1:4)     = [];

elseif mod(nsym,2) == 0
  p = nfold+1;
  p
  [u,s,v] = svd(squeeze(V(:,:,p)));
  pcs{p} = zeros(natom3,nmono3);
  for l=1:nsym
    pcs{p}((nmono3*(l-1)+1):(nmono3*(l-1)+nmono3),:) = T(p,l)*R^(l-1)*u;
  end
  pcs{p} = real(pcs{p}./sqrt(nsym));
  f{p} = sqrt(diag(s));

  %順番を低振動からに並び替える
  pcs{p} = pcs{p}(:,end:-1:1);
  f{p} = f{p}(end:-1:1);
end

pcs_cell = pcs;
f_cell = f;

pcs = [];
f = [];
for i = 1:nsym
  pcs = [pcs pcs_cell{i}];
  f = [f; f_cell{i}];
end
[f,I] = sort(f);
pcs = pcs(:,I);

B_sum = 0;
B_xyz = 0;
B_exp = 0;
ind = 0;

