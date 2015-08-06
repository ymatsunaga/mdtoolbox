function kai = calcstructuraloverlap(ref, trj, cutoff)
% calcstructuraloverlap
% calc structural overlap between ref and crd
% 
% example: 
% ref = readpdb('bln.pdb');
% trj = readdcd('bln.dcd');
% ref = ref./3.4;
% trj = trj./3.4;
% kai = calcstructuraloverlap(ref, trj, 0.2);
% 

%% setup
[nframe, natom3] = size(trj);
natom = natom3./3;
index = find(triu(true(natom), 3));

dmat1 = calcdistancematrix(ref);
dvec1 = dmat1(index);

%% calculation
kai = zeros(nframe, 1);

for iframe = 1:nframe
  dmat2 = calcdistancematrix(trj(iframe, :));
  dvec2 = dmat2(index);
  kai(iframe) = sum(abs(dvec1 - dvec2) < cutoff);
end

kai = 1 - kai * (2 ./ (natom.^2 - 5*natom + 6));







