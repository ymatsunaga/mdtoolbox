function [trj,box] = readfloattrj(natom, filename, index, isbox)
%% readfloattrj
% read float trajectory file
%
% function [trj, box, ititle] = readambertrj(natom, filename, index_atom, isbox)
%
% input: natom 原子数
%        filename ファイル名
%        index_atom 読み込む原子番号のリスト(省略可。省略された時には全原子を読む)
%        isbox ボックス付きのtrjファイルか否か。ブーリアン型(false or true)。(省略可。省略された時はボックスなし=false)
%
% output: trj (nframe x natom*3) トラジェクトリ each row containing coordinates in the order [x1 y1 z1 x2 y2 z2 ...]
%         box (nframe x 3) box
% 
% example:
% ボックス無しの場合
% natom = 3343;
% trj = readfloattrj(natom,'md.trj');
% ボックス有りの場合
% natom = 62475;
% [trj,box] = readfloattrj(natom,'md_with_box.trj',1:natom,true);
% 

natom3 = natom*3;

if ~exist('index', 'var') || isempty(index)
  index = 1:natom;
end

if ~exist('index', 'var')
  isbox = false;
end

if isbox
  nlimit = natom3 + 3;
else
  nlimit = natom3;
end

index3 = to3(index);

fid = fopen(filename, 'r');
trj = fread(fid, '*float');
fclose(fid);

nframe = length(trj) / nlimit;
trj = reshape(trj,nlimit,nframe)';

if isbox
  box = trj(:, end-2:end);
  trj(:, end-2:end) = [];
end

trj = trj(:, index3);

