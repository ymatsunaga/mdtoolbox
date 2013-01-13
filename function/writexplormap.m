function rc = writexplormap(filename,map,box)
%% writexplormap
% write xplor density format file
%
% function rc = writexplormap(filename,map,box)
%
% input: filename: filename
%        map (na x nb x nc): density data
%        box (3): box size
%
% output: rc 
%
% reference: http://www.scripps.edu/rc/softwaredocs/msi/xplor981/formats.html
% 
% example:
% map = randn(10,10,10);
% box = [10 10 10];
% writexplormap('tmp.xplor',map,box);
% 

%% check existing file
if exist(filename, 'file')
  filename_old = sprintf('%s.old', filename);
  display(sprintf('existing file %s is moved to %s', filename, filename_old));
  movefile(filename, filename_old);
end

%% set parameters
ang = 90;
[na,nb,nc] = size(map);

% amin = -na/2;
% bmin = -nb/2;
% cmin = -nc/2;
% amax = (na/2)-1;
% bmax = (nb/2)-1;
% cmax = (nc/2)-1;

amin = 1;
bmin = 1;
cmin = 1;
amax = na;
bmax = nb;
cmax = nc;

map_ave = mean(map(:));
map_std = std(map(:));

fid = fopen(filename,'w');

%% write 
% write header
fprintf(fid,'\n');
fprintf(fid,'      1 !NTITLE\n');
for i = 3:5
  fprintf(fid,' REMARKS\n');
end
fprintf(fid,'%8d',na,amin,amax,nb,bmin,bmax,nc,cmin,cmax);
fprintf(fid,'\n');
fprintf(fid,'%12.5E',box(1),box(2),box(3),ang,ang,ang);
fprintf(fid,'\n');
%fprintf(fid,'ZXY\n');
fprintf(fid,'ZYX\n');

% write density array
for ic = 1:nc
  map_2d = squeeze(map(:,:,ic));
  map_unroll = map_2d(:);
  nlength = length(map_unroll);
  fprintf(fid,'%8d\n',ic);
  for i = 1:6:nlength
    fprintf(fid,'%12.5E',map_unroll(i:min(i+5,nlength)));
    fprintf(fid,'\n');
  end
end

% write footer
fprintf(fid,'%8d\n',-9999);
fprintf(fid,'%12.4E %12.4E ',map_ave,map_std);
fprintf(fid,'\n');

fclose(fid);

rc = 0;

%   write(6,*)
%   write(6,*) '      1 !NTITLE'
%   write(6,*) 'REMARKS'
%   write(6,'(9I8)') ng3(1),-ng3(1)/2,ng3(1)/2-1,&
%                   &ng3(2),-ng3(2)/2,ng3(2)/2-1,&
%                   &ng3(3),-ng3(3)/2,ng3(3)/2-1
%   write(6,'(6E12.5)') box(1),box(2),box(3),ang,ang,ang
%   write(6,'(A3)') xyz

%   do k=1,ng3(3)
%     write(6,'(I8)') k-1-ng3(3)/2
%     write(6,'(6E12.5)') ((guv(i,j,k,nv),i=1,ng3(1)),j=1,ng3(2))
%   end do

%   write(6,'(I8)') -9999
%   write(6,'(2(E12.4,1X))') 0.0,1.0

