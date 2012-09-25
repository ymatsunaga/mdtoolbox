function rc = writevmdvector(filename,crd,pcs,scale,resolution,radius)
%% writevmdvector
% output tcl scripts for drawing vectors in VMD
% 
% function rc = writevmdvector(filename,crd,pcs)
% 
%

if nargin == 3
% for bpti ca
%  scale = 15.0;
%  resolution = 100;
%  radius = 0.08;

% for trap
  scale = 170.0;
  resolution = 100;
  radius = 0.2;

% for alanine
%  scale = 3.0;
%  resolution = 100;
%  radius = 0.1;

% for bpti all-atom
%  scale = 40.0;
%  resolution = 100;
%  radius = 0.06;
end

nmode = size(pcs,2);
for i=1:nmode
  outfilename = [filename,int2str(i),'.tcl'];
  fid = fopen(outfilename,'w');
  
  fprintf(fid,'# define vector drawing function.\n');
  fprintf(fid,'# contrary to the example from the users guide, the vector is\n');
  fprintf(fid,'# defined by its origin, a vector at this point, and a scaling factor.\n');
  fprintf(fid,'# the function returns a list of the graphics ids for easy deletion.\n');
  fprintf(fid,'proc vmd_draw_vector {mol cnt vec {scale 1.0} {res 10} {radius 0.2}} {\n');
  fprintf(fid,'    set vechalf [vecscale [expr $scale *0.5] $vec]\n');
  fprintf(fid,'    return [list \\\n');
  fprintf(fid,'      [graphics $mol cylinder [vecsub $cnt $vechalf] \\\n');
  fprintf(fid,'        [vecadd $cnt [vecscale 0.7 $vechalf]] \\\n');
  fprintf(fid,'        radius $radius resolution $res filled yes] \\\n');
  fprintf(fid,'      [graphics $mol cone [vecadd $cnt [vecscale 0.7 $vechalf]] \\\n');
  fprintf(fid,'        [vecadd $cnt $vechalf] radius [expr $radius * 1.7] \\\n');
  fprintf(fid,'        resolution $res]]\n');
  fprintf(fid,'}\n');
  fprintf(fid,'# alternative version of a vector drawing subroutine\n');
  fprintf(fid,'# give basepoint and direction\n');
  fprintf(fid,'proc vmd_draw_vector2 { mol pos val {scale 1.0} {res 10} {radius 0.2}} {\n');
  fprintf(fid,'    set cnt  [ vecadd $pos [ vecscale [expr 0.5 * $scale] $val ] ]\n');
  fprintf(fid,'    set rv [vmd_draw_vector $mol $cnt $val $scale $res $radius]\n');
  fprintf(fid,'    return $rv\n');
  fprintf(fid,'}\n');
  fprintf(fid,'\n');
  
  for j=1:(size(pcs,1)/3)
    fprintf(fid,'draw vector2 {%f %f %f} {%f %f %f} %f %d %f\n', ...
            crd(1,3*(j-1)+1),crd(1,3*(j-1)+2),crd(1,3*(j-1)+3), ...
            pcs(3*(j-1)+1,i),pcs(3*(j-1)+2,i),pcs(3*(j-1)+3,i), ...
            scale, resolution, radius);
  end
  fclose(fid);
end

rc = 0;
%draw vector2 {1 1 1} {-0.5 0.5 0.5} 4.0 100 0.05

