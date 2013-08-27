function d = periodic(x, center)
%
%

d = abs(x - center);
while d > 180
  d = d - 360;
end

