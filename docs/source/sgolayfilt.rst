.. sgolayfilt
.. highlight:: matlab

==============================================================
Savitzky-Golay Filter (sgolayfilt)
==============================================================

::

 f = 1;
 x = 0:(f/100):f*5;
 plot(x,sin(2*pi*f*x)+rand(size(x)))
 hold on
 plot(x,sgolayfilt(sin(2*pi*f*x)+rand(size(x)),2,21),'r-','LineWidth',2)



