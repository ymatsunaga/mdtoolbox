.. matlab
.. highlight:: matlab

グラフの保存と画像への出力(exportasなど)
==============================================================

MATLABはグラフの取扱いに優れている。プロットしたグラフは
::

 saveas(gcf,'plot.fig')

で保存することができる。ここで、 ``gcf`` は ``get current figure`` の略。

一旦保存したグラフは、後ほど
::

 open('plot.fig')

でオープンして、再利用することができる。

グラフを画像へ出力するには以下のように実行する。
::

 set(gcf,'PaperPositionMode','auto'); %output the figure in WYSWYG way
 print('-deps','plot.eps'); %monochrome EPS
 print('-depsc','plot.eps'); % color EPS
 print('-dtiff','plot.tiff'); %compressed TIFF
 print('-dtiffnocompression','plot.tiff'); %non-compressedTIFF
 print('-dpng','plot.png'); %PNG
 print('-depsc2','-adobecset','plot.eps'); % EPS editable with Adobe Illustrator

上記の一連の操作(``saveas``, ``set``, ``eps+png+tiff`` への出力)を行う関数
``exportas.m`` を用意してある。

以下のように実行する。
::

 plot(randn(1000,1));
 xlabel('step','FontSize',50);
 ylabel('value','FontSize',50);
 exportas('myplot');

