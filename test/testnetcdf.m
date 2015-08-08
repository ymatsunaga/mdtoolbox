%% usage:
%% >> test = testnetcdf;
%% >> result = run(test)
%%

classdef testnetcdf < matlab.unittest.TestCase

  methods(Test)
      
    function CheckConsistensyOfWriteReadOfTrj(testCase)
      natom = 1;
      nstep = 10;
      trj_expected = randn(nstep, natom*3);
      writenetcdf('tmp.nc', trj_expected);
      trj_actual = readnetcdf('tmp.nc');
      testCase.verifyEqual(trj_actual, trj_expected, 'RelTol', 10^(-7));
      delete 'tmp.nc'
    end
    
    function CheckConsistensyOfWriteReadOfBox(testCase)
      natom = 1;
      nstep = 10;
      trj_expected = randn(nstep, natom*3);
      box_expected = randn(nstep, 3);
      writenetcdf('tmp.nc', trj_expected, box_expected);
      [trj_actual, box_actual] = readnetcdf('tmp.nc');
      testCase.verifyEqual(box_actual, box_expected, 'RelTol', 10^(-14));
      delete 'tmp.nc'
    end
    
    function CheckSubsetReadByIndexAtom(testCase)
      natom = 3;
      nstep = 10;
      index_atom = [1 3];
      index_atom3 = to3(index_atom);
      trj_expected = randn(nstep, natom*3);
      writenetcdf('tmp.nc', trj_expected);
      trj_actual = readnetcdf('tmp.nc', index_atom);
      testCase.verifyEqual(trj_actual, trj_expected(:, index_atom3), 'RelTol', 10^(-7));
      delete 'tmp.nc'
    end
    
    function CheckSubsetReadByIndexTime(testCase)
      natom = 1;
      nstep = 10;
      index_time = 2:2:10;
      trj_expected = randn(nstep, natom*3);
      writenetcdf('tmp.nc', trj_expected);
      trj_actual = readnetcdf('tmp.nc', [], index_time);
      testCase.verifyEqual(trj_actual, trj_expected(index_time, :), 'RelTol', 10^(-7));
      delete 'tmp.nc'
    end
    
  end
  
end

