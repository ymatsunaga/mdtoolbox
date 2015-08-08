%% usage:
%% >> test = testdcd;
%% >> result = run(test)
%%

classdef testdcd < matlab.unittest.TestCase

  methods(Test)
      
    function CheckConsistensyOfWriteReadOfTrj(testCase)
      natom = 1;
      nstep = 10;
      trj_expected = randn(nstep, natom*3);
      writedcd('tmp.dcd', trj_expected);
      trj_actual = readdcd('tmp.dcd');
      testCase.verifyEqual(trj_actual, trj_expected, 'RelTol', 10^(-7));
      delete 'tmp.dcd'
    end
    
    function CheckConsistensyOfWriteReadOfBox(testCase)
      natom = 1;
      nstep = 10;
      trj_expected = randn(nstep, natom*3);
      box_expected = randn(nstep, 3);
      writedcd('tmp.dcd', trj_expected, box_expected);
      [trj_actual, box_actual] = readdcd('tmp.dcd');
      testCase.verifyEqual(box_actual, box_expected, 'RelTol', 10^(-14));
      delete 'tmp.dcd'
    end
    
    function CheckSubsetReadByIndexAtom(testCase)
      natom = 3;
      nstep = 10;
      index_atom = [1 3];
      index_atom3 = to3(index_atom);
      trj_expected = randn(nstep, natom*3);
      writedcd('tmp.dcd', trj_expected);
      trj_actual = readdcd('tmp.dcd', index_atom);
      testCase.verifyEqual(trj_actual, trj_expected(:, index_atom3), 'RelTol', 10^(-7));
      delete 'tmp.dcd'
    end
    
  end
  
end

