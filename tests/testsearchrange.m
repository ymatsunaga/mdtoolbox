%% usage:
%% >> test = testsearchrange;
%% >> result = run(test)
%%

classdef testsearchrange < matlab.unittest.TestCase

  methods(Test)
      
    function CompareWithExhaustiveNoPBC(testCase)
      natom = 16;
      crd = 16*rand(1, floor((natom^3)/3)*3);
      tic;
      pair_expected = searchrange_exhaustive(crd, [8.0 8.0 8.0], 5.0);
      toc
      tic;
      pair_actual = searchrange(crd, [8.0 8.0 8.0], 5.0);
      toc
      [~, id] = sort(pair_actual(:, 2));
      pair_actual = pair_actual(id, :);
      testCase.verifyEqual(pair_actual, pair_expected);
    end

    function CompareWithExhaustivePBC(testCase)
      natom = 16;
      crd = 16*rand(1, floor((natom^3)/3)*3);
      tic;
      pair_expected = searchrange_exhaustive(crd, [1.0 1.0 1.0], 5.0, [16.0 16.0 16.0]);
      toc
      tic;
      pair_actual = searchrange(crd, [1.0 1.0 1.0], 5.0, [16.0 16.0 16.0]);
      toc
      [~, id] = sort(pair_actual(:, 2));
      pair_actual = pair_actual(id, :);
      testCase.verifyEqual(pair_actual, pair_expected);
    end

  end
  
end

