%% usage:
%% >> test = testcalcpairlist;
%% >> result = run(test)
%%

classdef testcalcpairlist < matlab.unittest.TestCase

  methods(Test)
      
    function CompareWithExhaustiveNoPBC(testCase)
      natom = 16;
      crd = 16*rand(1, floor((natom^3)/3)*3);
      tic;
      pair_expected = calcpairlist_exhaustive(crd, 5.0);
      toc
      tic;
      pair_actual = calcpairlist(crd, 5.0);
      toc
      [~, id] = sort(pair_actual(:, 1));
      pair_actual = pair_actual(id, :);
      [~, id] = sort(pair_actual(:, 2));
      pair_actual = pair_actual(id, :);
      testCase.verifyEqual(pair_actual, pair_expected);
    end

    function CompareWithExhaustivePBC(testCase)
      natom = 16;
      crd = 16*rand(1, floor((natom^3)/3)*3);
      tic;
      pair_expected = calcpairlist_exhaustive(crd, 5.0, [16.0 16.0 16.0]);
      toc
      tic;
      pair_actual = calcpairlist(crd, 5.0, [16.0 16.0 16.0]);
      toc
      [~, id] = sort(pair_actual(:, 1));
      pair_actual = pair_actual(id, :);
      [~, id] = sort(pair_actual(:, 2));
      pair_actual = pair_actual(id, :);
      testCase.verifyEqual(pair_actual, pair_expected);
    end

  end
  
end


