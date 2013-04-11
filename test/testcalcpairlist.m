%% usage:
%% >> test = testsimulatepoissonprocess;
%% >> result = run(test)
%%

classdef testsimulatepoissonprocess < matlab.unittest.TestCase

  methods(Test)
      
    function AveragedNumberOfEventsShouldBeRateMutipliedByTimewidth(testCase)
      time_width = 10000;
      rate = 1;
      nevent = [];
      for i = 1:100
        t = simulatepoissonprocess(time_width, rate);
        nevent = [nevent; numel(t)];
      end
      nevent_expected = time_width*rate;
      nevent_actual = mean(nevent);
      testCase.verifyEqual(nevent_actual, nevent_expected, 'AbsTol', 10^(-7));
    end
    
  end
  
end

