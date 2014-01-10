%% usage:
%% >> test = testsimulatepoisson;
%% >> result = run(test)
%%

classdef testsimulatepoisson < matlab.unittest.TestCase

  methods(Test)
      
    function AveragedNumberOfEventsShouldBeTimewidthMutipliedByRate(testCase)
      time_width = 1;
      rate = 10000;
      nevent = [];
      for i = 1:100
        t = simulatepoisson(time_width, rate);
        nevent = [nevent; numel(t)];
      end
      nevent_expected = time_width*rate;
      nevent_actual = mean(nevent);
      testCase.verifyEqual(nevent_actual, nevent_expected, 'RelTol', 10^(-2));
    end
    
  end
  
end

