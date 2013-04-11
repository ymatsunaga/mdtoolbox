%% usage:
%% >> test = testsimulatenonhomopoisson;
%% >> result = run(test)
%%

classdef testsimulatenonhomopoisson < matlab.unittest.TestCase

  methods(Test)
      
    function AveragedNumberOfEventsShouldBeTimewidthMutipliedByRate(testCase)
      time_width = 2;
      rate = 10000;
      nevent = [];
      for i = 1:100
        [t_event, t] = simulatenonhomopoisson(time_width, [100*ones(1,1000) rate*ones(1,1000)]);
        nevent = [nevent; numel(t(t_event>1))];
      end
      nevent_expected = (time_width/2)*rate;
      nevent_actual = mean(nevent);
      testCase.verifyEqual(nevent_actual, nevent_expected, 'RelTol', 10^(-2));
    end
    
  end
  
end

