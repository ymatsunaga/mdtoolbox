%% usage:
%% >> test = testsuperimpose;
%% >> result = run(test)
%%

classdef testsuperimpose < matlab.unittest.TestCase

  methods(Test)
      
    function RMSDOfSameCrdShouldBeZero(testCase)
      natom = 100;
      ref = randn(1, natom*3);
      crd = ref;
      rmsd = superimpose(ref, crd);
      testCase.verifyEqual(rmsd, 0, 'AbsTol', 10^(-7));
    end
    
    function RMSDOfSameCrdButRotatedShouldBeZero(testCase)
      natom = 100;
      ref = randn(1, natom*3);
      ref = decenter(ref);
      theta = rand(1);
      rot = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
      crd = zeros(1, natom*3);
      for i = 1:natom
        x = ref((3*(i-1)+1):(3*i))';
        x = rot * x;
        crd((3*(i-1)+1):(3*i)) = x';
      end
      rmsd = superimpose(ref, crd);
      testCase.verifyEqual(rmsd, 0, 'AbsTol', 10^(-7));
    end
    
    function CheckEmptyArgument(testCase)
      natom = 100;
      ref = randn(1, natom*3);
      crd = ref;
      rmsd = superimpose(ref, crd, [], []);
      testCase.verifyEqual(rmsd, 0, 'AbsTol', 10^(-7));
    end
    
    function RMSDShouldBeEqualToExplicitRMSD(testCase)
      natom = 100;
      ref = randn(1, natom*3);
      crd = randn(1, natom*3);
      [rmsd, crd] = superimpose(ref, crd);
      rmsd_explicit = calcrmsd(ref, crd);
      testCase.verifyEqual(rmsd, rmsd_explicit, 'AbsTol', 10^(-7));
    end
    
  end
  
end

