%% this script carries out all unit tests in this folder
%% usage:
%% >> cd test/
%% >> testall
%% 

% import TestSuite Class
import matlab.unittest.TestSuite;

% create Suite from all test case files in current folder
suiteFolder = TestSuite.fromFolder(pwd);
result = run(suiteFolder)

