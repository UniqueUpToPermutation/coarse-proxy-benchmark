% Loads results from .mat file
% Results are in a table, first column is epsilon value, second column is
% the RBObject
data = load('result_table_BIE_more_samples.mat'); % For BIE example
%data = load('result_table_RTE.mat'); % For RTE example

result_table = data.result_table;

% Get the RBObject for epsilon = 1E-06
rb_object = result_table{7, 2};

% Plot error vs. epsilon, and output latex table
computeResults(result_table);

% View reduced basis, parameter is figure id
rb_object.viewReducedBasis(2);

% View a random result from the sample space, parameter is figure id
rb_object.viewRandomResult(3);

% View a specific result from the sample space
% First parameter is the index of the corresponding parameter in the
% parameter space, second parameter is figure id
rb_object.viewResult(1, 4);