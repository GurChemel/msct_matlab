function [Result] = Import_Loss_From_Csv(filename)

dataLines = [2, Inf];

%% Setup the Import Options
opts = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Walltime", "Step", "Value"];
opts.VariableTypes = ["double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
tbl = readtable(filename, opts);

%% Convert to output type
% Result.Walltime = tbl.Walltime;
% Result.Step = tbl.Step;
% Result.Value = tbl.Value;
Result = tbl.Value;
end