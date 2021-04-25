function [X, Y] = importAirfoil(filename, dataLines)
%IMPORTFILE Import data from a text file
%  [X, Y] = IMPORTFILE(FILENAME) reads data from text file FILENAME for
%  the default selection.  Returns the data as column vectors.
%
%  [X, Y] = IMPORTFILE(FILE, DATALINES) reads data for the specified row
%  interval(s) of text file FILENAME. Specify DATALINES as a positive
%  scalar integer or a N-by-2 array of positive scalar integers for
%  dis-contiguous row intervals.
%
%  Example:
%  [X, Y] = importfile("C:\Users\Slidgster\Desktop\School\Research\Code\morphing-af-calc\S826.dat", [1, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 04-Mar-2021 12:47:03

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [1, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = ["X", "Y", "Var3"];
opts.SelectedVariableNames = ["X", "Y"];
opts.VariableTypes = ["double", "double", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% Specify variable properties
opts = setvaropts(opts, "Var3", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Var3", "EmptyFieldRule", "auto");
opts = setvaropts(opts, "X", "TrimNonNumeric", true);
opts = setvaropts(opts, "X", "ThousandsSeparator", ",");

% Import the data
tbl = readtable(filename, opts);

%% Convert to output type
X = tbl.X;
Y = tbl.Y;
end