function mywritetable( data, cellName, outputPoisition)
%MYWRITETABLE write a floating point data to table
%	mywritetable( data, cellName, outputPoisition) output
%	the matrix data to a table whose title will be
%	the elements in cellName.
%	outputPosition specify the saved location
%   
%   data : should aligned in columns
%   cellName: A vector of strings contains the names 
%       ["name1","name2",...,"namen"]
%   outputPosition: rel or abs position are all acceptable.

G = double(data);
n = length(cellName);
path=outputPoisition;

G = array2table(G);
G.Properties.VariableNames(1:n) = cellName;
writetable(G,path);

end
