clear all; close all;

%Read CR6 TOA5 data files from Gill R3-50 and concatenate
dirin = 'C:\Campbellsci\LoggerNet\';
files = dir([dirin '*.dat']);
allTables = cell(1, length(files));

for jj = 1:length(files)
    allTables{jj} = readtable([dirin files(jj).name], 'NumHeaderLines', 4);
end

combinedTable = vertcat(allTables{:});
combinedTable.Properties.VariableNames = {'time_str','cr6_count','gill_count','gill_code','ux','uy','uz','temp','wspd'};

%GILL = table2struct(combinedTable);
GILL = table2struct(combinedTable,'ToScalar',true)

save([dirin '20250923_Gill_SOARS_Data.mat'],'GILL')