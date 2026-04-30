clear; close all;

%% Read CR6 TOA5 data files from Gill R3-50 and concatenate
dirin = 'H:\';
files = dir([dirin 'Gill_50Hz_*.dat']);

% Sort files by date
[~, sortIdx] = sort([files.datenum]);
files = files(sortIdx);

%%
start_date=datetime(2025,12,23,0,0,0);
end_date=datetime(2025,12,24,0,0,0);
%%
files_sameday=25:30;
allTables = cell(1, length(files_sameday));
fileDates = cell(1, length(files_sameday));
%%
for jj = 1:length(files_sameday)
    % Read the table
    tempTable = readtable([dirin files(files_sameday(jj)).name], 'NumHeaderLines', 4);

    % Assign variable names
    tempTable.Properties.VariableNames = {'time_str','cr6_count','gill_status_analog',...
        'gill_status_digital','ux','uy','uz','temp_gill','wspd','temp','rh'};

    % Convert time_str to datetime
    tempTable.time_datetime = datetime(tempTable.time_str, 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SS');

    % Get date range for this file
    startDate = datestr(min(tempTable.time_datetime), 'yyyy-mm-dd');
    endDate = datestr(max(tempTable.time_datetime), 'yyyy-mm-dd');

    fprintf('File %d: %s\n', jj, files(files_sameday(jj)).name);
    fprintf('  Start: %s\n', startDate);
    fprintf('  End: %s\n', endDate);
    fprintf('  Records: %d\n\n', height(tempTable));

    allTables{jj} = tempTable;
    fileDates{jj} = struct('filename', files(files_sameday(jj)).name, 'start', startDate, 'end', endDate);
end

% Combine all tables
combinedTable = vertcat(allTables{:});

%% Filter for one day only
combinedTable.time_datetime = datetime(combinedTable.time_str, 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SS');
oneday_idx = (combinedTable.time_datetime >= start_date) & ...
    (combinedTable.time_datetime < end_date);
filteredTable = combinedTable(oneday_idx, :);


% Keep only desired variables
filteredTable = filteredTable(:, {'time_str','ux','uy','uz','temp','rh'});

% Convert to structure
GILL = table2struct(filteredTable, 'ToScalar', true);

%% Save with date in filename
savename = [dirin '20251223_Gill_SOARS_Data.mat'];

save(savename, 'GILL')
fprintf('Data saved to: %s\n', savename);
fprintf('Total records saved: %d\n', length(GILL.time_str));