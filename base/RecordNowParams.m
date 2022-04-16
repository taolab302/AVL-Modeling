function RecordNowParams(params, now)

Folder = ['Parameter History\'];
if ~exist(Folder, 'dir')
    mkdir(Folder);
end

params_table = cell2table([fieldnames(params) struct2cell(params)], 'VariableNames', {'Parameter', 'Value'});
writetable(params_table, [Folder now '.csv']);