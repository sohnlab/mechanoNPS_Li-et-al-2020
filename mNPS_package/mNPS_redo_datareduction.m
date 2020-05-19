[xls_filename, ~] = uigetfile('*.xlsx','Choose spreadsheet'); 

opts = detectImportOptions(xls_filename,'NumHeaderLines',0);
opts.DataRange = 'A12'; % remove first few columns

T = readtable(xls_filename,opts); % import spreadsheet with cell data only

cells_index = table2array(T(:,2)); % reduced index of cells
cells_index = (cells_index-200)*20; % corrected for shift and downsample

file_bounds = find(diff(cells_index) < 0); % find where .mat files end
file_bounds = [0, file_bounds, length(cells_index)]; % add begin, end
num_files = length(file_bounds)-1; % number of files
fprintf('    > %d files found in spreadsheet\n',num_files),

data_from_all_files = [];

for file_index = 1:num_files
    cells_subset = cells_index(file_bounds(file_index)+1:file_bounds(file_index+1)); % per-file cell index
    
    prompt_msg = sprintf('Choose file number:  %d',file_index);
    [filename, path] = uigetfile('*.mat',prompt_msg);
    load(filename); % get first file
    
    data_from_one_file = zeros(2,length(cells_subset)*20000);
    
    % take data according to cell index
    for i = 1:length(cells_subset)
        data_subset = data(:,cells_subset(i):cells_subset(i)+20000-1);
        data_cat_start = 20000*(i-1)+1;
        data_cat_end = 20000*i;
        data_from_one_file(:,data_cat_start:data_cat_end) = data_subset;
    end
    
    data_from_all_files = [data_from_all_files, data_from_one_file];
    
end

data = data_from_all_files;
save_filename = strcat(xls_filename(1:end-4), 'mat');
save(save_filename,'data');