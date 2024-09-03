all_colors = load('D:\Roshan\Project\PythonCodes\Codes\Plotting\colorData.mat');
all_color_list = fieldnames(all_colors);

for i=1:size(all_color_list, 1)
    a=all_color_list(i);
    all_data.(a{:}) = othercolor(a{:});
end

filename = fullfile('D:\Roshan\Project\PythonCodes\Codes\Plotting', 'colorList.mat');
save(filename, 'all_data');