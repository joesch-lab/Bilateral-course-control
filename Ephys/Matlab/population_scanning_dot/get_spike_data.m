function [spike_data, mean_spike] = get_spike_data(rawtrace)
    %rawtrace is an array of size m X n, where m is number of samples and n
    %is the length of samples
    spike_data = zeros(size(rawtrace));
    mean_spike = zeros(size(rawtrace,1), 100);
    for i=1:size(rawtrace, 1)
        [pks, locs, width, prominence] = findpeaks(rawtrace(i, :),1:size(rawtrace(i, :), 1), ...
            'MinPeakProminence',0.04, 'MinPeakHeight', -0.5, 'MaxPeakWidth', 100, ...
            'MinPeakWidth', 20, 'MinPeakDistance', 50);
        spike_data(i, locs) = 1;
    
        all_spikes = zeros(size(pks, 2),100);
        for j=1:size(pks, 2)
            if (size(rawtrace, 2)-50>locs(j)) && (locs(j)>50)
                all_spikes(j, :) = rawtrace(i, locs(j)-49:locs(j)+50);
            end
        end
        mean_spike(i, :) = mean(all_spikes, 1);
    %     plot(mean(all_spikes, 1));
    end
end
