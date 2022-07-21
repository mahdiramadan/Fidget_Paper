% code written by Mahdi Ramadan to process a dataset of mouse visual cortex fluoresence data that can be downloaded 
% from this link: https://datadryad.org/stash/share/S-QAenoUCDeDypT64bgG6Vxz6HE0bQAFPomMx6eUW_k

% this code generates 1 plot: for each session, the average neuronal response to fidget behavior is plotted for all neurons. Fidget events occurr at index 100, and the sample
% rate is 100 indices/second. 

clear all; close all; clc;


% change to data directory
myDir = '/Users/mahdiramadan/Documents/Allen Analysis/Allen_Institute_Data/'; %gets directory
myFiles = dir(fullfile(myDir,'*.mat')); %gets all wav files in struct
all_lims = zeros(1,length(myFiles));
th_order = zeros(1,length(myFiles));




for k = 1:length(myFiles)

  % load data per session
  baseFileName = myFiles(k).name;
  fullFileName = fullfile(myDir, baseFileName);
  fprintf(1, 'Now reading %s\n', fullFileName);
  struct = load(fullFileName);
  lims = fullFileName(end-12:end-4);

  % Session ID
  all_lims(k) = string(lims);

  data = struct.arr;


  average_cell_response = squeeze(mean(data,2));

  av_cell_resp_centered = (average_cell_response - mean(average_cell_response(:,1:100),2))./std(average_cell_response(:,1:100),[],2);

   
  % order neurons by mean response to fidget ( each array is -2 seonds
  % relative to a fidget event to + 4 seconds after fidget, 50
  % samples/second)

    order = zeros(size(av_cell_resp_centered,1),1);
    [m,I] = min(mean(av_cell_resp_centered,2));
    
    for i = 1:size(av_cell_resp_centered,1)
        val = dtw(av_cell_resp_centered(I,:), av_cell_resp_centered(i,:));
        order(i,1) = val;
    end
    
    av_cell_resp_centered = sortrows([av_cell_resp_centered, order], -301);
    av_cell_resp_centered = av_cell_resp_centered(:,1:end-1);

    figure()
    imagesc(av_cell_resp_centered, [ -2 2])
    colorbar()
    set(gca,'fontsize',16)
    title(strcat(' Sorted Average Neuron Response to Fidget in Time ', ' ',lims))
    xlabel(' time index')
    ylabel( ' Neuron Number ')
    
end



