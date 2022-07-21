% code written by Mahdi Ramadan to process a dataset of mouse visual cortex fluoresence data that can be downloaded 
% from this link: https://datadryad.org/stash/share/S-QAenoUCDeDypT64bgG6Vxz6HE0bQAFPomMx6eUW_k

% this code generates 8 plots: 4 plots show the average neural response time aligened to fidget behaviors for each neuronal cluster type per session. Fidget events occurr 
% at index 100, and the frame rate is 50 indices/second. The next 4 plots show the average response of all neurons from all sessions to fidget for each cluster type. 



clear all; close all; clc;

% change to data directory
get = load('/Users/mahdiramadan/Documents/Allen Analysis/Allen_table.mat');
table = get.data;

myDir = '/Users/mahdiramadan/Documents/Allen Analysis/Allen_Institute_Data/'; %gets directory
myFiles = dir(fullfile(myDir,'*.mat')); 

all_neurons_1 = [];
all_neurons_2 = [];
all_neurons_3 = [];
all_neurons_4 = [];
    
for k = 1:length(myFiles)
  
  % load data   
  % ( each array is -2 seonds
  % relative to a fidget event to + 4 seconds after fidget, 50
  % samples/second)

  baseFileName = myFiles(k).name;
  fullFileName = fullfile(myDir, baseFileName);
  fprintf(1, 'Now reading %s\n', fullFileName);
  struct = load(fullFileName);
  % Session ID
  lims = fullFileName(end-12:end-4);


    
    data = struct.arr;
    
    cluster = load('Cluster_Values.mat');
    
    average_cell_response = squeeze(mean(data,2));
    
    av_cell_resp_centered = (average_cell_response - mean(average_cell_response(:,1:100),2))./std(average_cell_response(:,1:100),[],2);


    % Cluster Neurons 

    [IDX, C, SUMD] = kmeans(av_cell_resp_centered(:,101:300), 4, 'start', cluster.C(:,101:300));


    % Cluster IDS, 1 = depressed, 2 = phasic, 3 = neutral, 4 = active
    all_neurons_1 = vertcat( all_neurons_1, av_cell_resp_centered(IDX == 1,:));
    all_neurons_2 = vertcat( all_neurons_2, av_cell_resp_centered(IDX == 2,:));
    all_neurons_3 = vertcat( all_neurons_3, av_cell_resp_centered(IDX == 3,:));
    all_neurons_4 = vertcat( all_neurons_4, av_cell_resp_centered(IDX == 4,:));

    %  plot each cluster of neurons per sessions
    figure(1)
    imagesc( av_cell_resp_centered(IDX == 1,:) , [ -2 2])
    figure(2)
    imagesc( av_cell_resp_centered(IDX == 2,:) , [ -2 2])
    figure(3)
    imagesc( av_cell_resp_centered(IDX == 3,:) , [ -2 2])
    figure(4)
    imagesc( av_cell_resp_centered(IDX == 4,:) , [ -2 2])

    xlabel(' time index')
    ylabel( ' Neuron Number ')


end


% plot all neurons from all sessions from each cluster

figure()
imagesc(all_neurons_1(randperm(size(all_neurons_1, 1)), :), [ -2 2.5])
hold on

figure(2)
imagesc( all_neurons_2(randperm(size(all_neurons_2, 1)), :) , [ -2 2.5])
hold on
% colorbar()
figure(3)
imagesc( all_neurons_3(randperm(size(all_neurons_3, 1)), :), [  -2 2.5])
hold on
% colorbar()
figure(4)
imagesc(all_neurons_4(randperm(size(all_neurons_4, 1)), :), [  -2 2.5])
hold on
% colorbar()



