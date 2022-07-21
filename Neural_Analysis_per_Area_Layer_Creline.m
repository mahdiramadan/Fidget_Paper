% code written by Mahdi Ramadan to process a dataset of mouse visual cortex fluoresence data that can be downloaded 
% from this link: https://datadryad.org/stash/share/S-QAenoUCDeDypT64bgG6Vxz6HE0bQAFPomMx6eUW_k

% this code generates 3 plots: the first plot shows the distribution of neuronal cluster types across either cortical area, layer or cre-line. The user can switch between
% the 3 by setting the variable u equal to the strings associated with each category ( lines 39-43). The next plot shows the percentage of neurons that have 
% a fidget induced modulation that exceeds 2 stds from baseline for each category (areas, layers or cre-lines. The final plot shows the percentage of al neurons that have 
% a fidget induced modulation that exceeds 2 stds from baseline. 

clear all; close all; clc;

% change to data directory
get = load('/Users/mahdiramadan/Documents/Allen Analysis/Allen_table.mat');
table = get.data;

myDir = '/Users/mahdiramadan/Documents/Allen Analysis/Allen_Institute_Data'; %gets directory
myFiles = dir(fullfile(myDir,'*.mat')); %gets all wav files in struct

% initialize arrays
bin_size = 5;
time_series = ones(length(myFiles), 300);
raw = ones(length(myFiles), 300);
av_pre = ones(length(myFiles), 1);
av_post = ones(length(myFiles), 1);
zscore_2 = ones(length(myFiles), 300/bin_size);
trial_av = ones(length(myFiles), 1);
cluster = load('Cluster_Values.mat');

zscore_all = [];

depressed = ones(length(myFiles), 1);
bipolar = ones(length(myFiles), 1);
neutral = ones(length(myFiles), 1);
active = ones(length(myFiles), 1);
mis_class = ones(length(myFiles), 1);


% PICK either area,layer or creline for analysis by uncommenting

% u = {'Cux2', 'Emx1', 'Rbp4','Rorb','Scnn1a'};

% u = {'175', '275','350','375'};

u = {'VISl', 'VISpm','VISal','VISp'};

cell_type = struct;
d= ones(length(u),5);

for y= 1:length(u)

    c = 0;  

for k = 1:length(myFiles)
  
  % load data   
  % ( each array is -2 seonds
  % relative to a fidget event to + 4 seconds after fidget, 50
  % samples/second)

  baseFileName = myFiles(k).name;
  fullFileName = fullfile(myDir, baseFileName);
  fprintf(1, 'Now reading %s\n', fullFileName);
  struct = load(fullFileName);
  lims = fullFileName(end-12:end-4);


log = string(table(:,2)) == lims; 
in = find(log);
interest = table(in,6);


if string(interest) == string(u(y))

data = struct.arr;
    
average_cell_response = squeeze(mean(data,2));

av_cell_resp_centered = (average_cell_response - mean(average_cell_response(:,1:100),2))./std(average_cell_response(:,1:100),[],2);


% compute mean response per neuron during fidget

order = zeros(size(av_cell_resp_centered,1),1);
[m,I] = min(mean(av_cell_resp_centered,2));

for i = 1:size(av_cell_resp_centered,1)
    val = dtw(av_cell_resp_centered(I,:), av_cell_resp_centered(i,:));
    order(i,1) = val;
end

av_cell_resp_centered = sortrows([av_cell_resp_centered, order], -301);
av_cell_resp_centered = av_cell_resp_centered(:,1:end-1);

av_time = mean(abs(av_cell_resp_centered),1);
raw(c+1,:) = mean(av_cell_resp_centered,1);
time_series(c+1,:) = av_time;


% compute histogram for when neuronal firing is zscore > 2 
am = (average_cell_response - mean(average_cell_response(:,1:100),2))./std(average_cell_response(:,1:100),[],2);
counts = [];
for ib = 1 : bin_size: size(am,2) - bin_size + 1
    counts = horzcat(counts,sum(any(abs(am(:,ib: ib + bin_size-1)) > 2, 2)));
end
zscore_2(c+1,:) = (counts./size(av_cell_resp_centered,1))*100;

  
    
[IDX, C, SUMD] = kmeans(av_cell_resp_centered(:,101:300), 4, 'start', cluster.C(:,101:300));

cell_type.lims{k} =  lims;

cell_type.type{k} = IDX;

% seperate out neutral neurons that have large activations as mislabeled
% neurons

IDX((IDX == 2 | IDX == 1 | IDX == 4) & (max(abs(av_cell_resp_centered(:,101:300)),[],2) < 2)) = 0;

IDX((IDX == 3) & (max(abs(av_cell_resp_centered(:,101:300)),[],2) > 2)) = 5;
 
depressed(c+1,1) = (sum(IDX==1) / length(IDX(IDX ~= 0)))*100;
bipolar(c+1,1) = (sum(IDX==2) / length(IDX(IDX ~= 0)))*100;
neutral(c+1,1) = (sum(IDX==3) / length(IDX(IDX ~= 0)))*100;
active(c+1,1) = (sum(IDX==4) / length(IDX(IDX ~= 0)))*100;
mis_class(c+1,1) = (sum(IDX==5) / length(IDX(IDX ~= 0)))*100;

c = c + 1;


end

end



% plots distribution of clutser types per area,layer or cre 
figure()

ac = histc(active(1:c,1), [0:10:100]);

bi = histc(bipolar(1:c,1), [0:10:100]);

de = histc(depressed(1:c,1), [0:10:100]);

ne = histc(neutral(1:c,1), [0:10:100]);

bar3(0:10:100, [100*ac./(sum(ac)) 100*bi./(sum(bi)) 100*de./(sum(de)) 100*ne./(sum(ne))])

title( strcat(' Neuronal Reponse Distribution for - ', string(u(y))))
xlabel( 'Type')
ylabel( ' Percent of neurons ')
zlabel('Percent of sessions')
set(gca,'fontsize',16)
legend ('active', 'bipolar', 'depressed', 'neutral')


zscore_2 = mean(zscore_2(1:c,:));
zscore_all = vertcat(zscore_all, zscore_2);

% percent of neurons with zscore > 2 over time per area,layer or cre

figure()
bar(zscore_2, 'FaceColor', 'k')
title(strcat(' Neuron Percentage with abs(zscore) > 2 -',string(u(y))))
xlabel(' Time index')
ylabel( ' Percentage ')
ylim([0 60])
set(gca,'fontsize',16)


end

% percent of neurons with zscore > 2 over time for all conditions

figure()
bar(mean(zscore_all), 'FaceColor', 'k')
title(' Neuron Percentage with abs(zscore) > 2')
xlabel(' Time index')
ylabel( ' Percentage ')
ylim([0 60])
set(gca,'fontsize',16)


