% GraphStreamChangingPointDetection.m

%% To use all of this in the Ethereum database
% ----- Open Folder name
%nameFolder = 'ETH100_EdgeList_A/A';
nameFolder = 'ETH100_EdgeList_W/W';
%nameFolder = 'ETH100_AGD_EdgeList_A/AGDA';
%nameFolder = 'ETH100_AGD_EdgeList_W/AGDW';

% ----- Open Summary
tableSummary = readtable('ETH100summary.txt');
NNetworks = size(tableSummary,1);

%% For every network
for(kNet = 1:NNetworks)

%% Name Files
nameNetwork = char(table2array(tableSummary(kNet, 6)));
% ----- Path
PATH = [nameFolder nameNetwork '*'];

%PATH = 'RMIT_A_EdgeList/RMitA*'; 
%PATH = 'RMIT_W_EdgeList/RMitW*'; 
%PATH = 'RMIT_AGD_A_EdgeList/AGDRMitA*'; 
%PATH = 'RMIT_AGD_W_EdgeList/AGDRMitW*'; 

%% Previous
START_TIME = tic;
close all;
ZERO_TH = 1e-9;
rng(5);

%% Loading Files
% To split name of file...
strToSplitAux = strsplit(PATH, '/');
strToSplit = strToSplitAux{2}(1:(end-1));
% To open all files
fid_list = ls(PATH);
fid_list = strsplit(fid_list,{'\t','\n',' '}); fid_list(end)=[]; fileList = [];
for index_file=1:length(fid_list)
  fileList(index_file).name = fid_list{index_file};
  A = strsplit(fid_list{index_file},{strToSplit, '.'});
  fileList(index_file).order = str2num(A{2});
end;
[A,A] = sort([fileList.order]);
fileList = fileList(A);
NFiles = length(fileList);
disp(NFiles);

%% Main procedure ...
TotalEu = zeros(NFiles,1);
TotalKL = zeros(NFiles,1);
nNode = 46000;
%nNode = 1000;
toc(START_TIME);
  
readTime = 0;
for index_file=1:NFiles
  timeBeforeRead = tic;
  fid = fopen([fileList(index_file).name]);
  if (fid<0)
    continue;
  else
    edges = fscanf(fid,'%d %d %lf',[3 inf]);
    fclose(fid);
    readTime = readTime+(toc(timeBeforeRead));
    if (index_file==1)
      NPair = 10;
      NGroup = 25;
      rand_buffer = zeros(NGroup*2,NPair);
      % to generate nodes
      % left is smaller than right
      rand_index = randi(size(edges,2),[1,NGroup*NPair]);
      rand_index = sort(rand_index);
      rand_buffer = edges(1:2,rand_index);
      candidateCoding = (rand_buffer(1,:)-1)*nNode+rand_buffer(2,:)-1;
      node_freq = zeros([NGroup*NPair,NFiles]);
      Dist = zeros(NFiles,1); % Dist(i): diff(snap_i,snap_(i-1)). Dist(1)=0
      DistKL = zeros(NGroup,NFiles);
      EdgeCoding = edges(1,:)*nNode+edges(2,:)-nNode-1; % start from zero

      % to sample
      [sharedVals, idxsIntoA] = intersect(EdgeCoding,candidateCoding);
      [sharedVals, idxsToCandi] = intersect(candidateCoding,EdgeCoding);
      node_freq(idxsToCandi,index_file) = edges(3,idxsIntoA);
    else
      if(~isempty(edges)) % Added....
          EdgeCoding = edges(1,:)*nNode+edges(2,:)-nNode-1; % start from zero
          [sharedVals, idxsIntoA] = intersect(EdgeCoding,candidateCoding);
          [sharedVals, idxsToCandi] = intersect(candidateCoding,EdgeCoding);
          node_freq(idxsToCandi,index_file) = edges(3,idxsIntoA);
          % data smooth
          node_freq(find(node_freq(:,index_file)>1-ZERO_TH),index_file) = 1-ZERO_TH;
          node_freq(find(node_freq(:,index_file)<ZERO_TH),index_file) = ZERO_TH;
          % compare 
          Dist(index_file) = sum((node_freq(:,index_file)-node_freq(:,index_file-1)).^2);
          DistKL(:,index_file) = KLForIndependentJointSequence(node_freq(:,index_file)',node_freq(:,index_file-1)',NGroup);
      
          % if the sample too sparse, re-track
          if (length(find(node_freq(:,index_file)<=ZERO_TH))>0.25*NGroup*NPair)
             % re-track and resample
             rand_index = randi(size(edges,2),[1,NGroup*NPair]);
             rand_index = sort(rand_index);
             rand_buffer = edges(1:2,rand_index);
             candidateCoding = (rand_buffer(1,:)-1)*nNode+rand_buffer(2,:)-1;
             % to resample
             [sharedVals, idxsIntoA] = intersect(EdgeCoding,candidateCoding);
             [sharedVals, idxsToCandi] = intersect(candidateCoding,EdgeCoding);
             node_freq(idxsToCandi,index_file) = edges(3,idxsIntoA);
          end; % end of if (length(find(node_freq(:,index_file)<=ZERO_TH))>0.25*NGroup*NPair)
      end; % Added
    end; % end of elseif (index_file==1)
  end;
end;
% all files read

TotalEu = TotalEu+Dist;
TotalKL = TotalKL+median(DistKL,1)';
disp('KL+Eu 46k');
toc(START_TIME); 
pureTime = toc(START_TIME)-readTime;
disp(['Pure computation time: ',num2str(pureTime)]);
%return;

%% To Plot
figure('units','normalized','position',[.0,.0,1.0,.4]);
TotalEu = TotalEu+30;
%TotalKL(31) = TotalKL(31)+20;
%TotalKL(301) = TotalKL(301)+20;
if (NFiles>100)
  [ax,p1,p2] = plotyy((1:NFiles),TotalEu',(1:NFiles),TotalKL','plot','plot','r-O','b-+','FontSize',20);
  p1.LineWidth = 2;
  p2.LineWidth = 2;
  p2.LineStyle = '--';
  %[ax,p1,p2] = plotyy((1:NFiles),TotalEu',(1:NFiles),TotalKL'-10,'plot','plot','r-O','b-+','FontSize',20);
else
  [ax,p1,p2] = plotyy((1:NFiles),TotalEu',(1:NFiles),TotalKL','stem','stem','r-O','b-+','FontSize',20);
  p1.LineStyle = '-';
  p1.Marker = 'O';
  p1.LineWidth = 4;
  p2.LineStyle = '--';
  p2.Marker = '*';
  p2.LineWidth = 4;
end;
set(gca,'FontSize',20);
set(ax(1),'xlim',[0,NFiles]);
set(ax(2),'xlim',[0,NFiles]);
title(sprintf(['Window-averaged snapshots difference with ',num2str(NPair*NGroup),' node-pairs']));
xlabel(ax(1),'Window index','FontSize',20);
ylabel(ax(1),'Euclidean Distance','FontSize',25);
ylabel(ax(2),'KL Divergence','FontSize',25);
set(ax,'FontSize',30);
% Threshold Creation
firstFile = 1; lastFile = NFiles;
%medEu = median(TotalEu(firstFile:lastFile));
%mrEu = TotalEu((firstFile+1):end)-TotalEu(firstFile:(end-1));
%mrEu = sum(abs(mrEu))/length(mrEu);
%EuTh = medEu+2.66*mrEu;
nBootstrapSize = 1e6;
EuBoot = datasample(TotalEu(firstFile:lastFile),nBootstrapSize,'Replace',true);
EuTh = quantile(EuBoot,(1-0.05/6));

%EuTh = max(TotalEu(3:end))*0.3;
%KLTh = max(TotalKL(3:end))*0.2;
KLBoot = datasample(TotalKL,nBootstrapSize,'Replace',true);
KLTh = quantile(KLBoot,(1-0.05/6));

hold(ax(1),'on');
plot(ax(1),[1,NFiles],[EuTh,EuTh],'b','LineWidth',4);
hold(ax(2),'on');
plot(ax(2),[1,NFiles],[KLTh,KLTh],'r','LineWidth',4);

legend({'Euclidean dissimilarity',['Euclidean threshold = ',num2str(EuTh)],'KL disimilarity',['KL threshold = ',num2str(KLTh)]},'FontSize',25);
toc(START_TIME);

xlabh = get(gca,'XLabel');

%% ----- The Results ------
foundKL = find(TotalKL>KLTh);
foundEU = find(TotalEu>EuTh);
fprintf('\nAnomalies (KL disimilarity), indices: ');  disp(foundKL)
fprintf('Anomalies (Euclidean dissimilarity), indices: ');  disp(foundEU)
% Save Results
foundKL = foundKL';   foundEU = foundEU';
fidKL = fopen(['IDXTimeResul/' strToSplitAux{1} '_' nameNetwork '_YWangKL.txt'], 'w');
fprintf(fidKL, '%d\n', foundKL);
fclose(fidKL);
fidEU = fopen(['IDXTimeResul/' strToSplitAux{1} '_' nameNetwork '_YWangEU.txt'], 'w');
fprintf(fidEU, '%d\n', foundEU);
fclose(fidEU);

% -------- END FOR NETWORKS
end