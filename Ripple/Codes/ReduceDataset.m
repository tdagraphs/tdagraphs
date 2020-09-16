%----------------------------------------------------------------------------
%TAD: Topological Anomaly Detection in Dynamic Multilayer Blockchain Networks
%----------------------------------------------------------------------------

% Given a dynamic network with thousands of nodes, 
% this Matlab's code reduce network's size via maximum weight subgraph approximation.

%% To list files in the folder %%
nameFolder = 'PathToDataset' ;
lisFiles = dir([nameFolder '/*.txt']);
NFiles = length(lisFiles);
fprintf('Total of Files: %i  \n', NFiles);
NSetNodes = 90; 

%% Reduction %%
for kF = 1:NFiles
%for kF = 16:NFiles
    %kF = 1;
    fprintf(' Working on File number %i (out of %i), name: %s\n', kF, NFiles, lisFiles(kF).name);
    dataFile = dlmread([nameFolder '/' lisFiles(kF).name]); % Open File
    dataNoBug = dataFile(dataFile(:,4)<1e30,:); % Dataset without bugs
    % ----- Links clustering-strongest links
    cpyData = dataNoBug;
    matDiverseEdges = [];
    cDivEdges = [];
    while(size(cpyData,1)>0)
        bnry = (cpyData(:,1)==cpyData(1,1))&(cpyData(:,2)==cpyData(1,2)); % Where to find this Edges (logic)
        indFind = find(bnry); % Where to find this Edges (indices) 
        matDiverseEdges = [matDiverseEdges; [cpyData(1,1:2) length(indFind)]]; % How many repetitions
        auxCell = cell(1); 
        auxCell{1} = cpyData(indFind,:); % Rows corresponding to the edge
        cDivEdges = [cDivEdges auxCell]; % Saving previous rows
        cpyData = cpyData(~bnry,:); % Deleting these rows
        %fprintf('Strength of edges: %i,   Size left: %i \n', length(indFind), size(cpyData,1))
    end % End While
    fprintf('Strongest links-clustering have been found... \n')
    
    %%---- To choose the strongest connections ----%%
     [~, iSort] = sort(matDiverseEdges(:,3), 'descend');
     dataTrim = [];
     iNodes = 0;   % Number of different nodes added to dataTrim
     iRow = 1; % Row-matrix added
    while(iNodes<NSetNodes) 
        if(matDiverseEdges(iSort(iRow),1)~=matDiverseEdges(iSort(iRow),2)) % From to different node...
          dataTrim = [dataTrim; cDivEdges{iSort(iRow)}]; 
          iNodes = length(unique([dataTrim(:,1); dataTrim(:,2)]));
        end
        %fprintf(' iNodes: %i    iRow:%i \n', iNodes, iRow);
        iRow = iRow + 1; 
    end % End While
    fprintf('Strongest connections have been added... \n')
    
    %%---- Sort by timestamp and save it ----%%
    [~, kSort] = sort(dataTrim(:,3)); 
    dataTrimSort = dataTrim(kSort,:); 
    nameFolderOut = 'ReducedDataset';
    fileID = fopen([nameFolderOut '/TRIM_ORIGINAL/' lisFiles(kF).name(1:(end-4)) 'XXX.txt'], 'w');
    dataOut = dataTrimSort';
    %fprintf(fileID,'%d %d %d %-30.0f\n', dataOut); % but truncate it so it LOOKS like an int
    fprintf(fileID,'%d %d %d %-30.12f\n', dataOut); % but truncate it so it LOOKS like an int
    %dlmwrite([nameFolderOut '/' lisFiles(kF).name(1:(end-4)) 'XXX.txt'], dataTrimSort, 'delimiter', ' ', 'precision', '%ld')
    fprintf('File saved... \n')
    
    %%----- Save re-indexing -----%%
    lisNodes = unique([dataTrimSort(:,1); dataTrimSort(:,2)]);
    newIndexTrimSort = dataTrimSort;
    for i = 1:length(lisNodes)
        newIndexTrimSort(newIndexTrimSort(:,1)==lisNodes(i),1) = i;
        newIndexTrimSort(newIndexTrimSort(:,2)==lisNodes(i),2) = i;
    end
    fileID2 = fopen([nameFolderOut '/RI/RI' lisFiles(kF).name(1:(end-4)) 'XXX.txt'], 'w');
    dataOut2 = newIndexTrimSort';
    %fprintf(fileID2,'%d %d %d %-30.0f\n', dataOut2); % but truncate it so it LOOKS like an int
    fprintf(fileID2,'%d %d %d %-30.12f\n', dataOut2); % but truncate it so it LOOKS like an int
    dlmwrite([nameFolderOut '/INDEX/INDEX' lisFiles(kF).name(1:(end-4)) 'XXX.txt'], [(1:length(lisNodes))' lisNodes], 'delimiter', ' ', 'precision', '%d')
    fprintf('Files saved (RI and INDEX)... \n')
    
end % End For


