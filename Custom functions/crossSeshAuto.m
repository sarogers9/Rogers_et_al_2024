

function [longRegistered, coordinates] = crossSeshAuto(L,dA,locs,nSesh)
%this function aligns longitudinally registered datasets
%INPUTS:    L - the cell registration file from IDPS
%           d - a cell containing the cell trace files from each session in
%               rows
% 
%           nSesh - the number of sessions
%OUTPUTS:   longRegistered = a cell of matrices of longitudinally registered cell
%               traces from each session
%
%
%
%Written by Sophie A. Rogers, Corder Laboratory, University of Pennsylvania
%%  
    
    %extract information from longitudinal registration doc
    glo_ind = L(:,1);
    loc_ind = L(:,2);
    sesh = L(:,3);
    
    if ~isempty(locs)
        posX = table2array(locs(:,6));
        posY = table2array(locs(:,7));
        diameter = table2array(locs(:,9));
    end
    %create a list of local IDs of cells registered in all sessions
    cellIDs = cell(nSesh,1);

    vec = unique(glo_ind);
    cellIDs = [];
    for n = 1:length(vec)
        ind = find(glo_ind==vec(n));
        if length(ind)==nSesh
               cellIDs = [cellIDs; loc_ind(ind)'];
        end
    end
   
    for n=1:nSesh
        b = table2array(dA{n,1});
        b = b(:,2:end);
    if length(b(1,:))<cellIDs(end,n)
        cellIDs = cellIDs(1:end-1,:);
    end
    end
    
    
    
    %extract the traces of registered cells from each session
     newTraces = cell(nSesh,1);
    for n=1:nSesh
        a = unique(cellIDs(:,n))+1;
        b = dA{n,1};
        f = table2array(b);
        f = f(:,2:end);
        f(isnan(f)) = 0;
        newTraces{n,1} = f(:,a);
        
         if n == 5 & ~isempty(locs)
             coordinates = [posX(a) posY(a) diameter(a)];
         end
    end

    if isempty(locs)
        coordinates = [];
    end
    
    longRegistered = newTraces;
    disp(size(longRegistered{1,1}))
end
