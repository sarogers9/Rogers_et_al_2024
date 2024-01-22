function percAvgTraces = percBaselineAvg(seshMat,trialLength,timeBack,nTrials)
%this function takes finds the average trace with respect to baseline
%activity
% INPUTS:   seshMat - cell traces for the given set of trials (time in rows,
%                   cells in columns
%           trialLength - how many seconds per trial
%           timeBack - time back before stim onset to be considered
%               baseline
%           nTrials - number of trials
% OUTPUTS:  percAvgTraces - average traces over trials in terms of % baseline
%
%
%
%Written by Sophie A. Rogers, Corder Laboratory, University of Pennsylvania

%%
nCells = length(seshMat(1,:));
    for n=1:nCells
        a = [];
        b=[];
        for m=1:nTrials
            a(:,m) = seshMat(trialLength*(m-1)+1:trialLength*m,n);
        end
        b=mean(a,2);
        b=b.*100./mean(b(1:timeBack));
        percAvgTraces(:,n) = b;
    end