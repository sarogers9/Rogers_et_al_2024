function alignedTraces = aligntraces(traces,times,nSesh,dt,start,fin)
%this function aligns recorded neural activity to stimuli
%INPUTS =   traces - a cell with neural traces from each session in rows
%           times - a cell with stimulus times from each session in rows
%           nSesh - number of sessions
%           dt - sampling rate of video
%           start - amount of time prior to stimulus start to include
%           fin - amount of time after stimulus start to include
%OUTPUTS =  alignedTraces - a cell with matrices of aligned neural traces
%           concatenated in columns
%
%
%
%Written by Sophie A. Rogers, Corder Laboratory, University of Pennsylvania
%%  
    alignedTraces = cell(nSesh,1);
    for n=1:nSesh
        stimTimes = times{n,1};
        trials = length(stimTimes);
        
        seshTraces = traces{n,1};
        acts = zeros(round(length(seshTraces(:,1))/dt),length(seshTraces(1,:)));
        a=0;
        for m = 1:round(length(seshTraces(:,1))/dt)-1
            a=a+1;
            acts(m,:) = sum(seshTraces((a-1)*dt+1:a*dt,:))/dt;
        end
        
        a = [];
        for m=1:trials
            a = [a; acts(round(stimTimes(m))-start:round(stimTimes(m))+fin,:)];
        end
        
        alignedTraces{n,1} = a;
        
    end
