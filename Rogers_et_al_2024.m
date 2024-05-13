

%% In Vivo Imaging Analysis Pipeline for Rogers et al., 2024
%
%
%Each part depends on previous parts and will not run without one another
%unless specified.
%
%In the following sections you will:
% 1. Load the data of all animals in a given experimental group from
% storeRogers2024Data.m
% 2. Extract the traces of longitudinally registered cells from all files
% and align to user-defined stimulus times or times of interest
% 3. Plot spatial coordinates of longitudinally registered cells
% 4. Extract and plot average traces in terms of % baseline of all neurons 
% averaged over all trials in each session
% 5. Identify cells that are upregulated or downregulated in response to
% stimuli in a given session according to their average trace
% 6. Measure the number of sessions cells tend to maintain stimulus
% responsiveness
% 7. Calculate freezing encoding in each neuron in each session and compare
% within animals across groups
% 8. Load & plot TCA models, identifying session-dominant components and
% correlating to freezing
% 9. Calculate and plot normalized trial factors.
% 10. Calculate and plot normalized temporal factors.
% 11. Assess cumulative distribution of neuron weights
% 12. Identify Acq-Dominant, Ext1-Dominant, and Ext3-Dominant ensembles and
% their overlaps
% 13. Use a Fisher linear decoder to predict group based on average
% activity (z-score relative to Acquisition) of each ensemble over time
% 14. Extract average activity (z-score relative to Acquisition) of each
% cell in the ensemble during Extinction 1, Extinction 3
% 15. Plot example traces from specific ensembles.

% For aesthetic purposes, most data was exported from these sections in csv
% and txt files and imported in Prism, where most statistical tests were
% also done
%
%Written by Sophie A. Rogers, Corder Laboratory, University of Pennsylvania




%% 1. Load the data of all animals in a given experimental group this will take a minute or so


%the following structure contains all the raw data necessary to run the following code.
%the structure is broken down by animals (structures).
%within each animal field are subfields for sessions (structure), 
%cell registration table across sessions (table), 
%cell spatial properties (table), and tca output(table)

%sessions has 5 subfields (6 if broken videos), one for each session
%(habituation, acquisition, extinction 1, extinction 2, extinction 3). in
%each session subfield is calcium, a table of traces; stimulusTimes, a
%vector of stimulus delivery times; and freezing, a table of freezing


load('/Users/Rogers1/Desktop/rogers2024.mat')
%%
% metadata of experimental design.
trialStructure = table2array(rogers2024.trialStructure);
start = trialStructure(1); % trial start. defined here ashow much time back from stimulus times to begin collecting data. (seconds)
fin = trialStructure(2);   % how much time after stimulus start to stop collecting data. (seconds)
dt = trialStructure(3);    % sampling rate of recording (Hz)
nSesh = trialStructure(4); %number of sessions
trialLength = start+fin+1; %length of trial (seconds)
nTrials = trialStructure(5); %number of total trials
nAnimals = length(rogers2024.animals); %number of animals is the number of substructures
videoOffset = trialStructure(6);

%experimental groups
groupIDs{1}= rogers2024.groups(1).members;
groupIDs{2}= rogers2024.groups(2).members;
groupIDs{3}= rogers2024.groups(3).members;
groupIDs{4}= rogers2024.groups(4).members;
groupNames = {rogers2024.groups(1).name,rogers2024.groups(2).name,rogers2024.groups(3).name,rogers2024.groups(4).name};
nGroups = length(groupIDs); %number of sessions

%define trial event times
bl = rogers2024.stimuli(1).times;  %seconds
stim = rogers2024.stimuli(2).times; %seconds
trace = rogers2024.stimuli(3).times; %seconds
shock = rogers2024.stimuli(4).times; %seconds

stimuliNames = {rogers2024.stimuli(1).name,rogers2024.stimuli(2).name,rogers2024.stimuli(3).name,rogers2024.stimuli(4).name};

%assign trials to sessions
hab = rogers2024.sessions(1).trials; %trials
acq = rogers2024.sessions(2).trials; %trials
ext1 = rogers2024.sessions(3).trials; %trials
ext2 = rogers2024.sessions(4).trials; %trials
ext3 = rogers2024.sessions(5).trials; %trials
sessions = {hab, acq, ext1, ext2, ext3};
sessionNames = {rogers2024.sessions(1).name,rogers2024.sessions(2).name,rogers2024.sessions(3).name,rogers2024.sessions(4).name,rogers2024.sessions(5).name};


% 2. Extract traces of longitudinal cells and align to stimuli

%CUSTOM FUNCTIONS REQUIRED:
%   crossSeshAuto
%   aligntraces

%outputs
nCells = zeros(nAnimals,1);  %initialize array to store number of cells per animal
lr = cell(nAnimals,nSesh);   %initialize cell to store original activity matrices downsampled to 15Hz, for freezing decoding later
tensors = cell(nAnimals,1);  %initialize cell to store same information as poolMat but arranged in tensors of t seconds X c cells X T trials. Export for use in tensortools python kit
spatial = cell(nAnimals,1);  %initialize cell to store spatial coordinates of longitudinally registered cells
poolMat = [];

for n=1:nAnimals
    
    %initialize temp storage
    reRegistered={};
    dA = {};
    timesA = {};
    
    %loop through number of recordings and store raw traces
    nRec= length(rogers2024.animals(n).sessions);
    for m=1:nRec
        dA{m,1} = rogers2024.animals(n).sessions(m).calcium;
    end
    
    %call temporal regristration and spatial registration tables
    L = table2array(rogers2024.animals(n).cellreg);
    locs = rogers2024.animals(n).cellprops;
        
   %longitudinally register and spatially locate cells across 5 sessions
    [longRegistered, coordinates]= crossSeshAuto(L,dA,locs,nRec);
    
   %store spatial data
   spatial{n,1} = coordinates;
   
   %repair broken recordings in 3 animals
    if n == 12
        reRegistered{1,1} = [longRegistered{1,1}];
        reRegistered{2,1} = [longRegistered{2,1}];
        reRegistered{3,1} = [longRegistered{3,1}; longRegistered{4,1}];
        reRegistered{4,1} = [longRegistered{5,1}];
        reRegistered{5,1} = [longRegistered{6,1}];
        longRegistered = reRegistered;
    elseif n == 13
        x=longRegistered{1,1};
        y=longRegistered{2,1};
        reRegistered{1,1} = [x; x(end-16.8*dt/2+1:end,:); y(1:16.8*dt/2,:); y];
        reRegistered{2,1} = [longRegistered{3,1}];
        reRegistered{3,1} = [longRegistered{4,1}];
        reRegistered{4,1} = [longRegistered{5,1}];
        reRegistered{5,1} = [longRegistered{6,1}];
        longRegistered = reRegistered;
    elseif n == 24
        reRegistered{1,1} = [longRegistered{1,1}; longRegistered{2,1}];
        reRegistered{2,1} = [longRegistered{3,1}];
        reRegistered{3,1} = [longRegistered{4,1}];
        reRegistered{4,1} = [longRegistered{5,1}];
        reRegistered{5,1} = [longRegistered{6,1}];
        longRegistered = reRegistered;
    end
    
    %loop through sessions, downsample registered traces to 15Hz and store.
    %also store stimulus delivery times
     for m=1:nSesh
        x=resample(longRegistered{m,1},3,4);
        lr{n,m} = x;
        
        timesA{m,1} = rogers2024.animals(n).sessions(m).stimulusTimes;
     end
    
    %store number of cells
    nCells(n,1) = size(longRegistered{1,1},2);
    
    %align traces to stimulus times and store 1Hz data at the predefined
    %trial windows
    alignedTraces = aligntraces(longRegistered,timesA,nSesh,dt,start,fin);
    
    %concatenate aligned traces
    allAligned = [alignedTraces{1,1}; alignedTraces{2,1}; alignedTraces{3,1}; alignedTraces{4,1}; alignedTraces{5,1}];
    
    %transform into t x c x T tensor and store
    ten = reshape(allAligned,[trialLength, nCells(n), nTrials]);
    tensors{n,1} = ten;
    
    poolMat= [poolMat allAligned];
end


%% 3.Plot spatial coordinates of longitudinally registered cells. Fig 2B 

%enter animal number of interest; dataset 1 animal 1 is the representative image
animal = 1; 
propM39{1,1} = readtable('mm39hab-props.csv'); propM39{2,1} = readtable('mm39acq-props.csv'); propM39{3,1} = readtable('mm39ext1-props.csv'); propM39{4,1} = readtable('mm39ext2-props.csv'); propM39{5,1} = readtable('mm39ext3-props.csv');
data = spatial{animal};

for m=1:nSesh
    dat = propM39{m,1};
    posX{1,m} = table2array(dat(:,6));
    posY{1,m} = table2array(dat(:,7));
    diameter{1,m} = table2array(dat(:,9));
end


%plot cells at x coordinate by y coordinate at their size in pixels, scaled up to be visible
colors = {[.5 .5 .5], [1 0 0], [0.9290 0.6940 0.1250], [0 1 0], [0 1 1]};
figure
for m=1:nSesh
    hold on
    scatter(posX{1,m},posY{1,m},diameter{1,m}.*(20-3*(m-1)),colors{m},'filled')
end
hold on
scatter(data(:,1),data(:,2),data(:,3).*8,[0.4940 0.1840 0.5560]) 
%% 4. for fig 2D - heatmaps of average traces

%CUSTOM FUNCTIONS
% percBaselineAvg

%initialize storage cell average session activity of each animal by session 
%avPB = cell(nAnimals, nSesh); 
avPBt = cell(nAnimals, nSesh); 
%loop through animals
for n=1:nAnimals %can look at specific groups by replacing "1:nAnimals" with saline, responders, nonresponders, or non-shock
     
    %call t x c x T data from tensors and convert to time x cells matrix
    data = tensors{n};
    
    data = reshape(data,trialLength*nTrials,nCells(n));
    
    %initialize temp storage
    avPB = zeros(trialLength,nCells(n),nSesh);
    
    %loop through sessions
    for m=1:nSesh
        
        trials = sessions{m};
        
        %calculate average trial activity (% change from baseline) in each
        %session and store to pooled matrix across animals for plotting
        if m == 1
            act = percBaselineAvg(data(1:max(sessions{m})*trialLength,:),trialLength,start,length(trials));
        else
            act = percBaselineAvg(data(1 + max(sessions{m-1})*trialLength:max(sessions{m})*trialLength,:),trialLength,start,length(trials));
        end
        
        %store mean % baseline in time x cells x session tensor
        avPB(:,:,m) = act;
        
    end
    %store tensor for each animal
    avPBt{n} = avPB;
end

%pool cells from all animals within sessions
poolAvgAct = cell(nSesh,1);

%loop through animals to pool mean trial activity in a session
for n=1:nAnimals
    
    %call mean % baseline tensor from animal
    data1 = avPBt{n};
    
    %loop through sessions
    for m=1:nSesh
        
        %call a specific session
        data = data1(:,:,m);
        
        %store all average traces in a session from every animal
        pool = [poolAvgAct{m,1}, data];
        poolAvgAct{m,1} = pool;
    end
    
    %store pooled traces by session
    
end

%sort data to shock
acqAct = poolAvgAct{2,1};
[B,I] = sort(mean(acqAct(shock,:)));

%plot
titles = {'Hab','Acq','Ext1','Ext2','Ext3'};
figure
for m=1:nSesh
subplot(1,5,m)
    data2plot = poolAvgAct{m,1};
    heatmap(data2plot(:,I)','Colorlimits',[-1 200])
    title(titles{m})
    ylabel('Cells')
    xlabel('Time')
end
colormap(jet)

%% 5a. Identify cells that are upregulated or downregulated in response to stimuli in a given session according to their average trace for fig 2J-M

%CUSTOM FUNCTIONS
%permTest
tone = [11:35];
trace = [36:55];
shock = [56:58];
bl = {[1:10],[1:10],[1:10],[1:10]};
stimuli = {tone, trace, [tone trace], shock};
allCellsEver = cell(nAnimals, 8);
allCellsTrials = cell(nAnimals, 8);
numIt = 1000;

for a=1:nAnimals
    data = poolMat(:,sum(nCells(1:a-1))+1:sum(nCells(1:a)));
    
    disp(a)
    for m=1:5
        trialSet = sessions{m};
        
        x = zeros(60,size(data,2),length(trialSet));
        for t=1:length(trialSet)
            x(:,:,t) = data(60*(trialSet(t)-1)+1:60*(trialSet(t)),:);
        end
        
        data2 = squeeze(mean(x,3));
        
        for k=1:4
            stim = stimuli{k};
            bL = bl{k};
            
            
            responses = permTest(data2,bL,stim,numIt);
            
            sigCells{a,m,k} = responses;
        end
    end
end



%% 6. Measure the number of sessions cells tend to maintain stimulus responsiveness. For fig 2N-O

%CUSTOM FUNCTIONS
%permTest

%loop through animals to identify tone, trace, and shock responsive neurons
for a = 1:21%nAnimals
    %call average activity tensor 
    toTake = 0;
    if a>1
        toTake = sum(nCells(1:a-1));
    end
    %loop through stimuli after baseline
    for k = 1:length(rogers2024.stimuli)-1
            counter = zeros(5,nCells(a));
            for c=1:nCells(a)
                for m=1:5
                    counter(m,c) = ismember(c,sigCells{a,m,k}{1});
                    
                end
            end
            
            for m=1:5
                frac(a,m,k) = length(find(sum(counter)==m))./nCells(a);
            end
    end
         disp(a)
 end
    

%plot
for s=1:3
    labels{s} = strcat('Fraction  ',rogers2024.stimuli(s+1).name,'-responsive');
end

plotPersist(frac,groupIDs,groupNames,labels(1),labels)

%% 7. Measure the freezing encoding of single neurons. For fig 2P-Q, Supp Fig. 2B

%initialize storage matrix for the fraction of neurons and average freezing
%encoding
freezeCells = zeros(nSesh,nAnimals);
AUCsMean = zeros(nSesh,nAnimals);

delay = 30; %offset between miniscope and freezing recording (seconds)
dt = 15; %sampling rate of freezing (Hz)

%loop through animals
for g = 1:nAnimals
    
    %loop through sessions
    for h = 1:nSesh
        
        %initialize temporary storage counter and vector
        numCells = 0;
        AUCtot = [];
        
        data = lr{g,h};    % downsampled traces
        ac = zscore(data(delay*15:end,:)); %remove temporal offset of recordings, z-score over the session
        
        freezing = table2array(rogers2024.animals(g).sessions(h).freezing);      %loop through freezing
        freezing(freezing==100) = 1;  %turn into binary
        
        if length(freezing)>length(ac)  %make matrices the same length
            freezing = freezing(1:length(ac(:,1)));
        else
            ac = ac(1:length(freezing),:);
        end
        
        %initialize storage vectors for auROC and cell indices for auROC > .6
        AUCs = [];
        posCells = [];
        
        %loop through neuronsΩç
        for m=1:length(ac(1,:))
                act = ac(:,m);
                
                %create shuffle vector
                randN = randperm(round(length(freezing)));
                
                %pick random 50% of activity freezing values for training
                %set
                trainAct = act(randN(1:round(length(randN)/2)),:);
                trainFreeze = freezing(randN(1:round(length(randN)/2)),:);
                
                %pick random 50% of activity freezing values for test
                %set
                testAct = act(randN(round(length(randN)/2)+1:end),:);
                testFreeze = freezing(randN(round(length(randN)/2)+1:end),1);
                
                %train logistic regressor
                [B,DEV,STATS] = glmfit(trainAct,trainFreeze,'binomial','Link','logit');
                
                %test
                pred = glmval(B,testAct,'logit');
                
                %calculate auROC
                [X,Y,T,AUC] = perfcurve(testFreeze,pred,1, 'XCrit','FPR','YCrit','TPR');
                
                %record auROC
                AUCs = [AUCs AUC];
               
                if AUC > 0.6
                    
                    %count cells with auROC > 0.6
                    numCells = numCells+1;
                    
                    %record their index
                    posCells = [posCells m];

                end

        end
        
        %record average auROC of animal in session
        AUCsMean(h,g) = mean(AUCs);
        
        %record fraction of freezing encoding neurons
        freezeCells(h,g) = length(posCells)/nCells(g);
    end
    
    
    disp(g)
  
    
end


%plot
ylab1 = {'Mean freezing encoding'};
ylab2 = {'Fraction of Cells'};
ylim1 = [.5 .7];
ylim2 = [0 .7];

plot2metrics(AUCsMean,freezeCells,sessionNames,groupNames,groupIDs,ylab1,ylab2,ylim1,ylim2)

%% 8. Extract freezing and align to stimuli
%calculate tensors from Williams et al. TCA code (https://github.com/ahwillia/tensortools)
% load tensors 

%get trial by trial freezing
freezingInTrials = {};
dt = 15;
timesF{1,1} = rogers2024.animals(24).sessions(1).stimulusTimes-videoOffset;
timesF{2,1} = rogers2024.animals(24).sessions(2).stimulusTimes-videoOffset;
timesF{3,1} = rogers2024.animals(24).sessions(3).stimulusTimes-videoOffset; 
timesF{4,1} = rogers2024.animals(24).sessions(4).stimulusTimes-videoOffset; 
timesF{5,1} = rogers2024.animals(24).sessions(5).stimulusTimes-videoOffset; 

%create second by second % freezing
for n=1:nAnimals
    fr = cell(nSesh,1);
    for m=1:nSesh
        fr{m,1} = table2array(rogers2024.animals(n).sessions(m).freezing);
    end
    

    if n==5
        timesf = timesF;
        timesf{3,1} = [120; 120+round(2844/15); 234+round(2844/15); 331+round(2844/15); 441+round(2844/15); 543+round(2844/15)];
    else
        timesf = timesF;
    end

    fIT = aligntraces(fr,timesf,nSesh,dt,start,fin);
    freezingInTrials = [freezingInTrials fIT];
end


for n=1:nAnimals
    a=0;
    for m=1:nSesh
        fre = freezingInTrials{m,n};
        if ismember(m,[1 2])
            for t = 1:8
               a=a+1;
               freTot(a,n) = mean(fre(trialLength*(t-1)+1:trialLength*t));
            end
        else
            for t = 1:6
                a=a+1;
               freTot(a,n) = mean(fre(trialLength*(t-1)+1:trialLength*t));
            end
        end
    end
end

%% 9. Load & plot TCA models, identifying session-dominant components and correlating to freezing for fig. 3A,D-F

nDims = 5;

coefs = zeros(nAnimals,nDims);
%to plot representative TCAs, uncomment the figure file in the for-loop -
%the representative image in the paper is from dataset 1 animal 12
%below
for a=1:nAnimals
    tt = table2array(rogers2024.animals(a).tca);
    fr = freTot(:,a);
    
    timeFactor = tt(2:61,:);
    neuronFactor = tt(62:end-34,:);
    trialFactor = tt(end-33:end,:);
    
    %store neuron factors - used in Fig. 3G, 4-5
    neuFac{a,1} = neuronFactor;
    timFac{a,1} = timeFactor;
    triFac{a,1} = trialFactor;
    
    %Identify session-dominant factors - used in Fig. 3D-F, 4, 5
    for n=1:nSesh
        facLoads = mean(trialFactor(sessions{n},:));
        domFacs(a,n) = find(facLoads==max(facLoads));
    end
    
    %calculate strength of each factor in each session - used in Fig. D,E
    for n=1:nSesh
        facLoads = mean(trialFactor(sessions{n},:));
        
        for m=1:nDims
            strength(a,n,m) = facLoads(domFacs(a,m))/sum(facLoads(domFacs(a,:)));
        end
    end
    
    %calculate correlation of trial factor loadings with trial-by-trial
    %freezing over extinction with multiple linear regression. save stats
    [b, bint, r, rint,stats] = regress(fr, [trialFactor ones(nTrials,1)]);
    coefs(a,:) = b(1:5);
    R2(a,:) = stats(1); 
    fs(a,:) = stats(2); 
    ps(a,:) = stats(3); 

    %uncomment to plot representative figure
    %plotTCAmodel(timeFactor,neuronFactor,trialFactor,fr)
end

%% 10. Plot normalized trial factors. for Supp Fig. 4A-E

%initialize trial weight storage vectors
normTrialWeights = zeros(nTrials, nAnimals, nSesh);

for n=1:nAnimals
    
    %call trial weights from a given animal
    tWeights = triFac{n};
    
    %store normalized trial weights
    for m=1:nSesh
        normTrialWeights(:,n,m) = tWeights(:,domFacs(n,m))./max(tWeights(:,domFacs(n,m)));
    end
    
end

%plot
ys = [0 1];% {'Hab-Dom','Acq-Dom','Ext1-Dom','Ext2-Dom','Ext3-Dom'};
colors = {'k','b','r','g'};
figure
for n=1:nSesh
subplot(1,5,n)
    for m=1:nGroups
        data=normTrialWeights(:,rogers2024.groups(m).members,n);
        errorbar(1:length(data),mean(data,2),std(data')/sqrt(length(data)),colors{m},'Linewidth',1)
        hold on
    end
    ylim([ys])
    ylabel('Normalized trial weights')
    xlabel('Time')
    title(strcat(sessionNames{n},'-Dominant Component'))
    if n==5
        legend(groupNames)
    end
end


%% 11. Plot normalized temporal factors. for Supp Fig. 3F-J
%initialize trial weight storage vectors
normTrialWeights = zeros(nTrials, nAnimals, nSesh);

for n=1:nAnimals
    
    %call trial weights from a given animal
    tWeights = timFac{n};
    
    %store normalized trial weights
    for m=1:nSesh
        normTempWeights(:,n,m) = tWeights(:,domFacs(n,m))./max(tWeights(:,domFacs(n,m)));
    end
    
end

%plot
ys = [0 1];% {'Hab-Dom','Acq-Dom','Ext1-Dom','Ext2-Dom','Ext3-Dom'};
colors = {'k','b','r','g'};
figure
for n=1:nSesh
subplot(1,5,n)
    for m=1:nGroups
        data=normTempWeights(:,rogers2024.groups(m).members,n);
        errorbar(1:length(data),mean(data,2),std(data')/sqrt(length(data)),colors{m},'Linewidth',1)
        hold on
    end
    ylim([ys])
    ylabel('Normalized trial weights')
    xlabel('Time')
    title(strcat(rogers2024.sessions(n).name,'-Dominant Component'))
    if n==5
        legend(groupNames{1},groupNames{2},groupNames{3},groupNames{4})
    end
end


%% 12. Assess cumulative distribution of neuron weights. for Fig. 4G,H

% set threshold range
threshold = linspace(0,10,11);

%define storage vectors for fraction of cells included in Acq-, Ext1-, and
%Ext3-dominant ensembles
nComps = 3;
comps = [2 3 5];
cellsInComp = zeros(nAnimals, length(threshold),nComps);

for t = 1:length(threshold)
    aCells = [];
    e1Cells = [];
    e3Cells = [];
    
    for n = 1:nAnimals
        
        %call neuron factor weights from animal
        weights = neuFac{n};
        
        for c = 1:nComps
        %calculate fractions of neurons with w > threshold and store
            cellsInComp(n,t,c) = length(find(weights(:,domFacs(n,comps(c)))>threshold(t)))/nCells(n);
        end
        
        %record ensembles at each threshold
        if n == 1
            c = 0;
        else
            c = sum(nCells(1:n-1));
        end
    
        aCells = [aCells; find(weights(:,domFacs(n,2))>threshold(t))+c];
        e1Cells = [e1Cells; find(weights(:,domFacs(n,3))>threshold(t))+c];
        e3Cells = [e3Cells; find(weights(:,domFacs(n,5))>threshold(t))+c];

    end
    
    %collect all the cells in the Acq, Ext1, and Ext3 ensembles at the 
    %predetermined varying thresholds
    collection{t,1} = aCells;
    collection{t,2} = e1Cells;
    collection{t,3} = e3Cells;
end

%plot
colors = {'r','y','c'};
figure
for c=1:nComps
    data = cellsInComp(:,:,c);
    errorbar(threshold,mean(data),std(data)/sqrt(nAnimals-1),colors{c},'Linewidth',3)
    hold on
end
hold on
xline(1)
xlabel('Threshold weight')
ylabel('Fraction of neurons')
legend({'Acq-Dominant','Ext1-Dominant','Ext3-Dominant'})
title('Choosing dominant neurons')

%% 13. Identify Acq-Dominant, Ext1-Dominant, and Ext3-Dominant ensembles and their overlaps for Fig. 5B, 6A

%set threshTest to 0 if reproducing main figs. set threshTest to 1 if
%reproducing Supp Fig 6. note that null thresholds were not determined for
%non-shock mice
threshTest = 0;

%define groups
responders = rogers2024.groups(3).members;
nonresponders = rogers2024.groups(4).members;
rapid = rogers2024.groups(1).members;
slow = rogers2024.groups(2).members;
nonshock = rogers2024.groups(5).members;

if threshTest == 1
    threshes = table2array(readtable('fMeans.csv'));
    threshold = threshes([[15:21] [1:14]]);
    nAnimals = 21;
end

%create storage vectors for all Acq-, Ext1-, and Ext3-dominant cells
acqDomCells = [];
ext1DomCells = [];
ext3DomCells = [];

%create storage vectors for their overlaps
acqOnly = [];
ext1Only = [];
ext3Only = [];
acqExt1 = [];
acqExt3 = [];
ext1Ext3 = [];
acqExt1Ext3 = [];

%allocate vectors to assign cell indices to for each group
cR = [];
cNR = [];
cSS = [];
cSR = [];
cNS = [];
for a = 1:nAnimals
    if threshTest == 1
        th = threshold(a);
    else
        th = 1;
    end
    %number to add to convert within-animal cell index to pooled cell index
    if a == 1
        c = 0;
    else
        c = sum(nCells(1:a-1));
    end
    
    %load cell indices into their group
    if ismember(a,responders)
        cR = [cR [1+c:nCells(a)+c]];
    elseif ismember(a,nonresponders)
        cNR = [cNR [1+c:nCells(a)+c]];
    elseif ismember(a,rapid)
        cSR = [cSR [1+c:nCells(a)+c]];
    elseif ismember(a,slow)
        cSS = [cSS [1+c:nCells(a)+c]];
    elseif ismember(a,nonshock)
        cNS = [cNS [1+c:nCells(a)+c]];
    end
    
    %call neuron factor weights
    weights = neuFac{a};
    
    %dominant cells are defined by those weights > 1 in the component
    %dominating that session
    aCells = find(weights(:,domFacs(a,2))>th);
    e1Cells = find(weights(:,domFacs(a,3))>th);
    e3Cells = find(weights(:,domFacs(a,5))>th);

    overlap = aCells(ismember(aCells,e1Cells));
    
    %cells dominating all three sessions
    ae1e3 = overlap(ismember(overlap,e3Cells));
    
    %cells dominating all acq, ext1
    ae1 = overlap(~ismember(overlap,e3Cells));
    
    %cells dominating all acq, ext3
    overlap = aCells(ismember(aCells,e3Cells));
    
    ae3 = overlap(~ismember(overlap,e1Cells));
    
    %cells dominating all ext1, ext3
    overlap = e1Cells(ismember(e1Cells,e3Cells));
    mat(a,1) = length(overlap)./nCells(a);
    e1e3 = overlap(~ismember(overlap,aCells));
    
    %cells dominating only acq
    aOnly = aCells(~ismember(aCells,ae1e3));
    aOnly = aOnly(~ismember(aOnly,ae1));
    aOnly = aOnly(~ismember(aOnly,ae3));
    
    %cells dominating only ext1
    e1Only = e1Cells(~ismember(e1Cells,ae1e3));
    e1Only = e1Only(~ismember(e1Only,ae1));
    e1Only = e1Only(~ismember(e1Only,e1e3));
    
    %cells dominating only ext3
    e3Only = e3Cells(~ismember(e3Cells,ae1e3));
    e3Only = e3Only(~ismember(e3Only,ae1));
    e3Only = e3Only(~ismember(e3Only,e1e3));
    
    %store fraction of overlap for individual animals
    ensOverlap(a,1,1) = length(aOnly)/length(aCells);
    ensOverlap(a,2,1) = length(ae1)/length(aCells);
    ensOverlap(a,3,1) = length(ae3)/length(aCells);
    ensOverlap(a,4,1) = length(ae1e3)/length(aCells);
    
    ensOverlap(a,1,2) = length(e1Only)/length(e1Cells);
    ensOverlap(a,2,2) = length(ae1)/length(e1Cells);
    ensOverlap(a,3,2) = length(e1e3)/length(e1Cells);
    ensOverlap(a,4,2) = length(ae1e3)/length(e1Cells);
    
    ensOverlap(a,1,3) = length(e3Only)/length(e3Cells);
    ensOverlap(a,2,3) = length(ae3)/length(e3Cells);
    ensOverlap(a,3,3) = length(e1e3)/length(e3Cells);
    ensOverlap(a,4,3) = length(ae1e3)/length(e3Cells);
    
    %assign neurons indices for activity pooled across animals and store in
    %vectors
    
    ensAns(a).acqcells = aCells+c;
    ensAns(a).ext1cells = e1Cells+c;
    ensAns(a).ext3cells = e3Cells+c;
    ensAns(a).a = aOnly+c;
    ensAns(a).e1 = e1Only+c;
    ensAns(a).e3 = e3Only+c;
    ensAns(a).ae1 = ae1+c;
    ensAns(a).ae3 = ae3+c;
    ensAns(a).e1e3 = e1e3+c;
    ensAns(a).ae1e3 = ae1e3+c;
    
    acqDomCells = [acqDomCells; aCells + c];
    ext1DomCells = [ext1DomCells; e1Cells + c];
    ext3DomCells = [ext3DomCells; e3Cells + c];
    
    acqOnly = [acqOnly; aOnly+c];
    ext1Only = [ext1Only; e1Only+c];
    ext3Only = [ext3Only; e3Only+c];
    acqExt1 = [acqExt1; ae1+c];
    acqExt3 = [acqExt3; ae3+c];
    ext1Ext3 = [ext1Ext3; e1e3+c];
    acqExt1Ext3 = [acqExt1Ext3; ae1e3+c];
end


labels = {'Overlap with Acq-Dom Neurons','Overlap with Ext1-Dom Neurons','Overlap with Ext3-Dom Neurons'};
labels2{1} = ['Acq Only', 'Acq/Ext1', 'Acq/Ext3', 'Acq/Ext1/Ext3'];
labels2{2} = ['Ext1 Only', 'Acq/Ext1', 'Ext1/Ext3', 'Acq/Ext1/Ext3'];
labels2{3} = ['Ext3 Only', 'Acq/Ext3', 'Ext1/Ext3', 'Acq/Ext1/Ext3'];

figure
for n = 1:nComps
    for m = 1:nGroups
        subplot(3,4,m+4*(n-1))
        barWithError(ensOverlap(rogers2024.groups(m).members,:,n))
        ylabel('Fraction of Cells')
        xticklabels(labels2{n})
        ylabel(labels{n})
        ylim([0 .8])
        if n==1
        title(rogers2024.groups(m).name)
        end
    end
end



%% 14. Use a Fisher linear decoder to predict group based on average activity (z-score relative to Acquisition) of each ensemble over time for fig. 6B

%TO AVOID HAVING TO RE-RUN PREVIOUS SCRIPTS, LOAD THESE FILES WITH
%ENSEMBLES' CELL INDICES
%
% IF RUNNING THRESHTEST SKIP THIS SECTION
% 
% poolMat = [];
% 
% for n=1:nAnimals
%     poolMat = [poolMat reshape(tensors{n},60*34,nCells(n))];
% end

%define time windows
timeAcq = [trialLength*max(hab)+1:trialLength*max(acq)];
timeExt1 = [trialLength*max(acq)+1:trialLength*max(ext1)];
timeExt3 = [trialLength*max(ext2)+1:trialLength*max(ext3)];

%normalize activity to acquisition
poolMat = (poolMat-mean(poolMat(timeAcq,:)))./std(poolMat(timeAcq,:));



%define ensembles to loop through
ens = fieldnames(ensAns(1));



%create matrices of average ensembles acts over time in the session of
%interest
sess = 'Ex3'; %set to Acq, Ex1, or Ex3

if sess == 'Ex1'
    times = timeExt1;
elseif sess == 'Ex3'
    times = timeExt3;
elseif sess == 'Acq'
    times = timeAcq;
end

%define binary class vector
classes = [zeros(length(times),1); ones(length(times),1)];

ens = fieldnames(ensAns(1));
numEnsembles = length(ens);

ensActsR =[];
ensActsPR =[];
ensActsPS =[];

for l = 1:numEnsembles
    for a=1:nAnimals
        if ismember(a,rogers2024.groups(1).members)
            ensActsR = [ensActsR mean(poolMat(times,ensAns(a).(ens{l})),2)];
        elseif ismember(a,rogers2024.groups(4).members)
            ensActsPS = [ensActsPS mean(poolMat(times,ensAns(a).(ens{l})),2)];
        elseif ismember(a,rogers2024.groups(3).members)
            ensActsPR = [ensActsPR mean(poolMat(times,ensAns(a).(ens{l})),2)];
        end
    end
end

%number of iterations of decoder
numIt = 100;

%set of comparisons - res vs. nonres, res vs. sal, nonres vs. sal
comparisons = [0:2];

%store 100 iterations of the model (rows), for each ensemble plus a round
%where every ensemble is a predictor (columns), for each comparison
%(columns')
accuracy = zeros(numIt,numEnsembles+1,length(comparisons));
accuracyShuff = zeros(numIt,numEnsembles+1,length(comparisons));

for c = 1:length(comparisons)
for l = 1:numEnsembles
    if comparisons(c) == 0
        data = [ensActsPR(:,l); ensActsPS(:,l)]; %compare responders to nonresponders
    elseif comparisons(c) == 1
        data = [ensActsPR(:,l); ensActsR(:,l)];  %compare responders to saline
    else
        data = [ensActsPS(:,l); ensActsR(:,l)]; %compare nonresponders to saline
    end
    
    for n=1:100
    [X,Y,x,y] = split_data(data,classes,.5); %split_data is a custom function, described below that randomly selects test_size % of your data to test on and 1 - test_size to train on
    model = fitcdiscr(x,y); %fit discr fits a linear model to your training data & classes
    y_pred = predict(model, X); %predict applies your model to your test data to generate class predictions
    accuracy(n,l,c) = mean(y_pred == Y); %the accuracy of your model is the number of instances in which your predicted classes matched your a priori classes divided by total number of predictions

    shuffle = randperm(length(Y));
    Y = Y(shuffle);
    y_pred = predict(model, X); %predict applies your model to your test data to generate class predictions
    accuracyShuff(n,l,c) = mean(y_pred == Y); %the accuracy of your model is the number of instances in which your predicted classes matched your a priori classes divided by total number of predictions
    end
end
end


%plot accuracies and shuffles
titles = {'Acq-Dom','Ext1-Dom','Ext3-Dom','Acq Only', 'Ext1 Only', 'Ext3 Only', 'Acq/Ext1', 'Acq/Ext3', 'Ext1/Ext3', 'Acq/Ext1/Ext3','All ensembles'};

%loop through comparisons and create a figure for each one
for c=1:3 
    figure
    for l=1:11
    subplot(6,2,l)
    barWithError(accuracy(:,l,c))
    hold on
    barWithError(accuracyShuff(:,l,c))
    ylim([0 1])
    title(titles{l})
    ylabel('Accuracy')
    end
end

%% 15. Extract average activity (z-score relative to Acquisition) of each cell in the ensemble during Extinction 1, Extinction 3. for fig. 5D-F, fig. 6C-I, fig. 7B-D


%loop through animals
for n = 1:nAnimals
    
    %loop through all ensembles
    for l=1:10
        
        cellPop = ensAns(n).(ens{l});

        %calculate average activity of cells during Extinction 1 and
        %Extinction 3
        x =  [mean(poolMat(timeExt1,cellPop))' mean(poolMat(timeExt3,cellPop))'];
        
        %store activity matrices
        changedActivity2{l,n} = x;
       
    end
    
end

%pool cells
for n=1:4
    group = rogers2024.groups(n).members;
    for l=1:10
        pop=[];
        pgroup=[];
        for g=1:length(group)
            an = group(g);
            toAdd=0;
            if an >1
                toAdd= sum(nCells(1:n-1));
            end
            pop=[pop; changedActivity2{l,an}];
            pgroup = [pgroup; ensAns(an).(ens{l})];
        end
    popGroup{l,n} = pgroup;
    changedActivity{l,n} = pop;
    end
end

%plot
co = {[.5 .5 .5],'k','c','r'};
co2 = {'c','k',[.5 .5 .5],'r'};
figure
subplot(121)
for n=1:4
errorbar(1:length(poolMat),mean(poolMat(:,popGroup{1,n}),2),std(poolMat(:,popGroup{1,n})')./sqrt(length(popGroup{1,n})-1),'color',co{n})
hold on
end
title('Acquistion-dominant neurons during Extinction 1')
xlabel('Time (seconds)')
ylabel('z-score from Acq')
subplot(122)
for n=1:4
errorbar(1:length(poolMat),mean(poolMat(:,popGroup{4,n}),2),std(poolMat(:,popGroup{4,n})')./sqrt(length(popGroup{4,n})-1),'color',co2{n})
hold on
title('Extinction 3-dominant neurons during Extinction 3')
xlabel('Time (seconds)')
legend({'responders','slow','rapid','nonresponders'})
end

