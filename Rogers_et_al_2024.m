%% In Vivo Imaging Analysis Pipeline for Rogers et al., 2024
%
%
%Each part depends on previous parts and will not run without one another
%unless specified.
%
%In the following sections you will:
% 1. Load the data of all animals in a given experimental group
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
% 16. Plot average dFoF of non-intersecting Acq-dominant and Ext1/3-dominant
% ensembles and their ratios over seconds, trials
%
% For aesthetic purposes, most data was exported from these sections in csv
% and txt files and imported in Prism, where most statistical tests were
% also done with Prism software
%
%Written by Sophie A. Rogers, Corder Laboratory, University of Pennsylvania




%% 1. Load the data of all animals in a given experimental group

% VARIABLES: 
% dataset - Boolean variable that indicates which dataset to load. nAnimals
%   is the number of animals in the dataset you choose.
% d - a cell containing the recorded neural traces of each animal in
%   columns in each session in rows. If you have different numbers videos for 
%   each animal, create new cells "d2", "d3" etc.
% f - a similarly structured cell containing  data of each animal.
%   Freezing data is not analyzed in this original version of the pipeline.
% times - a similarly structured cell containing the times at which stimuli
%   were delivered to each animal in each session
% L - a cell with animals in rows where each element is the IDPS indexing 
%output of longitudinal registration
% props - a csv with index, spatial coordinates, and size of each cell in
% the Ext3 recording session

dataset =0;
if dataset==0 %saline
    nAnimals = 7;
elseif dataset == 1 %psilocybin
    nAnimals = 14;
elseif dataset == 2 %nonshock
    nAnimals = 3;
end

if dataset ==0
    d{1,1} =  readtable('mm43hab.csv');
    d{1,2} =  readtable('mm51hab.csv');
    d{1,3} =  readtable('mm55hab.csv');
    d{1,4} =  readtable('mm71hab.csv');
    d{1,5} =  readtable('mm56hab.csv');
    d{1,6} =  readtable('mm58hab.csv');
    d{1,7} =  readtable('mm57hab.csv');
    % 
    d{2,1} =  readtable('mm43acq.csv');
    d{2,2} =  readtable('mm51acq.csv');
    d{2,3} =  readtable('mm55acq.csv');
    d{2,4} =  readtable('mm71acq.csv');
    d{2,5} =  readtable('mm56acq.csv');
    d{2,6} =  readtable('mm58acq.csv');
    d{2,7} =  readtable('mm57acq.csv');

    d{3,1} =  readtable('mm43ext1.csv');
    d{3,2} =  readtable('mm51ext1.csv');
    d{3,3} =  readtable('mm55ext1.csv');
    d{3,4} =  readtable('mm71ext1.csv');
    d{3,5} =  readtable('mm56ext1.csv');
    d{3,6} =  readtable('mm58ext1.csv');
    d{3,7} =  readtable('mm57ext1.csv');

    d{4,1} =  readtable('mm43ext2.csv');
    d{4,2} =  readtable('mm51ext2.csv');
    d{4,3} =  readtable('mm55ext2.csv');
    d{4,4} =  readtable('mm71ext2.csv');
    d{4,5} =  readtable('mm56ext2.csv');
    d{4,6} =  readtable('mm58ext2.csv');
    d{4,7} =  readtable('mm57ext2.csv');

    d{5,1} =  readtable('mm43ext3.csv');
    d{5,2} =  readtable('mm51ext3.csv');
    d{5,3} =  readtable('mm55ext3.csv');
    d{5,4} =  readtable('mm71ext3.csv');
    d{5,5} =  readtable('mm56ext3.csv');
    d{5,6} =  readtable('mm58ext3.csv');
    d{5,7} =  readtable('mm57ext3.csv');

    f{1,1} =  table2array(readtable('fm43_hab.csv'));
    f{1,2} =  table2array(readtable('M51_hab.csv'));
    f{1,3} =  table2array(readtable('M55_hab.csv'));
    f{1,4} =  table2array(readtable('M71_hab.csv'));
    f{1,5} =  table2array(readtable('M56_hab.csv'));
    f{1,6} =  table2array(readtable('M58_hab.csv'));
    f{1,7} =  table2array(readtable('M57_hab.csv'));

    f{2,1} =  table2array(readtable('fm43_acq.csv'));
    f{2,2} =  table2array(readtable('M51_acq.csv'));
    f{2,3} =  table2array(readtable('M55_acq.csv'));
    f{2,4} =  table2array(readtable('M71_acq.csv'));
    f{2,5} =  table2array(readtable('M56_acq.csv'));
    f{2,6} =  table2array(readtable('M58_acq.csv'));
    f{2,7} =  table2array(readtable('M57_acq.csv'));

    f{3,1} =  table2array(readtable('fm43_ext1.csv'));
    f{3,2} =  table2array(readtable('M51_ext1.csv'));
    f{3,3} =  table2array(readtable('M55_ext1.csv'));
    f{3,4} =  table2array(readtable('M71_ext1.csv'));
    f{3,5} =  table2array(readtable('M56_ext1.csv'));
    f{3,6} =  table2array(readtable('M58_ext1.csv'));
    f{3,7} =  table2array(readtable('M57_ext1.csv'));

    f{4,1} =  table2array(readtable('fm43_ext2.csv'));
    f{4,2} =  table2array(readtable('M51_ext2.csv'));
    f{4,3} =  table2array(readtable('M55_ext2.csv'));
    f{4,4} =  table2array(readtable('M71_ext2.csv'));
    f{4,5} =  table2array(readtable('M56_ext2.csv'));
    f{4,6} =  table2array(readtable('M58_ext2.csv'));
    f{4,5} =  table2array(readtable('M56_ext2.csv'));
    f{4,6} =  table2array(readtable('M58_ext2.csv'));
    f{4,7} =  table2array(readtable('M57_ext2.csv'));

    f{5,1} =  table2array(readtable('fm43_ext3.csv'));
    f{5,2} =  table2array(readtable('M51_ext3.csv'));
    f{5,3} =  table2array(readtable('M55_ext3.csv'));
    f{5,4} =  table2array(readtable('M71_ext3.csv'));
    f{5,5} =  table2array(readtable('M56_ext3.csv'));
    f{5,6} =  table2array(readtable('M58_ext3.csv'));
    f{5,7} =  table2array(readtable('M57_ext3.csv'));

    %Create vectors containing the stimulus delivery times for each animal &
    %session.

    %Habituation and Acquisition = [150; 249; 364; 471; 572; 678; 794; 908]
    %Extinction sessions = [150; 264; 3trialLength; 471; 573; 682];
    %if you started the session late or early with respect to the miniscope
    %recording (that is 30s ± error), add that error into your stimulus times
    %for that session


    %Acq times
    times{1,1} = [150; 249; 364; 471; 572; 678; 794; 908];
    times{1,2} = [150; 249; 364; 471; 572; 678; 794; 908];
    times{1,3} = [150; 249; 364; 471; 572; 678; 794; 908];
    times{1,4} = [150; 249; 364; 471; 572; 678; 794; 908];
    times{1,5} = [150; 249; 364; 471; 572; 678; 794; 908];
    times{1,6} = [150; 249; 364; 471; 572; 678; 794; 908];
    times{1,7} = [150; 249; 364; 471; 572; 678; 794; 908];

    times{2,1} = [150; 249; 364; 471; 572; 678; 794; 908];
    times{2,2} = [150; 249; 364; 471; 572; 678; 794; 908];
    times{2,3} = [150; 249; 364; 471; 572; 678; 794; 908];
    times{2,4} = [150; 249; 364; 471; 572; 678; 794; 908];
    times{2,5} = [150; 249; 364; 471; 572; 678; 794; 908];
    times{2,6} = [150; 249; 364; 471; 572; 678; 794; 908];
    times{2,7} = [150; 249; 364; 471; 572; 678; 794; 908];

    times{3,1} = [150; 264; 361; 471; 573; 682];  
    times{3,2} = [150; 264; 361; 471; 573; 682]; 
    times{3,3} = [150; 264; 361; 471; 573; 682]; 
    times{3,4} = [150; 264; 361; 471; 573; 682]; 
    times{3,5} = [150; 249; 364; 471; 572; 678];
    times{3,6} = [150; 264; 361; 471; 573; 682]; 
    times{3,7} = [150; 264; 361; 471; 573; 682]; 

    times{4,1} = [124; 239; 341; 446; 556; 658]; 
    times{4,2} = [150; 264; 361; 471; 573; 682]; 
    times{4,3} = [150; 264; 361; 471; 573; 682]; 
    times{4,4} = [150; 264; 361; 471; 573; 682]; 
    times{4,5} = [150; 264; 361; 471; 573; 682]; 
    times{4,6} = [150; 264; 361; 471; 573; 682]; 
    times{4,7} = [150; 264; 361; 471; 573; 682]; 

    times{5,1} = [123; 238; 339; 445; 555; 657]; 
    times{5,2} = [150; 264; 361; 471; 573; 682]; 
    times{5,3} = [150; 264; 361; 471; 573; 682]; 
    times{5,4} = [150; 264; 361; 471; 573; 682]; 
    times{5,5} = [150; 264; 361; 471; 573; 682]; 
    times{5,6} = [150; 264; 361; 471; 573; 682]; 
    times{5,7} = [150; 264; 361; 471; 573; 682];


    L{1,1}= table2array(readtable('mm43long.csv'));
    L{2,1}= table2array(readtable('mm51long.csv'));
    L{3,1}= table2array(readtable('mm55long.csv'));
    L{4,1}= table2array(readtable('mm71long.csv'));
    L{5,1}= table2array(readtable('mm56long.csv'));
    L{6,1}= table2array(readtable('mm58long.csv'));
    L{7,1}= table2array(readtable('mm57long.csv'));

    props{1,1} = readtable('mm43ext3-props.csv');
    props{1,2} = readtable('mm51ext3-props.csv');
    props{1,3} = readtable('mm55ext3-props.csv');
    props{1,4} = readtable('mm71ext3-props.csv');
    props{1,5} = readtable('mm56ext3-props.csv');
    props{1,6} = readtable('mm58ext3-props.csv');
    props{1,7} = readtable('mm57ext3-props.csv');

elseif dataset == 1
    d{1,1} =  readtable('mm39hab.csv');
    d{1,2} =  readtable('mm41hab.csv');
    d{1,3} =  readtable('mm46hab.csv');
    d{1,4} =  readtable('mm48hab.csv');
    d{1,5} =  readtable('mm53hab.csv');
    d{1,6} =  readtable('mm68hab.csv');
    d{1,7} =  readtable('mm70hab.csv');
    d{1,8} =  readtable('mm59hab.csv');
    d{1,9} =  readtable('mm64hab.csv');
    d{1,10} =  readtable('mm65hab.csv');
    d{1,11} =  readtable('mm69hab.csv');
    d{1,12} =  readtable('mm61hab.csv');

    d{2,1} =  readtable('mm39acq.csv');
    d{2,2} =  readtable('mm41acq.csv');
    d{2,3} =  readtable('mm46acq.csv');
    d{2,4} =  readtable('mm48acq.csv');
    d{2,5} =  readtable('mm53acq.csv');
    d{2,6} =  readtable('mm68acq.csv');
    d{2,7} =  readtable('mm70acq.csv');
    d{2,8} =  readtable('mm59acq.csv');
    d{2,9} =  readtable('mm64acq.csv');
    d{2,10} =  readtable('mm65acq.csv');
    d{2,11} =  readtable('mm69acq.csv');
    d{2,12} =  readtable('mm61acq.csv');

    d{3,1} =  readtable('mm39ext1.csv');
    d{3,2} =  readtable('mm41ext1.csv');
    d{3,3} =  readtable('mm46ext1.csv');
    d{3,4} =  readtable('mm48ext1.csv');
    d{3,5} =  readtable('mm53ext1.csv');
    d{3,6} =  readtable('mm68ext1.csv');
    d{3,7} =  readtable('mm70ext1.csv');
    d{3,8} =  readtable('mm59ext1.csv');
    d{3,9} =  readtable('mm64ext1.csv');
    d{3,10} =  readtable('mm65ext1.csv');
    d{3,11} =  readtable('mm69ext1.csv');
    d{3,12} =  readtable('mm61ext1.csv');

    d{4,1} =  readtable('mm39ext2.csv');
    d{4,2} =  readtable('mm41ext2.csv');
    d{4,3} =  readtable('mm46ext2.csv');
    d{4,4} =  readtable('mm48ext2.csv');
    d{4,5} =  readtable('mm53ext2.csv');
    d{4,6} =  readtable('mm68ext2.csv');
    d{4,7} =  readtable('mm70ext2.csv');
    d{4,8} =  readtable('mm59ext2.csv');
    d{4,9} =  readtable('mm64ext2.csv');
    d{4,10} =  readtable('mm65ext2.csv');
    d{4,11} =  readtable('mm69ext2.csv');
    d{4,12} =  readtable('mm61ext2.csv');

    d{5,1} =  readtable('mm39ext3.csv');
    d{5,2} =  readtable('mm41ext3.csv');
    d{5,3} =  readtable('mm46ext3.csv');
    d{5,4} =  readtable('mm48ext3.csv');
    d{5,5} =  readtable('mm53ext3.csv');
    d{5,6} =  readtable('mm68ext3.csv');
    d{5,7} =  readtable('mm70ext3.csv');
    d{5,8} =  readtable('mm59ext3.csv');
    d{5,9} =  readtable('mm64ext3.csv');
    d{5,10} =  readtable('mm65ext3.csv');
    d{5,11} =  readtable('mm69ext3.csv');
    d{5,12} =  readtable('mm61ext3.csv');

    d2{1,1} =  readtable('mm45hab.csv');
    d2{1,2} =  readtable('mm49habA.csv');

    d2{2,1} =  readtable('mm45acq.csv');
    d2{2,2} =  readtable('mm49habB.csv');

    d2{3,1} =  readtable('mm45ext1a.csv');
    d2{3,2} =  readtable('mm49acq.csv');

    d2{4,1} =  readtable('mm45ext1b.csv');
    d2{4,2} =  readtable('mm49ext1.csv');

    d2{5,1} =  readtable('mm45ext2.csv');
    d2{5,2} =  readtable('mm49ext2.csv');

    d2{6,1} =  readtable('mm45ext3.csv');
    d2{6,2} =  readtable('mm49ext3.csv');

    %Create vectors containing the stimulus delivery times for each animal &
    %session.

    %Habituation and Acquisition = [150; 249; 364; 471; 572; 678; 794; 908]
    %Extinction sessions = [150; 264; 361; 471; 573; 682];
    %if you started the session late or early with respect to the miniscope
    %recording (that is 30s ± error), add that error into your stimulus times
    %for that session


    %Acq times
    times{1,1} = [150; 249; 364; 471; 572; 678; 794; 908]; %MM39
    times{1,2} = [150; 249; 364; 471; 572; 678; 794; 908];  %MM41
    times{1,3} = [150; 249; 364; 471; 572; 678; 794; 908];  %MM46
    times{1,4} = [160; 259; 374; 481; 582; 688; 804; 918];  %MM48
    times{1,5} = [150; 249; 367; 471; 572; 678; 794; 908];  %MM45
    times{1,6} = [150; 249; 364; 471; 572; 678; 794; 908];  %MM49
    times{1,7} = [150; 249; 364; 471; 572; 678; 794; 908];  %MM48
    times{1,8} = [150; 249; 367; 471; 572; 678; 794; 908];  %MM45
    times{1,9} = [150; 249; 364; 471; 572; 678; 794; 908];  %MM49
    times{1,10} = [150; 249; 364; 471; 572; 678; 794; 908];  %MM49
    times{1,11} = [150; 249; 364; 471; 572; 678; 794; 908];  %MM48
    times{1,12} = [150; 249; 367; 471; 572; 678; 794; 908];  %MM45
    times{1,13} = [150; 249; 364; 471; 572; 678; 794; 908];  %MM49
    times{1,14} = [150; 249; 364; 471; 572; 678; 794; 908];  %MM49

    times{2,1} = [150; 249; 364; 471; 572; 678; 794; 908];
    times{2,2} = [150; 249; 364; 471; 572; 678; 794; 908];
    times{2,3} = [150; 249; 364; 471; 572; 678; 794; 908];
    times{2,4} = [150; 249; 364; 471; 572; 678; 794; 908];
    times{2,5} = [150; 249; 364; 471; 572; 678; 794; 908];
    times{2,6} = [150; 249; 364; 471; 572; 678; 794; 908];
    times{2,7} = [150; 249; 364; 471; 572; 678; 794; 908];  
    times{2,8} = [150; 249; 367; 471; 572; 678; 794; 908];  
    times{2,9} = [150; 249; 364; 471; 572; 678; 794; 908];  
    times{2,10} = [150; 249; 364; 471; 572; 678; 794; 908];  %MM49
    times{2,11} = [150; 249; 364; 471; 572; 678; 794; 908];  %MM48
    times{2,12} = [150; 249; 367; 471; 572; 678; 794; 908];  %MM45
    times{2,13} = [150; 249; 364; 471; 572; 678; 794; 908];  %MM49
    times{2,14} = [150; 249; 364; 471; 572; 678; 794; 908];

    times{3,1} = [150; 264; 361; 471; 573; 682];  
    times{3,2} = [150; 264; 361; 471; 573; 682]; 
    times{3,3} = [150; 264; 361; 471; 573; 682]; 
    times{3,4} = [150; 264; 361; 471; 573; 682];  
    times{3,5} = [150; 343.8; 457.8; 554.8; 664.8; 764.8];
    times{3,6} = [150; 264; 361; 471; 573; 682]; 
    times{3,7} = [150; 264; 361; 471; 573; 682]; 
    times{3,8} = [150; 264; 361; 471; 573; 682]; 
    times{3,9} = [150; 264; 361; 471; 573; 682];
    times{3,10} = [150; 264; 361; 471; 573; 682]; 
    times{3,11} = [150; 264; 361; 471; 573; 682]; 
    times{3,12} = [150; 264; 361; 471; 573; 682]; 
    times{3,13} = [150; 264; 361; 471; 573; 682];
    times{3,14} = [150; 264; 361; 471; 573; 682];

    times{4,1} = [150; 264; 361; 471; 573; 682]; 
    times{4,2} = [150; 264; 361; 471; 573; 682]; 
    times{4,3} = [150; 264; 361; 471; 573; 682]; 
    times{4,4} = [150; 264; 361; 471; 573; 682]; 
    times{4,5} = [150; 264; 361; 471; 573; 682]; 
    times{4,6} = [150; 264; 361; 471; 573; 682]; 
    times{4,7} = [150; 264; 361; 471; 573; 682]; 
    times{4,8} = [150; 264; 361; 471; 573; 682]; 
    times{4,9} = [150; 264; 361; 471; 573; 682];
    times{4,10} = [150; 264; 361; 471; 573; 682]; 
    times{4,11} = [150; 264; 361; 471; 573; 682]; 
    times{4,12} = [150; 264; 361; 471; 573; 682]; 
    times{4,13} = [150; 264; 361; 471; 573; 682];
    times{4,14} = [150; 264; 361; 471; 573; 682];

    times{5,1} = [150; 264; 361; 471; 573; 682];
    times{5,2} = [150; 264; 361; 471; 573; 682]; 
    times{5,3} = [150; 264; 361; 471; 573; 682]; 
    times{5,4} = [150; 264; 361; 471; 573; 682]; 
    times{5,5} = [150; 264; 361; 471; 573; 682];
    times{5,6} = [150; 264; 361; 471; 573; 682];
    times{5,7} = [150; 264; 361; 471; 573; 682]; 
    times{5,8} = [150; 264; 361; 471; 573; 682]; 
    times{5,9} = [150; 264; 361; 471; 573; 682];
    times{5,10} = [150; 264; 361; 471; 573; 682]; 
    times{5,11} = [150; 264; 361; 471; 573; 682]; 
    times{5,12} = [150; 264; 361; 471; 573; 682]; 
    times{5,13} = [150; 264; 361; 471; 573; 682];
    times{5,14} = [150; 264; 361; 471; 573; 682];

    L{1,1}= table2array(readtable('mm39long.csv'));
    L{2,1}= table2array(readtable('mm41long.csv'));
    L{3,1}= table2array(readtable('mm46long.csv'));
    L{4,1}= table2array(readtable('mm48long.csv'));
    L{5,1}= table2array(readtable('mm45long.csv'));
    L{6,1}= table2array(readtable('mm49long.csv'));
    L{7,1}= table2array(readtable('mm53long.csv'));
    L{8,1}= table2array(readtable('mm68long.csv'));
    L{9,1}= table2array(readtable('mm70long.csv'));
    L{10,1}= table2array(readtable('mm59long.csv'));
    L{11,1}= table2array(readtable('mm64long.csv'));
    L{12,1}= table2array(readtable('mm65long.csv'));
    L{13,1}= table2array(readtable('mm69long.csv'));
    L{14,1}= table2array(readtable('mm61long.csv'));

    f{1,1} =  table2array(readtable('fm39_hab.csv'));
    f{1,2} =  table2array(readtable('fm41_hab.csv'));
    f{1,3} =  table2array(readtable('fm46_hab.csv'));
    f{1,4} =  table2array(readtable('fm48_hab.csv'));
    f{1,5} =  table2array(readtable('fm45_hab.csv'));
    f{1,6} =  table2array(readtable('fm49_hab.csv'));
    f{1,7} =  table2array(readtable('M53_hab.csv'));
    f{1,8} =  table2array(readtable('M68_hab.csv'));
    f{1,9} =  table2array(readtable('M70_hab.csv'));
    f{1,10} =  table2array(readtable('M59_hab.csv'));
    f{1,11} =  table2array(readtable('M64_hab.csv'));
    f{1,12} =  table2array(readtable('M65_hab.csv'));
    f{1,13} =  table2array(readtable('M69_hab.csv'));
    f{1,14} =  table2array(readtable('M61_hab.csv'));

    f{2,1} =  table2array(readtable('fm39_acq.csv'));
    f{2,2} =  table2array(readtable('fm41_acq.csv'));
    f{2,3} =  table2array(readtable('fm46_acq.csv'));
    f{2,4} =  table2array(readtable('fm48_acq.csv'));
    f{2,5} =  table2array(readtable('fm45_acq.csv'));
    f{2,6} =  table2array(readtable('fm49_acq.csv'));
    f{2,7} =  table2array(readtable('M53_acq.csv'));
    f{2,8} =  table2array(readtable('M68_acq.csv'));
    f{2,9} =  table2array(readtable('M70_acq.csv'));
    f{2,10} =  table2array(readtable('M59_acq.csv'));
    f{2,11} =  table2array(readtable('M64_acq.csv'));
    f{2,12} =  table2array(readtable('M65_acq.csv'));
    f{2,13} =  table2array(readtable('M69_acq.csv'));
    f{2,14} =  table2array(readtable('M61_acq.csv'));

    f{3,1} =  table2array(readtable('fm39_ext1.csv'));
    f{3,2} =  table2array(readtable('fm41_ext1.csv'));
    f{3,3} =  table2array(readtable('fm46_ext1.csv'));
    f{3,4} =  table2array(readtable('fm48_ext1.csv'));
    f{3,5} =  [table2array(readtable('fm45_ext1a.csv')); table2array(readtable('fm45_ext1b.csv'))];
    f{3,6} =  table2array(readtable('fm49_ext1.csv'));
    f{3,7} =  table2array(readtable('M53_ext1.csv'));
    f{3,8} =  table2array(readtable('M68_ext1.csv'));
    f{3,9} =  table2array(readtable('M70_ext1.csv'));
    f{3,10} =  table2array(readtable('M59_ext1.csv'));
    f{3,11} =  table2array(readtable('M64_ext1.csv'));
    f{3,12} =  table2array(readtable('M65_ext1.csv'));
    f{3,13} =  table2array(readtable('M69_ext1.csv'));
    f{3,14} =  table2array(readtable('M61_ext1.csv'));

    f{4,1} =  table2array(readtable('fm39_ext2.csv'));
    f{4,2} =  table2array(readtable('fm41_ext2.csv'));
    f{4,3} =  table2array(readtable('fm46_ext2.csv'));
    f{4,4} =  table2array(readtable('fm48_ext2.csv'));
    f{4,5} =  table2array(readtable('fm45_ext2.csv'));
    f{4,6} =  table2array(readtable('fm49_ext2.csv'));
    f{4,7} =  table2array(readtable('M53_ext2.csv'));
    f{4,8} =  table2array(readtable('M68_ext2.csv'));
    f{4,9} =  table2array(readtable('M70_ext2.csv'));
    f{4,10} =  table2array(readtable('M59_ext2.csv'));
    f{4,11} =  table2array(readtable('M64_ext2.csv'));
    f{4,12} =  table2array(readtable('M65_ext2.csv'));
    f{4,13} =  table2array(readtable('M69_ext2.csv'));
    f{4,14} =  table2array(readtable('M61_ext2.csv'));

    f{5,1} =  table2array(readtable('fm39_ext3.csv'));
    f{5,2} =  table2array(readtable('fm41_ext3.csv'));
    f{5,3} =  table2array(readtable('fm46_ext3.csv'));
    f{5,4} =  table2array(readtable('fm48_ext3.csv'));
    f{5,5} =  table2array(readtable('fm45_ext3.csv'));
    f{5,6} =  table2array(readtable('fm49_ext3.csv'));
    f{5,7} =  table2array(readtable('M53_ext3.csv'));
    f{5,8} =  table2array(readtable('M68_ext3.csv'));
    f{5,9} =  table2array(readtable('M70_ext3.csv'));
    f{5,10} =  table2array(readtable('M59_ext3.csv'));
    f{5,11} =  table2array(readtable('M64_ext3.csv'));
    f{5,12} =  table2array(readtable('M65_ext3.csv'));
    f{5,13} =  table2array(readtable('M69_ext3.csv'));
    f{5,14} =  table2array(readtable('M61_ext3.csv'));

    props{1,1} = readtable('mm39ext3-props.csv');
    props{1,2} = readtable('mm41ext3-props.csv');
    props{1,3} = readtable('mm46ext3-props.csv');
    props{1,4} = readtable('mm48ext3-props.csv');
    props{1,5} = readtable('mm45ext3-props.csv');
    props{1,6} = readtable('mm49ext3-props.csv');
    props{1,7} = readtable('mm53ext3-props.csv');
    props{1,8} = readtable('mm68ext3-props.csv');
    props{1,9} = readtable('mm70ext3-props.csv');
    props{1,10} = readtable('mm59ext3-props.csv');
    props{1,11} = readtable('mm64ext3-props.csv');
    props{1,12} = readtable('mm65ext3-props.csv');
    props{1,13} = readtable('mm69ext3-props.csv');
    props{1,14} = readtable('mm61ext3-props.csv');

    %responders and non-responders were determined in prism by identifying who
    %extinguished >50% of freezing
    responders = [1 3 7 8 10 11 12];
    nonresponders = [2 4 5 6 9 13 14];
else
    d{1,1} = readtable('M23_Hab_dFoFL.csv');
    d{1,2} = readtable('M27_Hab_dFoFL.csv');
    d2{1,1} = readtable('M28_Hab_dFoFa.csv');
    d2{2,1} = readtable('M28_Hab_dFoFb.csv');
    
    d{2,1} = readtable('M23_Acq_dFoFL.csv');
    d{2,2} = readtable('M27_Acq_dFoFL.csv');
    d2{3,1} = readtable('M28_Acq_dFoFL.csv');

    d{3,1} = readtable('M23_Ext1_dFoFL.csv');
    d{3,2} = readtable('M27_Ext1_dFoFL.csv');
    d2{4,1} = readtable('M28_Ext1_dFoFL.csv');

    d{4,1} = readtable('M23_Ext2_dFoFL.csv');
    d{4,2} = readtable('M27_Ext2_dFoFL.csv');
    d2{5,1} = readtable('M28_Ext2_dFoFL.csv');

    d{5,1} = readtable('M23_Ext3_dFoFL.csv');
    d{5,2} = readtable('M27_Ext3_dFoFL.csv');
    d2{6,1} = readtable('M28_Ext3_dFoFL.csv');

    times{1,1} = [150; 249; 364; 471; 572; 678; 794; 908];
    times{1,2} = [150; 249; 364; 471; 572; 678; 794; 908];
    times{1,3} = [150; 249; 364; 471; 572; 678; 794; 908];

    times{2,1} = [150; 249; 364; 471; 572; 678; 794; 908];
    times{2,2} = [150; 249; 364; 471; 572; 678; 794; 908];
    times{2,3} = [150; 249; 364; 471; 572; 678; 794; 908];

    times{3,1} = [150; 264; 361; 471; 573; 682]; 
    times{3,2} = [150; 264; 361; 471; 573; 682]; 
    times{3,3} = [150; 264; 361; 471; 573; 682]; 

    times{4,1} = [150; 264; 361; 471; 573; 682]; 
    times{4,2} = [150; 264; 361; 471; 573; 682]; 
    times{4,3} = [150; 264; 361; 471; 573; 682]; 


    times{5,1} = [150; 264; 361; 471; 573; 682]; 
    times{5,2} = [150; 265; 366; 472; 585; 687]; 
    times{5,3} = [150; 264; 361; 471; 573; 682]; 
    
    f{1,1} = table2array(readtable('mini23_hab.csv'));
    f{1,2} = table2array(readtable('mini27_hab.csv'));
    f{1,3} = table2array(readtable('mini28_hab.csv'));
    
    f{2,1} = table2array(readtable('mini23_acq.csv'));
    f{2,2} = table2array(readtable('mini27_acq.csv'));
    f{2,3} = table2array(readtable('mini28_acq.csv'));
    
    f{3,1} = table2array(readtable('mini23_ext1.csv'));
    f{3,2} = table2array(readtable('mini27_ext1.csv'));
    f{3,3} = table2array(readtable('mini28_ext1.csv'));
    
    f{4,1} = table2array(readtable('mini23_ext2.csv'));
    f{4,2} = table2array(readtable('mini27_ext2.csv'));
    f{4,3} = table2array(readtable('mini28_ext2.csv'));
    
    f{5,1} = table2array(readtable('mini23_ext3.csv'));
    f{5,2} = table2array(readtable('mini27_ext3.csv'));
    f{5,3} = table2array(readtable('mini28_ext3.csv'));
    
    L{1,1}= table2array(readtable('mm23longitude3.csv'));
    L{2,1}= table2array(readtable('mm27longitude.csv'));
    L{3,1}= table2array(readtable('mm28longitude.csv'));

end


%% 2. Extract traces of longitudinal cells and align to stimuli

%CUSTOM FUNCTIONS REQUIRED:
%   crossSeshAuto
%   aligntraces

%variables
start = 10; %how much time back from stimulus times to begin collecting data. (seconds)
fin = 69;   % how much time after stimulus start to stop collecting data. (seconds)
dt = 20;    % sampling rate of recording (Hz)
nSesh = 5;  % number of recordings
trialLength = start+fin+1; %length of trial (seconds)
nTrials = 34; %number of total trials


%outputs
poolMat = []; %initialize matrix to store vertically concatenated neural activity across trials and sessions, summed within seconds
nCells = zeros(nAnimals,1);  %initialize array to store number of cells per animal
lr = cell(nAnimals,nSesh);   %initialize cell to store original activity matrices downsampled to 15Hz, for freezing decoding later
tensors = cell(nAnimals,1);  %initialize cell to store same information as poolMat but arranged in tensors of t seconds X c cells X T trials. Export for use in tensortools python kit
spatial = cell(nAnimals,1);  %initialize cell to store spatial coordinates of longitudinally registered cells
cR = []; %initialize storage vector for cell indices belonging to psilocybin responders
cNR = []; %initialize storage vector for cell indices belonging to psilocybin responders

for n=1:nAnimals
    
    %initialize temp storage
    reRegistered={};
    dA = {};
    timesA = {};
    
    %define special cases for broken recordings
    if dataset == 1
        if n==5
            nSesh = 6;
            dA = d2(:,1);
        elseif n==6
            nSesh = 6;
            dA = d2(:,2);
        elseif n>6
            nSesh = 5;
            for m=1:nSesh
                dA{m,1} = d{m,n-2};
            end
        else
            nSesh = 5;
            for m=1:nSesh
                dA{m,1} = d{m,n};
            end
        end
    elseif dataset == 2
        if n==3
            nSesh =6;
            dA = d2(:,1);
        else
            for m=1:nSesh
                dA{m,1} = d{m,n};
            end
        end
    else
        nSesh = 5;
        for m=1:nSesh
            dA{m,1} = d{m,n};
        end
    end
    
    %extract spatial data 
    if dataset < 2
        locs = props{1,n};
    else
        locs = [];
    end
    
   %longitudinally register and spatially locate cells across 5 sessions
    [longRegistered, coordinates]= crossSeshAuto(L{n,1},dA,locs,nSesh);
    
   %store spatial data
    spatial{n,1} = coordinates;
   
   %stitch information from broken recordings in 2 animals
    if dataset == 1 & n == 5
        reRegistered{1,1} = [longRegistered{1,1}];
        reRegistered{2,1} = [longRegistered{2,1}];
        reRegistered{3,1} = [longRegistered{3,1}; longRegistered{4,1}];
        reRegistered{4,1} = [longRegistered{5,1}];
        reRegistered{5,1} = [longRegistered{6,1}];
        longRegistered = reRegistered;
    elseif dataset == 1 & n == 6
        x=longRegistered{1,1};
        y=longRegistered{2,1};
        reRegistered{1,1} = [x; x(end-16.8*dt/2+1:end,:); y(1:16.8*dt/2,:); y];
        reRegistered{2,1} = [longRegistered{3,1}];
        reRegistered{3,1} = [longRegistered{4,1}];
        reRegistered{4,1} = [longRegistered{5,1}];
        reRegistered{5,1} = [longRegistered{6,1}];
        longRegistered = reRegistered;
    elseif dataset == 2 & n == 3
        reRegistered{1,1} = [longRegistered{1,1}; longRegistered{2,1}];
        reRegistered{2,1} = [longRegistered{3,1}];
        reRegistered{3,1} = [longRegistered{4,1}];
        reRegistered{4,1} = [longRegistered{5,1}];
        reRegistered{5,1} = [longRegistered{6,1}];
        longRegistered = reRegistered;
    end
    
    %loop through sessions, downsample registered traces to 15Hz and store
    nSesh = 5;
     for m=1:nSesh
         timesA{m,1} = times{m,n};
        x=resample(longRegistered{m,1},3,4);
        lr{n,m} = x;
     end
    
    %count and store number of cells
    dims = size(longRegistered{1,1});
    nCells(n,1) = dims(2);
    
    %align traces to stimulus times and store 1Hz data at the predefined
    %trial windows
    alignedTraces = aligntraces(longRegistered,timesA,nSesh,dt,start,fin);
    
    %concatenate aligned traces
    allAligned = [alignedTraces{1,1}; alignedTraces{2,1}; alignedTraces{3,1}; alignedTraces{4,1}; alignedTraces{5,1}];
    
    %store session concatenated traces
    poolMat = [poolMat allAligned];
    
    %transform into t x c x T tensor and store
    ten = reshape(allAligned,[trialLength, nCells(n), nTrials]);
    tensors{n,1} = ten;
    
    if dataset == 1
        if ismember(n,responders)
            if n==1
                cR = [cR 1:nCells(n)];
            else
                cR = [cR [1+sum(nCells(1:n-1)):sum(nCells(n))]];
            end
        else
            cNR = [cNR [1+sum(nCells(1:n-1)):sum(nCells(n))]];
        end
    end
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

%define trial event times
bl = [1:10];  %seconds
stim = [11:35]; %seconds
trace = [36:55]; %seconds
shock = [56:58]; %seconds

%assign trials to sessions
hab = [1:8]; %trials
acq = [9:16]; %trials
ext1 = [17:22]; %trials
ext2 = [23:28]; %trials
ext3 = [29:34]; %trials
sessions = {hab, acq, ext1, ext2, ext3};

%define time window prior to stimulus for statistical comparison
timeBack = 10; %seconds
nTrials = 34;

%initialize storage cell average session activity of each animal by session 
avPB = cell(nAnimals, nSesh); 

%loop through animals
for n=1:nAnimals
     
    %call t x c x T data from tensors 
    data = tensors{n};
    data = reshape(data,trialLength*nTrials,nCells(n));
    
    for m=1:nSesh
        
        trials = sessions{m};
        
        %calculate average trial activity (% change from baseline) in each
        %session and store to pooled matrix across animals for plotting
        if m == 1
            act = percBaselineAvg(data(1:max(sessions{m})*trialLength,:),trialLength,timeBack,length(trials));
        else
            act = percBaselineAvg(data(1 + max(sessions{m-1})*trialLength:max(sessions{m})*trialLength,:),trialLength,timeBack,length(trials));
        end
        
        %store
        avPB{n,m} = act;

    end
end

%pool cells from all animals within sessions
poolAvgAct = cell(nSesh,1);

for m=1:nSesh
    
    %initialize temp storage
    pool = [];
    
    %store all average traces in a session from every animal
    for n=1:nAnimals
        pool = [pool avPB{n,m}];
    end
    
    %store pooled traces by session
    poolAvgAct{m,1} = pool;
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

%% 5. Identify cells that are upregulated or downregulated in response to stimuli in a given session according to their average trace for fig 2G-J

%CUSTOM FUNCTIONS
%isSig
%barWithError

%initialize storage cells for indices of responsive neurons in each animal
%in each session
stimCells = cell(nAnimals,nSesh);
traceCells = cell(nAnimals,nSesh);
shockCells = cell(nAnimals,nSesh);

%initialize storage matrices for indices of up- and downregulated neurons in each animal
%in each session
fracStimUp = zeros(nSesh, nAnimals);
fracStimDo = zeros(nSesh, nAnimals);
fracTraceUp = zeros(nSesh, nAnimals);
fracTraceDo = zeros(nSesh, nAnimals);
fracShockUp = zeros(nSesh, nAnimals);
fracShockDo = zeros(nSesh, nAnimals);

%initialize storage matrices for indices of up- and downregulated neurons in each animal
%in each session
fracStimTrace = zeros(nSesh, nAnimals);

%loop through animals and sessions
for n=1:nAnimals
    for m=1:nSesh
        %test if activity during stim diff from baseline (wilcoxon), store
        %cell indices
        stimCells{n,m} = isSig(avPB{n,m},bl,stim);
        
        %test if activity during trace diff from baseline (wilcoxon)
        traceCells{n,m} = isSig(avPB{n,m},bl,trace);
        
        %test if activity during shock diff from baseline (wilcoxon)
        shockCells{n,m} = isSig(avPB{n,m},bl,shock);
        
        
        
        
        %calculate fractions of modulated cells - stim
        data = stimCells{n,m};
        
        up1 = data{1};
        fracStimUp(m,n) = length(up1)/nCells(n);
        down1 = data{2};
        fracStimDo(m,n) = length(down1)/nCells(n);
        
        %calculate fractions of modulated cells - trace
        data = traceCells{n,m};
        
        up = data{1};
        fracTraceUp(m,n) = length(up)/nCells(n);
        down = data{2};
        fracTraceDo(m,n) = length(down)/nCells(n);
        
        %measure overlap of stim and trace cells
        fracStimTrace(m,n) = mean(ismember(up1,up));
        
        %calculate fractions of modulated cells - shock
        data = shockCells{n,m};
        
        up = data{1};
        fracShockUp(m,n) = length(up)/nCells(n);
        down = data{2};
        fracShockDo(m,n) = length(down)/nCells(n);
        
        
    end
end

%get rid of nans if a mouse had zero stim responsive neurons in denominator
fracStimTrace(isnan(fracStimTrace)) = 0;

if dataset == 1
    %responders
    figure
    subplot(221)
    barWithError(fracStimUp(:,responders)')
    ylabel('Fraction of Cells')
    xticklabels(titles)
    title('Stimulus-Responsive Cells')
    
    subplot(222)
    barWithError(fracTraceUp(:,responders)')
    ylabel('Fraction of Cells')
    xticklabels(titles)
    title('Trace-Responsive Cells')
    
    subplot(223)
    barWithError(fracShockUp(:,responders)')
    ylabel('Fraction of Cells')
    xticklabels(titles)
    title('Shock-Responsive Cells')
    
    subplot(224)
    barWithError(fracStimTrace(:,responders)')
    ylabel('Fraction of Cells')
    xticklabels(titles)
    title('Stimulus-and-Trace-Responsive Cells')
    
    % %nonresponders
    figure
    subplot(221)
    barWithError(fracStimUp(:,nonresponders)')
    ylabel('Fraction of Cells')
    xticklabels(titles)
    title('Stimulus-Responsive Cells')
    
    subplot(222)
    barWithError(fracTraceUp(:,nonresponders)')
    ylabel('Fraction of Cells')
    xticklabels(titles)
    title('Trace-Responsive Cells')
    
    subplot(223)
    barWithError(fracShockUp(:,nonresponders)')
    ylabel('Fraction of Cells')
    xticklabels(titles)
    title('Shock-Responsive Cells')
    
    subplot(224)
    barWithError(fracStimTrace(:,nonresponders)')
    ylabel('Fraction of Cells')
    xticklabels(titles)
    title('Stimulus-and-Trace-Responsive Cells')


else

    figure
    subplot(221)
    barWithError(fracStimUp')
    ylabel('Fraction of Cells')
    xticklabels(titles)
    title('Stimulus-Responsive Cells')

    subplot(222)
    barWithError(fracTraceUp')
    ylabel('Fraction of Cells')
    xticklabels(titles)
    title('Trace-Responsive Cells')

    subplot(223)
    barWithError(fracShockUp')
    ylabel('Fraction of Cells')
    xticklabels(titles)
    title('Shock-Responsive Cells')

    subplot(224)
    barWithError(fracStimTrace')
    ylabel('Fraction of Cells')
    xticklabels(titles)
    title('Stimulus-and-Trace-Responsive Cells')

end
%% 6. Measure the number of sessions cells tend to maintain stimulus responsiveness. For fig 2L-K

%CUSTOM FUNCTIONS
%barWithError

%initialize storage matrices for number of session neurons are responsive
fracPersistStim = zeros(nSesh,nAnimals);
fracPersistTrace = zeros(nSesh,nAnimals);

whichSeshS = zeros(nSesh-1,nAnimals);
whichSeshT = zeros(nSesh-1,nAnimals);

%loop through animals
for n=1:nAnimals
    
    %initialize temp storage vector
    persistStim = zeros(nCells(n),1);
    persistTrace = zeros(nCells(n),1);
    
    wS = zeros(5,nCells(n));
    wT = zeros(5,nCells(n));
    
    %loop through cells
    for c = 1:nCells(n)
        
        %initialize counters for stimulus and trace
        numSessions = 0;
        numSessionsT = 0;
        
        %loop through sessions
         for m = 1:nSesh
             
             %call indices of stim and trace responsive neurons
            data = stimCells{n,m};
            data1 = traceCells{n,m};
                
                %if this  cell is in the list of stim or trace responsive
                %neurons in that session, add one to counter
               if ismember(c,data{1})
                   numSessions = numSessions+1;
                   wS(numSessions,c) = m;
               end

               if ismember(c,data1{1})
                   numSessionsT = numSessionsT+1;
                   wT(numSessionsT,c) = m;
               end
         end
     
     %save final counter values
     persistStim(c,1) = numSessions;
     persistTrace(c,1) = numSessionsT;
     
    end
     
    %get rid of cells that were never responsive
    persistStim(persistStim==0) = [];
    persistTrace(persistTrace==0) = [];
    
    for m=1:nSesh-1
        combo = [m; m+1];
        for c=1:length(wS(1,:))
            if sum(wS(1:2,c) == combo) == 2
                whichSeshS(m,n) = whichSeshS(m,n)+1;
            end
            if sum(wT(1:2,c) == combo) == 2
                whichSeshT(m,n) = whichSeshT(m,n)+1;
            end
        end
    end
    
    %store the fraction of stim and trace responsive neurons that last for
    %1, 2, 3, 4, and 5 sessions
    for m = 1:5
        fracPersistStim(m,n) = length(find(persistStim==m))/length(persistStim);
        fracPersistTrace(m,n) = length(find(persistTrace==m))/length(persistTrace);
    end
    
    
end

whichSeshS = whichSeshS./nCells(n);
whichSeshT = whichSeshT./nCells(n);

%get rid of nans if an animal had no trace active neurons
fracPersistTrace(isnan(fracPersistTrace)) = 0;

if dataset == 1
    %responders
    figure
    subplot(121)
    barWithError(fracPersistStim(:,responders)')
    ylabel('Fraction of Cells')
    title('Persistency of Stimulus-Responsive Cells')
    
    subplot(122)
    barWithError(fracPersistTrace(:,responders)')
    ylabel('Fraction of Cells')
    title('Persistency of Trace-Responsive Cells')

    %nonresponders
    figure
    subplot(121)
    barWithError(fracPersistStim(:,nonresponders)')
    ylabel('Fraction of Cells')
    title('Persistency of Stimulus-Responsive Cells')
    
    subplot(122)
    barWithError(fracPersistTrace(:,nonresponders)')
    ylabel('Fraction of Cells')
    title('Persistency of Trace-Responsive Cells')
    
    figure
    subplot(121)
    errorbar(1:4,mean(whichSeshS(:,responders),2),std(whichSeshS(:,responders)')./sqrt(7))
    hold on
    errorbar(1:4,mean(whichSeshT(:,responders),2),std(whichSeshT(:,responders)')./sqrt(7))
    ylabel('Fraction of Cells')
    legend({'Stim','Trace'})
    subplot(221)
    errorbar(1:4,mean(whichSeshS(:,nonresponders),2),std(whichSeshS(:,nonresponders)')./sqrt(7))
    hold on
    errorbar(1:4,mean(whichSeshT(:,nonresponders),2),std(whichSeshT(:,nonresponders)')./sqrt(7))
    ylabel('Fraction of Cells')
    legend({'Stim','Trace'})
else
    figure
    subplot(121)
    barWithError(fracPersistStim')
    ylabel('Fraction of Cells')
    title('Persistency of Stimulus-Responsive Cells')

    subplot(122)
    barWithError(fracPersistTrace')
    ylabel('Fraction of Cells')
    title('Persistency of Trace-Responsive Cells')
end


%% 7. Calculate freezing encoding in each neuron in each session and compare within animals across groups for fig 2M-N

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
        
        freezing = f{h,g};      %loop through freezing
        freezing(freezing==100) = 1;  %turn into binary
        
        if length(freezing)>length(ac)  %make matrices the same length
            freezing = freezing(1:length(ac(:,1)));
        else
            ac = ac(1:length(freezing),:);
        end
        
        %initialize storage vectors for auROC and cell indices for auROC > .6
        AUCs = [];
        posCells = [];
        
        %loop through neurons
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

%plot freezing encoding stats
if dataset == 1
    %responders
    figure
    subplot(121)
    barWithError(AUCsMean(:,responders)')
    ylabel('Mean freezing encoding')
    title('Freezing encoding in RSC neurons (Fig. 2M)')
    xticklabels(titles)
    
    subplot(122)
    barWithError(freezeCells(:,responders)')
    ylabel('Fraction of Cells')
    title('Freezing encoding neurons')
    xticklabels(titles)
    
    %non-responders
    figure
    subplot(121)
    barWithError(AUCsMean(:,nonresponders)')
    ylabel('Mean freezing encoding')
    title('Freezing encoding in RSC neurons (Fig. 2M)')
    xticklabels(titles)
    
    subplot(122)
    barWithError(freezeCells(:,responders)')
    ylabel('Fraction of Cells')
    title('Freezing encoding neurons')
    xticklabels(titles)
    
else
    figure
    subplot(121)
    barWithError(AUCsMean')
    ylabel('Mean freezing encoding')
    title('Freezing encoding in RSC neurons (Fig. 2M)')
    xticklabels(titles)

    subplot(122)
    barWithError(freezeCells')
    ylabel('Fraction of Cells')
    title('Freezing encoding neurons')
    xticklabels(titles)
end



%plot representative freezing encoding cells - image in paper is dataset 0
%mouse 7 session 5
timeWin = [1000:1500];

figure
[B,I] = sort(AUCs(posCells));
figure
subplot(211)
heatmap(freezing(timeWin)')
ylabel('Freezing')
colormap(parula)
subplot(212)
heatmap(zscore(ac(timeWin,posCells(I)))','Colorlimits',[-2 2])
ylabel('Freezing encoding cells')
colormap(jet)

%% 8. Load & plot TCA models, identifying session-dominant components and correlating to freezing for fig. 3A,D-F

%calculate tensors from Williams et al. TCA code (https://github.com/ahwillia/tensortools)
% load tensors 
start = 10; %how much time back from stimulus times to begin collecting data. (seconds)
fin = 49;   % how much time after stimulus start to stop collecting data. (seconds)
dt = 20;    % sampling rate of recording (Hz)
nSesh = 5;  % number of recordings
trialLength = start+fin+1; %length of trial (seconds)
nTrials = 34; %number of total trials

hab = [1:8]; %trials
acq = [9:16]; %trials
ext1 = [17:22]; %trials
ext2 = [23:28]; %trials
ext3 = [29:34]; %trials

if dataset == 1
    nAnimals = 14;
    tds{1,1} = table2array(readtable('td39.csv'));
    tds{2,1} = table2array(readtable('td41.csv'));
    tds{3,1} = table2array(readtable('td46.csv'));
    tds{4,1} = table2array(readtable('td48.csv'));
    tds{5,1} = table2array(readtable('td45.csv'));
    tds{6,1} = table2array(readtable('td49.csv'));
    tds{7,1} = table2array(readtable('td53.csv'));
    tds{8,1} = table2array(readtable('td68.csv'));
    tds{9,1} = table2array(readtable('td70.csv'));
    tds{10,1} = table2array(readtable('td59.csv'));
    tds{11,1} = table2array(readtable('td64.csv'));
    tds{12,1} = table2array(readtable('td65.csv'));
    tds{13,1} = table2array(readtable('td69.csv'))
    tds{14,1} = table2array(readtable('td61.csv'));
elseif dataset == 0
    tds{1,1} = table2array(readtable('td43.csv'));
    tds{2,1} = table2array(readtable('td51.csv'));
    tds{3,1} = table2array(readtable('td55.csv'));
    tds{4,1} = table2array(readtable('td71.csv'));
    tds{5,1} = table2array(readtable('td56.csv'));
    tds{6,1} = table2array(readtable('td58.csv'));
    tds{7,1} = table2array(readtable('td57.csv'));
elseif dataset == 2
    tds{1,1} = table2array(readtable('td23.csv'));
    tds{2,1} = table2array(readtable('td27.csv'));
    tds{3,1} = table2array(readtable('td28.csv'));
end

%get trial by trial freezing
freezingInTrials = {};
dt = 15;
timesF{1,1} =[120; 219; 334; 441; 542; 648; 764; 878];
timesF{2,1} =[120; 219; 334; 441; 542; 648; 764; 878];
timesF{3,1} = [120; 234; 331; 441; 543; 652]; 
timesF{4,1} = [120; 234; 331; 441; 543; 652]; 
timesF{5,1} = [120; 234; 331; 441; 543; 652]; 

%create second by second % freezing
for n=1:nAnimals
    fr = f(:,n);
    if dataset == 0
        if n==5
            timesf = timesF;
            timesf{3,1} = [120; 120+round(2844/15); 234+round(2844/15); 331+round(2844/15); 441+round(2844/15); 543+round(2844/15)];
        else
            timesf = timesF;
        end
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

%%

nDims = 5;

coefs = zeros(nAnimals,nDims);
sessions = {hab, acq, ext1, ext2, ext3};
%to plot representative TCAs, uncomment the figure file in the for-loop -
%the representative image in the paper is from dataset 1 animal 12
%below
for a=1:nAnimals
    tt = tds{a,1};
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

    %plot representative figure
    [B,I] = sort(neuronFactor(:,1));
%     figure
%     for n=1:nDims
%         neus = neuronFactor(:,n);
% 
%         subplot(nDims,3,1+(n-1)*3)
%         plot(timeFactor(:,n),'Linewidth',3) %for relative time loadings, timeFactor(:,n)/sum(timeFactor(:,n))
%         ylabel(n)
%         ylim([0 max(timeFactor(:,n))+.2])
%         xticks([0 10 35 55])
%         ylim([0 max(timeFactor(:,n)+.2)])
% 
%         subplot(nDims,3,2+(n-1)*3)
%         bar(neus(I),'c')
% 
%         subplot(nDims,3,3+(n-1)*3)
%         scatter([1:nTrials],trialFactor(:,n),50,fr,'filled')
%         hold on
%         xline(8)
%         hold on
%         xline(16)
%         hold on
%         xline(22)
%         hold on
%         xline(28)
%
%
%     end
end

%figure; bar([mean(R2(responders)) mean(R2(nonresponders))]); hold on; errorbar([1:2],[mean(R2(responders)) mean(R2(nonresponders))], [std(R2(responders))./sqrt(7) std(R2(nonresponders))./sqrt(7)])
%figure; bar([mean(ps(responders)) mean(ps(nonresponders))]); hold on; errorbar([1:2],[mean(ps(responders)) mean(ps(nonresponders))], [std(ps(responders))./sqrt(7) std(ps(nonresponders))./sqrt(7)])
% figure; bar([mean(R2)]); hold on; errorbar(1:5,[mean(R2)], [std(R2)./sqrt(7)])
% figure; bar([mean(R2(responders,:))]); hold on; errorbar(1:5,[mean(R2(responders,:))], [std(R2(responders,:))./sqrt(7)])
% figure; bar([mean(R2(nonresponders,:))]); hold on; errorbar(1:5,[mean(R2(nonresponders,:))], [std(R2(nonresponders,:))./sqrt(7)])

coefOrdered = zeros(nAnimals,nDims);
coefs = coefs(:,1:5);
for a=1:nAnimals
for n=1:nDims
    coefOrdered(a,n) = coefs(a,domFacs(a,n));
end
end
% figure
% errorbar(1:5,mean(coefs(responders,:)),std(coefs(responders,:))./sqrt(7))
% hold on
% errorbar(1:5,mean(coefs(nonresponders,:)),std(coefs(nonresponders,:))./sqrt(7))
%% 9. Plot normalized trial factors. for Supp Fig. 3A-E

%initialize trial weight storage vectors
ht = zeros(nTrials, nAnimals);
at = zeros(nTrials, nAnimals);
e1t = zeros(nTrials, nAnimals);
e2t = zeros(nTrials, nAnimals);
e3t = zeros(nTrials, nAnimals);

for n=1:nAnimals
    
    %call trial weights from a given animal
    tWeights = triFac{n};
    
    %store normalized trial weights
    ht(:,n) = tWeights(:,domFacs(n,1))./max(tWeights(:,domFacs(n,1)));
    at(:,n) = tWeights(:,domFacs(n,2))./max(tWeights(:,domFacs(n,2)));
    e1t(:,n) = tWeights(:,domFacs(n,3))./max(tWeights(:,domFacs(n,3)));
    e2t(:,n) = tWeights(:,domFacs(n,4))./max(tWeights(:,domFacs(n,4)));
    e3t(:,n) = tWeights(:,domFacs(n,5))./max(tWeights(:,domFacs(n,5)));
end

ys = [0 1];
if dataset == 1
    figure
    subplot(321)
    plot(ht(:,responders))
    hold on
    errorbar(1:length(ht(:,1)),mean(ht(:,responders),2),std(ht(:,responders)')/sqrt(length(responders)),'Linewidth',3)
    ylim([ys])
    ylabel('Normalized trial weights')
    xlabel('Time')
    title('Hab')
    subplot(322)
    plot(at(:,responders))
    hold on
    errorbar(1:length(at(:,1)),mean(at(:,responders),2),std(at(:,responders)')/sqrt(length(responders)),'Linewidth',3)
    ylim([ys])
    ylabel('Normalized trial weights')
    xlabel('Time')
    title('Acq')
    subplot(323)
    plot(e1t(:,responders))
    hold on
    errorbar(1:length(e1t(:,1)),mean(e1t(:,responders),2),std(e1t(:,responders)')/sqrt(length(responders)),'Linewidth',3)
    ylim([ys])
    ylabel('Normalized trial weights')
    xlabel('Time')
    title('Ext1')
    subplot(324)
    plot(e2t(:,responders))
    hold on
    errorbar(1:length(e2t(:,1)),mean(e2t(:,responders),2),std(e2t(:,responders)')/sqrt(length(responders)),'Linewidth',3)
    ylim([ys])
    ylabel('Normalized trial weights')
    xlabel('Time')
    title('Ext2')
    subplot(325)
    plot(e3t(:,responders))
    hold on
    errorbar(1:length(e3t(:,1)),mean(e3t(:,responders),2),std(e3t(:,responders)')/sqrt(length(responders)),'Linewidth',3)
    ylim([ys])
    ylabel('Normalized trial weights')
    xlabel('Time')
    title('Ext3')

    figure
    subplot(321)
    plot(ht(:,nonresponders))
    hold on
    errorbar(1:length(ht(:,1)),mean(ht(:,nonresponders),2),std(ht(:,nonresponders)')/sqrt(length(responders)),'Linewidth',3)
    ylim([ys])
    ylabel('Normalized trial weights')
    xlabel('Time')
    title('Hab')
    subplot(322)
    plot(at(:,nonresponders))
    hold on
    errorbar(1:length(at(:,1)),mean(at(:,nonresponders),2),std(at(:,nonresponders)')/sqrt(length(responders)),'Linewidth',3)
    ylim([ys])
    ylabel('Normalized trial weights')
    xlabel('Time')
    title('Acq')
    subplot(323)
    plot(e1t(:,nonresponders))
    hold on
    errorbar(1:length(e1t(:,1)),mean(e1t(:,nonresponders),2),std(e1t(:,nonresponders)')/sqrt(length(responders)),'Linewidth',3)
    ylim([ys])
    ylabel('Normalized trial weights')
    xlabel('Time')
    title('Ext1')
    subplot(324)
    plot(e2t(:,nonresponders))
    hold on
    errorbar(1:length(e2t(:,1)),mean(e2t(:,nonresponders),2),std(e2t(:,nonresponders)')/sqrt(length(responders)),'Linewidth',3)
    ylim([ys])
    ylabel('Normalized trial weights')
    xlabel('Time')
    title('Ext2')
    subplot(325)
    plot(e3t(:,nonresponders))
    hold on
    errorbar(1:length(e3t(:,1)),mean(e3t(:,nonresponders),2),std(e3t(:,nonresponders)')/sqrt(length(responders)),'Linewidth',3)
    ylim([ys])
    ylabel('Normalized trial weights')
    xlabel('Time')
    title('Ext3')

else
    figure
    subplot(321)
    plot(ht)
    hold on
    errorbar(1:length(ht(:,1)),mean(ht,2),std(ht')/sqrt(length(responders)),'Linewidth',3)
    ylim([0 1])
    ylabel('Normalized trial weights')
    xlabel('Time')
    title('Hab')
    subplot(322)
    plot(at)
    hold on
    errorbar(1:length(at(:,1)),mean(at,2),std(at')/sqrt(length(responders)),'Linewidth',3)
    ylim([0 1])
    ylabel('Normalized trial weights')
    xlabel('Time')
    title('Acq')
    subplot(323)
    plot(e1t)
    hold on
    errorbar(1:length(e1t(:,1)),mean(e1t,2),std(e1t')/sqrt(length(responders)),'Linewidth',3)
    ylim([0 1])
    ylabel('Normalized trial weights')
    xlabel('Time')
    title('Ext1')
    subplot(324)
    plot(e2t)
    hold on
    errorbar(1:length(e2t(:,1)),mean(e2t,2),std(e2t')/sqrt(length(responders)),'Linewidth',3)
    ylim([0 1])
    ylabel('Normalized temporal weights')
    xlabel('Time')
    title('Ext2')
    subplot(325)
    plot(e3t)
    hold on
    errorbar(1:length(e3t(:,1)),mean(e3t,2),std(e3t')/sqrt(length(responders)),'Linewidth',3)
    ylim([0 1])
    ylabel('Normalized trial weights')
    xlabel('Time')
    title('Ext3')

end

%% 10. Plot normalized temporal factors. for Supp Fig. 3F-J
ht = zeros(trialLength, nAnimals);
at = zeros(trialLength, nAnimals);
e1t = zeros(trialLength, nAnimals);
e2t = zeros(trialLength, nAnimals);
e3t = zeros(trialLength, nAnimals);

for n=1:nAnimals
    
    tWeights = timFac{n};
    ht(:,n) = tWeights(:,domFacs(n,1))./max(tWeights(:,domFacs(n,1)));
    at(:,n) = tWeights(:,domFacs(n,2))./max(tWeights(:,domFacs(n,2)));
    e1t(:,n) = tWeights(:,domFacs(n,3))./max(tWeights(:,domFacs(n,3)));
    e2t(:,n) = tWeights(:,domFacs(n,4))./max(tWeights(:,domFacs(n,4)));
    e3t(:,n) = tWeights(:,domFacs(n,5))./max(tWeights(:,domFacs(n,5)));
end

ys = [0 1];
if dataset == 1
    figure
    subplot(321)
    plot(ht(:,responders))
    hold on
    errorbar(1:length(ht(:,1)),mean(ht(:,responders),2),std(ht(:,responders)')/sqrt(length(responders)),'Linewidth',3)
    ylim([ys])
    ylabel('Normalized temporal weights')
    xlabel('Time')
    title('Hab')
    subplot(322)
    plot(at(:,responders))
    hold on
    errorbar(1:length(at(:,1)),mean(at(:,responders),2),std(at(:,responders)')/sqrt(length(responders)),'Linewidth',3)
    ylim([ys])
    ylabel('Normalized temporal weights')
    xlabel('Time')
    title('Acq')
    subplot(323)
    plot(e1t(:,responders))
    hold on
    errorbar(1:length(e1t(:,1)),mean(e1t(:,responders),2),std(e1t(:,responders)')/sqrt(length(responders)),'Linewidth',3)
    ylim([ys])
    ylabel('Normalized temporal weights')
    xlabel('Time')
    title('Ext1')
    subplot(324)
    plot(e2t(:,responders))
    hold on
    errorbar(1:length(e2t(:,1)),mean(e2t(:,responders),2),std(e2t(:,responders)')/sqrt(length(responders)),'Linewidth',3)
    ylim([ys])
    ylabel('Normalized temporal weights')
    xlabel('Time')
    title('Ext2')
    subplot(325)
    plot(e3t(:,responders))
    hold on
    errorbar(1:length(e3t(:,1)),mean(e3t(:,responders),2),std(e3t(:,responders)')/sqrt(length(responders)),'Linewidth',3)
    ylim([ys])
    ylabel('Normalized temporal weights')
    xlabel('Time')
    title('Ext3')

    figure
    subplot(321)
    plot(ht(:,nonresponders))
    hold on
    errorbar(1:length(ht(:,1)),mean(ht(:,nonresponders),2),std(ht(:,nonresponders)')/sqrt(length(responders)),'Linewidth',3)
    ylim([ys])
    ylabel('Normalized temporal weights')
    xlabel('Time')
    title('Hab')
    subplot(322)
    plot(at(:,nonresponders))
    hold on
    errorbar(1:length(at(:,1)),mean(at(:,nonresponders),2),std(at(:,nonresponders)')/sqrt(length(responders)),'Linewidth',3)
    ylim([ys])
    ylabel('Normalized temporal weights')
    xlabel('Time')
    title('Acq')
    subplot(323)
    plot(e1t(:,nonresponders))
    hold on
    errorbar(1:length(e1t(:,1)),mean(e1t(:,nonresponders),2),std(e1t(:,nonresponders)')/sqrt(length(responders)),'Linewidth',3)
    ylim([ys])
    ylabel('Normalized temporal weights')
    xlabel('Time')
    title('Ext1')
    subplot(324)
    plot(e2t(:,nonresponders))
    hold on
    errorbar(1:length(e2t(:,1)),mean(e2t(:,nonresponders),2),std(e2t(:,nonresponders)')/sqrt(length(responders)),'Linewidth',3)
    ylim([ys])
    ylabel('Normalized temporal weights')
    xlabel('Time')
    title('Ext2')
    subplot(325)
    plot(e3t(:,nonresponders))
    hold on
    errorbar(1:length(e3t(:,1)),mean(e3t(:,nonresponders),2),std(e3t(:,nonresponders)')/sqrt(length(responders)),'Linewidth',3)
    ylim([ys])
    ylabel('Normalized temporal weights')
    xlabel('Time')
    title('Ext3')

else
    figure
    subplot(321)
    plot(ht)
    hold on
    errorbar(1:length(ht(:,1)),mean(ht,2),std(ht')/sqrt(length(responders)),'Linewidth',3)
    ylim([0 1])
    ylabel('Normalized temporal weights')
    xlabel('Time')
    title('Hab')
    subplot(322)
    plot(at)
    hold on
    errorbar(1:length(at(:,1)),mean(at,2),std(at')/sqrt(length(responders)),'Linewidth',3)
    ylim([0 1])
    ylabel('Normalized temporal weights')
    xlabel('Time')
    title('Acq')
    subplot(323)
    plot(e1t)
    hold on
    errorbar(1:length(e1t(:,1)),mean(e1t,2),std(e1t')/sqrt(length(responders)),'Linewidth',3)
    ylim([0 1])
    ylabel('Normalized temporal weights')
    xlabel('Time')
    title('Ext1')
    subplot(324)
    plot(e2t)
    hold on
    errorbar(1:length(e2t(:,1)),mean(e2t,2),std(e2t')/sqrt(length(responders)),'Linewidth',3)
    ylim([0 1])
    ylabel('Normalized temporal weights')
    xlabel('Time')
    title('Ext2')
    subplot(325)
    plot(e3t)
    hold on
    errorbar(1:length(e3t(:,1)),mean(e3t,2),std(e3t')/sqrt(length(responders)),'Linewidth',3)
    ylim([0 1])
    ylabel('Normalized temporal weights')
    xlabel('Time')
    title('Ext3')

end
%% 11. Assess cumulative distribution of neuron weights. for Fig. 3G,H

% set threshold range
threshold = linspace(0,10,11);

%define storage vectors for fraction of cells included in Acq-, Ext1-, and
%Ext3-dominant ensembles
a = zeros(nAnimals, length(threshold));
e1 = zeros(nAnimals, length(threshold));
e3 = zeros(nAnimals, length(threshold));

for t = 1:length(threshold)
    aCells = [];
    e1Cells = [];
    e3Cells = [];
    
    for n = 1:nAnimals
        
        %call neuron factor weights from animal
        weights = neuFac{n};
        
        %calculate fractions of neurons with w > threshold and store
        a(n,t) = length(find(weights(:,domFacs(n,2))>threshold(t)))/nCells(n);
        e1(n,t) = length(find(weights(:,domFacs(n,3))>threshold(t)))/nCells(n);
        e3(n,t) = length(find(weights(:,domFacs(n,5))>threshold(t)))/nCells(n);
        
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
figure
errorbar(threshold,mean(a),std(a)/sqrt(nAnimals-1),'r')
hold on
errorbar(threshold,mean(e1),std(e1)/sqrt(nAnimals-1),'y')
hold on
errorbar(threshold,mean(e3),std(e3)/sqrt(nAnimals-1),'c')
hold on
xline(1)
xlabel('Threshold weight')
ylabel('Number of neurons')
title('Choosing dominant neurons')
%% 12. Identify Acq-Dominant, Ext1-Dominant, and Ext3-Dominant ensembles and their overlaps for Fig. 4B

%set threshTest to 0 if reproducing main figs. set threshTest to 1 if
%reproducing Supp Fig 6. note that null thresholds were not determined for
%non-shock mice

threshTest = 1;

if threshTest == 1
    threshes = table2array(readtable('fMeans.csv'));
    threshes = threshes(2:22);

    if dataset == 0 
        threshold = threshes(15:21);
    elseif dataset == 1
        threshold = threshes(1:14);
    end
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

cR = [];
CNR = [];
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
    
    if ismember(a,responders)
        cR = [cR [1+c:nCells(a)+c]];
    else
        cNR = [cNR [1+c:nCells(a)+c]];
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
    ae3 = overlap(~ismember(overlap,ae1e3));
    
    %cells dominating all ext1, ext3
    overlap = e1Cells(ismember(e1Cells,e3Cells));
    e1e3 = overlap(~ismember(overlap,ae1e3));
    
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
    acqDomOverlap(1,a) = length(aOnly)/length(aCells);
    acqDomOverlap(2,a) = length(ae1)/length(aCells);
    acqDomOverlap(3,a) = length(ae3)/length(aCells);
    acqDomOverlap(4,a) = length(ae1e3)/length(aCells);
    
    ext1DomOverlap(1,a) = length(e1Only)/length(e1Cells);
    ext1DomOverlap(2,a) = length(ae1)/length(e1Cells);
    ext1DomOverlap(3,a) = length(e1e3)/length(e1Cells);
    ext1DomOverlap(4,a) = length(ae1e3)/length(e1Cells);
    
    ext3DomOverlap(1,a) = length(e3Only)/length(e3Cells);
    ext3DomOverlap(2,a) = length(ae3)/length(e3Cells);
    ext3DomOverlap(3,a) = length(e1e3)/length(e3Cells);
    ext3DomOverlap(4,a) = length(ae1e3)/length(e3Cells);
    
    %assign neurons indices for activity pooled across animals and store in
    %vectors
    acqDomCells = [acqDomCells; aCells + c];
    ext1DomCells = [ext1DomCells; e1Cells + c];
    ext3DomCells = [ext3DomCells; e3Cells + c];
    
    acqOnly = [acqOnly; aOnly+c];
    ext1Only = [ext1Only; e1Only+c];
    ext3Only = [ext3Only; e3Only+c];
    acqExt1 = [acqExt1; ae1+c];
    acqExt3 = [acqExt3; ae3+c];
    ext1Ext3 = [ext1Ext3; ae3+c];
    acqExt1Ext3 = [acqExt1Ext3; ae1e3+c];
end

if dataset == 1
    aR = acqOnly(ismember(acqOnly,cR));
    e1R = ext1Only(ismember(ext1Only,cR));
    e3R = ext3Only(ismember(ext3Only,cR));
    ae1R = acqExt1(ismember(acqExt1,cR));
    ae3R = acqExt3(ismember(acqExt3,cR));
    e1e3R = ext1Ext3(ismember(ext1Ext3,cR));
    ae1e3R = ext1Ext3(ismember(ext1Ext3,cR));
    
    acqDomR = acqDomCells(ismember(acqDomCells,cR));
    ext1DomR = acqDomCells(ismember(ext1DomCells,cR));
    ext3DomR = acqDomCells(ismember(ext3DomCells,cR));

    aNR = acqOnly(ismember(acqOnly,cNR));
    e1NR = ext1Only(ismember(ext1Only,cNR));
    e3NR = ext3Only(ismember(ext3Only,cNR));
    ae1NR = acqExt1(ismember(acqExt1,cNR));
    ae3NR = acqExt3(ismember(acqExt3,cNR));
    e1e3NR = ext1Ext3(ismember(ext1Ext3,cNR));
    ae1e3NR = ext1Ext3(ismember(ext1Ext3,cNR));
    
    acqDomNR = acqDomCells(ismember(acqDomCells,cNR));
    ext1DomNR = ext1DomCells(ismember(ext1DomCells,cNR));
    ext3DomNR = ext3DomCells(ismember(ext3DomCells,cNR));
end

% plot overlap fractions

if dataset == 1
    %responders
    figure
    subplot(131)
    barWithError(acqDomOverlap(:,responders)')
    ylabel('Fraction of Cells')
    xticklabels({'Acq Only', 'Acq/Ext1', 'Acq/Ext3', 'Acq/Ext1/Ext3'})
    title('Overlap with Acq-Dominant Neurons')
    
    subplot(132)
    barWithError(ext1DomOverlap(:,responders)')
    ylabel('Fraction of Cells')
    xticklabels({'Ext1 Only', 'Acq/Ext1', 'Ext1/Ext3', 'Acq/Ext1/Ext3'})
    title('Overlap with Ext1-Dominant Neurons')
    
    subplot(133)
    barWithError(ext3DomOverlap(:,responders)')
    ylabel('Fraction of Cells')
    xticklabels({'Ext3 Only', 'Acq/Ext3', 'Ext1/Ext3', 'Acq/Ext1/Ext3'})
    title('Overlap with Ext3-Dominant Neurons')
    
    
    %nonresponders
    figure
    subplot(131)
    barWithError(acqDomOverlap(:,nonresponders)')
    ylabel('Fraction of Cells')
    xticklabels({'Acq Only', 'Acq/Ext1', 'Acq/Ext3', 'Acq/Ext1/Ext3'})
    title('Overlap with Acq-Dominant Neurons')
    
    subplot(132)
    barWithError(ext1DomOverlap(:,nonresponders)')
    ylabel('Fraction of Cells')
    xticklabels({'Ext1 Only', 'Acq/Ext1', 'Ext1/Ext3', 'Acq/Ext1/Ext3'})
    title('Overlap with Ext1-Dominant Neurons')
    
    subplot(133)
    barWithError(ext3DomOverlap(:,nonresponders)')
    ylabel('Fraction of Cells')
    xticklabels({'Ext3 Only', 'Acq/Ext3', 'Ext1/Ext3', 'Acq/Ext1/Ext3'})
    title('Overlap with Ext3-Dominant Neurons')
else
    figure
    subplot(131)
    barWithError(acqDomOverlap')
    ylabel('Fraction of Cells')
    xticklabels({'Acq Only', 'Acq/Ext1', 'Acq/Ext3', 'Acq/Ext1/Ext3'})
    title('Overlap with Acq-Dominant Neurons')

    subplot(132)
    barWithError(ext1DomOverlap')
    ylabel('Fraction of Cells')
    xticklabels({'Ext1 Only', 'Acq/Ext1', 'Ext1/Ext3', 'Acq/Ext1/Ext3'})
    title('Overlap with Ext1-Dominant Neurons')

    subplot(133)
    barWithError(ext3DomOverlap')
    ylabel('Fraction of Cells')
    xticklabels({'Ext3 Only', 'Acq/Ext3', 'Ext1/Ext3', 'Acq/Ext1/Ext3'})
    title('Overlap with Ext3-Dominant Neurons')
end


%% 13. Use a Fisher linear decoder to predict group based on average activity (z-score relative to Acquisition) of each ensemble over time for fig. 4C

%TO AVOID HAVING TO RE-RUN PREVIOUS SCRIPTS, LOAD THESE FILES WITH
%ENSEMBLES' CELL INDICES
%
% IF RUNNING THRESHTEST SKIP THIS SECTION

trialLength = 60;
hab = [1:8];
acq = [9:16];
ext1 = [17:22];
ext2 = [23:28];
ext3 = [29:34];

%saline ensembles
aS = table2array(readtable('acqOnlyS.txt'));
e1S = table2array(readtable('ext1OnlyS.txt'));
e3S = table2array(readtable('ext3OnlyS.txt'));
ae1S = table2array(readtable('acqExt1S.txt'));
ae3S = table2array(readtable('acqExt3S.txt'));
e1e3S = table2array(readtable('ext1Ext3S.txt'));
ae1e3S = table2array(readtable('acqExt1Ext3S.txt'));

%responders ensembles
aR = table2array(readtable('acqOnlyR.txt'));
e1R = table2array(readtable('ext1OnlyR.txt'));
e3R = table2array(readtable('ext3OnlyR.txt'));
ae1R = table2array(readtable('acqExt1R.txt'));
ae3R = table2array(readtable('acqExt3R.txt'));
e1e3R = table2array(readtable('ext1Ext3R.txt'));
ae1e3R = table2array(readtable('acqExt1Ext3R.txt'));

%nonresponders ensembles
aNR =table2array(readtable('acqOnlyNR.txt'));
e1NR = table2array(readtable('ext1OnlyNR.txt'));
e3NR = table2array(readtable('ext3OnlyNR.txt'));
ae1NR = table2array(readtable('acqExt1NR.txt'));
ae3NR = table2array(readtable('acqExt3NR.txt'));
e1e3NR = table2array(readtable('ext1Ext3NR.txt'));
ae1e3NR = table2array(readtable('acqExt1Ext3NR.txt'));

%nonshock ensembles
aN =table2array(readtable('acqOnlyN.txt'));
e1N = table2array(readtable('ext1OnlyN.txt'));
e3N = table2array(readtable('ext3OnlyN.txt'));
ae1N = table2array(readtable('acqExt1N.txt'));
ae3N = table2array(readtable('acqExt3N.txt'));
e1e3N = table2array(readtable('ext1Ext3N.txt'));
ae1e3N = table2array(readtable('acqExt1Ext3N.txt'));

%load pooled data for both datasets
poolMatPsil = table2array(readtable('totApoolPsil.csv'));
poolMatSal = table2array(readtable('totApoolSal.csv'));
poolMatNonshock = table2array(readtable('totApoolNonshock.csv'));

%define time windows
timeAcq = [trialLength*max(hab)+1:trialLength*max(acq)];
timeExt1 = [trialLength*max(acq)+1:trialLength*max(ext1)];
timeExt3 = [trialLength*max(ext2)+1:trialLength*max(ext3)];


%zscore with respect to acq
poolMatPsil = (poolMatPsil-mean(poolMatPsil(timeAcq,:)))./std(poolMatPsil(timeAcq,:));
poolMatSal = (poolMatSal-mean(poolMatSal(timeAcq,:)))./std(poolMatSal(timeAcq,:));
poolMatNonshock = (poolMatNonshock-mean(poolMatNonshock(timeAcq,:)))./std(poolMatNonshock(timeAcq,:));

%ensemble sets for responders, nonresponders, and saline mice. Acq-Only,
%Ext1-Only, Ext3-Only, Acq/Ext1, Acq/Ext3, Ext1/Ext3, Acq/Ext1/Ext3,
%Acq-Dom, Ext1-Dom, Ext3-Dom
resEns = {aR, e1R, e3R, ae1R, ae3R, e1e3R, ae1e3R, [aR; ae1R; ae3R; ae1e3R], [e1R; ae1R; e1e3R; ae1e3R], [e3R;  ae3R; e1e3R; ae1e3R]};
nonresEns = {aNR, e1NR, e3NR, ae1NR, ae3NR, e1e3NR, ae1e3NR, [aNR; ae1NR; ae3NR; ae1e3NR], [e1NR; ae1NR; e1e3NR; ae1e3NR], [e3NR;  ae3NR; e1e3NR; ae1e3NR]};
salEns = {aS, e1S, e3S, ae1S, ae3S, e1e3S, ae1e3S,[aS; ae1S; ae3S; ae1e3S], [e1S; ae1S; e1e3S; ae1e3S], [e3S;  ae3S; e1e3S; ae1e3S]};
nsEns = {aN, e1N, e3N, ae1N, ae3N, e1e3N, ae1e3N,[aN; ae1N; ae3N; ae1e3N], [e1N; ae1N; e1e3N; ae1e3N], [e3N;  ae3N; e1e3N; ae1e3N]};

numEnsembles = 10;


for g = 1:4
    if g==1
        data = resEns;
    elseif g==2
        data = nonresEns;
    elseif g == 3
        data = salEns;
    else
        data = nsEns;
    end
    for n=1:7
        l(n,g) = length(data{n});
    end
end


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

%create matrices of average ensembles acts over time in the session of
%interest
ensActsR =[];
ensActsNR =[];
ensActsS =[];
for l = 1:numEnsembles
    ensActsR = [ensActsR mean(poolMatPsil(times,resEns{l}),2)];
    ensActsNR = [ensActsNR mean(poolMatPsil(times,nonresEns{l}),2)];
    ensActsS = [ensActsS mean(poolMatPsil(times,salEns{l}),2)];
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
for l = 1:numEnsembles+1
    if l < 10 %for each individual ensemble
        if comparisons(c) == 0
            data = [ensActsR(:,l); ensActsNR(:,l)]; %compare responders to nonresponders
        elseif comparisons(c) == 1
            data = [ensActsR(:,l); ensActsS(:,l)];  %compare responders to saline
        else
            data = [ensActsNR(:,l); ensActsS(:,l)]; %compare nonresponders to saline
        end
    else %with every ensemble as a predictor
        if comparisons(c) == 0
            data = [ensActsR(:,1:7); ensActsNR(:,1:7)]; %compare responders to nonresponders
        elseif comparisons(c) == 1
            data = [ensActsR(:,1:7); ensActsS(:,1:7)]; %compare responders to saline
        else
            data = [ensActsNR(:,1:7); ensActsS(:,1:7)]; %compare nonresponders to saline
        end
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
titles = {'Acq Only', 'Ext1 Only', 'Ext3 Only', 'Acq/Ext1', 'Acq/Ext3', 'Ext1/Ext3', 'Acq/Ext1/Ext3','Acq-Dom','Ext1-Dom','Ext3-Dom','All ensembles'};

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

%% 14. Extract average activity (z-score relative to Acquisition) of each cell in the ensemble during Extinction 1, Extinction 3. for fig 4D-J, fig. 5B,C,D

%returns activity during Ext1 and Ext3, z-scored with respect to Acq

if threshTest == 1 %(Supp Fig 6 D-H)
    resEns = {acqDomR, ext1DomR, ext3DomR};
    nonresEns = {acqDomNR, ext1DomNR, ext3DomNR};
    salEns = {acqDomCells, ext1DomCells, ext3DomCells};

    nGroups = 3;
    
    %load pooled data for both datasets
    poolMatPsil = table2array(readtable('totApoolPsil.csv'));
    poolMatSal = table2array(readtable('totApoolSal.csv'));
    poolMatNonshock = table2array(readtable('totApoolNonshock.csv'));

    %define time windows
    timeAcq = [trialLength*max(hab)+1:trialLength*max(acq)];
    timeExt1 = [trialLength*max(acq)+1:trialLength*max(ext1)];
    timeExt3 = [trialLength*max(ext2)+1:trialLength*max(ext3)];

    %zscore with respect to acq
    poolMatPsil = (poolMatPsil-mean(poolMatPsil(timeAcq,:)))./std(poolMatPsil(timeAcq,:));
    poolMatSal = (poolMatSal-mean(poolMatSal(timeAcq,:)))./std(poolMatSal(timeAcq,:));
    
else
    nGroups = 4;
end
%initialize storage cell for average changed activity matrices of each ensemble, in
%each group
changedActivity = cell(length(resEns),nGroups);

%loop through animals
for n = 1:nGroups
    
    %loop through all ensembles
    for l=1:length(resEns)
        
        %determine which group
        if n == 1
        
            cellPop = resEns{l};
        
            %calculate average activity of cells during Extinction 1 and
            %Extinction 3
            x =  [mean(poolMatPsil(trialLength*max(acq)+1:trialLength*max(ext1),cellPop))'  mean(poolMatPsil(trialLength*max(ext2)+1:trialLength*max(ext3),cellPop))'];
            
        elseif n == 2
            
            cellPop = nonresEns{l};
        
            %calculate average activity of cells during Extinction 1 and
            %Extinction 3
            x =  [mean(poolMatPsil(trialLength*max(acq)+1:trialLength*max(ext1),cellPop))'  mean(poolMatPsil(trialLength*max(ext2)+1:trialLength*max(ext3),cellPop))'];
        elseif n==3
            cellPop = salEns{l};
        
            %calculate average activity of cells during Extinction 1 and
            %Extinction 3
            x =  [mean(poolMatSal(trialLength*max(acq)+1:trialLength*max(ext1),cellPop))'  mean(poolMatSal(trialLength*max(ext2)+1:trialLength*max(ext3),cellPop))'];
        else
            cellPop = nsEns{l};
            x=[mean(poolMatNonshock(trialLength*max(acq)+1:trialLength*max(ext1),cellPop))'  mean(poolMatNonshock(trialLength*max(ext2)+1:trialLength*max(ext3),cellPop))'];
        end
        
        %store activity matrices
        changedActivity{l,n} = x;
        
    end
    
end

%to calculate results with the fixed threshold variation (Supp Fig 6 C)

% nEns = 3;
% nThr = length(threshes);
% 
% for n=1:nEns
%     for t = 1:nThr
%        
%        cellPop = collectionR{t,n};
%        x =  [mean(poolMatPsil(trialLength*max(acq)+1:trialLength*max(ext1),cellPop))'  mean(poolMatPsil(trialLength*max(ext2)+1:trialLength*max(ext3),cellPop))'];
%        domMat1(t,1:3,n) = [mean(x(:,1)) std(x(:,1)) length(x(:,1))];
%        domMat3(t,1:3,n) = [mean(x(:,2)) std(x(:,2)) length(x(:,2))];
%        
%        cellPop = collectionNR{t,n};
%        x =  [mean(poolMatPsil(trialLength*max(acq)+1:trialLength*max(ext1),cellPop))'  mean(poolMatPsil(trialLength*max(ext2)+1:trialLength*max(ext3),cellPop))'];
%        domMat1(t,4:6,n) = [mean(x(:,1)) std(x(:,1)) length(x(:,1))];
%        domMat3(t,4:6,n) = [mean(x(:,2)) std(x(:,2)) length(x(:,2))];
%        
%        cellPop = collection{t,n};
%        x =  [mean(poolMatSal(trialLength*max(acq)+1:trialLength*max(ext1),cellPop))'  mean(poolMatSal(trialLength*max(ext2)+1:trialLength*max(ext3),cellPop))'];
%        domMat1(t,7:9,n) = [mean(x(:,1)) std(x(:,1)) length(x(:,1))];
%        domMat3(t,7:9,n) = [mean(x(:,2)) std(x(:,2)) length(x(:,2))];
%        
%    end
% end

%% 15. Plot example traces from specific ensembles. for fig. 5A
cellPop = aS; %for ext3-dom example, plug in "e3S", "e3R","e3NR" for the cell pop variable. for acq-dom example, plug in "ae1e3S", "ae1e3R","ae1e3NR" for the cell pop variable. 

cellPopR = aR;
cellPopNR = aNR;


figure
subplot(321)
for n=1:15
    plot(poolMatSal(times+1:(trialLength*(2*length(hab)+length(ext1))),cellPop(n))+10*n)
    hold on
end
ylabel('Saline')
title('Ext1')
subplot(322)
for n=1:15
    plot(poolMatSal(trialLength*(2*length(hab)+2*length(ext1))+1:(trialLength*(2*length(hab)+3*length(ext1))),cellPop(n))+10*n)
    hold on
end
title('Ext3')
subplot(323)
for n=1:15
    plot(poolMatPsil(trialLength*(2*length(hab))+1:(trialLength*(2*length(hab)+length(ext1))),cellPopNR(n))+10*n)
    hold on
end
ylabel('Non-Responders')
title('Ext1')
subplot(324)
for n=1:15
    plot(poolMatPsil(trialLength*(2*length(hab)+2*length(ext1))+1:(trialLength*(2*length(hab)+3*length(ext1))),cellPopNR(n))+10*n)
    hold on
end
title('Ext3')
subplot(325)
for n=1:15
    plot(poolMatPsil(trialLength*(2*length(hab))+1:(trialLength*(2*length(hab)+length(ext1))),cellPopR(n))+10*n)
    hold on
end
ylabel('Responders')
title('Ext1')
subplot(326)
for n=1:15
    plot(poolMatPsil(trialLength*(2*length(hab)+2*length(ext1))+1:(trialLength*(2*length(hab)+3*length(ext1))),cellPopR(n))+10*n)
    hold on
end
title('Ext3')
%export data into visualizer


%% 16. Plot average dFoF of non-intersecting Acq-dominant and Ext1/3-dominant ensembles and their ratios over seconds, trials 

poolMatPsil = table2array(readtable('totApoolPsil.csv'));
poolMatSal = table2array(readtable('totApoolSal.csv'));

time = timeExt3; %set time to desired session (timeAcq, timeExt1, timeExt3)

x=mean(poolMatPsil(time,[aR; ae1R; ae1R; ae1e3R]),2)./mean(poolMatPsil(time,[e1R; e3R; e1e3R]),2);

a=mean(poolMatPsil(time,[aNR; ae1NR; ae1NR; ae1e3NR]),2)./mean(poolMatPsil(time,[e1NR; e3NR; e1e3NR]),2);

d = mean(poolMatSal(time,[aS; ae1S; ae1S; ae1e3S]),2)./mean(poolMatSal(time,[e1S; e3S; e1e3S]),2);

if length(time) == length(timeAcq)
    nTrials = 8;
    y=[mean(x(1:trialLength)) mean(x(61:trialLength*2)) mean(x(trialLength*2+1:trialLength*3)) mean(x(trialLength*3+1:trialLength*4)) mean(x(trialLength*4+1:trialLength*5)) mean(x(trialLength*5+1:trialLength*6)) mean(x(trialLength*6+1:trialLength*7)) mean(x(trialLength*7+1:trialLength*8))];
    z=[std(x(1:trialLength)) std(x(61:trialLength*2)) std(x(trialLength*2+1:trialLength*3)) std(x(trialLength*3+1:trialLength*4)) std(x(trialLength*4+1:trialLength*5)) std(x(trialLength*5+1:trialLength*6)) std(x(trialLength*6+1:trialLength*7)) std(x(trialLength*7+1:trialLength*8))];

    b=[mean(a(1:trialLength)) mean(a(61:trialLength*2)) mean(a(trialLength*2+1:trialLength*3)) mean(a(trialLength*3+1:trialLength*4)) mean(a(trialLength*4+1:trialLength*5)) mean(a(trialLength*5+1:trialLength*6)) mean(a(trialLength*6+1:trialLength*7)) mean(a(trialLength*7+1:trialLength*8))];
    c=[std(a(1:trialLength)) std(a(61:trialLength*2)) std(a(trialLength*2+1:trialLength*3)) std(a(trialLength*3+1:trialLength*4)) std(a(trialLength*4+1:trialLength*5)) std(a(trialLength*5+1:trialLength*6)) std(a(trialLength*6+1:trialLength*7)) std(a(trialLength*7+1:trialLength*8))];
    
    e = [mean(d(1:trialLength)) mean(d(61:trialLength*2)) mean(d(trialLength*2+1:trialLength*3)) mean(d(trialLength*3+1:trialLength*4)) mean(d(trialLength*4+1:trialLength*5)) mean(d(trialLength*5+1:trialLength*6)) mean(d(trialLength*6+1:trialLength*7)) mean(d(trialLength*7+1:trialLength*8))];
    f = [std(d(1:trialLength)) std(d(61:trialLength*2)) std(d(trialLength*2+1:trialLength*3)) std(d(trialLength*3+1:trialLength*4)) std(d(trialLength*4+1:trialLength*5)) std(d(trialLength*5+1:trialLength*6)) std(d(trialLength*6+1:trialLength*7)) std(d(trialLength*7+1:trialLength*8))];
else
    nTrials = 6;
    y=[mean(x(1:trialLength)) mean(x(61:trialLength*2)) mean(x(trialLength*2+1:trialLength*3)) mean(x(trialLength*3+1:trialLength*4)) mean(x(trialLength*4+1:trialLength*5)) mean(x(trialLength*5+1:trialLength*6)) ];
    z=[std(x(1:trialLength)) std(x(61:trialLength*2)) std(x(trialLength*2+1:trialLength*3)) std(x(trialLength*3+1:trialLength*4)) std(x(trialLength*4+1:trialLength*5)) std(x(trialLength*5+1:trialLength*6)) ];

    b=[mean(a(1:trialLength)) mean(a(61:trialLength*2)) mean(a(trialLength*2+1:trialLength*3)) mean(a(trialLength*3+1:trialLength*4)) mean(a(trialLength*4+1:trialLength*5)) mean(a(trialLength*5+1:trialLength*6)) ];
    c=[std(a(1:trialLength)) std(a(61:trialLength*2)) std(a(trialLength*2+1:trialLength*3)) std(a(trialLength*3+1:trialLength*4)) std(a(trialLength*4+1:trialLength*5)) std(a(trialLength*5+1:trialLength*6)) ];
    
    e = [mean(d(1:trialLength)) mean(d(61:trialLength*2)) mean(d(trialLength*2+1:trialLength*3)) mean(d(trialLength*3+1:trialLength*4)) mean(d(trialLength*4+1:trialLength*5)) mean(d(trialLength*5+1:trialLength*6)) ];
    f = [std(d(1:trialLength)) std(d(61:trialLength*2)) std(d(trialLength*2+1:trialLength*3)) std(d(trialLength*3+1:trialLength*4)) std(d(trialLength*4+1:trialLength*5)) std(d(trialLength*5+1:trialLength*6)) ];
end

figure
subplot(131)
errorbar(1:length(time),mean(poolMatSal(time,[aS; ae1S; ae1S; ae1e3S]),2),std(poolMatSal(time,[aS; ae1S; ae3S; ae1e3S])')./sqrt(length([aR; ae1R; ae3R; ae1e3R])-1),'r')
hold on
errorbar(1:length(time),mean(poolMatSal(time,[e1S; e3S; e1e3S]),2),std(poolMatSal(time,[e1S; e3S; e1e3S])')./sqrt(length([e1S; e3S; e1e3S])-1),'b')
hold on
if length(time) == length(timeAcq)
for n=1:8
    xline((n-1)*trialLength+56)
end
end
ylabel('Responders, avg dFoF')
legend({'Acq-Dom', 'Ext Only-Dom'})
subplot(132)
plot(d)
ylabel('Activity ratio (Acq/Ext)')
subplot(133)
bar(e)
hold on
errorbar(1:nTrials,e,f./sqrt(59))
ylabel('Activity ratio within trials')

% figure
% subplot(231)
% errorbar(1:length(time),mean(poolMatPsil(time,[aR; ae1R; ae1R; ae1e3R]),2),std(poolMatPsil(time,[aR; ae1R; ae3R; ae1e3R])')./sqrt(length([aR; ae1R; ae3R; ae1e3R])-1),'r')
% hold on
% errorbar(1:length(time),mean(poolMatPsil(time,[e1R; e3R; e1e3R]),2),std(poolMatPsil(time,[e1R; e3R; e1e3R])')./sqrt(length([e1R; e3R; e1e3R])-1),'b')
% hold on
% if length(time) == length(timeAcq)
% for n=1:8
%     xline((n-1)*trialLength+56)
% end
% end
% ylabel('Responders, avg dFoF')
% legend({'Acq-Dom', 'Ext Only-Dom'})
% subplot(232)
% plot(x)
% ylabel('Activity ratio (Acq/Ext)')
% subplot(233)
% bar(y)
% hold on
% errorbar(1:nTrials,y,z./sqrt(59))
% ylabel('Activity ratio within trials')
% 
% subplot(234)
% errorbar(1:length(time),mean(poolMatPsil(time,[aNR; ae1NR; ae3NR; ae1e3NR]),2),std(poolMatPsil(time,[aNR; ae1NR; ae3NR; ae1e3NR])')./sqrt(length([aNR; ae1NR; ae3NR; ae1e3NR])-1),'r')
% hold on
% errorbar(1:length(time),mean(poolMatPsil(time,[e1NR; e3NR; e1e3NR]),2),std(poolMatPsil(time,[e1NR; e3NR; e1e3NR])')./sqrt(length([e1NR; e3NR; e1e3NR])-1),'b')
% hold on
% if length(time) == length(timeAcq)
% for n=1:8
%     xline((n-1)*trialLength+56)
% end
% end
% ylabel('Non-Responders, avg dFoF')
% legend({'Acq-Dom', 'Ext Only-Dom'})
% subplot(235)
% plot(a)
% ylabel('Activity ratio (Acq/Ext)')
% subplot(236)
% bar(b)
% hold on
% errorbar(1:nTrials,b,c./sqrt(59))
% ylabel('Activity ratio within trials')


