%%%%%%%%%%%%%%%%%%%%%%%%%%% overnight_CE_v2.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Runs "CONSTRUCT_CE_v4.m" on a loop the numbers of times specified. 
%%%% Output from "CONSTRUCT_CE_v4.m" is placed in a folder called "LG_EICs",
%%%% which must exist for this script to work.

%%%%%%% USER INPUT %%%%%%%%%
numRuns = 500;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

backhome = cd('./LG_ECIs');
for i = 1:numRuns    
    cd(backhome)
    run('CONSTRUCT_CE_v4_4.m')
    thisCV = sprintf('%6.6f',CV_score);
    backhome = cd('./LG_ECIs');
    movefile('../LG_ECIs.txt',['LG_ECIs.' thisCV '.txt']);
end
cd(backhome)
toc