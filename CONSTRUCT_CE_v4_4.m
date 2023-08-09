%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTRUCT_CE.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A script to hopefully find the optimal set of clusters that best
% describes and predicts the training data provided in matrix 
% "surf_normalized_counts" and either "surf_energy" or "mf_resids". 

% Must have CV_calc_vX.m in the active directory.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format short g

% If you'd like to use the mean field residuals as your fitting energies,
% please set "mf_flag" to "1". Else, set to "0". Default is "0"
                
mf_flag = 0;

% If you'd like to use the adsorption energy to fit the energies, set
% "ads_flag" to "1". Else, set to "0". Default is "0"

ads_flag = 1;

% Would you like to add clusters 1 at a time of both 1 and 2 at a time?
% Adding 2 at a time increases the total number of gradent calculation 
% loops by M choose 2...so it scales as M^2 ('M squared') where M is the 
% total number of clusters available in "INTERACTIONS.txt". 
% Select 1 at a time (1) or 1 AND 2 at a time (2)

at_a_time = 1;

% User defined percentage (fraction) of total number of known energies that become
% validation set during CV calculation (for a constant leave-x-out, set "a"
% to x/N)

a = 0.60;

% If a cluster is underrepresented in the known structures, it can be
% removed from consideration to speed up the algorithm. Set "rep" as the
% minimum number of structures that have to have a cluster for that cluster 
% to be considered 

rep = 5;

% If you'd like to specify that certain clusters be added at the start (no
% guarantee the algorithm won't remove them, mind you) add them here. Note:
% if you don't want ANY starting clusters, just delete the numbers and
% leave an empty set...DO NOT comment this out.

start_clusters = [1	2];%	11	16	28	29	34	51	60	74	78	83	107	114	161	186	207	372	404	416	533	608	614	620	622	623	656	718	728	770];
len_start_cl = round(length(start_clusters)+5,1);
start_protect_clusters = [1 2];

% v4.4 added feature: mandatory clusters. These clusters will never be
% removed by the algorithm.

mandatory_clusters = [1 2];

% User defined fraction of clusters to delete from start_clusters (must be
% less than 1)

start_fraction_to_remove = 0.8*rand; 

% User defined initial size of *randomly* generated clusters

max_stsz = len_start_cl - round(start_fraction_to_remove*size(start_clusters,2));
min_stsz = 3;

start_sz = round(min_stsz +(max_stsz-min_stsz)*rand);

if start_sz < min_stsz
    start_sz = min_stsz;
end
% Specify a CV lowering tolerance, making it so additions/removals only 
% occur if the CV lowers by AT LEAST this amount (0.005 eV is a good start)

tol = 0.00005;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
unrep = find(~any(surf_normalized_counts,1)); 
zero_clust=find(sum(surf_normalized_counts,1)==0);

%%%% Find all unique structures and amongst duplicates, select the
%%%% structure corresponding to the lowest energy inputted.
unrep = find(~any(surf_normalized_counts,1)); 
zero_clust=find(sum(surf_normalized_counts,1)==0);
%%%% Find all unique structures and amongst duplicates, select the
%%%% structure corresponding to the lowest energy inputted.
uCovs = unique(coverages);
sMat = surf_normalized_counts;
indMat = 1:size(surf_normalized_counts,1);
sEn = surf_energy;
totoss = [];
for ii = uCovs'
    sub_mat= sMat(abs(coverages-ii)<0.001,:);
    sub_en = sEn(abs(coverages-ii)<0.001);
    sub_ind = indMat(abs(coverages-ii)<0.001);
    [uSubMat,jja,jjb] = unique(sub_mat,'rows','stable');
    sz_dup = numel(jjb)-numel(jja);
    if sz_dup == 0
        continue
    end
    [count,~,idxcount] = histcounts(jjb,numel(jja));
    count(count>1) = 2;
    kk = 0;
    for jj = 1:numel(count)
        if count(jj)>1
            count(jj) = count(jj) +kk;
            kk = kk + 1;
        end
    end
    sub_thing = min(count(count>1));
    max_thing = max(count);
%     count = count - sub_thing + 2;
    what = count(idxcount);
    for jj = sub_thing:max_thing
        test = find(what==jj);
        [~, iik] = min(sub_en(test));
        totoss = [totoss sub_ind(test(test~=test(iik)))]; %#ok<AGROW>  
    end
 
    
end
tokeep = setdiff(indMat,totoss);    
    
MAT = surf_normalized_counts(tokeep,:);
N=size(MAT,1);                                                             % total number of known structures
M=size(MAT,2);                                                             % total number of clusters
mod_coverages = coverages(tokeep);                                         % grab the associated coverages


deficient = (all(MAT(:,start_clusters)==0,1));
deficient_clusters = start_clusters(deficient);

if any(deficient)
    start_clusters(deficient) = [];
    fprintf('There are no non-zero terms in the counts matrix for the following cluster(s)\n ');
    disp(deficient_clusters)
    fprintf('This(these) clusters will be removed.\n')
    
end
N=size(MAT,1);                                                              % total number of known (unique) structures
candidate_clusters = find(sum(counts > 0,1)>=rep);                                                            
M=size(MAT,2);                                                             % total number of clusters
if isempty(start_clusters)
    nonzero_M = union(candidate_clusters, start_clusters)';
else
    nonzero_M = union(candidate_clusters, start_clusters);
end
% Check to ensure each cluster included actually has representation in the
% MAT matrix, remove and alert user

% Grab appropriate energies
if mf_flag == 1
    energy = mf_resids;
else
    energy = surf_energy(tokeep);
end                       

if ads_flag == 1 && mf_flag == 0
    energy = ads_energy(tokeep);
    
    MAT = MAT./mod_coverages;
    MAT(~isfinite(MAT)) = 0;
else
    energy = surf_energy(tokeep);
end
               
% create random vector of interactions that will be used as a starting CE
% in addition to specified starting CE
starting_flag = 0;
too_big = 0;
too_small = 0;
while starting_flag ==0
    
    out_CE = setdiff(nonzero_M,start_clusters);                 % which clusters aren't part of the starting clusters?
    rsCE = randperm(size(out_CE,2),start_sz);                   % A randon permutation to use as indices
    test_CE = out_CE(rsCE);
    sz_CE = size(test_CE,2);
    
    sz_start_clusters = size(start_clusters,2);                             % How many clusters in the start_clusters vector?
    start_remove = round(start_fraction_to_remove*sz_start_clusters);       % How many to remove?
    rrCE = randperm(sz_start_clusters,start_remove);                        % a random permutation to use as indices
    start_clusters(rrCE) = [];                                              % Remove those indices from the start clusters
    start_clusters = union(start_clusters,start_protect_clusters);                % add back the protected clusters in any have been removed.
    CE = sort([test_CE start_clusters]);
    CE = union(CE,mandatory_clusters);                                      % Mandatory clusters added back in case they were removed.
    
    
    test_mat = MAT(:,CE);
    if size(CE,2) >= rank(MAT)
        error('There are not enough known structures for this many starting clusters. Reduce the amount or add more structures and try again.\n')
    end    
    [CV_score,std_dev] = CV_calc_v3(test_mat,energy,a);        %calulate the intial CV score
    if CV_score < 100
        starting_flag = 1;
    else
        too_big = too_big + 1;
    end
    if too_big > 5 && too_small == 0
        fprintf('Starting size potentially too big. Reducing starting size by 1\n')
        start_sz = start_sz - 1;
        too_big = 0;
        if start_sz == 1
            too_small = 1;
        end
    elseif too_big > 5 && too_small == 1
        fprintf('Oops now we are too small. Lets try adding instead...\n');
        start_sz = start_sz + 1;
        too_big = 0;
    end
        
end
fprintf('----------------------------------------------------------------\n')
fprintf('----------------------------------------------------------------\n')
fprintf('\nThe starting Cluster Expansion (CE) is:\n')
disp(CE)
fprintf('Its CV score is: %8.6g eV/site\n',CV_score)
fprintf('Its standard deviation is: %8.6g eV/site\n',std_dev)
fprintf('----------------------------------------------------------------\n')
fprintf('----------------------------------------------------------------\n')
% Run through every cluster and either adding it or subtracting it from the
% current CE. Keep a record of 

flag = 1;
AA = zeros(M.^2,4);        % Column 1 and 2 are clusters ii and jj. Col 3 is used to track what kind of addition/removal was performed. Col 4 is the CV score
newl = 0;
nonzero_M = setdiff(nonzero_M,mandatory_clusters);      % don't loop over the mandatory clusters
if at_a_time == 2
    totloops = 1/2*size(nonzero_M,2)*(size(nonzero_M,2)-1);
else
    totloops = size(nonzero_M,2);
end

while flag < 2
    fprintf('Calculating the gradient of the current CE...\n(This can take a while)\n')
    loop = 0;
    for kk = nonzero_M
        loop = loop + 1;
        AA(loop,1) = kk;
        if any(zero_clust==kk)
            AA(loop,4) = 1000;
        elseif any(CE == kk)
            AA(loop,3) = -1;                       % This is used to track which clusters were added or removed along the 1 thru M cluster additions/removals
            test_CE = setdiff(CE,kk);
            out_CE = setdiff(nonzero_M,test_CE);
            test_mat = MAT(:,test_CE);
            [AA(loop,4),~] = CV_calc_v3(test_mat,energy,a);
                       
        else
            AA(loop,3) = 1;
            test_CE = sort([CE kk]);
            out_CE = setdiff(nonzero_M,test_CE);
            test_mat = MAT(:,test_CE);
            [AA(loop,4),~] = CV_calc_v3(test_mat,energy,a);
            
        end     
        if mod(loop,round(totloops/20)) == 0                    % Display a percent complete
            pcent = round(round(loop/totloops/.05)*0.05*100);
            if pcent >= 100
                pcent = 99.99;
            end
            fprintf('%2g%%...',pcent);
            newl = newl + 1;
        end
        if mod(newl,10) < 0.001 && mod(loop,round(totloops/20)) == 0
            fprintf('\n');
        end
    end
    if at_a_time == 2
        tt = 0;
        for kk = nonzero_M
            tt = tt + 1;
            for jj = nonzero_M(tt+1:end)
                loop = loop + 1;
                AA(loop,[1 2]) = [kk jj];    
                if any(zero_clust==kk) || any(zero_clust==jj)
                    AA(loop,4) = 1000;        
                elseif any(CE == kk) && any(CE == jj)
                    AA(loop,3) = -3;                                               % This is used to track which clusters were added or removed along the 1 thru Mchoose2 cluster additions/removals
                    test_CE = setdiff(CE,[kk jj]);
                    out_CE = setdiff(nonzero_M,test_CE);
                    test_mat = MAT(:,test_CE);
                    [AA(loop,4),~] = CV_calc_v3(test_mat,energy,a);                  
                elseif any(CE == kk) && ~any(CE == jj)
                    AA(loop,3) = -2;
                    test_CE = setdiff(CE,kk);
                    test_CE = [test_CE jj];                                         %#ok<AGROW>
                    out_CE = setdiff(nonzero_M,test_CE);
                    test_mat = MAT(:,test_CE);
                    [AA(loop,4),~] = CV_calc_v3(test_mat,energy,a);
                elseif ~any(CE == kk) && any(CE == jj)
                    AA(loop,3) = 2;
                    test_CE = setdiff(CE,jj);
                    test_CE = [test_CE kk];                                         %#ok<AGROW>
                    out_CE = setdiff(nonzero_M,test_CE);
                    test_mat = MAT(:,test_CE);
                    [AA(loop,4),~] = CV_calc_v3(test_mat,energy,a);
                else
                    AA(loop,3) = 3;               
                    test_CE = [CE kk jj];                                     
                    out_CE = setdiff(nonzero_M,test_CE);
                    test_mat = MAT(:,test_CE);
                    [AA(loop,4),~] = CV_calc_v3(test_mat,energy,a);
                end    
                if mod(loop,round(totloops/20)) == 0                    % Display percent complete
                    pcent = round(round(loop/totloops/.05)*0.05*100);
                    if pcent >= 100
                        pcent = 99.99;
                    end
                    fprintf('%2g%%...',pcent);
                    newl = newl + 1;
                end
                if mod(newl,10) < 0.001 && mod(loop,round(totloops/20)) == 0
                    fprintf('\n');
                end
            end
        end
    end
    fprintf('...done!\n\n')
    fprintf('----------------------------------------------------------------\n')
    fprintf('Adding and removing clusters along the gradient\nuntil CV score no longer decreases..\n\n')
    AA = AA(any(AA ~= 0,2),:);  % Remove all unused rows
    sorted_AA = sortrows(AA,4);
    sz_AA = size(AA,1);
  
    new_CE = CE;
    attempt = 0;
    for kk = 1:sz_AA   
        ii = sorted_AA(kk,1);
        jj = sorted_AA(kk,2);
        if at_a_time == 2
            if (kk ~= 1 && any(any(sorted_AA(1:kk-1,1:2) == ii)))==1 || (kk ~= 1 && any(any(sorted_AA(1:kk-1,1:2) == jj)))==1
                 sorted_AA(kk,[1,2]) = [-1 -1];
                continue
            end
        else
            if (kk ~= 1 && any(any(sorted_AA(1:kk-1,1:2) == ii)))==1
                 sorted_AA(kk,[1,2]) = [-1 -1];
                continue
            end
        end
        if sorted_AA(kk,3) == -3
            test_CE = setdiff(CE,[ii jj]);
            out_CE = setdiff(nonzero_M,test_CE);
            test_mat = MAT(:,test_CE);
            [new_CV_score,new_std_dev] = CV_calc_v3(test_mat,energy,a);

        elseif sorted_AA(kk,3) == -2
            test_CE = setdiff(CE,ii);
            test_CE = [test_CE jj];                        %#ok<AGROW>
            out_CE = setdiff(nonzero_M,test_CE);
            test_mat = MAT(:,test_CE);
            [new_CV_score,new_std_dev] = CV_calc_v3(test_mat,energy,a);

        elseif sorted_AA(kk,3) == -1
            test_CE = setdiff(CE,ii);
            out_CE = setdiff(nonzero_M,test_CE);
            test_mat = MAT(:,test_CE);
            [new_CV_score,new_std_dev] = CV_calc_v3(test_mat,energy,a);

        elseif sorted_AA(kk,3) == 1
            test_CE = sort([CE ii]);
            out_CE = setdiff(nonzero_M,test_CE);
            test_mat = MAT(:,test_CE);
            [new_CV_score,new_std_dev] = CV_calc_v3(test_mat,energy,a);

        elseif sorted_AA(kk,3) == 2
            test_CE = setdiff(CE,jj);
            test_CE = [test_CE ii];                        %#ok<AGROW>
            out_CE = setdiff(nonzero_M,test_CE);
            test_mat = MAT(:,test_CE);
            [new_CV_score,new_std_dev] = CV_calc_v3(test_mat,energy,a);

        elseif sorted_AA(kk,3) == 3
            test_CE = [CE ii jj];          
            out_CE = setdiff(nonzero_M,test_CE);
            test_mat = MAT(:,test_CE);
            [new_CV_score,new_std_dev] = CV_calc_v3(test_mat,energy,a);

        end

        if CV_score - new_CV_score >= tol           
            flag = 0;
            CE = sort(test_CE);
            CV_score = new_CV_score;
            std_dev = new_std_dev; 

            if sorted_AA(kk,3) == -3
                fprintf('\nCluster %d and %d removed!\n',ii,jj)
            elseif sorted_AA(kk,3) == -2
                fprintf('\nCluster %d removed and %d added!\n',ii,jj)
            elseif sorted_AA(kk,3) == -1
                fprintf('\nCluster %d removed!\n',ii)
            elseif sorted_AA(kk,3) == 1
                fprintf('\nCluster %d added!\n',ii)
            elseif sorted_AA(kk,3) == 2   
                fprintf('\nCluster %d added and %d removed!\n',ii,jj)
            elseif sorted_AA(kk,3) == 3
                fprintf('\nCluster %d and %d added!\n',ii,jj)
            end           
            fprintf('\nThe new CE is:\n')
            disp(CE)
            fprintf('Its CV score is: %8.6g eV/site\n',CV_score)
            fprintf('Its standard deviation is: %8.6g eV/site\n',std_dev)
            fprintf('----------------------------------------------------------------\n')       
        elseif CV_score - new_CV_score < tol    
           
            if flag < 2  
                if attempt == 4 
                    flag = flag + 1;
                    if flag < 2
                        fprintf('\nReached a minimum along this gradient!\n')
                        fprintf('----------------------------------------------------------------\n')
                        fprintf('----------------------------------------------------------------\n')
                    end
                    break
                end
                attempt = attempt + 1;
                continue
            elseif flag >= 2
                continue
            end
        end        
    end
    if kk == sz_AA   
        fprintf('\nAll possible additions/removals along this gradient exhausted!\n')
        fprintf('----------------------------------------------------------------\n')
    end
end
interactions = INTERACTIONS(CE,:); 

final_mat = MAT(:,CE);
final_coeffs = (final_mat'*final_mat)\(final_mat'*energy);
final_en = final_mat*final_coeffs;
final_resids = energy - final_en;
RMSE = sqrt(mean(final_resids.^2));

%%% Write the output to file "LG_ECIs.txt"
dlmwrite('LG_ECIs.txt',CE,'delimiter','\t')
dlmwrite('LG_ECIs.txt',final_coeffs,'delimiter','\t','-append')
fID = fopen('LG_ECIs.txt','a');
fprintf(fID,'%8.6f eV/site',CV_score);
fclose(fID);
fprintf('----------------------------------------------------------------\n')
fprintf('----------------------------------------------------------------\n')
fprintf('\nThe algorithm has found a local minimum!\nNo further cluster additions or removals lower the CV score\n\n')
fprintf('\nThe final CE is:\n')
disp(CE)
fprintf('Its CV score is: %8.6g eV/site\n',CV_score)
fprintf('Its standard deviation is: %8.6g eV/site\n',std_dev)
fprintf('The RMSE of the final fit is: %8.6g eV/site\n',RMSE)
fprintf('The LG ECI for this CE are:\n')
disp(final_coeffs)
%fprintf('The corresponding residuals are:\n')
%display(final_resids)
fprintf('This CE corresponds to the following interactions:\n')
disp(interactions)

clearvars -except CV_score ads_energy part_en orig_surf_normalized_counts atat_names newest_configs residuals mf_resids counts surf_normalized_counts coverages surf_energy sMax INTERACTIONS

