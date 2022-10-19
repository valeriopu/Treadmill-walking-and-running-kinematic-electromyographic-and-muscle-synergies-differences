function sorted_synergy_struct = cosinesimilarity_sorting_TR(muscle_synergy_struct,n_samples,n_muscles,patient)

% Inputs: 
% - muscle_synergy_struct = structure containing TR muscle synergies
% - n_samples = length of gait cycle
% - n_muscles = number of muscles considered
% - patient = string vector containing every subject code

% 'cosinesimilarity_sorting_TR' is a Matlab function that sorts the extracted
% muscle synergies (both motor modules and motor primitives) all over a
% population, selecting the same number of synergies for each and every
% subject.

% This script can be applied only when the chosen number of optimal
% muscle synergies is equal to 4 and only to this specific population.
% In fact, after the sorting step, every subject has been individually
% evaluated to check the final muscle synergies sorting; 
% in a few patients cosine similarity didn't provide a correct result and
% modules and primitives have been sorted manually. 

% Sorting method is based on cosine similarity (CS) between motor modules: 
% firstly, it's chosen a reference patient (P0007), whose muscle synergies
% are sorted as in literature, and then CS is performed between a given
% motor module of the reference and every motor module of another subject
% to choose the best match. This procedure is performed within the entire
% population. 

n_opt_mode_R = 4;

%% Definition of the reference subject (P0007)

%C
sorted_synergy_struct.P0007.C = zeros(n_opt_mode_R,n_samples);
sorted_synergy_struct.P0007.C(1,:) = muscle_synergy_struct.P0007(4).C(1,:);
sorted_synergy_struct.P0007.C(2,:) = muscle_synergy_struct.P0007(4).C(2,:);
sorted_synergy_struct.P0007.C(3,:) = muscle_synergy_struct.P0007(4).C(4,:);
sorted_synergy_struct.P0007.C(4,:) = muscle_synergy_struct.P0007(4).C(3,:);
%W
sorted_synergy_struct.P0007.W = zeros(n_muscles,n_opt_mode_R);
sorted_synergy_struct.P0007.W(:,1) = muscle_synergy_struct.P0007(4).W(:,1);
sorted_synergy_struct.P0007.W(:,2) = muscle_synergy_struct.P0007(4).W(:,2);
sorted_synergy_struct.P0007.W(:,3) = muscle_synergy_struct.P0007(4).W(:,4);
sorted_synergy_struct.P0007.W(:,4) = muscle_synergy_struct.P0007(4).W(:,3);

%Swapping muscle synergies 3-4 in order to perform cosine similarity
%following the 1-2-4-3 sorting
sorted_synergy_struct.P0007.C([3,4],:) = sorted_synergy_struct.P0007.C([4,3],:);
sorted_synergy_struct.P0007.W(:,[3,4]) = sorted_synergy_struct.P0007.W(:,[4,3]);


%% Cosine Similarity Implementation

for k=1:length(patient) %for every subject    
    for j=1:4   %for every synergy of k-th subject
        sj = sorted_synergy_struct.P0007.W(:,j);    %reference
        for i=1:4   %for j-th synergy of k-th subject
            si = muscle_synergy_struct.(patient{k})(4).W(:,i);
            CS(j,i) = dot(sj,si)/(norm(sj)*norm(si));
            if k == 7
                CS_7(j,i) = CS(j,i);
            end
        end
    [maxCS(j),indCS(j)] = max(CS(j,:));
    end
    
    % Sorting the 1° synergy of k-th subject    
    sorted_synergy_struct.(patient{k}).C(1,:) = muscle_synergy_struct.(patient{k})(4).C(indCS(1),:);
    sorted_synergy_struct.(patient{k}).W(:,1) = muscle_synergy_struct.(patient{k})(4).W(:,indCS(1));
    
    % Sorting the 2°-3°-4° synergies of k-th subject, avoiding possible
    % repetitions
    
    %Syn2
    if indCS(2) ~= indCS(1)
        sorted_synergy_struct.(patient{k}).C(2,:) = muscle_synergy_struct.(patient{k})(4).C(indCS(2),:);
        sorted_synergy_struct.(patient{k}).W(:,2) = muscle_synergy_struct.(patient{k})(4).W(:,indCS(2));
    else        
        CS(2,indCS(2)) = 0;                     
        [maxCS(2),indCS(2)] = max(CS(2,:));
        sorted_synergy_struct.(patient{k}).C(2,:) = muscle_synergy_struct.(patient{k})(4).C(indCS(2),:);
        sorted_synergy_struct.(patient{k}).W(:,2) = muscle_synergy_struct.(patient{k})(4).W(:,indCS(2));
    end

    %Syn3
    if indCS(3) ~= indCS(1) && indCS(3) ~= indCS(2)
        sorted_synergy_struct.(patient{k}).C(3,:) = muscle_synergy_struct.(patient{k})(4).C(indCS(3),:);
        sorted_synergy_struct.(patient{k}).W(:,3) = muscle_synergy_struct.(patient{k})(4).W(:,indCS(3));
    elseif indCS(3) == indCS(1)        
        CS(3,indCS(3)) = 0;                     
        [maxCS(3),indCS(3)] = max(CS(3,:));
        if indCS(3) ~= indCS(2)
            sorted_synergy_struct.(patient{k}).C(3,:) = muscle_synergy_struct.(patient{k})(4).C(indCS(3),:);
            sorted_synergy_struct.(patient{k}).W(:,3) = muscle_synergy_struct.(patient{k})(4).W(:,indCS(3));
        else            
            CS(3,indCS(3)) = 0;
            [maxCS(3),indCS(3)] = max(CS(3,:)); 
            sorted_synergy_struct.(patient{k}).C(3,:) = muscle_synergy_struct.(patient{k})(4).C(indCS(3),:);
            sorted_synergy_struct.(patient{k}).W(:,3) = muscle_synergy_struct.(patient{k})(4).W(:,indCS(3));
        end
    elseif indCS(3) == indCS(2)        
        CS(3,indCS(3)) = 0;                     
        [maxCS(3),indCS(3)] = max(CS(3,:));
        if indCS(3) ~= indCS(1)
            sorted_synergy_struct.(patient{k}).C(3,:) = muscle_synergy_struct.(patient{k})(4).C(indCS(3),:);
            sorted_synergy_struct.(patient{k}).W(:,3) = muscle_synergy_struct.(patient{k})(4).W(:,indCS(3));
        else            
            CS(3,indCS(3)) = 0;
            [maxCS(3),indCS(3)] = max(CS(3,:)); 
            sorted_synergy_struct.(patient{k}).C(3,:) = muscle_synergy_struct.(patient{k})(4).C(indCS(3),:);
            sorted_synergy_struct.(patient{k}).W(:,3) = muscle_synergy_struct.(patient{k})(4).W(:,indCS(3));
        end
    end

    %Syn4
    if indCS(4) ~= indCS(1) && indCS(4) ~= indCS(2) && indCS(4) ~= indCS(3)
        sorted_synergy_struct.(patient{k}).C(4,:) = muscle_synergy_struct.(patient{k})(4).C(indCS(4),:);
        sorted_synergy_struct.(patient{k}).W(:,4) = muscle_synergy_struct.(patient{k})(4).W(:,indCS(4));
    else
        sort_ind = unique(indCS);
        for index=1:length(sort_ind)
            a = diff(sort_ind);
            if max(a) == 1
                indCS(4) = 4;
            else
                [m,z] = max(a);
                indCS(4) = z+1;
            end
        end
        sorted_synergy_struct.(patient{k}).C(4,:) = muscle_synergy_struct.(patient{k})(4).C(indCS(4),:);
        sorted_synergy_struct.(patient{k}).W(:,4) = muscle_synergy_struct.(patient{k})(4).W(:,indCS(4));
    end 
end

% Swapping back muscle synergies 3-4
for k=1:length(patient)
    sorted_synergy_struct.(patient{k}).C([3,4],:) = sorted_synergy_struct.(patient{k}).C([4,3],:);
    sorted_synergy_struct.(patient{k}).W(:,[3,4]) = sorted_synergy_struct.(patient{k}).W(:,[4,3]);
end

%% Manual sorting of critical subjects

% P0022 - Running
sorted_synergy_struct.P0022.C([3,4],:) = sorted_synergy_struct.P0022.C([4,3],:);
sorted_synergy_struct.P0022.W(:,[3,4]) = sorted_synergy_struct.P0022.W(:,[4,3]);
% P0027 - Running 
sorted_synergy_struct.P0027.C([2,4],:) = sorted_synergy_struct.P0027.C([4,2],:);
sorted_synergy_struct.P0027.C([1,4],:) = sorted_synergy_struct.P0027.C([4,1],:);
sorted_synergy_struct.P0027.C([3,4],:) = sorted_synergy_struct.P0027.C([4,3],:);
sorted_synergy_struct.P0027.W(:,[2,4]) = sorted_synergy_struct.P0027.W(:,[4,2]);
sorted_synergy_struct.P0027.W(:,[1,4]) = sorted_synergy_struct.P0027.W(:,[4,1]);
sorted_synergy_struct.P0027.W(:,[3,4]) = sorted_synergy_struct.P0027.W(:,[4,3]);
% P0030 - Running 
sorted_synergy_struct.P0030.C([3,4],:) = sorted_synergy_struct.P0030.C([4,3],:);
sorted_synergy_struct.P0030.W(:,[3,4]) = sorted_synergy_struct.P0030.W(:,[4,3]);    

end


