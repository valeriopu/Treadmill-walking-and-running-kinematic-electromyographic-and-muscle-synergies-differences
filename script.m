%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                Extraction and Analysis of Muscle Synergies              %
%                  during Walking and Running on treadmill                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear variables 
close all
clc

%% LOADING sEMG SIGNALS

disp ('Loading sEMG signals...');

% Define Output folder:
outputFolderName = 'Results';

if ~exist(outputFolderName,'dir')
    mkdir(outputFolderName);
end
withtol
currDir = pwd;
outputDir = strcat(currDir,'/',outputFolderName);

% Load the .MAT file containing the sEMG signals:
% -----------------------------------------------
% RAW_EMG ---> structure that contains sEMG signals from 14-channels from
% 30 patients during walking and running on treadmill
% RAW_EMG.RAW_EMG_P000x_TR_01 --> sEMG signal for patient x during running
% RAW_EMG.RAW_EMG_P000x_TW_01 --> sEMG signal for patient x during walking
%
% CYCLE_TIMES --> structure with 60 fields: Each field is structured as Nx2
% (N is the number of gait cycles acquired for each trial). 
% The first column contains the touchdown incremental times in seconds 
% The second column contains the duration of each stance phase in seconds 

disp('Load the file with raw sEMG signals');
[file1,path] = uigetfile('*.mat','Select the file with raw sEMG signals');
load([path file1])
disp('Load the file with cycle times');
[file2,path] = uigetfile('*.mat','Select the file with cycle times');
load([path file2])

disp('   Done')

% Load the Functions' Library:
mySyn   = myMuscleSynergiesLibrary;

%% DATA for PRE-PROCESSING

% In  this section the aim is to define a new structure, called EMG, which
% contains the 13-channels sEMG signals related to the desired gait cycles.
% In particular, because of the initial and final 10 gait cycles linked to
% adaptation phenomenon on treadmill, and because of the choice to study
% muscle synergies for a smaller number of gait cycles, the selected cycles
% are 30 (15th - 45th).
%
% Moreover, the structure TIME will contain, for each walking and running
% tasks, the samples related to the beginning of given gate cycle and to
% the end of stance phase.
% It's important to highlight that in CYCLE_TIMES, the second column refers
% to stance phase duration in seconds (for example, for the patient 1 in
% running task, stance lasts 0.2680 s); in TIME, the second column refers
% to ending stance phase sample (for example, for the the patient 1 in
% running task, stance ends at the sample 24312, so after 550 samples after
% the beginning of that gait cycle, at the sample 23682). It's the same
% information reported in a different way.

% NOTE: because the selected cycles are 30 (15th - 45th), in the first four 
% cycle, 31 rows of the general field CYCLE_TIMES.CYCLE_TIMES_P0xx are 
% considered to have the duration of the last cycle. 

global task
global patient
global muscle_code

task = {'TR';'TW'};

patient = {'P0001';'P0002';'P0003';'P0004';'P0005';'P0006';'P0007';'P0008';'P0009';'P0010';...
           'P0011';'P0012';'P0013';'P0014';'P0015';'P0016';'P0017';'P0018';'P0019';'P0020';...
           'P0021';'P0022';'P0023';'P0024';'P0025';'P0026';'P0027';'P0028';'P0029';'P0030';};

muscle_code = {'ME','MA','FL','RF','VM','VL','ST','BF','TA','PL','GM','GL','SOL'};       

fs = 2000;                                              % fs = 2000 Hz sampling frequency
T = length(RAW_EMG.RAW_EMG_P0001_TR_01(:,1))/fs;        % acquisition time interval
n_pat = 30;                                             % number of patients
n_cycles = 30;                                          % number of cycles
num_muscles = length(muscle_code);                      % number of muscles

% Reorganizing CYCLE_TIMES and RAW_EMG structures

%CYCLE_TIMES structure
for t=1:length(task)
    for p=1:length(patient)
        if t==1
            str = ['CYCLE_TIMES_',patient{p},'_TR','_01'];
            oldField = join(str);
            newField = patient{p};
            [CYCLE_TIME.TR.(newField)] = CYCLE_TIMES.(oldField);
            %CYCLE_TIMES = rmfield(CYCLE_TIMES,oldField); 
        else
            str = ['CYCLE_TIMES_',patient{p},'_TW','_01'];
            oldField = join(str);
            newField = patient{p};
            [CYCLE_TIME.TW.(newField)] = CYCLE_TIMES.(oldField);
            %CYCLE_TIMES = rmfield(CYCLE_TIMES,oldField); 
        end
    end
end

% RAW_EMG structure
for t=1:length(task)
    for p=1:length(patient)
        if t==1
            str = ['RAW_EMG_',patient{p},'_TR','_01'];
            oldField = join(str);
            newField = patient{p};
            [RAW_SEMG.TR.(newField)] = RAW_EMG.(oldField);
            %RAW_EMG = rmfield(RAW_EMG,oldField); 
        else
            str = ['RAW_EMG_',patient{p},'_TW','_01'];
            oldField = join(str);
            newField = patient{p};
            [RAW_SEMG.TW.(newField)] = RAW_EMG.(oldField);
            %RAW_EMG = rmfield(RAW_EMG,oldField); 
        end
    end
end

% Definition of TIME structure

for t=1:length(task)
    for p=1:length(patient)
        TIME.(task{t}).(patient{p})(:,1) = round(CYCLE_TIME.(task{t}).(patient{p})(15:45,1)*fs);
        TIME.(task{t}).(patient{p})(:,2) = round(CYCLE_TIME.(task{t}).(patient{p})(15:45,2)*fs + TIME.(task{t}).(patient{p})(:,1));
    end
end

% Processing time values of every first column of the general
% RAW_SEMG.Tx.Px with a sample rate of 2000 Hz, the signal has been sampled 
% every 0.0005 seconds. 
% Some of the time values in the first column of the general RAW_SEMG.Tx.Px 
% need to be rounded to 4 decimal digits to have values every 0.0005 s.

for t=1:length(task)
    for p=1:length(patient)
        RAW_SEMG.(task{t}).(patient{p})(:,1) = round(RAW_SEMG.(task{t}).(patient{p})(:,1),4);
    end
end

% Definition of EMG structure

for t=1:length(task)
    for p=1:length(patient)
        init_s = find(RAW_SEMG.(task{t}).(patient{p})(:,1) == (TIME.(task{t}).(patient{p})(1,1)/fs));
        end_s = find(RAW_SEMG.(task{t}).(patient{p})(:,1) == (TIME.(task{t}).(patient{p})(end,1)/fs));
        EMG.(task{t}).(patient{p}) = RAW_SEMG.(task{t}).(patient{p})(init_s:end_s,:);
    end
end
        
%% GAIT PARAMETERS

% In this section the aim is to calculate swing and stance times and
% cadence during treadmill walking and running for every patient. To save
% these informations a new structure will be created. This structure will
% have three fields (stance, swing, cadence). Each of these fields contanis
% for every patient the corresponding values during walking and running
% (i.e. GAIT_parameters.stance.TR.P00xx).

for t = 1:length(task)
    for p = 1:length(patient)
        % Stance field
        GAIT_parameters.stance.(task{t}).(patient{p})= (TIME.(task{t}).(patient{p})(1:end-1,2)-TIME.(task{t}).(patient{p})(1:end-1,1))/fs;
        % Swing field
        GAIT_parameters.swing.(task{t}).(patient{p})= (TIME.(task{t}).(patient{p})(2:end,1)-TIME.(task{t}).(patient{p})(1:end-1,2))/fs;
        % Cadence field (steps/min)
        GAIT_parameters.cadence.(task{t}).(patient{p})= 2*n_cycles*60/((TIME.(task{t}).(patient{p})(end,1)-TIME.(task{t}).(patient{p})(1,1))/fs);
    end
end

fprintf('\n\nPress a key to calculate the gait parameters...');
pause

% TIME ANALYSIS

% In this section the aim is studying the gait parameters by a
% statistical point of view. For each gait parameter (stance, swing, cadence)
% walking task and running task (30 patients) will be compared representing
% boxplots and doing a t-test (p-value < 0.001 will be considered for statistical
% significant difference). 

% Mean value of stance and swing for each patient
stance_run = zeros(n_pat,1);    % stance
stance_walk = zeros(n_pat,1);
swing_run = zeros(n_pat,1);     % swing
swing_walk = zeros(n_pat,1);
cadence_run = zeros(n_pat,1);   % cadence
cadence_walk = zeros(n_pat,1);

for p = 1:length(patient)
    
    % Running task
    stance_run(p,1) = mean(GAIT_parameters.stance.TR.(patient{p}));  % stance
    swing_run(p,1) = mean(GAIT_parameters.swing.TR.(patient{p}));    % swing
    cadence_run(p,1) = GAIT_parameters.cadence.TR.(patient{p});      % cadence

    % walking task
    stance_walk(p,1) = mean(GAIT_parameters.stance.TW.(patient{p}));  % stance
    swing_walk(p,1) = mean(GAIT_parameters.swing.TW.(patient{p}));    % swing
    cadence_walk(p,1) = GAIT_parameters.cadence.TW.(patient{p});      % cadence
        
end

fprintf('\n\nPress a key for t-test during running and walking task for gait parameters...\n');
pause

% t-test between running and walking 
% (paired samples because they refers to the same patients in two different tasks)

% stance
[H_stance,p_stance] = ttest(stance_walk,stance_run,'alpha',0.001);
% swing
[H_swing,p_swing] = ttest(swing_walk,swing_run,'alpha',0.001);
% cadence
[H_cadence,p_cadence] = ttest(cadence_walk,cadence_run,'alpha',0.001);

if H_stance == 1
    fprintf(['Running stance and walking stance are statistically different (p = ',num2str(p_stance),')\n']);
end
if H_swing == 1
    fprintf(['Running swing and walking swing are statistically different (p = ',num2str(p_swing),')\n']);
end
if H_cadence == 1
    fprintf(['Running cadence and walking cadence are statistically different(p = ',num2str(p_cadence),')']);
end

fprintf('\n\nPress a key to see the boxplots of the gait parameters...')
pause

fprintf('\n')
% BOXPLOTS

% Boxplot stance
figure,
T1 = table([stance_run*1000;stance_walk*1000], categorical([repmat("Running",n_pat,1);repmat("Walking",n_pat,1)]),'VariableNames',["Stance","model"]);
T1.idx = grp2idx(T1.model); % convert categories into group indices
T1.idx(1:30) = 2;
T1.idx(31:60) = 1;
figure(1),
boxchart(T1.idx, T1.Stance); % group by index
hold on
% overlay the scatter plots
for n=1:max(unique(T1.idx))
    hs = scatter(ones(sum(T1.idx==n),1) + n-1, T1.Stance(T1.idx == n),"filled",'jitter','on','JitterAmount',0.1);
    %hs.MarkerEdgeColor = 'k';
    hs.MarkerFaceColor = '#EDB120';
    hs.MarkerFaceAlpha = 0.5;
end
set(gca,"XTick", unique(T1.idx),"XTickLabel",['Walking';'Running']),
%annotation('textbox', [0.6, 0.8, 0.1, 0.1], 'String', "Walking vs Running: p < 0.001"),
title('Stance [ms]');

% Boxplot swing
T1 = table([swing_run*1000;swing_walk*1000], categorical([repmat("Running",n_pat,1);repmat("Walking",n_pat,1)]),'VariableNames',["Swing","model"]);
T1.idx = grp2idx(T1.model); % convert categories into group indices
T1.idx(1:30) = 2;
T1.idx(31:60) = 1;
figure,
boxchart(T1.idx, T1.Swing); % group by index
hold on
% overlay the scatter plots
for n=1:max(unique(T1.idx))
    hs = scatter(ones(sum(T1.idx==n),1) + n-1, T1.Swing(T1.idx == n),"filled",'jitter','on','JitterAmount',0.1);
    %hs.MarkerEdgeColor = 'k';
    hs.MarkerFaceColor = '#EDB120';
    hs.MarkerFaceAlpha = 0.5;
end
set(gca,"XTick", unique(T1.idx),"XTickLabel",['Walking';'Running']),title('Swing [ms]');

% Boxplot cadence
T1 = table([cadence_run;cadence_walk], categorical([repmat("Running",n_pat,1);repmat("Walking",n_pat,1)]),'VariableNames',["Cadence","model"]);
T1.idx = grp2idx(T1.model); % convert categories into group indices
T1.idx(1:30) = 2;
T1.idx(31:60) = 1;
figure,
boxchart(T1.idx, T1.Cadence); % group by index
hold on
% overlay the scatter plots
for n=1:max(unique(T1.idx))
    hs = scatter(ones(sum(T1.idx==n),1) + n-1, T1.Cadence(T1.idx == n),"filled",'jitter','on','JitterAmount',0.1);
    %hs.MarkerEdgeColor = 'k';
    hs.MarkerFaceColor = '#EDB120';
    hs.MarkerFaceAlpha = 0.5;
end
set(gca,"XTick", unique(T1.idx),"XTickLabel",['Walking';'Running']),title('Cadence [steps/min]');

%% SIGNAL PROCESSING

fprintf('\n\nPress a key to proceed with the EMG signal processing...\n')
pause

fNy=fs/2;

% High-pass filter parameters
Wp_high = 35/fNy;  % Passband frequency (normalized from 0 to 1)
Ws_high = 25/fNy;  %  Stopband frequency (normalized from 0 to 1)
Rp_high = 1;  % Passband ripple (dB)
Rs_high = 60;  % Stopband attenuation (dB)
Order_high=8;

% Buttord
[bh,ah] = butter(Order_high,Wp_high,'high');
% freqz(bh,ah)

% Low-pass filter (Passband frequency 12 Hz and Order 5) parameters
Wp_low = 12/fNy;
Ws_low = 20/fNy;
Rp_low = 1;
Rs_low = 60;
Order_low=5;
[bl,al] = butter(Order_low,Wp_low);
% freqz(bl,al)

for t=1:length(task)
    for p=1:length(patient)
        s = RAW_SEMG.(task{t}).(patient{p});
        init_s = find(s(:,1) == (TIME.(task{t}).(patient{p})(1,1)/fs));
        end_s = find(s(:,1) == (TIME.(task{t}).(patient{p})(end,1)/fs));
        EMG_CYCLES.(task{t}).(patient{p}) = s(init_s:end_s,:);
    end
end

% HIGH PASS FILTER

% In this section the aim is to define a new structure, called EMG_FPH, which
% contains the 13-channels EMG signals high-pass filtered

for t = 1:length(task)
    for p = 1:length(patient)
        EMG_HPF.(task{t}).(patient{p})=filtfilt(bh,ah,EMG_CYCLES.(task{t}).(patient{p})(:,2:14));
        EMG_HPF.(task{t}).(patient{p})=[EMG_CYCLES.(task{t}).(patient{p})(:,1),EMG_HPF.(task{t}).(patient{p})];
    end
end

% RECTIFICATION

% In this section the aim is to define a new structure, called EMG_rect, which
% contains the 13-channels EMG original signals rectified.

for t = 1:length(task)
    for p = 1:length(patient)
        EMG_rect.(task{t}).(patient{p})=abs(EMG_HPF.(task{t}).(patient{p})(:,2:14));
        EMG_rect.(task{t}).(patient{p})=[EMG_CYCLES.(task{t}).(patient{p})(:,1),EMG_rect.(task{t}).(patient{p})];
    end 
end 

% LOW-PASS FILTER

% In this section the aim is to define a new structure, called EMG_LPF, which
% contains the EMG signals rectified, filtered by a low-pass filter 


for t = 1:length(task)
    for p = 1:length(patient)
        EMG_LPF.(task{t}).(patient{p})=filtfilt(bl,al,EMG_rect.(task{t}).(patient{p})(:,2:14));
        EMG_LPF.(task{t}).(patient{p})=[EMG_CYCLES.(task{t}).(patient{p})(:,1),EMG_LPF.(task{t}).(patient{p})];     
    end 
end 

% NEGATIVE VALUES TO ZEROS

for t = 1:length(task)
    for p = 1:length(patient)
        x = EMG_LPF.(task{t}).(patient{p});
        [m,n] = size(x);
        for i = 1:m
            for j = 1:n
                if x(i,j) < 0
                    x(i,j) = 0;
                end
            end
        end
        EMG_LPF.(task{t}).(patient{p}) = x;
    end
end

% TIME NORMALIZATION (LP-FILTERED EMG SIGNALS)

% In this section the aim is to time-normalize, for each patient end for each
% task, stance and swing phases to 200 points (100 samples assigned to
% stance and 100 samples assigned to swing) and to take a certain number of
% cycles (30). Furthermore a lower number of points allows lower computational
% time and cost.
% This way the interpretation of the results are indipendent from the
% absolute duration of gait events (which can be different in different
% patients) and, also, the division in two macro-phases makes the
% the temporal contribution of synergies easier to understand. 

for t = 1:length(task)
    for p = 1:length(patient)
        emgmat = EMG_LPF.(task{t}).(patient{p});
        timemat = TIME.(task{t}).(patient{p});
        EMG_LPF_SEGM.(task{t}).(patient{p}) = emg_segment(emgmat,timemat,fs,n_cycles);
    end
end

% AMPLITUDE NORMALIZATION
% In this section the aim is to normalize the envelopes in amplitude with
% respect to its global maximum.

for t =1:length(task)
    for p = 1:length(patient)
        EMG_LPF_SEGM.(task{t}).(patient{p}) = EMG_LPF_SEGM.(task{t}).(patient{p})./max(EMG_LPF_SEGM.(task{t}).(patient{p}),[],1);
    end
end


%% Plot of 30 gait cycle EMG envelopes of a patient muscle (running task)

global color_label
global color_label_1
color_label = {'#0072BD';'#A2142F'};
color_label_1 = {'#B0E0E6';'#FF6347'};

% Mean of envelope of each muscle within the 30 gait cycles for every
% patient for both tasks. No dependence from patients.
for t = 1:length(task)
    for m = 1:length(muscle_code)
        EMG_mean.(task{t}).(muscle_code{m}) = [];
        EMG_std_err.(task{t}).(muscle_code{m}) = [];
        for  p = 1: length(patient)
            x = [];
            for i = 1:n_cycles
                x = [x EMG_LPF_SEGM.(task{t}).(patient{p})((i-1)*200+1:i*200,m)];
            end
            x_std_err = std(x,0,2)/sqrt(n_cycles);
            x_mean = mean(x,2);
            EMG_mean.(task{t}).(muscle_code{m}) = [EMG_mean.(task{t}).(muscle_code{m}) x_mean];
            EMG_std_err.(task{t}).(muscle_code{m}) = [EMG_std_err.(task{t}).(muscle_code{m}) x_std_err];
        end
        EMG_mean.(task{t}).(muscle_code{m}) = EMG_mean.(task{t}).(muscle_code{m})';
        EMG_std_err.(task{t}).(muscle_code{m}) = EMG_std_err.(task{t}).(muscle_code{m})';
    end
end

% Plot of 30 gait cyle EMG envelopes of a patient muscle

fprintf('\nPress a key to plot the 30 gait cyle EMG envelopes of a patient muscle, their mean and standard deviation...\n')
pause
fprintf('Patients: 30\n')
pat = input('Write the patient number you want to consider: ');
fprintf('Muscles: ')
for m = 1:length(muscle_code)
    fprintf((muscle_code{m}))
    fprintf(' ')
end
muscle = input('\nWrite the muscle number you want to plot: ');

for t = 1:length(task)
figure,
for m = muscle
    subplot(2,1,1)
    for p = pat
        for i=1:n_cycles
            plot(EMG_LPF_SEGM.(task{t}).(patient{p})((i-1)*200+1:i*200,m),'LineWidth',0.8);
            hold on
            xline(100,'k'),xlabel('Stance                     Swing','fontweight','bold'),
            ylabel(muscle_code{m}),xlim([1 200]),ylim([0 1]),
            if t == 1
                subtitle({'Gait cycles envelopes - Running';patient{p}},'fontweight','bold','fontsize',14)    
            else
                subtitle({'Gait cycles envelopes - Walking';patient{p}},'fontweight','bold','fontsize',14) 
            end

        end
        subplot(2,1,2)
        x_mean = EMG_mean.(task{t}).(muscle_code{m})(p,:);
        x_std = EMG_std_err.(task{t}).(muscle_code{m})(p,:);
        errorbar((1:1:length(x_mean)),x_mean,x_std*sqrt(30),'color',color_label_1{t},'LineWidth',2,'CapSize',0),hold on,
        plot((1:1:length(x_mean)),x_mean,'color',color_label{t},'LineWidth',3), hold off, ylabel(muscle_code{m}),
        xline(100,'k'),xlabel('Stance                     Swing','fontweight','bold'),legend('std deviation','average')
        ylim([0 1])
        if t == 1
                subtitle({'Gait cycles envelopes - Running';patient{p}},'fontweight','bold','fontsize',14)    
            else
                subtitle({'Gait cycles envelopes - Walking';patient{p}},'fontweight','bold','fontsize',14) 
        end
    end
end
end

%% EMG analysis
% The aim of this section is to evaluate differences between running and
% walking tasks. 

% Structure to apply the t-test for every of each patient for both tasks

for t = 1:length(task)
    for p = 1:length(patient)
        for m = 1:length(muscle_code)
            x = EMG_LPF_SEGM.(task{t}).(patient{p})(:,m);
            EMG_TEST.(task{t}).(patient{p}).(muscle_code{m})= [];
            for n = 1:n_cycles
                EMG_TEST.(task{t}).(patient{p}).(muscle_code{m}) = [EMG_TEST.(task{t}).(patient{p}).(muscle_code{m}) ; x((n-1)*200+1:200*n)'];
            end
        end
    end
end

% In order to compare the two tasks a t-test is performed for every patient
% muscle  considering 30 gait cycles (p value = 0.001). In output
% the statistical significance percentage is calculated as ratio between
% the sum of non null hypothesis and the total number of hypothesis.
% For every muscle the mean and standard error of statistical significance
% percentage among the 30 patients are calculated.

fprintf('\n\nPress a key to compare EMG envelopes among the two tasks...')
pause

% t_test on envelopes

fprintf('\n\nt-test on envelopes:\n')

for m = 1:length(muscle_code)
    for p = 1:length(patient)
        H.(muscle_code{m})(p,:) = ttest(EMG_TEST.TW.(patient{p}).(muscle_code{m}),EMG_TEST.TR.(patient{p}).(muscle_code{m}),'alpha',0.001);
    end
end

fprintf('\n')

for m= 1:length(muscle_code)
    H_emg_tot.(muscle_code{m}) = mean(sum(H.(muscle_code{m}),2)*100/200);
    H_emg_tot_std.(muscle_code{m}) = std(sum(H.(muscle_code{m}),2)*100/200)/sqrt(length(patient));
    fprintf(['Muscle ',(muscle_code{m}),': statistical significance = ',...
        num2str(H_emg_tot.(muscle_code{m})),' ± ',num2str(H_emg_tot_std.(muscle_code{m})),' (mean ± std error)'])
    disp(' % of the gait cycle')
end

fprintf('\n')

for m = 1:length(muscle_code)
    H_emg_stance_tot.(muscle_code{m}) = mean(sum(H.(muscle_code{m})(:,1:100),2));
    H_emg_stance_tot_std.(muscle_code{m}) = std(sum(H.(muscle_code{m})(:,1:100),2))/sqrt(length(patient));
    fprintf(['Muscle ',(muscle_code{m}),': statistical significance = ',...
        num2str(H_emg_stance_tot.(muscle_code{m})),' ± ',num2str(H_emg_stance_tot_std.(muscle_code{m})),' (mean ± std error)'])
    disp(' % of the stance phase')
end

fprintf('\n')

for m = 1:length(muscle_code)
    H_emg_swing_tot.(muscle_code{m}) = mean(sum(H.(muscle_code{m})(:,101:200),2));
    H_emg_swing_tot_std.(muscle_code{m}) = std(sum(H.(muscle_code{m})(:,1:100),2))/sqrt(length(patient));
    fprintf(['Muscle ',(muscle_code{m}),': statistical significance = ',...
        num2str(H_emg_swing_tot.(muscle_code{m})),' ± ',num2str(H_emg_swing_tot_std.(muscle_code{m})),' (mean ± std error)'])
    disp(' % of the swing phase')
end
            

%% EXTRACTING MUSCLE SYNERGY WEIGHTS & ACTIVATIONS

% The aim of this section is to perform the muscle synergies extraction
% through the NNMF algorithm, implemented by the functions in the library
% 'myMuscleSynergiesLibrary'.

% disp('Extracting Muscle Synergy Weights & Activations...');

% Number of synergies to compute:
% -------------------------------
minSyn = 1;    % Minimum number of synergies   
maxSyn = 8;    % Maximum number of synergies
cycles = 30;   % Number of cycles for each subgroup

if maxSyn > num_muscles
    maxSyn = num_muscles;
end

n = minSyn:1:maxSyn;                 % Number of synergies to extract
num_subgroups = n_cycles/cycles;     % We concatenate 10 adjacent gait cycles

n_samples = 200;

% The following rows are commented because of big computational time due to 
% the synergy extraction. The results (MuscleSynergy_TR.mat, 
% MuscleSynergy_TW_alen.mat) are uploaded with the script. 

% for subgroups = 1:round(num_subgroups) 
%     for t = 1:length(task)
%         for p = 1:length(patient)
%             x = EMG_LPF_SEGM.(task{t}).(patient{p})(1+(subgroups-1)*n_samples*cycles:subgroups*n_samples*cycles,:)';
%             Synergy.(task{t}).(patient{p})= mySyn.SynergiesExtraction(x,n,n_samples);
%         end
%     end
%     
%     MuscleSynergy_TR{subgroups,1} = Synergy.TR;
%     MuscleSynergy_TW{subgroups,1} = Synergy.TW;
% end
% 
% fprintf('\n\nPress a key to save the results...\n\n')
% pause

% Save results structure (structGait) in output folder:
% -----------------------------------------------------
% cd(outputDir);
% save ('MuscleSynergy_TR_alen.mat','MuscleSynergy_TR')
% save ('MuscleSynergy_TW_alen.mat','MuscleSynergy_TW')
% cd(currDir);

%% VAF Curves
% VAF is a struct with 2 fields (TR and TW) for the 2 tasks, whose dimensions are
% (n_pat x maxSyn), so in this case (30x8), and which contains the VAF
% values for every number of synergies (performed by the NNMF algorithm to
% extract the muscle synergies), for each subject.

% To choose the optimal number of synergies for each subject the Threshold
% method on the VAF curves was performed: it was selected as n_opt the first
% value of synergy which overcomes this threshold.

fprintf('\n\nPress a key to load the synergy matrices...')
pause

fprintf('\n\nPress a key to plot the Variance Accounted For (VAF) curves...')
pause

% Uncomment to upload the results
load MuscleSynergy_TR_alen.mat
load MuscleSynergy_TW_alen.mat

VAF_th = 90; % Wright here the VAF threshold level (%) to choose the optimal number of synergies

VAF.TR = [];
VAF.TW = [];

for p = 1:length(patient)
    VAF.TR = [VAF.TR; MuscleSynergy_TR{1, 1}.(patient{p}).VAF];
    VAF.TW = [VAF.TW; MuscleSynergy_TW{1, 1}.(patient{p}).VAF];
end

figure()
for i=1:n_pat
    subplot(6,5,i)
    plot(VAF.TR(i,:),'-o')
    yline(VAF_th,'r'),
    xlabel('Number of synergies'),ylabel('VAF'),
    title('Subject ',num2str(i))
end
sgtitle('Running')

figure()
for i=1:n_pat
    subplot(6,5,i)
    plot(VAF.TW(i,:),'-o')
    yline(VAF_th,'r'),
    xlabel('Number of synergies'),ylabel('VAF'),
    title('Subject ',num2str(i))
end
sgtitle('Walking')

n_opt = zeros(n_pat-1,2);
%n_opt = zeros(30,2) % first column = n_opt of running, second column = n_opt for walking of the same subject
for i = 1:n_pat
    n_opt(i,1) = min(find(VAF.TR(i,:)>=VAF_th));
    n_opt(i,2) = min(find(VAF.TW(i,:)>=VAF_th));
end


%% Evaluation of number of synergies for a comparison between the two tasks

fprintf('\n\nPress a key to visualize the optimal number of synergies...')
pause

% For a comparison between walking and running synergies, it is necessary
% to select a certain number of synergies, rather the same for both tasks.
% For this purpose the mode has been selected.
% Most recurrent value in n_opt for walk and run
n_opt_mode_W = mode(n_opt(:,2));
n_opt_mode_R = mode(n_opt(:,1));

fprintf(['\nWalking: most recurrent optimal number of synergies = ',num2str(n_opt_mode_W)])
fprintf(['\nRunning: most recurrent optimal number of synergies = ',num2str(n_opt_mode_R),'\n'])

% For a further check the mean and standard deviation of the optimal
% numbers of synergies for the 30 patients have been calculated for both
% tasks.
n_opt_mean = mean(n_opt,1); % mean optimal number of synergies 
n_opt_std = std(n_opt,0,1); % standard deviation optimal number of synergies 
err_std = n_opt_std/sqrt(length(patient));     % standard error

fprintf(['\nWalking: optimal number of synergies = ',num2str(n_opt_mean(2)),'±',num2str(err_std(2)),' (mean ± std error)'])
fprintf(['\nRunning: optimal number of synergies = ',num2str(n_opt_mean(1)),'±',num2str(err_std(1)),' (mean ± std error)\n'])

VAF_4_R = [];
VAF_4_W = [];

for p=1:length(patient)
    VAF_4_R = [VAF_4_R MuscleSynergy_TR{1, 1}.(patient{p})(4).VAF];
    VAF_4_W = [VAF_4_W MuscleSynergy_TW{1, 1}.(patient{p})(4).VAF];
end

VAF_4_R_mean = mean(VAF_4_R);
VAF_4_W_mean = mean(VAF_4_W);

VAF_4_R_std = std(VAF_4_R)/sqrt(length(patient));
VAF_4_W_std = std(VAF_4_W)/sqrt(length(patient));

fprintf(['\nWalking: VAF (',num2str(n_opt_mode_W),' synergies) = ',num2str(VAF_4_W_mean),'±',num2str(VAF_4_W_std)])
disp(' % (mean ± std error)')
fprintf(['Running: VAF (',num2str(n_opt_mode_R),' synergies) = ',num2str(VAF_4_R_mean),'±',num2str(VAF_4_R_std)])
disp(' % (mean ± std error)')

%% Local VAF
% It is verified that the number of synergies chosen for a specific patient
% is accurated enough for each muscle. In particularly, if the local VAF
% was lower than a preestablished threshold ('VAF_muscle_th'), then the
% optimal number of synergies for that patient was increased by 1.

% fprintf('\n\nPress a key to verify the local accuracy of the optimal number of synergies chosen...\n\n')
% pause
% 
% VAF_muscle_th = 75;
% for i = 1:n_pat
%     if i<=9
%         isopt = eval(['find([MuscleSynergy_TR{1, 1}.P000',num2str(i),'(n_opt(',num2str(i),',1)).VAF_muscle]>=VAF_muscle_th)']);
%     else
%         isopt = eval(['find([MuscleSynergy_TR{1, 1}.P00',num2str(i),'(n_opt(',num2str(i),',1)).VAF_muscle]>=VAF_muscle_th)']);
%     end
%     if length(isopt)==num_muscles
%         fprintf('\nThe local VAF of each muscle for the running task is higher than %d perc --> the chosen n_opt is enough accurated!\n',VAF_muscle_th);
%     else
%         fprintf('\nThe local VAF of some muscles for the running task is lower than %d perc --> the chosen n_opt is not enough accurated!\n',VAF_muscle_th);
%         %n_opt(i,1) = n_opt(i,1)+1;
%     end
% end
% 
% for i = 1:n_pat
%     if i<=9
%         isopt = eval(['find([MuscleSynergy_TW{1, 1}.P000',num2str(i),'(n_opt(',num2str(i),',1)).VAF_muscle]>=VAF_muscle_th)']);
%     else
%         isopt = eval(['find([MuscleSynergy_TW{1, 1}.P00',num2str(i),'(n_opt(',num2str(i),',1)).VAF_muscle]>=VAF_muscle_th)']);
%     end
%     if length(isopt)==num_muscles
%         fprintf('\nThe local VAF of each muscle for the walking task is higher than %d perc --> the chosen n_opt is enough accurated!\n',VAF_muscle_th);
%     else
%         fprintf('\nThe local VAF of some muscles for the walking task is lower than %d perc --> the chosen n_opt is not enough accurated!\n',VAF_muscle_th);
%         %n_opt(i,2) = n_opt(i,2)+1;
%     end
% end

fprintf('\n\nPress a key to compare synergies among the two tasks...')
pause

%% MUSCLE SYNERGIES ORDERING - COSINESIMILARITY_SORTING FUNCTION

MS_TR = cosinesimilarity_sorting_TR(MuscleSynergy_TR{1, 1},n_samples,length(muscle_code),patient);
MS_TW = cosinesimilarity_sorting_TW(MuscleSynergy_TW{1, 1},n_samples,length(muscle_code),patient);

%% COMPARING MUSCLE PRIMITIVES AND MUSCLE MODULES IN DIFFERENT TASKS

% The aim of this section is to evaluate muscle synergies all over the
% population, after imposing as optimal number of synergies for every and
% each subject the statistical mode of n_opt across the population itself.
% This way it's possible to compare the muscle synergies between the 2
% tasks, through average values, standard deviations (or std errors) and
% statistical differences (if present).

%Re-organization of data
global syn_label
syn_label = {'S1';'S2';'S3';'S4'};

for t=1:length(task)
    for s=1:length(syn_label)
        weights_s = [];
        primitives_s = [];
        for p=1:length(patient) 
            if t==1                
                current_C = MS_TR.(patient{p}).C(s,:);
                current_W = MS_TR.(patient{p}).W(:,s);
            else
                current_C = MS_TW.(patient{p}).C(s,:);
                current_W = MS_TW.(patient{p}).W(:,s);
            end 
            primitives_s = [primitives_s; current_C];   % a subject in each row
            weights_s = [weights_s, current_W];         % a subject in each column
        end
        primitives.(task{t}).(syn_label{s}) = primitives_s;
        weights.(task{t}).(syn_label{s}) = weights_s;
    end
end

% Averaging across the population for each task
for t=1:length(task)
    for s=1:length(syn_label)
        primitives_mean.(task{t})(s,:) = mean(primitives.(task{t}).(syn_label{s}));
        weights_mean.(task{t})(:,s) = mean(weights.(task{t}).(syn_label{s}),2);
        primitives_stderr.(task{t})(s,:) = std(primitives.(task{t}).(syn_label{s}))/sqrt(length(patient));
        weights_stderr.(task{t})(:,s) = std(weights.(task{t}).(syn_label{s}),0,2)/sqrt(length(patient));
    end
end

%% Similarities between sinergies among the two tasks

% t-test between primitives 

fprintf('\n\nt-test on primitives:\n')

for s = 1:length(syn_label)
    [H_primitives.gait.(syn_label{s}),p_primitives.gait.(syn_label{s})] = ttest(primitives.TR.(syn_label{s}),primitives.TW.(syn_label{s}),'alpha',0.001);
    H_primitives_tot.(syn_label{s}) = sum(H_primitives.gait.(syn_label{s}))*100/200; % percentage
    fprintf(['Synergy ',num2str(s),': statistical significance in the ',num2str(H_primitives_tot.(syn_label{s}))])
    disp(' % of the gait cycle')

end

fprintf('\n')

for s = 1:length(syn_label)
    [H_primitives.stance.(syn_label{s}),p_primitives.stance.(syn_label{s})] = ttest(primitives.TR.(syn_label{s})(:,1:100),primitives.TW.(syn_label{s})(:,1:100),'alpha',0.001);
    H_primitives_stance_tot.(syn_label{s}) = sum(H_primitives.stance.(syn_label{s}))*100/100; % percentage
    fprintf(['Synergy ',num2str(s),': statistical significance in the ',num2str(H_primitives_tot.(syn_label{s}))])
    disp(' % of the stance phase')

end

fprintf('\n')

for s = 1:length(syn_label)
    [H_primitives.swing.(syn_label{s}),p_primitives.swing.(syn_label{s})] = ttest(primitives.TR.(syn_label{s})(:,101:end),primitives.TW.(syn_label{s})(:,101:end),'alpha',0.001);
    H_primitives_swing_tot.(syn_label{s}) = sum(H_primitives.swing.(syn_label{s}))*100/100; % percentage
    fprintf(['Synergy ',num2str(s),': statistical significance in the ',num2str(H_primitives_tot.(syn_label{s}))])
    disp(' % of the swing phase')
end

% CS
% fprintf('\nCosine similarity on primitives\n')
% for s=1:length(syn_label)    
%     CS_primitives.(syn_label{s}) = dot(primitives_mean.TR(s,:),primitives_mean.TW(s,:))/(norm(primitives_mean.TR(s,:))*norm(primitives_mean.TW(s,:)));
%     fprintf(['Synergy',num2str(s),': CS = ',num2str(CS_primitives.(syn_label{s})),' %\n'])
% end

%% Similarities between weights

% t_test

% fprintf('\nt-test on weights:\n')
% 
% for s = 1:length(syn_label)
%     [H_weights.(syn_label{s}),p_weights.(syn_label{s})] = ttest(weights.TR.(syn_label{s})',weights.TW.(syn_label{s})','alpha',0.001);
%     H_weights_tot.(syn_label{s}) = sum(H_weights.(syn_label{s})); % percentage
% end

% CS

fprintf('\nCosine similarity on weights\n')

for s=1:length(syn_label)    
    CS_weights.(syn_label{s}) = dot(weights_mean.TR(:,s),weights_mean.TW(:,s))/(norm(weights_mean.TR(:,s))*norm(weights_mean.TW(:,s)));
    fprintf(['Synergy',num2str(s),': CS = ',num2str(CS_weights.(syn_label{s})*100)])
    disp(' %')
end

%% HeatMap

fprintf('\nPress a key to visualize the synergy pimitives heatmaps...')
pause

hm_run = [];
hm_walk = [];
for s = 1:length(syn_label)
    hm_run = [hm_run ; primitives.TR.(syn_label{s})];
    hm_walk = [hm_walk ; primitives.TW.(syn_label{s})];
end

figure,
heatmap(hm_run),grid off,
colorMap = [linspace(1,0,256)' linspace(1,0.4470,256)' linspace(1,0.7410,256)'];
colormap(colorMap),title('Running'),
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
xlabel('Stance                                  Swing'),
figure,
heatmap(hm_walk),grid off,
colorMap = [linspace(1,0.6350,256)' linspace(1,0.0780,256)' linspace(1,0.1840,256)'];
colormap(colorMap),title('Walking'),
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
xlabel('Stance                                  Swing'),

%% Plot of average muscle synergies for both running and walking task

fprintf('\n\nPress a key to visualize synergy primitives and weights...\n')
pause

color_label = {'#0072BD';'#A2142F'};

global syn_name 
syn_name = {'Weight                               ';...
            'acceptance                               ';...
            'Propulsion                               ';...
            'Early swing                               ';...
            'Late swing                               '};

for t=1:length(task)
    figure()
    for s=1:length(syn_label)
        subplot(4,2,2*s-1)
        sub = primitives_mean.(task{t})(s,:);
        errorbar((1:1:n_samples),sub,primitives_stderr.(task{t})(s,:)*sqrt(length(patient)),'color',color_label_1{t},'LineWidth',1.5,'CapSize',0);
        hold on       
        plot((1:1:n_samples),sub,'color',color_label{t},'LineWidth',3)
        xline(100,'k'),xlim([1 200]),
        yl = ylim;
        ylim([0 yl(2)]),
        if s==1
            title('Motor primitives')
         elseif s==4
         xlabel('Stance                                Swing','fontweight','bold')
        end
        subplot(4,2,2*s)
        sub = weights_mean.(task{t})(:,s);
        bar(reordercats(categorical(muscle_code),muscle_code),sub,'FaceColor',color_label{t}),
        hold on
        x = [];
        y = [];
        for m = 1:length(muscle_code)
            x = [x; m*ones(30,1)];
            y = [y; weights.(task{t}).(syn_label{s})(m,:)'];
        end
        % overlay the scatter plots
        for n=1:length(muscle_code)
            hs = scatter(ones(sum(x==n),1) + n-1, y(x == n),"filled",'jitter','on','JitterAmount',0.1);
            %hs.MarkerEdgeColor = 'k';
            hs.MarkerFaceColor = (color_label_1{t});
            hs.MarkerFaceAlpha = 0.5;
        end
        if s == 1
            ylabel({syn_name{s},syn_name{s+1}},'Rotation',0,'fontweight','bold','fontsize',16)
        else
            ylabel(syn_name{s+1},'Rotation',0,'fontweight','bold','fontsize',16)
        end
        if s==1
            title('Motor modules')
        end
    end
    if t==1
        sgtitle('Running - Average Muscle Synergies across the population','fontweight','bold')
    else
        sgtitle('Walking - Average Muscle Synergies across the population','fontweight','bold')
    end
end


%% 5 - MUSCLE SYNERGIES PLOT (for a single subject)

fprintf('\n\nPress a key to visualize synergy primitives and weights of a single subject...\n')
pause
p = input('Write the number of the patient: '); % number of the subject considered

n_opt_4 = 4*ones(30,2);

for t=1:length(task)
    figure()
    for s=1:length(syn_label)
        subplot(4,2,2*s-1)
        if t==1
            sub = MS_TR.(patient{p});
        else
            sub = MS_TW.(patient{p});
        end        
        plot(sub.C(s,:),'LineWidth',1,'color',color_label{t}),
        xline(100,'k'),
        if s==1
            title('Motor primitives')
         elseif s==4
         xlabel('Stance                     Swing')
        end
        subplot(4,2,2*s)
        bar(reordercats(categorical(muscle_code),muscle_code),sub.W(:,s),'FaceColor',color_label{t})
        if s == 1
            ylabel({syn_name{s},syn_name{s+1}},'Rotation',0,'fontweight','bold','fontsize',16)
        else
            ylabel(syn_name{s+1},'Rotation',0,'fontweight','bold','fontsize',16)
        end
        if s==1
            title('Motor modules')
        end
    end
    if t==1
        sgtitle(['Running ',(patient{p})],'fontweight','bold')
    else
        sgtitle(['Walking ',(patient{p})],'fontweight','bold')
    end
end
