%% OMEGA: Script for automatic preprocessing
% 
% ** OMEGA_BSTdb_v1_guio **
%
%%% 1) Import MEG recordings (resting state)
% 2) Compute sources
% 3) PSD on sensors (in all freqs)
% 4) PSD on sensors (in bands)
% 5) PSD on sources (in bands)
%
% Guiomar Niso, 6 Apr 2015 (v0)
% Guiomar Niso, 30 May 2016 (v1)

clc; clear;

%% ==== PARAMETERS ================================================

% 1) MEG datasets storage
% mydirMEG = '/meg/meg2/omega/OMEGA_SCANS/';
% 2) Dir to save progress report
% mydirBST = '/meg/meg2/omega/OMEGA_BSTdb/';


% Frequency bands of interest
freq_bands = {'delta', '2, 4', 'mean';
              'theta', '5, 7', 'mean'; 
              'alpha', '8, 12', 'mean'; 
              'beta', '15, 29', 'mean'; 
              'gamma1', '30, 59', 'mean'; 
              'gamma2', '60, 90', 'mean';
              'ripple1', '90, 250', 'mean';
              'ripple2', '250, 500', 'mean'};

% Window length and overlap for PSD Welch method
win_length = 2; % sec
win_overlap = 50; % percentage


datefield = 4; % Field containing the DATE in the recordings name
runfield = 5; % Field containing the RUN in the recordings name


summary_noise = {};
summary_head = {};
summary_sources = {};
summary_psd_sensors_all = {};
summary_psd_sensors_all_relative = {};
summary_psd_sensors_bands = {};
summary_psd_sensors_bands_relative = {};
summary_psd_sources_bands = {};
summary_psd_sources_bands_relative = {};
summary_default = {};
summary_default_relative = {};

% =========================================================================

%% Prepare MEG files

sSubjects = bst_get('ProtocolSubjects');
SubjectNames = {sSubjects.Subject.Name};



for iSubject = 2:numel(SubjectNames) % 1: Group analysis

try  
%% 0) SELECT RECORDINGS

% For Brainstorm
sFiles0 = [];

% Start a new report
bst_report('Start', sFiles0);

% Process: Select file names with tag: SUBJECT NAME
sFilesMEG = bst_process('CallProcess', 'process_select_files_data', ...
    sFiles0, [], ...
    'tag', '', ...
    'subjectname', SubjectNames{iSubject}, ...
    'condition', '');

if isempty(sFilesMEG), continue; end

% Process: Select file names with tag: high
sFilesMEG = bst_process('CallProcess', 'process_select_tag', ...
    sFilesMEG, [], ...
    'tag', 'high', ...
    'search', 1, ... % 1: Filename, 2: Comments
    'select', 1);  % Select only the files with the tag

% ==== SELECT RESTING AND NOISE ====

% Process: Select file names with tag: resting
sFilesRESTING = bst_process('CallProcess', 'process_select_tag', ...
    sFilesMEG, [], ...
    'tag', 'rest', ...
    'search', 1, ... % 1: Filename, 2: Comments
    'select', 1);  % Select only the files with the tag

% Process: Select file names with tag: noise
sFilesNOISE = bst_process('CallProcess', 'process_select_tag', ...
    sFilesMEG, [], ...
    'tag', 'emptyroom', ...
    'search', 1, ... % 1: Filename, 2: Comments
    'select', 1);  % Select only the files with the tag


% Load subject study
SubjectFile = bst_get('Subject', SubjectNames{iSubject});
SubjectStudy = bst_get('StudyWithSubject', SubjectFile.FileName);

    
%% ==== 1) Compute noise covariance ===============================

for iN = 1:numel(sFilesNOISE)
    
% If already computed Ncov, skip
indN = strcmp([SubjectStudy.Condition],sFilesNOISE(iN).Condition);
sFilesNcov = SubjectStudy(indN).NoiseCov;

if isempty(sFilesNcov)

    % Process: Compute noise covariance
    sFilesNcov = bst_process('CallProcess', 'process_noisecov', ...
        sFilesNOISE(iN), [], ...
        'baseline', [], ... % All file
        'target', 1, ...
        'dcoffset', 1, ...  % Block by block, to avoid effects of slow shifts in data
        'method', 1, ...  % Full noise covariance matrix
        'copycond', 0, ... 
        'copysubj', 0);

    % ** ONLY copy 'noisecov_full.mat' to directories with same date **********
    Omg_copy_noisecov( sFilesNcov, sFilesRESTING )
    
    % Update summary report
    summary_noise{end+1} = sFilesNOISE(iN).FileName; %%%%%%
end
end


%% ANALYZE RESTING FILES

for iR = 1:numel(sFilesRESTING)

disp(sFilesRESTING(iR).Condition)

% Check session number (with file date) and Run
[SESSION,RUN] = Omg_check_session( sFilesRESTING(iR), SubjectNames{iSubject} );
    
% ==== Check if SSP have been reviewed
isReviewed = regexp(sFilesRESTING(iR).Comment, 'clean', 'once'); %%%%
if isempty(isReviewed), continue; end


%% ==== 2) Compute head model =====================================
 
% If already Computed head model, skip
indH = strcmp([SubjectStudy.Condition],sFilesRESTING(iR).Condition);
sFilesHM = SubjectStudy(indH).HeadModel;

if isempty(sFilesHM)
    
    % Process: Compute head model
    % Overlapping spheres 
    sFilesHM = bst_process('CallProcess', 'process_headmodel', ...
        sFilesRESTING(iR), [], ...
        'comment', '', ...
        'sourcespace', 1, ...
        'meg', 3, ...  % Overlapping spheres
        'eeg', 1, ...  % 
        'ecog', 1, ...  % 
        'seeg', 1, ...
        'openmeeg', struct(...
             'BemFiles', {{}}, ...
             'BemNames', {{'Scalp', 'Skull', 'Brain'}}, ...
             'BemCond', [1, 0.0125, 1], ...
             'BemSelect', [1, 1, 1], ...
             'isAdjoint', 0, ...
             'isAdaptative', 1, ...
             'isSplit', 0, ...
             'SplitLength', 4000));
      
      % Update summary report
      summary_head{end+1} = sFilesRESTING(iR).FileName; %%%%%%

end


%% ==== 3) Compute sources ========================================

% If already Computed sources , skip
sFilesSRC = {};
    
[sStudy, iStudy, iResult] = bst_get('ResultsForDataFile', sFilesRESTING(iR).FileName);
for indS=1:size(iResult)
    sFilesSRC{end+1} = sStudy.Result(iResult(indS)).FileName;
end

if isempty(sFilesSRC)
    
    % Process: Compute sources *********************************** BEST??!
    sFilesSRC = bst_process('CallProcess', 'process_inverse', ...
        sFilesRESTING(iR), [], ... 
        'comment', '', ...
        'method', 2, ...  % dSPM
        'wmne', struct(...
             'NoiseCov', [], ...
             'InverseMethod', 'wmne', ...
             'ChannelTypes', {{}}, ...
             'SNR', 3, ...
             'diagnoise', 0, ...
             'SourceOrient', {{'free'}}, ...
             'loose', 0.2, ...
             'depth', 1, ...
             'weightexp', 0.5, ...
             'weightlimit', 10, ...
             'regnoise', 1, ...
             'magreg', 0.1, ...
             'gradreg', 0.1, ...
             'eegreg', 0.1, ...
             'ecogreg', 0.1, ...
             'seegreg', 0.1, ...
             'fMRI', [], ...
             'fMRIthresh', [], ...
             'fMRIoff', 0.1, ...
             'pca', 1), ...
        'sensortypes', 'MEG', ...
        'output', 1);  % Kernel only: shared

      % Update summary report
      summary_sources{end+1} = sFilesRESTING(iR).FileName; %%%%%%

end


%% ==== 4) PSD on sensors (in all freqs) ==========================

% Process: Select file names with tag: SUBJECT NAME
sFilesPSD = bst_process('CallProcess', 'process_select_files_timefreq', ...
    sFiles0, [], ...
    'subjectname', SubjectNames{iSubject}, ...
    'condition', sFilesRESTING(iR).Condition);


% Process: Select file comments with tag: XX
sFilesPSD_sens_all = bst_process('CallProcess', 'process_select_tag', ...
    sFilesPSD, [], ...
    'tag', [SESSION,' run-',RUN,' PSD sensors_all_total'], ...
    'search', 2, ... % 1: Filename, 2: Comments
    'select', 1);  % Select only the files with the tag

% Process: Select file comments with tag: XX
sFilesPSD_sens_all_relative = bst_process('CallProcess', 'process_select_tag', ...
    sFilesPSD, [], ...
    'tag', [SESSION,' run-',RUN,' PSD sensors_all_relative'], ...
    'search', 2, ... % 1: Filename, 2: Comments
    'select', 1);  % Select only the files with the tag

if isempty(sFilesPSD_sens_all)

    % Process: Power spectrum density (Welch)
    sFilesPSD_sens_all = bst_process('CallProcess', 'process_psd', ...
        sFilesRESTING(iR), [], ... %%%%%%
        'timewindow', [], ...
        'win_length', win_length, ...
        'win_overlap', win_overlap, ...
        'sensortypes', 'MEG, EEG', ...
        'edit', struct(...
             'Comment', 'Avg,Power', ...
             'TimeBands', [], ...
             'Freqs', [], ...
             'ClusterFuncTime', 'none', ...
             'Measure', 'power', ...
             'Output', 'all', ...
             'SaveKernel', 0));

    % Process: Set comment: sensor_all
    sFilesPSD_sens_all = bst_process('CallProcess', 'process_set_comment', ...
        sFilesPSD_sens_all, [], ...
        'tag', [SESSION,' run-',RUN,' PSD sensors_all_total'], ...
        'isindex', 1);

    % Add healthy/patient tags to PSD files
    Omg_check_rename(sFilesPSD_sens_all);

    % Update summary report
    summary_psd_sensors_all{end+1} = sFilesRESTING(iR).FileName; %%%%%%

end

if isempty(sFilesPSD_sens_all_relative)
     
    % Process: Spectral flattening
    sFilesPSD_sens_all_relative = bst_process('CallProcess', 'process_tf_norm', ...
        sFilesPSD_sens_all, [], ...
        'normalize', 'relative', ...  % Relative power (divide by total power)
        'overwrite', 0);

    % Process: Set comment: sensor_all_relative
    sFilesPSD_sens_all_relative = bst_process('CallProcess', 'process_set_comment', ...
        sFilesPSD_sens_all_relative, [], ...
        'tag', [SESSION,' run-',RUN,' PSD sensors_all_relative'], ...
        'isindex', 1);
    
    % Add healthy/patient tags to PSD files
    Omg_check_rename(sFilesPSD_sens_all_relative);

    % Update summary report
    summary_psd_sensors_all_relative{end+1} = sFilesPSD_sens_all.FileName; %%%%%%

end


%% ==== 5) PSD on sensors (in bands) ==============================

% Process: Select file comments with tag: XX
sFilesPSD_sens_bands = bst_process('CallProcess', 'process_select_tag', ...
    sFilesPSD, [], ...
    'tag', [SESSION,' run-',RUN,' PSD sensors_bands_total'], ...
    'search', 2, ... % 1: Filename, 2: Comments
    'select', 1);  % Select only the files with the tag

% Process: Select file comments with tag: XX
sFilesPSD_sens_bands_relative = bst_process('CallProcess', 'process_select_tag', ...
    sFilesPSD, [], ...
    'tag', [SESSION,' run-',RUN,' PSD sensors_bands_relative'], ...
    'search', 2, ... % 1: Filename, 2: Comments
    'select', 1);  % Select only the files with the tag

if isempty(sFilesPSD_sens_bands)

    % Process: Power spectrum density (Welch)
    sFilesPSD_sens_bands = bst_process('CallProcess', 'process_psd', ...
        sFilesRESTING(iR), [], ... %%%%
        'timewindow', [], ...
        'win_length', win_length, ...
        'win_overlap', win_overlap, ...
        'sensortypes', 'MEG, EEG', ...
        'edit', struct(...
             'Comment', 'Avg,Power,FreqBands', ...
             'TimeBands', [], ...
             'Freqs', {freq_bands}, ...
             'ClusterFuncTime', 'none', ...
             'Measure', 'power', ...
             'Output', 'all', ...
             'SaveKernel', 0));

    % Process: Set comment: sensor_bands
    sFilesPSD_sens_bands = bst_process('CallProcess', 'process_set_comment', ...
        sFilesPSD_sens_bands, [], ...
        'tag', [SESSION,' run-',RUN,' PSD sensors_bands_total'], ...
        'isindex', 1);

    % Add healthy/patient tags to PSD files
    Omg_check_rename(sFilesPSD_sens_bands);
    
    % Update summary report
    summary_psd_sensors_bands{end+1} = sFilesRESTING(iR).FileName; %%%%%%

end

if isempty(sFilesPSD_sens_bands_relative)

    % Process: Spectral flattening
    sFilesPSD_sens_bands_relative = bst_process('CallProcess', 'process_tf_norm', ...
        sFilesPSD_sens_bands, [], ...
        'normalize', 'relative', ...  % Relative power (divide by total power)
        'overwrite', 0);

    % Process: Set comment: sensor_bands_relative
    sFilesPSD_sens_bands_relative = bst_process('CallProcess', 'process_set_comment', ...
        sFilesPSD_sens_bands_relative, [], ...
        'tag', [SESSION,' run-',RUN,' PSD sensors_bands_relative'], ...
        'isindex', 1);

    % Add healthy/patient tags to PSD files
    Omg_check_rename(sFilesPSD_sens_bands_relative);
    
    % Update summary report
    summary_psd_sensors_bands_relative{end+1} = sFilesPSD_sens_bands.FileName; %%%%%%

end


%% ==== 6) PSD on sources (in bands) ==============================

% Process: Select file comments with tag: XX
sFilesPSD_sources_bands = bst_process('CallProcess', 'process_select_tag', ...
    sFilesPSD, [], ...
    'tag', [SESSION,' run-',RUN,' PSD sources_bands_total'], ...
    'search', 2, ... % 1: Filename, 2: Comments
    'select', 1);  % Select only the files with the tag

% Process: Select file comments with tag: XX
sFilesPSD_sources_bands_relative = bst_process('CallProcess', 'process_select_tag', ...
    sFilesPSD, [], ...
    'tag', [SESSION,' run-',RUN,' PSD sources_bands_relative'], ...
    'search', 2, ... % 1: Filename, 2: Comments
    'select', 1);  % Select only the files with the tag

if isempty(sFilesPSD_sources_bands)
    
    if isempty(sFilesSRC), continue; end

    % Process: Power spectrum density (Welch)
    sFilesPSD_sources_bands = bst_process('CallProcess', 'process_psd', ...
        sFilesSRC, [], ...
        'timewindow', [], ...
        'win_length', win_length, ...
        'win_overlap', win_overlap, ...
        'clusters', [], ...
        'scoutfunc', 1, ...
        'edit', struct(...
             'Comment', 'Avg,Power,FreqBands', ...
             'TimeBands', [], ...
             'Freqs', {freq_bands}, ...
             'ClusterFuncTime', 'none', ...
             'Measure', 'power', ...
             'Output', 'all', ...
             'SaveKernel', 0));

    % Process: Set comment: sources_bands
    sFilesPSD_sources_bands = bst_process('CallProcess', 'process_set_comment', ...
        sFilesPSD_sources_bands, [], ...
        'tag', [SESSION,' run-',RUN,' PSD sources_bands_total'], ...
        'isindex', 1);

    % Add healthy/patient tags to PSD files
    Omg_check_rename(sFilesPSD_sources_bands);
    
    % Update summary report
    summary_psd_sources_bands{end+1} = sFilesSRC.FileName; %%%%%%
end

if isempty(sFilesPSD_sources_bands_relative)

    % Process: Spectral flattening
    sFilesPSD_sources_bands_relative = bst_process('CallProcess', 'process_tf_norm', ...
        sFilesPSD_sources_bands, [], ...
        'normalize', 'relative', ...  % Relative power (divide by total power)
        'overwrite', 0);

    % Process: Set comment: sources_bands_relative
    sFilesPSD_sources_bands_relative = bst_process('CallProcess', 'process_set_comment', ...
        sFilesPSD_sources_bands_relative, [], ...
        'tag', [SESSION,' run-',RUN,' PSD sources_bands_relative'], ...
        'isindex', 1);

    % Add healthy/patient tags to PSD files
    Omg_check_rename(sFilesPSD_sources_bands_relative);

    % Update summary report
    summary_psd_sources_bands_relative{end+1} = sFilesPSD_sources_bands.FileName; %%%%%%
end


%% ==== 7) PSD on scouts (in bands) =======================================

% MORE EFFICIENT WITH SCRIPT


%% PROJECT TO DEFAULT TEMPLATE

% ==== TOTAL PSD
% Process: Select file names with tag: SUBJECT NAME
sFilesPSD_sources_bandsDEF = bst_process('CallProcess', 'process_select_files_timefreq', ...
    sFiles0, [], ...
    'tag', [SubjectNames{iSubject},'/',SESSION,' run-',RUN,' PSD sources_bands_total'], ...
    'subjectname', 'Group_analysis', ...
    'condition', '');

if isempty(sFilesPSD_sources_bandsDEF)

    % Process: Project on default anatomy: surface
    sFilesPSD_sources_bandsDEF = bst_process('CallProcess', 'process_project_sources', ...
        sFilesPSD_sources_bands, [], ...
        'headmodeltype', 'surface');  % Cortex surface

    % Process: Spatial smoothing (3.00)
    bst_process('CallProcess', 'process_ssmooth_surfstat', ...
        sFilesPSD_sources_bandsDEF, [], ...
        'fwhm',      3, ...
        'overwrite', 1);

    % Update summary report
    summary_default{end+1} = sFilesPSD_sources_bands.FileName; %%%%%%
end

% ==== RELATIVE PSD
% Process: Select file names with tag: SUBJECT NAME
sFilesPSD_sources_bands_relativeDEF = bst_process('CallProcess', 'process_select_files_timefreq', ...
    sFiles0, [], ...
    'tag', [SubjectNames{iSubject},'/',SESSION,' run-',RUN,' PSD sources_bands_relative'], ...
    'subjectname', 'Group_analysis', ...
    'condition', '');

if isempty(sFilesPSD_sources_bands_relativeDEF)

    % Process: Project on default anatomy: surface
    sFilesPSD_sources_bands_relativeDEF = bst_process('CallProcess', 'process_project_sources', ...
        sFilesPSD_sources_bands_relative, [], ...
        'headmodeltype', 'surface');  % Cortex surface

    % Process: Spatial smoothing (3.00)
    bst_process('CallProcess', 'process_ssmooth_surfstat', ...
        sFilesPSD_sources_bands_relativeDEF, [], ...
        'fwhm',      3, ...
        'overwrite', 1);

    % Update summary report
    summary_default_relative{end+1} = sFilesPSD_sources_bands_relative.FileName; %%%%%%
end


%% FIND MAX FREQ AND AMPLITUD OF PEAK

for iF = 1:size(freq_bands,1)-2
    
% Process: Select file names with tag: SUBJECT NAME
sFilesF = bst_process('CallProcess', 'process_select_files_timefreq', ...
    sFiles0, [], ...
    'tag', [SESSION,' run-',RUN,' PSD_relative Freq ',freq_bands{iF,1}], ...
    'subjectname', SubjectNames{iSubject}, ...
    'condition', '');

if isempty(sFilesF)
 
    % Process: Find maximum amplitude:  X-XXHz, frequency
    sFilesF = bst_process('CallProcess', 'process_extract_maxfreq', ...
        sFilesPSD_sens_all_relative, [], ...
        'freqrange',   str2num(freq_bands{iF,2}), ...
        'method',      'max', ...  % Maximum value  (positive peak)
        'output',      'frequency', ...  % Frequency of the peak  (for each signal separately)
        'overwrite',   0);

    % Process: Set comment: sources_bands_relative
    sFilesF = bst_process('CallProcess', 'process_set_comment', ...
        sFilesF, [], ...
        'tag', [SESSION,' run-',RUN,' PSD_relative Freq ',freq_bands{iF,1}], ...
        'isindex', 1);

    % Add healthy/patient tags to PSD files
    Omg_check_rename(sFilesF);
end

% Process: Select file names with tag: SUBJECT NAME
sFilesA = bst_process('CallProcess', 'process_select_files_timefreq', ...
    sFiles0, [], ...
    'tag', [SESSION,' run-',RUN,' PSD_relative Ampl ',freq_bands{iF,1}], ...
    'subjectname', SubjectNames{iSubject}, ...
    'condition', '');

if isempty(sFilesA)

    sFilesA = bst_process('CallProcess', 'process_extract_maxfreq', ...
        sFilesPSD_sens_all_relative, [], ...
        'freqrange',   str2num(freq_bands{iF,2}), ...
        'method',      'max', ...  % Maximum value  (positive peak)
        'output',      'amplitude', ...  % Peak amplitude  (for each signal separately)
        'overwrite',   0);

    % Process: Set comment: sources_bands_relative
    sFilesA = bst_process('CallProcess', 'process_set_comment', ...
        sFilesA, [], ...
        'tag', [SESSION,' run-',RUN,' PSD_relative Ampl ',freq_bands{iF,1}], ...
        'isindex', 1);
    
    % Add healthy/patient tags to PSD files
    Omg_check_rename(sFilesA);
end
end
end

catch
 disp ('*Error');% ERROR
end
end

% Save and display report
ReportFile = bst_report('Save', sFiles0);
bst_report('Open', ReportFile);

% Add healthy/patient tags to PSD files
Omg_check_rename([]); 

% Summary report
disp(' %%%% ==== SUMMARY ==== %%%%');
disp('* summary_noise');
disp(summary_noise');
disp('* summary_head');
disp(summary_head');
disp('* summary_sources');
disp(summary_sources');
disp('* summary_psd_sensors_all');
disp(summary_psd_sensors_all');
disp('* summary_psd_sensors_all_relative');
disp(summary_psd_sensors_all_relative');
disp('* summary_psd_sensors_bands');
disp(summary_psd_sensors_bands');
disp('* summary_psd_sensors_bands_relative');
disp(summary_psd_sensors_bands_relative');
disp('* summary_psd_sources_bands');
disp(summary_psd_sources_bands');
disp('* summary_psd_sources_bands_relative');
disp(summary_psd_sources_bands_relative');
disp('* summary_default');
disp(summary_default');
disp('* summary_default_relative');
disp(summary_default_relative');

clear




