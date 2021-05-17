%% OMEGA: Script for automatic preprocessing
% 
% 1) Import MEG + Convert to continuous (CTF) + Refine registration
% 2) PSD on sensosrs (pre)
% 3) Notch filter (60Hz and harmonics) + High bandpass (0.3 Hz)
% 4) Detect and remove blinks and cardiac
% 5) PSD on sensosrs (post) *skip for now
% 6) Mark noisy segments
% 7) SSP for low frequencies (saccades) 1-7Hz
% 8) SSP for high frequencies (muscle) 40-400Hz
% Adds HC/PA label (healthy/patient) to every subject file
% Saves Report and Exports subject to .zip
%
% This script needs function: create_html_report.m and Omg_check_rename.m
%
% Guiomar Niso, 26 May 2016 (v1)
% Guiomar Niso, 6 May 2015 (v0)

clc; clear;

%% ==== PARAMETERS ================================================

% 1) MEG datasets storage
mydirMEG = 'C:\Users\Maria\OneDrive - Universidad Politécnica de Madrid\TFM\Pruebas_sujetos\';
% 2) Dir to save progress report
mydirBST = 'C:\Users\Maria\OneDrive - Universidad Politécnica de Madrid\TFM\Reports';

% OMEGA protocol
omegaProt = 'CAMCAN_BSTdb_v1';

% -------------------------------------------------------------------------

% Frequencies to filter with the noth (power line 60Hz and harmonics)
freqs_notch = [50:50:450];

% Filters
highpass = 0.3;
lowpass = 0; % 0: no filter

% Window length and overlap for PSD Welch method
win_length = 4; % sec
win_overlap = 50; % percentage

% =========================================================================

%% Prepare MEG files

%SubjectNames = {'sub-CC110033'};
SubjectNames = {'sub-CC110037'};

for iSubject=1:numel(SubjectNames)
       
%% ==== 1) Import MEG files =======================================

% For Brainstorm
sFiles0 = [];
% Start a new report
bst_report('Start', sFiles0);

% SELECCIONAR SUJETO
%
% Process: Select file names with tag: SUBJECT NAME
sFilesMEG = bst_process('CallProcess', 'process_select_files_data', ...
    sFiles0, [], ...
    'tag', '', ...
    'subjectname', SubjectNames{iSubject}, ...
    'condition', '');

% Process: Select file names with tag: rest
sFilesR = bst_process('CallProcess', 'process_select_tag', ...
    sFilesMEG, [], ...
    'tag', 'rest', ...
    'search', 1, ... % 1: Filename, 2: Comments
    'select', 1);  % Select only the files with the tag

if isempty(sFilesMEG), return; end

% ** Process: Snapshot: Sensors/MRI registration
bst_process('CallProcess', 'process_snapshot', ...
    sFilesR, [], ...
    'target', 1, ...  % Sensors/MRI registration
    'modality', 1, ...% MEG (All)
    'orient', 1, ...  % left
    'time', 0, ...
    'contact_time', [0, 0.1], ...
    'contact_nimage', 12, ...
    'threshold', 30, ...
    'comment', '');


%% ==== 2) PSD on sensors (before filtering) ======================

% Process: Power spectrum density (Welch)
sFilesPSDpre = bst_process('CallProcess', 'process_psd', ...
    sFilesR, [], ...
    'timewindow', [], ...
    'win_length', win_length, ...
    'win_overlap', win_overlap, ...
    'sensortypes', 'MEG, EEG', ...
    'edit', struct(...
         'Comment', 'Power', ...
         'TimeBands', [], ...
         'Freqs', [], ...
         'ClusterFuncTime', 'none', ...
         'Measure', 'power', ...
         'Output', 'all', ...
         'SaveKernel', 0));

% ** Process: Snapshot: Frequency spectrum
bst_process('CallProcess', 'process_snapshot', ...
    sFilesPSDpre, [], ...
    'target', 10, ...  % Frequency spectrum
    'modality', 1, ...  % MEG (All)
    'orient', 1, ...  % left
    'time', 0, ...
    'contact_time', [0, 0.1], ...
    'contact_nimage', 12, ...
    'threshold', 30, ...
    'comment', '');


%% ==== 3)  Notch filter + High pass (0.3 Hz) =====================

% Process: Notch filter: 60Hz 120Hz 180Hz 240Hz 300Hz 360Hz 420Hz 480Hz 540Hz 600Hz
sFilesNotch = bst_process('CallProcess', 'process_notch', ...
    sFilesR, [], ...
    'freqlist', freqs_notch, ...
    'sensortypes', 'MEG, EEG', ...
    'read_all', 0); 

% Process: High-pass:0.3Hz
sFilesNotchHigh = bst_process('CallProcess', 'process_bandpass', ...
    sFilesNotch, [], ...
    'highpass', highpass, ...
    'lowpass', lowpass, ...
    'mirror', 1, ...
    'sensortypes', 'MEG, EEG', ...
    'read_all', 0);

% Delete intermediate files (Notch) 
for iRun=1:numel(sFilesNotch)
    % Process: Delete data files
    bst_process('CallProcess', 'process_delete', ...
        sFilesNotch(iRun).FileName, [], ...
        'target', 2);  % Delete conditions
end

%% ==== 4) SSP Detect and remove blinks and cardiac ===============

% ==== SELECT RESTING AND NOISE ====

% Process: Select file names with tag: rest
sFilesRESTING = bst_process('CallProcess', 'process_select_tag', ...
    sFilesNotchHigh, [], ...
    'tag', 'rest', ...
    'search', 1, ...
    'select', 1);  % Select only the files with the tag

% Process: Select file names with tag: emptyroom
% sFilesNOISE = bst_process('CallProcess', 'process_select_tag', ...
%     sFilesMEG, [], ...
%     'tag', 'emptyroom', ...
%     'search', 1, ...
%     'select', 1);  % Select only the files with the tag

% ==== DETECT ====
for iRun=1:numel(sFilesRESTING)
    
    % Read the channel file
    ChannelMat = in_bst_channel(sFilesRESTING(iRun).ChannelFile);

    % Look for ECG channel
    iChannelECG = channel_find(ChannelMat.Channel, 'ECG');
%     if isempty(iChannelECG)
%         iChannelECG = channel_find(ChannelMat.Channel, 'EEG057');
%     end

    % Look for EOG channel
    iChannelVEOG = channel_find(ChannelMat.Channel, 'EOG062');
%     if isempty(iChannelVEOG)
%         iChannelVEOG = channel_find(ChannelMat.Channel, 'EEG058');
%     end
%     if isempty(iChannelVEOG)
%         iChannelVEOG = channel_find(ChannelMat.Channel, 'EOG');
%     end

    % heog = 'HEOG'; 
    % heog = 'EEG059';

    % Process: Detect heartbeats
    if ~isempty(iChannelECG)
        bst_process('CallProcess', 'process_evt_detect_ecg', ...
            sFilesRESTING(iRun), [], ...
            'channelname', ChannelMat.Channel(iChannelECG).Name, ...
            'timewindow', [], ...
            'eventname', 'cardiac');
    end

    % Process: Detect eye blinks
    if ~isempty(iChannelVEOG)
        bst_process('CallProcess', 'process_evt_detect_eog', ...
            sFilesRESTING(iRun), [], ...
            'channelname', ChannelMat.Channel(iChannelVEOG).Name, ...
            'timewindow', [], ...
            'eventname', 'blink');
    end
end

% ==== REMOVE SIMULTANEOUS ====

% Process: Remove simultaneous
bst_process('CallProcess', 'process_evt_remove_simult', ...
    sFilesRESTING, [], ...
    'remove', 'cardiac', ...
    'target', 'blink', ...
    'dt', 0.25, ...
    'rename', 0);

% ==== SSP (force 1st component) ====

% Process: SSP ECG: cardiac
bst_process('CallProcess', 'process_ssp_ecg', ...
    sFilesRESTING, [], ...
    'eventname', 'cardiac', ...
    'sensortypes', 'MEG', ...
    'usessp', 1, ...
    'select', 1);

% Process: SSP EOG: blink
bst_process('CallProcess', 'process_ssp_eog', ...
    sFilesRESTING, [], ...
    'eventname', 'blink', ...
    'sensortypes', 'MEG', ...
    'usessp', 1, ...
    'select', 1);

% % Process: ICA components: Infomax for EEG
% bst_process('CallProcess', 'process_ica', ...
%     sFilesRESTING, [], ...
%     'timewindow',   [], ...
%     'eventname',    '', ...
%     'eventtime',    [-0.2, 0.2], ...
%     'bandpass',     [0, 0], ...
%     'nicacomp',     20, ...
%     'sensortypes',  'EEG', ...
%     'usessp',       0, ...
%     'ignorebad',    1, ...
%     'saveerp',      0, ...
%     'method',       1);  % Infomax:    EEGLAB / RunICA

% ** Process: Snapshot: SSP projectors
bst_process('CallProcess', 'process_snapshot', ...
    sFilesRESTING, [], ...
    'target', 2, ...  % SSP projectors
    'modality', 1, ...  % MEG (All)
    'orient', 1, ...  % left
    'time', 0, ...
    'contact_time', [0, 0.1], ...
    'contact_nimage', 12, ...
    'threshold', 30, ...
    'comment', '');


%% ==== 5) PSD on sensors (after filtering) ======================

% Process: Power spectrum density (Welch)
sFilesPSDpost = bst_process('CallProcess', 'process_psd', ...
    sFilesRESTING, [], ...
    'timewindow', [], ...
    'win_length', win_length, ...
    'win_overlap', win_overlap, ...
    'sensortypes', 'MEG, EEG', ...
    'edit', struct(...
         'Comment', 'Power', ...
         'TimeBands', [], ...
         'Freqs', [], ...
         'ClusterFuncTime', 'none', ...
         'Measure', 'power', ...
         'Output', 'all', ...
         'SaveKernel', 0));
  
% ** Process: Snapshot: Frequency spectrum
bst_process('CallProcess', 'process_snapshot', ...
    sFilesPSDpost, [], ...
    'target', 10, ...  % Frequency spectrum
    'modality', 1, ...  % MEG (All)
    'orient', 1, ...  % left
    'time', 0, ...
    'contact_time', [0, 0.1], ...
    'contact_nimage', 12, ...
    'threshold', 30, ...
    'comment', '');


%% ==== 6) Mark noisy segments ====================================

% Check if they have HLC (continuous head localizer, in the questionnaire or easier
% on the channels). If they do: Detect movements ********

% % Process: Detect other artifacts
% bst_process('CallProcess', 'process_evt_detect_badsegment', ...
%     sFilesRESTING, [], ...
%     'timewindow', [], ...
%     'sensortypes', 'MEG, EEG', ...
%     'threshold', 3, ...  % 3
%     'isLowFreq', 1, ...
%     'isHighFreq', 1);


%% ==== 7) SSP for low frequencies (saccades) 1-7Hz ===============

% % Process: SSP
% bst_process('CallProcess', 'process_ssp', ...
%     sFilesRESTING, [], ...
%     'timewindow',  [], ...
%     'eventname',   '1-7Hz', ...
%     'eventtime',   [], ...
%     'bandpass',    [1, 7], ...
%     'sensortypes', 'MEG', ...
%     'usessp',      1, ...
%     'saveerp',     0, ...
%     'method',      1, ...  % PCA: One component per sensor
%     'select',      1);


%% ==== 8) SSP for high frequencies (muscle) 40-400Hz =============

% % Process: SSP
% bst_process('CallProcess', 'process_ssp', ...
%     sFilesRESTING, [], ...
%     'timewindow',  [], ...
%     'eventname',   '', ...
%     'eventtime',   [], ...
%     'bandpass',    [40, 400], ...
%     'sensortypes', 'MEG', ...
%     'usessp',      1, ...
%     'saveerp',     0, ...
%     'method',      1, ...  % PCA: One component per sensor
%     'select',      1);



%% SAVE REPORT ====================================================

% Save and display report
ReportFile = bst_report('Save', sFiles0);
bst_report('Open', ReportFile);
bst_report('Export', ReportFile, mydirBST);

end