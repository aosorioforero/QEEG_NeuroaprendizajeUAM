function DMS = f_PowerAnalysisSupplementary(DMS,str_Name,str_Path)
% f_PowerAnalysisSupplementary
%   This function performs the absolute power analysis as described in 
% Osorio-Forero et al., 2021 in the supplementary materials. In short, 
% we filtered the data in the different bands of interest (Delta: 1-4 Hz; 
% Theta: 4-8 Hz; % Alpha: 8-13 Hz; and Beta 13-30 Hz). In fixed windows 
% of 2 s. Then we used the fuction "bandpower" in the signal processing 
% toolbox in of Matlab. These results are saved under the feature
% WelchPower.
%   In addition, the original method developed in the f_PowerAnalysis
% intended to mantain the energy contain at different frequency bands using
% the same number of cycles per band. However, fixed windows are commonly
% used in spectral analysis and comparisons across subjects. Therefore, we
% also performed the same analysis as in the f_PowerAnalysis for absolute
% powers but this time we used a fixed window of 2 s.
%
%
% Inputs:
%   DMS: structure containing the features of Multiple Subjects.
%   str_Name: String containing the name of the file to load (string).
%   str_Path: String containing the the path of the file (string).
%   Remember: the file to load should be the 
%
% Outputs
%   DMS: structure containing the Data of Multiple Subjects including the
% power analysis. i.e. inside the field named by each subject, it will
% contain the fields: WelchPower and FixedWindow; and inside of each one 
% of them, there will be the fields of m_Delta, m_Theta, m_Alpha, m_Beta as
% previously described.
%
% See also
% f_CoherenceAnalysis f_PowerAnalysis


%% Loading files
clc
disp('Loading file to analyze...')
if nargin<2 
    [str_Name, str_Path] = uigetfile('*.mat', 'Select the subject to analyze');   
end
if contains(str_Name,'C')
    str_Subject = str_Name(strfind(str_Name,'C'):strfind(str_Name,'C')+2);
else contains(str_Name,'P')
    str_Subject = str_Name(strfind(str_Name,'P'):strfind(str_Name,'P')+2);
end
st_File = matfile(fullfile(str_Path,str_Name));
if nargin<1
    disp('Loading DMS file...')
    [str_DMSFile, str_DMSPath] = uigetfile('*.mat', 'Select the DMS file');
    DMS = load(fullfile(str_DMSPath,str_DMSFile));
    DMS = DMS.DMS;
end

%% Getting started

m_EEG = st_File.data;
s_Fs =  250; % For this study, we used a sampling frequency of 250 Hz.
             % Change here if necessary.

v_BandNames = {'Delta','Theta','Alpha','Beta'};
m_Bands =  [1,4;    % Delta
            4,8;    % Theta
            8,13;   % Alpha
            13,30;  % Beta
            1,30];  % Total
        
v_WindowSize = [4       % Delta
                1       % Theta
                .12     % Alpha
                0.08];  % Beta

m_WelchPower = [];
m_FixedWindow = [];

%% Start analysis
disp('Starting power Analysis...')
for idxBand = 1:size(m_Bands,1)-1
    %% Setup filters
    disp([v_BandNames{idxBand},' power...'])
    N      = 100;  % Order
    Fstop1 = m_Bands(idxBand,1);   % First Stopband Frequency
    Fstop2 = m_Bands(idxBand,2);   % Second Stopband Frequency

    bpFiltBand = designfilt('bandpassiir','FilterOrder',N, ...
         'StopbandFrequency1',Fstop1,'StopbandFrequency2',Fstop2, ...
         'SampleRate',s_Fs,'DesignMethod','cheby2');
    
    for idxCh = 1:size(m_EEG,1);
        
        v_FilterBand = filtfilt(bpFiltBand,double(m_EEG(idxCh,:)));
        %% For Welch power
        % Divide in different windows sizes
        m_Temp = bandpower(v_FilterBand,s_Fs,[m_Bands(idxBand,1) m_Bands(idxBand,2)]);
        m_WelchPower(idxCh,idxBand)= mean(m_Temp);     
        
        %% For Fixed window power
        % Divide in windows of 2s
        m_FixedWindow(idxCh,idxBand)= mean(sum(reshape(v_FilterBand,s_Fs*2,[]).^2)/v_WindowSize(idxBand));
    
    end
end
%% Load from DMS
if isfield(DMS,'m_DeltaP')
    
    % Welch Power
    m_DeltaWP = DMS.(char(str_Subject)).WelchPower.m_Delta;
    m_ThetaWP = DMS.(char(str_Subject)).WelchPower.m_Theta;
    m_AlphaWP = DMS.(char(str_Subject)).WelchPower.m_Alpha;
    m_BetaWP = DMS.(char(str_Subject)).WelchPower.m_Beta;
    
    m_DeltaWP(:,end+1) = m_WelchPower(:,1);
    m_ThetaWP(:,end+1) = m_WelchPower(:,2);
    m_AlphaWP(:,end+1) = m_WelchPower(:,3);
    m_BetaWP(:,end+1)  = m_WelchPower(:,4);
    
    % Fixed window Power
    m_DeltaFW = DMS.(char(str_Subject)).FixedWindow.m_Delta;
    m_ThetaFW = DMS.(char(str_Subject)).FixedWindow.m_Theta;
    m_AlphaFW = DMS.(char(str_Subject)).FixedWindow.m_Alpha;
    m_BetaFW  = DMS.(char(str_Subject)).FixedWindow.m_Beta;
    
    m_DeltaFW(:,end+1) = m_FixedWindow(:,1);
    m_ThetaFW(:,end+1) = m_FixedWindow(:,2);
    m_AlphaFW(:,end+1) = m_FixedWindow(:,3);
    m_BetaFW(:,end+1)  = m_FixedWindow(:,4);
    
else % In case is the first subject
    
    % Welch Power
    m_DeltaWP = m_WelchPower(:,1);
    m_ThetaWP = m_WelchPower(:,2);
    m_AlphaWP = m_WelchPower(:,3);
    m_BetaWP  = m_WelchPower(:,4);

    % Fixed Window Power
    m_DeltaFW = m_FixedWindow(:,1);
    m_ThetaFW = m_FixedWindow(:,2);
    m_AlphaFW = m_FixedWindow(:,3);
    m_BetaFW  = m_FixedWindow(:,4);
    
end

%% Saving in DMS
% Welch Power
DMS(1).(char(str_Subject)).WelchPower.m_Delta = m_DeltaWP;
DMS(1).(char(str_Subject)).WelchPower.m_Theta = m_ThetaWP;
DMS(1).(char(str_Subject)).WelchPower.m_Alpha = m_AlphaWP;
DMS(1).(char(str_Subject)).WelchPower.m_Beta  = m_BetaWP;

% Fixed Window Power
DMS(1).(char(str_Subject)).FixedWindow.m_Delta = m_DeltaFW;
DMS(1).(char(str_Subject)).FixedWindow.m_Theta = m_ThetaFW;
DMS(1).(char(str_Subject)).FixedWindow.m_Alpha = m_AlphaFW;
DMS(1).(char(str_Subject)).FixedWindow.m_Beta  = m_BetaFW;

%% Saving single subject
if nargin<1    
    save(fullfile(str_DMSPath,str_DMSFile),'DMS','-v7.3')
end

end