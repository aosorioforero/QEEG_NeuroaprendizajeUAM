function DMS = f_PowerAnalysis(DMS,str_Name,str_Path)
% f_PowerAnalysis
%   This function performs the absolute and relative power analysis as
% described in Osorio-Forero et al., 2021. In short, we filtered the data
% in the different bands of interest (Delta: 1-4 Hz; Theta: 4-8 Hz;
% Alpha: 8-13 Hz; and Beta 13-30 Hz). In windows ofT = 4 cycles (T = 4/f), 
% where f being the lowest frequency per band. Then we used the summation
% of the squared values divided by the size of the window.
%   For the Relative power, we used the same procedure but we divided each
% window by the power value between 1-30 Hz.
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
%   power analysis. i.e. inside the field named by each subject, it will
%   contain the fields: AbsPow and RelPow; inside each one of them, there
%   will be the fields of m_Delta, m_Theta, m_Alpha, m_Beta which contains
%   the information of absolute and relative power per band, respectively.
%       Finally, as described in the methods, we saved the standard
%	deviation of these activities in the field PowVar.
%
% See also
% f_CoherenceAnalysis f_PowerAnalysisSupplementary


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

m_AbsPowers = [];
m_RelPowers = [];
m_PowVar = [];
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
    
    Fstop1 = m_Bands(end,1);   % First Stopband Frequency
    Fstop2 = m_Bands(end,2);   % Second Stopband Frequency

    bpFiltTotal = designfilt('bandpassiir','FilterOrder',N, ...
         'StopbandFrequency1',Fstop1,'StopbandFrequency2',Fstop2, ...
         'SampleRate',s_Fs,'DesignMethod','cheby2');
    
    for idxCh = 1:size(m_EEG,1);
        v_FilterBand = filtfilt(bpFiltBand,double(m_EEG(idxCh,:)));
        v_FilterTotal = filtfilt(bpFiltTotal,double(m_EEG(idxCh,:)));
        %% For absolute power
        % Divide in different windows sizes
        m_Temp = reshape(v_FilterBand,s_Fs*v_WindowSize(idxBand),[]);        
        m_AbsPowers(idxCh,idxBand)= mean(sum(m_Temp.^2)/v_WindowSize(idxBand));     
        m_PowVar(idxCh,idxBand)= std(sum(m_Temp.^2)/v_WindowSize(idxBand)); 
        %% For relative power
        % Divide in windows of 2s
        v_Band = sum(reshape(v_FilterBand,s_Fs*2,[]).^2)/v_WindowSize(idxBand); 
        v_Total = sum(reshape(v_FilterTotal,s_Fs*2,[]).^2)/v_WindowSize(idxBand);
        m_RelPowers(idxCh,idxBand)= mean(v_Band./v_Total);
    end
end
%% Load from DMS
if isfield(DMS,'m_DeltaP')
    % Abs Power
    m_DeltaP = DMS.(char(str_Subject)).AbsPow.m_Delta;
    m_ThetaP = DMS.(char(str_Subject)).AbsPow.m_Theta;
    m_AlphaP = DMS.(char(str_Subject)).AbsPow.m_Alpha;
    m_BetaP = DMS.(char(str_Subject)).AbsPow.m_Beta;
    
    m_DeltaP(:,end+1) = m_AbsPowers(:,1);
    m_ThetaP(:,end+1) = m_AbsPowers(:,2);
    m_AlphaP(:,end+1) = m_AbsPowers(:,3);
    m_BetaP(:,end+1)  = m_AbsPowers(:,4);
    
    % Power variability
    m_DeltaPV = DMS.(char(str_Subject)).PowVar.m_Delta;
    m_ThetaPV = DMS.(char(str_Subject)).PowVar.m_Theta;
    m_AlphaPV = DMS.(char(str_Subject)).PowVar.m_Alpha;
    m_BetaPV = DMS.(char(str_Subject)).PowVar.m_Beta;
    
    m_DeltaPV(:,end+1) = m_PowVar(:,1);
    m_ThetaPV(:,end+1) = m_PowVar(:,2);
    m_AlphaPV(:,end+1) = m_PowVar(:,3);
    m_BetaPV(:,end+1)  = m_PowVar(:,4);
    
    % Rel Power
    m_DeltaRP = DMS.(char(str_Subject)).RelPow.m_Delta;
    m_ThetaRP = DMS.(char(str_Subject)).RelPow.m_Theta;
    m_AlphaRP = DMS.(char(str_Subject)).RelPow.m_Alpha;
    m_BetaRP  = DMS.(char(str_Subject)).RelPow.m_Beta;
    
    m_DeltaRP(:,end+1) = m_RelPowers(:,1);
    m_ThetaRP(:,end+1) = m_RelPowers(:,2);
    m_AlphaRP(:,end+1) = m_RelPowers(:,3);
    m_BetaRP(:,end+1)  = m_RelPowers(:,4);
    
else % In case is the first subject
    
    % Abs Power
    m_DeltaP = m_AbsPowers(:,1);
    m_ThetaP = m_AbsPowers(:,2);
    m_AlphaP = m_AbsPowers(:,3);
    m_BetaP  = m_AbsPowers(:,4);

    % Power Variability
    m_DeltaPV = m_PowVar(:,1);
    m_ThetaPV = m_PowVar(:,2);
    m_AlphaPV = m_PowVar(:,3);
    m_BetaPV  = m_PowVar(:,4);
    
    % Rel Power
    m_DeltaRP = m_RelPowers(:,1);
    m_ThetaRP = m_RelPowers(:,2);
    m_AlphaRP = m_RelPowers(:,3);
    m_BetaRP  = m_RelPowers(:,4);
    
end

%% Saving in DMS
% Abs Power
DMS(1).(char(str_Subject)).AbsPow.m_Delta = m_DeltaP;
DMS(1).(char(str_Subject)).AbsPow.m_Theta = m_ThetaP;
DMS(1).(char(str_Subject)).AbsPow.m_Alpha = m_AlphaP;
DMS(1).(char(str_Subject)).AbsPow.m_Beta  = m_BetaP;

% Abs Power
DMS(1).(char(str_Subject)).PowVar.m_Delta = m_DeltaPV;
DMS(1).(char(str_Subject)).PowVar.m_Theta = m_ThetaPV;
DMS(1).(char(str_Subject)).PowVar.m_Alpha = m_AlphaPV;
DMS(1).(char(str_Subject)).PowVar.m_Beta  = m_BetaPV;

% Rel Power
DMS(1).(char(str_Subject)).RelPow.m_Delta = m_DeltaRP;
DMS(1).(char(str_Subject)).RelPow.m_Theta = m_ThetaRP;
DMS(1).(char(str_Subject)).RelPow.m_Alpha = m_AlphaRP;
DMS(1).(char(str_Subject)).RelPow.m_Beta  = m_BetaRP;

%% Saving single subject
if nargin<1    
    save(fullfile(str_DMSPath,str_DMSFile),'DMS','-v7.3')
end

end