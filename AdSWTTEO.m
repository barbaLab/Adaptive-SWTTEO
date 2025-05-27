function [ts,p2pamp,pmin,pW,E] = AdSWTTEO(data, fs, varargin)
%
% New adaptive threshold SWTTEO, same procedure as normal SWTTEO but the
% threshold on the energy is computed locally.
%
%   Usage:
%       [ts,p2pamp,pmin,pW,E] = AdSWTTEO(data, fs);  % using default parameters
%       
%       % Using name-value pairs for custom parameters
%       [ts,p2pamp,pmin,pW,E] = AdSWTTEO(data, fs, 'MultCoeff', 8, 'RefrTime', 1);
%
%   Input parameters:
%       data:   Input single-channel signal (1xN or Nx1 vector)
%       fs:     Sampling frequency in Hz
%       
%       Optional parameters (name-value pairs or structure):
%           wavLevel:   Wavelet decomposition level (default: 2)
%           waveName:   Wavelet type (default: 'sym5')
%           winType:    Function handle for smoothing window (default: @hamming)
%           smoothN:    Size of the smoothing operator in ms (default: 1, set to 0 to turn off)
%           winPars:    Optional parameters for the smoothing window (default: {'symmetric'})
%           RefrTime:   Refractory time in ms (default: 0.5)
%           MultCoeff:  Multiplication coefficient for SWTTEO thresholding (default: 6)
%           Polarity:   Spike polarity: -1 for negative, 1 for positive (default: -1)
%           PeakDur:    Max peak duration (pulse lifetime) in ms (default: 2.3)
%           medWdw:     Window length for moving quantile in seconds (default: 0.003)
%           k:          Number of samples for TEO (default: 3)
%           overshoot:  Secondary window length for boundary control in seconds (default: 0.0007)
%
%   Output parameters:
%       p2pamp: Peak-to-peak amplitude  
%       pmin: Amplitude in volt at spike's positions  
%       ts: Timestamps of the detected spikes stored in 1xNspks vector  
%       pW: Pulse width of detected spikes  
%       E: Energy values from the SWTTEO algorithm
%
%   Description:
%       AdSWTTEO(data,pars) computes the location of action potential in noisy
%       MEA measurements. This method is based on the work of Lieb et al. 2017
%       "A Stationary Wavelet Transform and a Time-Frequency Based Spike
%       Detection Algorithm for Extracellular Recorded Data". The algorithm
%       therein was further improved by using a windowed quantile and a peak
%       refinement strategy. Modify 'MultCoeff' parameter to increase/decrease the
%       threshold. A larger value results in fewer detected spikes but also the
%       number of false positives decrease. Decreasing this factor makes it more
%       sensitive to detect spikes.
%
%   References:
%       tbd.
%
%
%   Author: Tommaso Lambresa 19/06/2024

% Default parameters
defaults.wavLevel   = 2;              % Wavelet decomposition level
defaults.waveName   = 'sym5';         % Wavelet type
defaults.winType    = @hamming;       % Function handle for smoothing window
defaults.smoothN    = 1;              % [ms] Size of the smoothing operator. Set to 0 to turn off
defaults.winPars    = {'symmetric'};  % Optional parameters for the smoothing window
defaults.RefrTime   = 0.5;            % [ms] Refractory time
defaults.MultCoeff  = 6;              % Multiplication coefficient for SWTTEO thresholding
defaults.Polarity   = -1;             % {-1; 1} Spike polarity: -1 for negative, 1 for positive
defaults.PeakDur    = 2.3;            % [ms] Max peak duration (pulse lifetime)
defaults.medWdw     = 0.003;          % [s] Window length for moving quantile
defaults.k          = 3;              % Number of samples for TEO
defaults.overshoot  = 0.0007;         % [s] Secondary window length for boundary control

% Create input parser
p = inputParser;

addRequired(p, 'data', @(x) isrow(x) || iscolumn(x));
addRequired(p, 'fs', @(x) isnumeric(x) && isscalar(x) && x > 0);

% Add parameters
addParameter(p, 'wavLevel', defaults.wavLevel, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'waveName', defaults.waveName, @ischar);
addParameter(p, 'winType', defaults.winType, @(x) isa(x, 'function_handle'));
addParameter(p, 'smoothN', defaults.smoothN, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'winPars', defaults.winPars, @iscell);
addParameter(p, 'RefrTime', defaults.RefrTime, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'MultCoeff', defaults.MultCoeff, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'Polarity', defaults.Polarity, @(x) isnumeric(x) && isscalar(x) && (x == -1 || x == 1));
addParameter(p, 'PeakDur', defaults.PeakDur, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'medWdw', defaults.medWdw, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'k', defaults.k, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'overshoot', defaults.overshoot, @(x) isnumeric(x) && isscalar(x) && x > 0);

% Parse inputs
parse(p, data, fs, varargin{:});

data = p.Results.data;
fs = p.Results.fs;
pars = rmfield(p.Results, {'data', 'fs'});

pars.smoothN = round(pars.smoothN*1e-3*fs);

% Defining TEO operator
TEO = @(x,k) (x.^2 - TEOshift(x,[-k, 0]).*TEOshift(x,[k, 0]));

L = length(data);
data = data(:);     % ensure input is column

% Zero padding if the L is not divisible by a power of two
% TODO use nextpow2 here
pow = 2^pars.wavLevel;
if rem(L,pow) > 0
    Lok = ceil(L/pow)*pow;
    Ldiff = Lok - L;
    data = [data; zeros(Ldiff,1)];
end

%vectorized version:
[lo_D,~] = wfilters(pars.waveName);
out_ = zeros(size(data));
ss = data;
for k=1:pars.wavLevel
    %Extension
    lf = length(lo_D);
    ss = extendswt(ss,lf);
    %convolution
    swa = conv(ss,lo_D','valid');
    swa = swa(2:end,:); %even number of filter coeffcients
    %apply teo to swt output
    temp = abs(TEO(swa,pars.k));
    %smoothing
    if pars.smoothN
        wind = window(pars.winType,pars.smoothN,pars.winPars{:});
        temp2 = conv(temp,wind','same');
    else
        temp2 = temp;
    end
    out_ = out_ + temp2;
    %dyadic upscaling of filter coefficients
    lo_D = dyadup(lo_D,0,1);
    %updates
    ss = swa;
end
clear('ss');

% Standard detection
lambda_swtteo   = moving_quantile(out_,.99,round(pars.medWdw*fs));
lambda_data      =  pars.MultCoeff*median(abs(data));
data_th = zeros(size(data));
data_th(out_>lambda_swtteo) = pars.Polarity .* data(out_>lambda_swtteo);
minTime = 1e-3*pars.RefrTime; % parameter in milliseconds
[ts_tmp,~] = faster_findpeaks(abs(data_th),round(minTime*fs));

%% GET MAX VALUE AROUND THE DETECTION
h = 1;
for i = 1:length(ts_tmp)
    window_size = round(0.001 * fs);
    interval = max(1, ts_tmp(i) - window_size):min(length(data), ts_tmp(i) + window_size);
    signal_window = data(interval);
    [pmin_tmp, locs] = max(pars.Polarity * signal_window);
    pmin_tmp = pars.Polarity * pmin_tmp;
    % Check the boundaries of the main window
    if locs == 1 || locs == length(signal_window)
        % Define secondary, smaller interval
        overshoot = pars.overshoot;
        secondary_size = round(overshoot * fs); % Half the size of the main window
        if locs == 1
            secondary_interval = max(1, interval(1) - secondary_size):interval(1);
            % Search for the maximum in the secondary interval
            secondary_window = data(secondary_interval);
            [pmin_tmp_secondary, locs_secondary] = max(pars.Polarity * secondary_window);
            pmin_tmp_secondary = pars.Polarity * pmin_tmp_secondary;
            % Update peak info if the secondary search improves the result
            if pars.Polarity < 0
                if pmin_tmp_secondary < pmin_tmp
                    locs = locs-(length(secondary_window)-locs_secondary);
                    pmin_tmp = pmin_tmp_secondary;
                end
            else
                if pmin_tmp_secondary > pmin_tmp
                    locs = locs-(length(secondary_window)-locs_secondary);
                    pmin_tmp = pmin_tmp_secondary;
                end
            end
        else
            secondary_interval = interval(end):min(length(data), interval(end) + secondary_size);
            % Search for the maximum in the secondary interval
            secondary_window = data(secondary_interval);
            [pmin_tmp_secondary, locs_secondary] = max(pars.Polarity * secondary_window);
            pmin_tmp_secondary = pars.Polarity * pmin_tmp_secondary;
            % Update peak info if the secondary search improves the result
            if pars.Polarity < 0
                if pmin_tmp_secondary < pmin_tmp
                    locs = locs+locs_secondary;
                    pmin_tmp = pmin_tmp_secondary;
                end
            else
                if pmin_tmp_secondary > pmin_tmp
                    locs = locs+locs_secondary;
                    pmin_tmp = pmin_tmp_secondary;
                end
            end
        end
    end
    if pars.Polarity < 0
        if pmin_tmp < pars.Polarity * lambda_data
            ts(h) = interval(1) + locs - 1;
            pmin(h) = pmin_tmp;
            h = h + 1;
        end
    else
        if pmin_tmp > pars.Polarity * lambda_data
            ts(h) = interval(1) + locs - 1;
            pmin(h) = pmin_tmp;
            h = h + 1;
        end
    end
end
if h > 1
    [ts, ia] = unique(ts, 'stable');
    pmin = pmin(ia);
    %% DOUBLE CHECK FOR REFRACTORY PERIOD
    period = diff(ts);
    candidates = find(period<round(minTime*fs));
    if ~isempty(candidates)
        if abs(pmin(candidates))>abs(pmin(candidates+1))
            ts(candidates+1)=[];
            pmin(candidates+1) = [];
        else
            ts(candidates)=[];
            pmin(candidates) = [];
        end
    end
    E = out_(ts);
    %% GET PEAK-TO-PEAK VALUES
    PLP = round(pars.PeakDur*1e-3*fs); % from ms to samples
    tloc = repmat(ts,2*PLP+1,1) + (-PLP:PLP).';
    tloc(tloc < 1) = 1;
    tloc(tloc > numel(data)) = numel(data);
    [pmax,~] = max(data(tloc));
    p2pamp = pmax + pmin;
    %% Get peak width and exlude pw > pars.PeakDur
    tlocmin = flipud( repmat(ts,PLP+1,1) + (-PLP:0).');
    tlocmin(tlocmin < 1) = 1;
    tlocmin(tlocmin > numel(data)) = numel(data);
    tlocmax = repmat(ts,PLP,1) + (1:PLP).';
    tlocmax(tlocmax < 1) = 1;
    tlocmax(tlocmax > numel(data)) = numel(data);
    for ii=1:size(tlocmin,2)
        [thispeak, ~] = faster_findpeaks(data(tlocmin(:,ii)),[],[],1);
        if isempty(thispeak)
            thispeak = PLP+1;
        end
        Imax1(ii) = 1-thispeak;
        [thispeak,~] = faster_findpeaks(data(tlocmax(:,ii)),[],[],1);
        if isempty(thispeak)
            thispeak = PLP;
        end
        Imax2(ii) = thispeak;
    end
    pW = Imax2 - Imax1;
    %% EXCLUDE VALUES
    pw_ex = pW>PLP;
    pm_ex = pmax<=0;
    ex = pw_ex | pm_ex;
    ts(ex) = [];
    p2pamp(ex) = [];
    pmax(ex) = [];
    pmin(ex) = [];
    E(ex) = [];
    pW(ex) = [];
else
    ts = [];
    p2pamp = [];
    pmax = [];
    pmin = [];
    E = [];
    pW = [];
end



%% FUNCTIONS
function y = extendswt(x,lf)
%EXTENDSWT extends the signal periodically at the boundaries
[r,c] = size(x);
y = zeros(r+lf,c);
y(1:lf/2,:) = x(end-lf/2+1:end,:);
y(lf/2+1:lf/2+r,:) = x;
y(end-lf/2+1:end,:) = x(1:lf/2,:);

function X = TEOshift(Y,k)
%circshift without the boundary behaviour...
colshift = k(1);
rowshift = k(2);
temp  = circshift(Y,k);
if colshift < 0
    temp(end+colshift+1:end,:) = flipud(Y(end+colshift+1:end,:));
elseif colshift > 0
    temp(1:1+colshift-1,:) = flipud(Y(1:1+colshift-1,:));
else
end
if rowshift<0
    temp(:,end+rowshift+1:end) = fliplr(Y(:,end+rowshift+1:end));
elseif rowshift>0
    temp(:,1:1+rowshift-1) = fliplr(Y(:,1:1+rowshift-1));
else
end
X = temp;

function quant_X = moving_quantile(X, p, n)
%% moving_quantile - calculate moving quantile
%   Input parameters:
%      X - 1-dimensional array
%      p - Quantile to calculate. Must be between 0 and 1 (both included).
%      n - Width of the moving window
%
% Written by Tommaso Lambresa 2024

quant_X = zeros(size(X));

% Padding if signal length is not divisible for n
remainder = mod(length(X), n);
num_windows = floor(length(X) / n)+1;
X_ax_pad = [X; zeros(n - remainder,1)];

% Reshape the signal into a matrix where each column is a window
X_matrix = reshape(X_ax_pad', n, num_windows);

% Compute the quantile for each column
quant = quantile(X_matrix, p, 1);

% Check if the p-th quantile is the maximum value, if so replace with the last
% but one
max_Xwin = max(X_matrix,[],1);
critical_values = find(max_Xwin == quant);
% Replace the last value with the last but one
if ~isempty(critical_values)
    data2sort = X_matrix(:,critical_values);
    sorted_win = sort(data2sort);
    quant(critical_values)=sorted_win(end-1,:);
end
quant_vector = repmat(quant, n, 1);  % Repeat each quantile value for each sample in the window
quant_vector = quant_vector(:);  % Convert to a row vector
quant_X = quant_vector(1:length(X));


function [locs, pks]=faster_findpeaks(x,minpeakdist,minpeakh,npeaks)
% Alternative to the findpeaks function.  This thing runs much much faster.
%
% x is a vector input (generally a timecourse)
% minpeakdist is the minimum desired distance between peaks (optional, defaults to 1)
% minpeakh is the minimum height of a peak (optional)
% npeaks is the number of peak to detect in order of appearance
%
% modified by Tommaso Lambresa 2024, introduced the "npeaks" input

switch nargin
    case 1
        minpeakdist = 1;
        minpeakh = -inf;
        npeaks = [];
    case 2
        minpeakh = -inf;
        npeaks = [];
    case 3
        npeaks = [];
    case 4
        ...
    otherwise
        error("not ienough input arguments");
end
if isempty(minpeakdist)
    minpeakdist = 1;
end
if isempty(minpeakh)
    minpeakh = -inf;
end
if size(x,2)==1, x=x'; end
% Find all maxima and ties
locs=find(x(2:end-1)>=x(1:end-2) & x(2:end-1)>=x(3:end))+1;
locs(x(locs)<=minpeakh)=[];
if minpeakdist>1
    while 1
        del=diff(locs)<minpeakdist;
        if ~any(del), break; end
        pks=x(locs);
        [~, mins]=min([pks(del) ; pks([false del])]);
        deln=find(del);
        deln=[deln(mins==1) deln(mins==2)+1];
        locs(deln)=[];
    end
end
if ~isempty(npeaks) && ~isempty(locs)
    locs = locs(1:npeaks);
end
if nargout>1
    pks=x(locs);
end



