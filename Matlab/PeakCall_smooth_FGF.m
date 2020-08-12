% -------------------------------------------------------------------------
%   This script accompanies the manuscript                                 
%   Simon et al., (2020) Developmental Cell                                
%   Repository available on https://github.com/therealkatlab               
%   Please consult READ_ME for more information                            
% ------------------------------------------------------------------------
% 
% 
% ------------------------------------------------------------------------
%                   Peak call in FGF treated time-lapse movies                
% ------------------------------------------------------------------------

%% Load in file from R that contains CN measurements to analyse
% Specify variable columns in data set fron R and change to array format

Time = 8 ;
SpotID = 17 ;
MtUniqueID = 20 ;
Value = 2 ;
data = embryosCNFGF ;

tTime = data(:, Time) ;
tSpotID = data(:, SpotID) ;
tMtUniqueID = data(:, MtUniqueID) ;
tValue = data(:, Value) ;

aTime = table2array(tTime) ;
aSpotID = table2array(tSpotID) ;
aMtUniqueID = table2array(tMtUniqueID) ;
aValue = table2array(tValue) ;

ndata = data(1,:) ; % make dummy data to add smoothing to
ndata.smoothed = 0 ;

%% For each unique cell calculate peaks

MtTracks = unique(aMtUniqueID) ; % Unique cells

tracks = length(MtTracks) ; 

% c = 0 ; % counter

pks = 0 ;
locs = 0 ;
w = 0 ;
p = 0 ;
ID = "NA" ;
TimeM = 0 ;

allpks = table(pks, locs, w, p, ID, TimeM) ; 
% Define variables and make dummy table


for i = 1:length(MtTracks) ;
   % c = c +i ; 
    iMtTrack = MtTracks([i]) ; % Each unique cell
    idata = data(data.MtUniqueID == iMtTrack, :) ; % Subset data on cell
    values = idata(:, 2) ; % Extract C:N values
    timeM = idata(:, 26) ; % Extract Time in min
    CN = table2array(values) ;
    if length(CN) > 4
    smoothed = smoothdata(CN, 'gaussian', 3) ;
    [pks, locs, w, p] = findpeaks(smoothed, 'MinPeakProminence', 0.15) ; 
    % Find peaks, location, width (at half prominence), prominence.
    l = length(pks) ;
    tID = idata(1:l, 20) ; % Get MtUniqueID that corresponds to these peaks
    ID = table2array(tID) ;
    T = table(pks, locs, w, p) ; 
    T.ID = [ID] ; % Bind ID with peak info
    tptime = timeM(locs, :) ; 
    ptime = table2array(tptime) ;
    T.TimeM = [ptime] ; 
    
    idata.smoothed = [smoothed] ;
    ndata = [ndata; idata] ;
    
    allpks = [allpks; T] ; % Concatenate with other cell peaks
    else
        
        
    end

end

ndata([1],:) = []; % Remove initial dummy row
allpks([1],:) = []; % Remove initial dummy row

writetable(allpks, 'smoothpks_FGF.csv')
writetable(ndata, 'smoothdat_FGF.csv')
