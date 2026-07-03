    % Notice: 
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility 
%   for its use
%
% File name: SetUpL22_14vLF_128RyLns.m.m - Example of scanline imaging with focused transmits at a single focal depth.
%
% Description: 
%   Sequence programming file for L22-14vLF Linear array, using 128 ray lines
%   (focus transmits) and 128 receive acquisitions. Of the 128 transmit
%   channels, the active transmit aperture is limited based on user-entered
%   transmit focus and f-number. All 128 receive channels are active for
%   each acquisition. This script uses 4X sampling with A/D sample rate of 
%   62.5 MHz for a 15.625 MHz processing center frequency.  Transmit is at
%   17.8 MHz and receive bandpass filter has been shifted to 18 MHz center 
%   frequency, 13.9MHz -3 dB bandwidth to support the 12 MHz bandwidth of 
%   the L22-14v (18MHz center frequency, 67% bandwidth). Processing is 
%   asynchronous with respect to acquisition.
%
% Last update:
% 12/07/2015 - modified for SW 3.0

clear all
close all
%% HAR load parameters

pushV=55;
trackVolt='55';
flagHighFs =2;
nAng=3; 
if nAng>1
    PRF=16;
else
    PRF=10;
end;
PRFcom=PRF/nAng; % 10 kHZ PRF
flag=1;

if flag ==1
     SWEI_push_loc = repmat([30 50 70 100],1,1) ;%[72 56 64];% 30 45 55 80  % interms of element number; 64 is the center element
else
    SWEI_push_loc = [22];%[30 95];% 30 45 55 80  % interms of element number; 64 is the center element
end

filename_short = 'L7-4pwHSWEI_';
phantom='Cylinder_';


% pulse sequence                       
SVns.npre   =round(60/nAng);                                % Number of reference (pre) tracking lines before first ARF push
arfpulslen=1200;
arfpulse_max=max(arfpulslen);
SVns.relaxed=round(465/nAng);  % add extra time points to see sw propagation
SVns.track_btwn_num=round(0/nAng);
SVns.na=SVns.npre+SVns.track_btwn_num+SVns.relaxed;
naSWI=SVns.na;
pushFreq= 3.968; % in Mhz; even number of 250 MHZ clock cycle
trackFreq= 6.097;% in Mhz; even number of 250 MHZ clock cycle
pushDuration_all=arfpulslen*1e-6/pushFreq;
%% HARF PUsh and track parameters
nAnglesB=61;
angRanB=30; % [deg], -angRanB/2 to angRanB/2 will be used
nAngles=nAng;
angleRange=8; % [deg], -1 to 1 will be 2

sweiFrames    = numel(SWEI_push_loc);
kApod         = 0; % 1= kaiser apodization
SVFNum        = 2.5;                              % F Number for the SV Seq [can be changed in the GUI]                          % Focal depth for SV Seq in mm [can be changed in the GUI]
SVfdepth      = 32;                             % depth of the Push Sequence
P.startDepthMm      = 2;
P.endDepthMm        = 45; 

% bmode
if (nAnglesB> 1) 
    dthetaB = (angRanB*pi/180)/(nAnglesB-1); 
    P.startAngleB = -angRanB*pi/180/2; 
else
    dthetaB = 0; 
    P.startAngleB=0;
end
% SWEI
if (nAngles > 1) 
    dtheta = (angleRange*pi/180)/(nAngles-1); 
    P.startAngle = -angleRange*pi/180/2; 
else
    dtheta = 0; 
    P.startAngle=0;
end
P.dtheta=dtheta;
P.dthetaB=dthetaB;
%% HARF post processing control
ApertureGrowth  = 1;
figs            = 1;
maxDisp         = 80E-6;
kern            = 3;
bpf             = 0;
ReceiveFNum     = 3; % this dynamic receive parameters
%% Common parameters
BmodeFrames         = 4;
maxVoltage          = 65;                               % Maximum High Voltage for the transducer
c                   = 1540;                             % Speed of Sound
numRays             = 128;
FocalDepthMm        = SVfdepth;    % Bmode and HARF has same Focal Depth
%% Define system parameters.
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;    % number of receive channels.
Resource.Parameters.speedOfSound = c;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
Resource.Parameters.connector =1;
%%
% HIFU % The Resource.HIFU.externalHifuPwr parameter must be specified in a
% script using TPC Profile 5 with the HIFU option, to inform the system
% that the script intends to use the external power supply.  This is also
% to make sure that the script was explicitly written by the user for this
% purpose, and to prevent a script intended only for an Extended Transmit
% system from accidentally being used on the HIFU system.
Resource.HIFU.externalHifuPwr = 1;

% HIFU % The string value assigned to the variable below is used to set the
% port ID for the virtual serial port used to control the external HIFU
% power supply.  The port ID was assigned by the Windows OS when it
% installed the SW driver for the power supply; the value assigned here may
% have to be modified to match.  To find the value to use, open the Windows
% Device Manager and select the serial/ COM port heading.  If you have
% installed the driver for the external power supply, and it is connected
% to the host computer and turned on, you should see it listed along with
% the COM port ID number assigned to it.
Resource.HIFU.extPwrComPortID = 'COM5';
%% Specify Trans structure array.
Trans.name = 'L7-4';
Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);  % L7-4 transducer is 'known' transducer so we can use computeTrans.
% note nominal center frequency from computeTrans is 6.25 MHz
Trans.maxHighVoltage = maxVoltage;  % set maximum high voltage limit for pulser supply.

wvlngthMm = c/(Trans.frequency*1E3);
lambda=c/(Trans.frequency*1e6);
SVfdwl = SVfdepth/(wvlngthMm);

txNumElARF = round((SVfdwl/SVFNum)/Trans.spacing/2);  % no. of elements in 1/2 aperture.
for ii=1:numel(SWEI_push_loc)
    if SWEI_push_loc(ii) < 64
        if (SWEI_push_loc(ii)-txNumElARF) < 0            
            SWEI_push_loc(ii)=txNumElARF+1;
        end
    else
        if (SWEI_push_loc(ii)+txNumElARF) > 128  