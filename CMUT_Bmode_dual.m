clear all; clc;

P.startDepth = 8e-3; % start depth (m)
P.endDepth = 60e-3;  % end depth (m)
P.numRays = 256;   % number of lines per image
P.focus = 42e-3; %20.0e-3;   % focal depth (m)
P.fps = 5; %Frame Rate
P.prf = round((P.numRays*P.fps)); % Hz -- 1.05 is to give a bit of room for transferring data and beamforming

P.fnum = 2;  % larger F-number == broader focus, more depth of field, lower pressure
             % smaller F-number == tighter focus, shallow depth of field, higher pressure

usrSParam.maxVol = 50;%32;%67;             % what max, approximate Vpp is desired?
usrSParam.freqTx = 4.5;            % desired transmit frequency in MHz
usrSParam.freqRx = 4.5;% 3.5;%14;             % desired receive frequency in MHz
usrSParam.dcBias = -50;            % desired dc bias
usrSParam.folderModifiers = 'canDelete';    % any descriptive save info you want in the foldername

%usrElement.pitch = 97.4280*10^-6; % set CMUT parameters (m) (97.4280*10^-6 for UWB)
usrElement.pitch = 86.6E-6; %For dual CMUT
usrElement.kerf = 0e-6; % set CMUT parameters -- should be 0

usrTXParam.numCyc = 0.5; %Number of full (bi-phasic) cycles applied
usrTXParam.eqPulse = 0; % 0 if trailing or leading pulse is uneccessary
usrTXParam.percentOn = 1; % 0.67 is good for sine, 1 is good for square
usrTXParam.polarityFirstHalf = 1; % % Polarity of the first half cycle, -1 for negative, 1 for positive
usrTXParam.freq = usrSParam.freqTx; % pullin data over to this struct

% Define system parameters.
Resource.Parameters.numTransmit = 256;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 256;    % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;            %Error, warning, and status messages are displayed. (default)
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.Connector = [1,2]; % 1 (left), 2 (right) connectors activated -- [1,2] for both
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'Dual-256-UWB-Probe-Bernuro17';
Trans.units = 'wavelengths'; % mm or wavelengths;string % Explicit declaration avoids warning message when selected by default
Trans.id = 1877; % L12-5 38mm (cable with hole drilled) -- id will change eventually
Trans.connType = 1; % =1 means Connector is standard HDI-format connector; =2 means "break-out board" and needs more scrupt
Trans.frequency = usrSParam.freqTx;  % center frequency in megahertz;double
Trans.type = 0;  % 0=Lin(y=z=0), 1=CurvedLin(y=0), 2=2D(z=0);double
Trans.numelements = Resource.Parameters.numTransmit;  % number of transducer elements;double
Trans.elementWidth = (usrElement.pitch-usrElement.kerf)/(Resource.Parameters.speedOfSound/Trans.frequency/1E6); 
% (spacing-kerf)-> wavelengths width in mm or wavelengths (spacingkerf);double
Trans.spacing = usrElement.pitch/(Resource.Parameters.speedOfSound/Trans.frequency/1E6);  % element spacing in wavelengths; double
Trans.ElementPos = zeros(Trans.numelements,4); % position(wvlngths),angle(rad) (x,y,z,alpha); [numelex4 double]
Trans.ElementPos(:,1) = (-(Trans.numelements/2 - 0.5):1:(Trans.numelements/2 - 0.5))*Trans.spacing;
Theta = (-pi/2:pi/100:pi/2); % from pg 16 in prog manual
X = Trans.elementWidth*pi*sin(Theta); % from pg 16 in prog manual
Trans.ElementSens = abs(cos(Theta).*(sin(X)./X)); % from pg 16 in prog manual
Trans.ElementSens(isnan(Trans.ElementSens)) = 1; % fixes divide by zero
Trans.Bandwidth = [3,10];%[0.7,1.3];%[0.003,300]; %  % high and low BW (MHz) relative to center frequency ; set to be absurdly large numbers to disable this filter no matter the tx freq
Trans.maxHighVoltage = usrSParam.maxVol;  % max. high voltage limit (optional).; double
Trans.impedance = 1E6;

%For UWB device ************************
%Trans.Connector = [80,177,112,145,73,184,105,152,72,185,104,153,65,192,97,160,66,191,98,159,71,186,103,154,74,183,106,151,79,178,111,146,78,179,110,147,75,182,107,150,70,187,102,155,67,190,99,158,68,189,100,157,69,188,101,156,76,181,108,149,77,180,109,148,93,140,85,132,92,139,86,133,91,138,87,134,90,137,88,135,89,144,81,136,96,143,82,129,95,142,83,130,94,141,84,131,122,175,114,167,123,174,115,166,126,171,118,163,127,170,119,162,128,169,120,161,125,172,117,164,124,173,116,165,121,176,113,168,40,241,48,249,37,244,45,252,36,245,44,253,33,248,41,256,34,247,42,255,35,246,43,254,38,243,46,251,39,242,47,250,3,212,13,222,2,211,14,223,1,210,15,224,8,209,16,217,7,216,9,218,6,215,10,219,5,214,11,220,4,213,12,221,20,237,52,205,21,236,53,204,28,229,60,197,29,228,61,196,30,227,62,195,27,230,59,198,22,235,54,203,19,238,51,206,18,239,50,207,23,234,55,202,26,231,58,199,31,226,63,194,32,225,64,193,25,232,57,200,24,233,56,201,17,240,49,208]';


%For dual-mode device *************************
Trans.Connector = [177,80,145,112,184,73,152,105,185,72,153,104,192,65,160,97,191,66,159,98,186,71,154,103,183,74,151,106,178,79,146,111,179,78,147,110,182,75,150,107,187,70,155,102,190,67,158,99,189,68,157,100,188,69,156,101,181,76,149,108,180,77,148,109,140,93,132,85,139,92,133,86,138,91,134,87,137,90,135,88,144,89,136,81,143,96,129,82,142,95,130,83,141,94,131,84,175,122,167,114,174,123,166,115,171,126,163,118,170,127,162,119,169,128,161,120,172,125,164,117,173,124,165,116,176,121,168,113,241,40,249,48,244,37,252,45,245,36,253,44,248,33,256,41,247,34,255,42,246,35,254,43,243,38,251,46,242,39,250,47,212,3,222,13,211,2,223,14,210,1,224,15,209,8,217,16,216,7,218,9,215,6,219,10,214,5,220,11,213,4,221,12,237,20,205,52,236,21,204,53,229,28,197,60,228,29,196,61,227,30,195,62,230,27,198,59,235,22,203,54,238,19,206,51,239,18,207,50,234,23,202,55,231,26,199,58,226,31,194,63,225,32,193,64,232,25,200,57,233,24,201,56,240,17,208,49]';

Trans.ConnectorES = Trans.Connector;
%not working elements for (center edge) PPC17 array
%rem_PPS1_notWorking = [19	22	27	40	44	46	49	63	67	71	73	75	79	83	87	97	98	99	102	174	178	182	186	190	235];
%rem_PPC17_notWorking = [7 17 29 43 45 53 54 66 68 70 72 74 76 78 80 82 87 94 110 114 118  122  124 126 127 129 144  153 157 161 241 250 256];
rem_PPC17_notWorking = [16,39,53,65,66,67,69,73,78,79,80,102,105,113,118,119,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,146,148,149,156,157,169,170,171,172,173,174,175,176,178,179,180,181,182,183,186,188,189,191,196,201,204,208,209,228,229,244,245,246]';

elOffEnable_TX = 0; % 0 means they will transmit, 1 means they will not transmit
elOffEnable_RX = 0; % 0 means they will receive, 1 means they will not receive

elOn = 1:256;
tx_elOff = [rem_PPC17_notWorking];%;%rem_PPC17_notWorking]; 
rx_elOff = [rem_PPC17_notWorking];%, rem_PPC17_notWorking];
notOffElements = setdiff(elOn,rx_elOff);

Trans.HVMux = struct('highVoltageRails', 100, ...
                        'logicRail', 10.5, ...
                        'clock', 5, ...
                        'clockInvert', 0, ...
                        'polarity', 0, ... 
                        'latchInvert', 0, ...
                        'Aperture', Trans.Connector, ...
                        'VDASAperture', uint8([2; 0; 0]));

P.lambda = Resource.Parameters.speedOfSound / (Trans.frequency*1e6);
P.startDepthWvl = P.startDepth/P.lambda;   % Acquisition depth in wavelengths
P.endDepthWvl = P.endDepth/P.lambda;   % This should preferrably be a multiple of 128 samples.

% Specify PData structure array.
PData(1).PDelta = [Trans.spacing, 0, 0.5];
PData(1).Size(1) = ceil((P.endDepthWvl-P.startDepthWvl)/PData(1).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
PData(1).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepthWvl]; % x,y,z of upper lft crnr.

if P.numRays > 1
    P.tx_positions = linspace(PData(1).Origin(1),-PData(1).Origin(1),P.numRays);
else
    P.tx_positions = 0;
end



for i = 1:P.numRays
    PData(1).Region(i).Shape = struct('Name','Rectangle',...
                                      'Position', [P.tx_positions(i), 0, P.startDepthWvl],...
                                      'width',1, ... % 1
                                      'height',PData(1).Size(1));
end


% No PData.Region specified, so a default Region for the entire PData array will be created by computeRegions.

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';

%%%
maxAcqLength = ceil(sqrt(P.endDepthWvl^2 + ((Trans.numelements-1)*Trans.spacing)^2));  
P.acqsamp = 4*round(maxAcqLength-P.startDepthWvl)*2;
P.acqsamp = P.acqsamp+128-mod(P.acqsamp,128);
Resource.RcvBuffer(1).rowsPerFrame = P.numRays*P.acqsamp+128; % total buffer length must be greater than num samples needed . . . see documentation
%%%

% Resource.RcvBuffer(1).rowsPerFrame = 4096 * P.numRays;   % this size allows for maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 10;       % 100 frames used for RF cineloop.
Resource.InterBuffer(1).numFrames = 1;  % one intermediate buffer needed.
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = 'L11-5 Focused';
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = 20;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);

% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [usrTXParam.freq,usrTXParam.percentOn,2*usrTXParam.numCyc,usrTXParam.polarityFirstHalf];


% num_active = round(Trans.numelements*P.apod);
% num_off = Trans.numelements-num_active;
% apod_vec = [zeros(1,floor(num_off/2)) ones(1,num_active) zeros(1,ceil(num_off/2))];
% disp(length(apod_vec));

TX = repmat(struct('waveform',1,...
                   'Origin',[0,0,0],...
                   'focus',P.focus/P.lambda,... % make me 0!
                   'Steer',[0,0,0],...
                   'aperture', 1, ...
                   'Apod',ones(1,Trans.numelements)),1,P.numRays);
% %                    'Apod',apod_vec),1,P.numRays);

D = P.focus / P.fnum; % aperture width (m)
for i = 1:P.numRays
    TX(i).Origin(1) = P.tx_positions(i);    
    TX(i).Delay = computeTXDelays(TX(i));
    
    tmp = TX(i).Origin(1)*P.lambda; % center of tx in m
    window = [tmp-(D/2), tmp+(D/2)];
    apod_idx = (Trans.ElementPos(:,1)*P.lambda) >= window(1);
    apod_idx = apod_idx & ((Trans.ElementPos(:,1)*P.lambda)<=window(2));
    
    %TX(i).Apod = double(apod_idx');
    %TX(i).Apod = 1';
% %     startN = Trans.numelements/2+ceil(P.tx_positions(i) - Resource.Parameters.numTransmit/2);
% %     endN = Trans.numelements/2+floor(P.tx_positions(i) + Resource.Parameters.numTransmit/2);
% %     
% %     TX(i).Apod(startN:endN) = hamming(length(startN:endN))'
    if elOffEnable_TX
        TX(i).Apod(tx_elOff) = 0;
        TX(i).Delay(tx_elOff) = 0;
    end
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [450,590,700,800,875,925,950,975]; % --type your TGC values here :)
TGC.rangeMax = P.endDepthWvl;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays -
%   endDepth - add additional acquisition depth to account for some channels
%              having longer path lengths.

% maxAcqLength -- distance between edge element and opposite corner of the pixel data
maxAcqLength = ceil(sqrt(P.endDepthWvl^2 + ((Trans.numelements-1)*Trans.spacing)^2));  

LPF = [+0.01100 +0.00720 -0.00420 -0.02020 -0.03340 -0.03370...
 -0.01280 +0.03110 +0.09140 +0.15380 +0.20080 +0.21830];
 
BPF = [+0.00064 -0.00095 -0.00247 +0.00259 +0.00421 -0.00238 -0.00037 ...
 -0.00253 -0.01376 +0.00998 +0.03030 -0.01151 -0.02365 +0.00098 ...
 -0.03329 +0.01553 +0.13797 -0.02243 -0.24423 +0.01065 +0.28943];

Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P.startDepthWvl, ...
                        'endDepth', maxAcqLength,...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'aperture', 1, ...,
                        'LowPassCoef',LPF,...
                        'InputFilter',BPF,...
                        'callMediaFunc', 1),1,Resource.RcvBuffer(1).numFrames*P.numRays);                                                                                  
                    
% - Set event specific Receive attributes.
for i = 1:length(Receive)
    % Setting frame number
    Receive(i).framenum = floor((i-1)/P.numRays)+1;
    % Setting acquisition number
    Receive(i).acqNum = mod(i-1,P.numRays)+1;
    %Receive(i).Apod = TX(Receive(i).acqNum).Apod;
        if elOffEnable_RX
           Receive(i).Apod(rx_elOff) = 0;
        end
end

% Specify Recon structure arrays.
Recon = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame', -1, ...     % use most recently transferred frame
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...  % auto-increment ImageBuffer each recon
               'RINums', 1:P.numRays);

% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'replaceIntensity', ...
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 1),1,Resource.RcvBuffer.numFrames * P.numRays);

for i = 1:length(ReconInfo)
    ReconInfo(i).txnum = mod(i-1,P.numRays)+1;
    ReconInfo(i).rcvnum = i;
    ReconInfo(i).regionnum = mod(i-1,P.numRays)+1;
end
               
% Specify Process structure array.
pers = 20;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',1.0,...            % pgain is image processing gain
                         'reject',2,...      % reject level
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',1};

% Specify SeqControl structure arrays.
nsc = 1;
SeqControl(nsc).command = 'jump'; % jump back to start.
SeqControl(nsc).argument = 1; %in microseconds
JUMP = nsc;
nsc = nsc + 1;
SeqControl(nsc).command = 'timeToNextAcq';  % time between frames
SeqControl(nsc).argument = (1/P.prf)*1e6;  % microseconds
TTNA = nsc;
nsc = nsc + 1;
SeqControl(nsc).command = 'returnToMatlab';
RTM = nsc;
nsc = nsc + 1; % nsc is count of SeqControl objects
SeqControl(nsc).command = 'triggerOut';
TRIG_OUT = nsc;
nsc = nsc + 1;

n = 1; % n is count of Events

% Acquire all frames defined in RcvBuffer
for i = 1:Resource.RcvBuffer(1).numFrames
    
    % Transmit/receive all lines for a single frame
    for j = 1:P.numRays   
        Event(n).info = 'Full aperture.';
        Event(n).tx = j;
        Event(n).rcv = (i-1)*P.numRays + j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = [TRIG_OUT,TTNA];           
        n = n+1;
    end

% % %For Hydrophone measurements
% % 	Event(n).info = 'Aquisition.';
% % 	Event(n).tx = 1;         % use 1st TX structure.
% % 	Event(n).rcv = i;      % use ith Rcv structure of frame.
% % 	Event(n).recon = 0;      % no reconstruction.
% % 	Event(n).process = 0;    % no processing
% % 	Event(n).seqControl = [TRIG_OUT,TTNA]; % set wait time and transfer data
% % % % 	   SeqControl(nsc).command = 'transferToHost';
% % % % 	   nsc = nsc + 1;
% %       n = n+1;
    
    % Add a transfer to the last Event of the frame. 
    Event(n-1).seqControl = [TRIG_OUT,TTNA,nsc];
       SeqControl(nsc).command = 'transferToHost';
       nsc = nsc + 1;
    
    % Recon the frame we just acquired
    Event(n).info = 'Reconstruct';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 1;
    if floor(i/5) == i/5     % Exit to Matlab every 5th frame
        Event(n).seqControl = RTM;
    else
        Event(n).seqControl = 0;
    end
    n = n+1;
end

Event(n).info = 'Jump back to first event';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = JUMP;

% % User specified UI Control Elements
% % - Sensitivity Cutoff
% UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
%                   'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
%                   'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
% UI(1).Callback = text2cell('%SensCutoffCallback');
% 
% % - Range Change
% MinMaxVal = [64,300,P.endDepth]; % default unit is wavelength
% AxesUnit = 'wls';
% if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
%     if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
%         AxesUnit = 'mm';
%         MinMaxVal = MinMaxVal * (Resource.Parameters.speedOfSound/1000/Trans.frequency);
%     end
% end
% UI(2).Control = {'UserA1','Style','VsSlider','Label',['Range (',AxesUnit,')'],...
%                  'SliderMinMaxVal',MinMaxVal,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f'};
% UI(2).Callback = text2cell('%RangeChangeCallback');
% 
% % Specify factor for converting sequenceRate to frameRate.
% frameRateFactor = 5;

% Save all the structures to a .mat file.
% save('MatFiles/L11-5Flash');

filename = 'MatFiles/cmut_bmode_calibration.mat';
save(filename);     
return

% % **** Callback routines to be converted by text2cell function. ****
% %SensCutoffCallback - Sensitivity cutoff change
% ReconL = evalin('base', 'Recon');
% for i = 1:size(ReconL,2)
%     ReconL(i).senscutoff = UIValue;
% end
% assignin('base','Recon',ReconL);
% Control = evalin('base','Control');
% Control.Command = 'update&Run';
% Control.Parameters = {'Recon'};
% assignin('base','Control', Control);
% return
% %SensCutoffCallback
% 
% %RangeChangeCallback - Range change
% simMode = evalin('base','Resource.Parameters.simulateMode');
% % No range change if in simulate mode 2.
% if simMode == 2
%     set(hObject,'Value',evalin('base','P.endDepth'));
%     return
% end
% Trans = evalin('base','Trans');
% Resource = evalin('base','Resource');
% scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
% 
% P = evalin('base','P');
% P.endDepth = UIValue;
% if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
%     if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
%         P.endDepth = UIValue*scaleToWvl;
%     end
% end
% assignin('base','P',P);
% 
% evalin('base','PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));');
% evalin('base','PData(1).Region = computeRegions(PData(1));');
% evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');
% Receive = evalin('base', 'Receive');
% maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
% for i = 1:size(Receive,2)
%     Receive(i).endDepth = maxAcqLength;
% end
% assignin('base','Receive',Receive);
% evalin('base','TGC.rangeMax = P.endDepth;');
% evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
% Control = evalin('base','Control');
% Control.Command = 'update&Run';
% Control.Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','Receive','TGC','Recon'};
% assignin('base','Control', Control);
% assignin('base', 'action', 'displayChange');
% return
% %RangeChangeCallback

