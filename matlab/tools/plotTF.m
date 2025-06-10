function varargout = plotTF(t,x,NFFT,sp, scaling, varargin)
% [f, psdX] = plotTF(t,x,NFFT)   - NFFT is optional

if nargin<5; scaling=1; end
if nargin<4 || isempty(sp);   sp = [subplot(2,1,1) subplot(2,1,2)]; end
if nargin<3 || isempty(NFFT); NFFT = 512; end
if size(x,2)<size(x,1), x=x'; end

dt = t(2)-t(1);

[XX, f] = VincePSD(x,NFFT,dt, scaling);

psdX = zeros(size(XX,1),size(XX,2));
for i=1:size(x,1)
  psdX(:,i) = XX(:,i,i); psdX(1,i) = nan;
end
    
axes(sp(1))
if isempty(varargin)
    plot(t,x), hold on, xlabel ('Time [s]')
else
    plot(t,x, varargin), hold on, xlabel ('Time [s]')
end

axes(sp(2))
semilogy(f,psdX); hold on, xlabel ('Frequency [Hz]')

if nargout
    varargout{1} = f;
    varargout{2} = psdX;
end

function [PSD, FREQ, OMEGA] = VincePSD (x, NFFT, DT, scaling)
% [PSD FREQ] = VincePSD (x, NFFT, DT)
%
% This routine computes the psd matrix of signals x
% signals in x are divided in smaller segments with length NFFT
% These segments can't overlap
% For each segment, a periodogram is estimated, using a Hanning window
% The psd are estiamted as the mean periodogram and returend in PSD
% An estimation of the corresponding frequencies can be obtained if the
% time step DT is given
% NOTE : x must be arranged as x(signal number, sampling) in such a way
% that the time is running along the columns.
%
% SEE ALSO: VincePSD_ (PSD is such that integral over f returns the
% variance).
%
NS= size(x,1);  % Number of signals
N = size(x,2);  % Number of time steps in x

if NS>N; x=x';NS= size(x,1);N = size(x,2);end

%trim nan's at beginning of signal and nan's at end of signal
for i=1:NS
    x(i,isnan(x(i,:)))= nanmean(x(i,:));
end


if nargin==3, scaling=1; end



Nblocs = floor(N/NFFT);

PSD=zeros(NFFT/2,NS,NS);

% Build Hanning window
t = (0:NFFT-1)/(NFFT-1);
hann = sin(pi*t).^2;
hann = hanning(NFFT)';
hann = rectwin(NFFT)';
W = sum(hann.^2);

% Loop along the blocs
for i=1:Nblocs
    for s1 = 1:NS
        % Extract segment for signal 1 and window it
        xx1 = x(s1, (i-1)*NFFT+1:i*NFFT);
        xx1 = xx1.* hann;
        XX1 = fft(xx1 - mean(xx1));
        for s2 = 1:NS
            % Extract segment for signal 2 and window it
            xx2 = x(s2, (i-1)*NFFT+1:i*NFFT);
            xx2 = xx2.* hann;
            XX2 = fft(xx2 - mean(xx2));
            periodogram = XX1.*conj(XX2);

            PSD(:,s1,s2) = PSD(:,s1,s2) + ( periodogram(1:NFFT/2) )';
        end
    end
end
PSD = PSD / Nblocs / W;

if nargin >= 3
    DF = 1 / (NFFT*DT);
    FREQ = (0:NFFT/2-1)*DF;
    if scaling == 1
        PSD = PSD / NFFT / DF * 2; % scale PSD such that intergral of PSD over FREQ returns variance
    end
    if scaling == 2
        PSD = PSD / NFFT / DF * 2 / 4/pi; % scale PSD such that twice intergral of PSD over OMEGA (OMEGA=0,..inf) returns variance
    end
else
    FREQ = 0;
end
OMEGA = 2*pi*FREQ;