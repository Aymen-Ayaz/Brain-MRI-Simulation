function x = fftc(x,n,dim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% One-dimensional inverse Fourier transform
%
% Christopher W. Roy 2018-12-04
% fetal.xcmr@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fctr = size(x,1)*size(x,2);
fctr=1;
x = (1/sqrt(fctr))*fftshift(fft(ifftshift(x),n,dim));