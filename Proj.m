function [y, a]=Proj(x, a, npred, pos)
% LPREDICT estimates the values of a data set before/after the observed set.
%
% LPREDICT uses linear prediction to extrapolate data, typically a
% timeseries. Note that this is not the same as linear extrapolation. 
% A window of autocorrelation coefficients is moved beyond the data
% limits to extrapolate the data. For a discussion, see Press et. al. [1].
%
% The required coefficients are derived from a call to LPC in MATLAB's
% Signal Processing Toolbox
%
% Example:
% y=LPREDICT(x, np, npred, pos)
% [y, a]=LPREDICT(x, np, npred, pos)
%      x:       the input data series as a column vector or a matrix 
%                   with series organized in columns
%      np:      the number of predictor coefficients to use (>=2)
%      npred:   the number of data values to return in the output
%      pos:     a string 'pre' or 'post' (default: post)
%                   This determines whether extrapolation occurs before or
%                   after the observed series x.
%
%      y:       the output, appropriately sequenced for concatenation with
%                   input x
%      a:       the coefficients returned by LPC (organized in rows).
%                   These can be used to check the quality/stability of the
%                   fit to the observed data as described in the
%                   documenation to the LPC function.
%
% The output y is given by:
%       y(k) = -a(2)*y(k-1) - a(3)*y(k-2) - ... - a(np)*y(k-np)
%                where y(n) => x(end-n) for n<=0
% 
% Note that sum(-a(2:end))is always less than unity. The output will
% therefore approach zero as npred increases. This may be a problem if x
% has a large DC offset. Subtract the the column mean(s) of x from x on
% input and add them to the output column(s) to restore DC. For a more
% accurate DC correction, see [1].
%
% To pre-pend data, the input sequence is reversed and the output is
% similarly reversed before being returned. The output may always be
% vertically concatenated with the input to extend the timeseries e.g:
%       k=(1:100)';
%       x=exp(-k/100).*sin(k/5);
%       x=[lpredict(x, 5, 100, 'pre'); x; lpredict(x, 5, 100, 'post')];
% 
% 
% See also LPC
%
% References:
% [1] Press et al. (1992) Numerical Recipes in C. (Ed. 2, Section 13.6).
%
% Toolboxes Required: Signal Processing
%
% Revisions:    10.07 renamed to avoid filename clash with System ID
%                     Toolbox
%                     DC correction help text corrected.
%
% -------------------------------------------------------------------------
% Author: Malcolm Lidierth 10/07
% Copyright � The Author & King's College London 2007
% -------------------------------------------------------------------------


% ---------------Argument checks---------------
if nargin<3
    error('Not enough input arguments');
end

if nargin<4
    pos='post';
end
cols=size(x,2);
if cols>1
    y=zeros(npred, cols);
    
    for k=1:size(x,2)
        [y(:,k) ]=Proj(x(:,k), a, npred, pos);
    end
    return

if nargin==4 && strcmpi(pos,'pre')
    x=x(end:-1:1);
end


cc=-a(2:end);

y=zeros(npred,1);
y(1)=cc*x(end:-1:end-np+1);
for k=2:min(np,npred)
    y(k)=cc*[y(k-1:-1:1); x(end:-1:end-np+k)];
end
for k=np+1:npred
    y(k)=cc*y(k-1:-1:k-np);
end
if nargin==4 && strcmpi(pos,'pre')
    y=y(end:-1:1);
end

return
end
