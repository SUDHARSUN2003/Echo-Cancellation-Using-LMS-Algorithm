function [A,E] = lms(x,d,mu,nord,a0)

%LMS	Adaptive filtering using the Widrow-Hoff LMS algorithm.

%---

%USAGE [A,E] = lms(x,d,mu,nord,a0)

%

%           x    : input data to the adaptive filter.

%           d    : desired output

%           mu   : adaptive filtering update (step-size) parameter

%           nord : number of filter coefficients

%           a0   : (optional) initial guess for FIR filter 

%		   coefficients - a row vector.  If a0 is omitted

%		   then a0=0 is assumed.

%

%     The output matrix A contains filter coefficients.

%        - The n'th row contains the filter coefficients at time n

%        - The m'th column contains the m'th filter coeff vs. time.

%        - The output vector E contains the error sequence versus time.

%

%  see also NLMS and RLS

%

%---------------------------------------------------------------

% copyright 1996, by M.H. Hayes.  For use with the book 

% "Statistical Digital Signal Processing and Modeling"

% (John Wiley & Sons, 1996).

%---------------------------------------------------------------



X=convm(x,nord);

[M,N] = size(X);

if nargin < 5,   a0 = zeros(1,N);   end

a0   = a0(:).';

E(1) = d(1) - a0*X(1,:).'; 

A(1,:) = a0 + mu*E(1)*conj(X(1,:));

if M>1

for k=2:M;

    E(k) = d(k) - A(k-1,:)*X(k,:).';

    A(k,:) = A(k-1,:) + mu*E(k)*conj(X(k,:));

    end;

end;



