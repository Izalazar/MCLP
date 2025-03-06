function [a_vt,h_e] = MCLP(s,p,ep1,ep2,alpha)
% Description
%   This function computes vocal tract filter (VTF) coefficients to the
%   voice inverse filtering using the Maximum correntropy criterion-based 
%   linear prediction method (MCLP).

% Inputs
%   s           : Voice signal frame [samples]
%   p           : Order of VTF
%   ep1, ep2    : Algorithm thresholds
%   alpha       : Pre-emphasis coefficient
%
% Outputs
%   a_vt        : VTF coefficients
%   h_e         : Resulting correntropy weighting function

% References
% [1] Zalazar, I. A., Alzamendi, G. A., Zañartu, M., & Schlotthauer, G. (2024).
%     Maximum Correntropy Linear Prediction for Voice Inverse Filtering: 
%     Theoretical Framework and Practical Implementation. IEEE/ACM Transactions 
%     on Audio, Speech, and Language Processing.

% Author
% Iván Ariel Zalazar izalazar@ingenieria.uner.edu.ar. 

s_v1=filter([1 -alpha],1,s);%pre-emphasis
s_v1=s_v1/max(abs(s_v1));%normalization

d=[s_v1(2:end);zeros(p,1)];

N=length(s_v1);

X=zeros(N+p-1,p);

for k=1:1:p
    X(k:N+k-1,k)=s_v1;%Data matrix
end

w=lpc(s_v1,p);%inicialization by LPC
w=w(2:end)';
e=d-X*w;% LP error
sigma=1.06*min([std(e),(quantile(e,0.75)-quantile(e,0.25))/1.34])*length(e)^(-1/5);% kernel size inicialization by silverman's rule.

while true
    h_e=exp(-(e.^2)/(2*sigma^2));
    H=diag(h_e);
    R=(X')*H*X;%weighted autocorrelation matrix
    r=(X')*H*d;%weighted autocorrelation vector

    w_old=w;%old coefficents
    w=R\r;%compute new coefficents
    e=d-X*w;% update LP error

    if norm(w-w_old,2)<ep1
        sigma=1.06*min([std(e),(quantile(e,0.75)-quantile(e,0.25))/1.34])*length(e)^(-1/5);% update kernel size by Silverman's rule.
    end

    if norm(w-w_old,2)<ep2
        break;%stop do-while
    end
end
h_e=exp(-(e.^2)/(2*sigma^2));
a_vt=[1;-w]';%vocal tract filter

end