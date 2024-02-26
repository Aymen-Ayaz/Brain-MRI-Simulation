function S = steady_state_signal(M0,T1,T2,TR,TE,alpha,TYPE,TI,T2p)
%
%S = steady_state_signal(M0,T1,T2,TR,alpha,TYPE)
%
%Compute MR steady state signal amplitude
%
%Input
%       T1,T2 : relaxation times [ms]
%       M0    : equilib. longitudinal magnetization
%       TR    : repetition time [ms]
%       FA    : flip angle [degree]
%       TYPE  : string, describing the signal type, one of
%               'Ernst' : Ernst Equation
%               'FID'   : The FID or S+ signal in unspoiled steady state seqs.
%               'ECHO'  : The ECHO or S- signal in unspoiled steady state seqs.


c  = cos(alpha);
s  = sin(alpha);
t  = tan(alpha/2);
E1 = exp(-TR./T1);
E2 = exp(-TR./T2);
EE = exp(-TE./T2);
E3 = 2.*exp(-TI./T1);


switch lower(TYPE)
    case 'ernst'
        S = M0.*s.*(1-E1)./(1-E1.*c);
    case 'fid'
        S = M0.*t.*(1-r.*(E1-c));
    case 'echo'
        S = M0.*t.*(1-r.*(1-E1.*c));
    case 'spinecho'
        S = M0.*(1-E1).*EE;
    case 'ir'
        S = M0.*(1-E3+E1).*EE;
    otherwise
        error('unkown signal type')
        
end