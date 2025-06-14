function [H] = FCTTransfertMDDL (M,K,C,f)
% [H] = FCTTransfertMDDL (M,K,C,f)
% 
% Calcule la fonction de Transfert d'un système M DDL
%
% see ALSO: FCTTransfertMDDL

Nfr = length(f); N = size(M,1);
H=zeros(N,N,Nfr);
for i=1:Nfr
    omega = 2*pi*f(i);    
    H(1:N,1:N,i)= inv(-M*omega^2 + sqrt(-1)*omega*C + K);
end