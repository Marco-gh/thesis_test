function [A_bar,B_bar,C_bar,T] = detectabilty_decomposition(O,A,B,C)
%DETECTABILTY_DECOMPOSITION Summary of this function goes here
%   Detailed explanation goes here
rO = rank(O);
n = size(A,1);

if rO == n
    A_bar = A;
    B_bar = B;
    C_bar = C;

    T = eye(n);
    return
end

U = null(O);
D = null(U');
T = [D,U];

A_bar = T'*A*T;
B_bar = T'*B;
C_bar = C*T;
end