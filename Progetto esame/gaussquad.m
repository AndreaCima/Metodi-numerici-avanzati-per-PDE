function [x, w] = gaussquad(q)
% nodi e pesi di quadratura di Gauss su [0 1]
B = ( 1:(q-1) )./ sqrt( 4*( 1:(q-1) ).^2 -1 );
[V, D] = eig( diag(B, -1) + diag(B, 1) );
x = ( diag(D)+1 )/2;
w = ( V(1, :).*V(1, :) )';

end