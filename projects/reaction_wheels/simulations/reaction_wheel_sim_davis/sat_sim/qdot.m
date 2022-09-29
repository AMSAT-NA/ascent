function q = qdot( q, omega )
%  QDOT Differential of a quaternian in a rotating frame.

% time evolution equations of the Euler parameters
%	see Fukushima_2008_AJ_135_2298.pdf

q = [ -q(2) -q(3) -q(4);
       q(1) -q(4)  q(3);
       q(4) +q(1) -q(2);
      -q(3) +q(2) +q(1) ] * omega / 2;
q = q';
end
