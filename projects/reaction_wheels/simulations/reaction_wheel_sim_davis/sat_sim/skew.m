function s = skew(x)
%SKEW Skew matrix of a vector.

	s = [ 0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];
end



% Unused miscellanious fuctions follow:

% approximation of tanh(S)
%
function y = fastTanh(x)
	y = x - 1/3.*x.^3 + 2/15.*x.^5;
end


% FAST_ACOS fast approximation of y = acos(x)
%
function y = fast_acos(x)
    a = sqrt(2+2*x);
    b = sqrt(2-2*x);
    c = sqrt(2-a);
    y = (8*c-b)/3;
end

