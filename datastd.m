function [Ystd] = datastd(Y)

% Data Standardization


[T, ~] = size(Y);
Ystd = Y - repmat(mean(Y), T, 1);
Ystd = Ystd./repmat(std(Y,0,1), T, 1);

end

