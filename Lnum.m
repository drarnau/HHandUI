function [ latex ] = Lnum( num )
%Lnum transforms a decimal number into a Latex friednly format: 1.23 is
%1&23
%   num numerical
%   latex string

% Get numbers
x = sscanf(strrep(num2str(num,8),'.',''),'%1d');

if length(x) == 3
    if num < 10
        latex = [num2str(x(1)) '&' num2str(x(2)) num2str(x(3))];
    else
         latex = [num2str(x(1)) num2str(x(2))  '&' num2str(x(3)) '0'];
    end
elseif length(x) == 2
    if num < 10
        latex = [num2str(x(1)) '&' num2str(x(2)) '0'];
    else
         latex = [num2str(x(1)) num2str(x(2))  '&00'];
    end
elseif length(x) == 1
    latex = [num2str(x(1)) '&00'];    
else
   latex = [num2str(x(1)) num2str(x(2))  '&' num2str(x(3)) num2str(x(4))]; 
end
    
% Put all togheter
    
    


end

