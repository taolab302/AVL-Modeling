function [I, I_start, I_end] = set_constant_current_sequence(time_length, I_val, I_blank)

if nargin < 2
    I_val = 10;
    I_blank = 100000;
elseif nargin < 3
    I_blank = 100000;
end

I_start = 1 + I_blank;
I_end = time_length-I_blank;

I = zeros(1, time_length);
I(I_start:I_end) = I_val;
