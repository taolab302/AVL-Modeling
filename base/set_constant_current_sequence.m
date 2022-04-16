function [I, I_start, I_end] = set_constant_current_sequence(time_length, I_val)

if nargin < 2
    I_val = 10;
end
I_start = 100001;
I_end = time_length-100000;

I = zeros(1, time_length);
I(I_start:I_end) = I_val;