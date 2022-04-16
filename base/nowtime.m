function now = nowtime()
time = num2str(clock);
s = regexp(time, '\s+', 'split');
s = char(s);
for i = 2:6 
    if s(i, 2) == ' ' || s(i, 2) == '.'
        s(i, 2) = s(i, 1);
        s(i, 1) =  '0';
    end
end
now = char(join(string(s)));
now(isspace(now)) = [];
now = now(3:14);
end

