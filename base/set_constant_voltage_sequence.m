function V = set_constant_voltage_sequence(V_max, time_length, V_holding)

V = zeros(1, time_length) + V_holding;
V(5001:end-10000) = V_max;
end
