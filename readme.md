The AVL model was developed using MATLAB R2018a (MathWorks). You can run the script `main.m` to plot the modeling figures in the paper. 

***

`PlotWholeCellVoltage.m` runs a current clamp simulation of the AVL model, and analyzes the temporal evolution of the normalized conductance of single channel. It uses `ode15s` function to run the simulation. The current clamp stimulation is generated by  the function `set_constant_current_sequence.m`. Each current-clamp step is set in 2 pA increment within a range of -4 pA ~ 12 pA. 

Input:

```matlab
strains: "wt", "exp-2(lf)", "exp-2(gf)", "shl-1" or "exp-2;shl-1"

time: length of the simulation (in milliseconds), larger than 2000.
```

`PlotWholeCellCurrent.m` runs voltage clamp simulations of the wild-type, *exp-2(lf)* and *exp-2(gf)* AVL model. It uses `ode15s` function to run the simulation. The voltage clamp stimulation is generated by the function`set_constant_voltage_sequence.m`. The holding potential of the voltage-clamp steps is -60 mV and each step is set in 10 mV increment within a range of -120 mV ~ 70 mV. 

Input:

	time: time length of the simulation (in milliseconds), larger than 1500.

`PlotTailPeakIV.m` calculates the normalized I-V relationships of the repolarization-activated currents for the wild-type, exp-2(lf) and exp-2(gf) AVL model, based on the voltage clamp simulation results. The I-V curve is fitted with Boltzmann function.

Input:

```
time: time length of the simulation (in milliseconds), same to the time used in PlotWholeCellCurrent.m.
```

***

`AVLModel.m` includes all of the ordinary differential equations of AVL model.

`WholeCellCurrent.m` includes all of the ordinary differential equations of the channel models in AVL model.

`AVLParameters.m` sets the parameters and initial value of the wild-type and mutant AVL model. You can adjust the parameters and the initial value here.

`set_constant_current_sequence.m` generates a current array that start at 1 s and end at the last 1 s.

`set_constant_voltage_sequence.m` generates a voltage array that start at 0.5 s and end at the last 1 s.

`RecordNowParams.m` records the parameters in a simulation. The parameters are saved as a CSV file.

`nowtime.m` generates a unique filename based on the current time.

