To run the simulator, double click on simulator.exe
The Matlab Runtim may be required to run the application. If it is not already installed, it should automatically download from the internet.

This default configuration uses 500 particles and inter-particle collisions and effusion are dissabled initially.
Enable collisions and effusion with the corresponding check boxes.
When effusion is enabled, the simulation may hang for a short time whenever a prticle escapes the box. This is normal, as the program must re-calculate a list of combinations of the particles (500 choose 2 combinations)
Use the slider at the bottom to control the simulation speed (set to maximum by default)
Use the buttons on the top toolbar to orbit, pan, and zoom the camera on the 3D view.

A histogram of the current particle velocities is shown. The particles are creaded with a random normal distribution, which will only change when particle collisions are enabled. 
There are also plots of the number of particles and average pressure over time. These will not change unless particle collisions are enabled.

Because the particles are visualized as markers and not 3D spheres, their marker size does not neccisarily correlate to their actual size used in the simulation. (These markers are proportional to the camera zoom, not distance from the camera.)