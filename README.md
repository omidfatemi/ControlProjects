# Rectilinear motion in a space	
in this code, the aim is to design a kalman filter that estimate a motion in 2d space, like a penalty kick, during a certain time limit. 
Main assumptions are :
1. Time limit is 0.2 seconds which is the time required for the goalkeeper to detect ball approximate position when it reaches the goal
2. Initial velocity of the ball is within two fixed values, [20,40]m/s.
3. Penalty kick is always inside the goal
