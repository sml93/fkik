## This is a repository for fkik for BLASTER UAV ##  
05 Jun 22  
The fkik algorithm will aid to allow BLASTER UAV to optimise the (position) placement of each of its relative joints, i.e. from the origin (inertial) frame to the POC frame in 3-dimension (3D).  
  
The origin fkik3D_single.py is written to optimise per step iteration.  
  
To work on refreshing fkik3D_single.py to allow for better fkik as a framework instead which will allow the algorithm to be called in realtime and repetitively.  

12 Jun 22  
The recent commit is corrected for optimizing over a range of deb_x and deb_y continuously.  
  
To insert the function for IK and FKIK in subsequently.  


18 Jun 22  
Adding an edition for hkf:  
 - hybrid_kinematics_force.py  
 - main_hkf.py  

22 Jun 22:  
Function for IK in hybrid_kinematics_force.py is edited for the inclusion for standoff distance, d_so. (line 177)  
However, to look into the expression. (line 209)  
Figure out an equation that relates angles and distances to forces.  

Fixed:  
updated objective function to optimize altitude and the trajectory of the water.
 - as this affects the force at the POC in a whole.  

updated constraints and bounds to:  
 - nozzle angle, x/y/z-water with their respective bounds.  


TODO:  
Look into pub/sub position to ros so that this will enable easy sharing of data, and potentially as position reference.  
