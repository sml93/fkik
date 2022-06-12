## This is a repository for fkik for BLASTER UAV ##  
05 Jun 22  
The fkik algorithm will aid to allow BLASTER UAV to optimise the (position) placement of each of its relative joints, i.e. from the origin (inertial) frame to the POC frame in 3-dimension (3D).  
  
The origin fkik3D_single.py is written to optimise per step iteration.  
  
To work on refreshing fkik3D_single.py to allow for better fkik as a framework instead which will allow the algorithm to be called in realtime and repetitively.  

12 Jun 22  
The recent commit is corrected for optimizing over a range of deb_x and deb_y continuously.  
  
To insert the function for IK and FKIK in subsequently.  
  
