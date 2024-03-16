import numpy as np
import matplotlib.pyplot as plt
import math 

#This code will find the final temperature distribution of a 2D square plate undergoing unsteady heat diffusion 
#the boundary conditions are such that the bottom row of the plate is held steady at 1 C where C = Celsius 
#The leftmost and rightmost column of the plate are both held at 0 C 
#and the top row of the plate is similarly held at 0 C 

#step 1
dom_size = 1; #defines both the length and width of the 2D square plate (i.e. domain) as being 1 meter



#step 2
# here the mesh of a 2d square plate (i.e. domain) will be defined
npoints = 20; #defines both the number of columns and rows for a 2d square plate
h = dom_size/(npoints-1); #defines the length of the edge of each cell in the mesh
#the value of h corresponds to delta x and delta y of the 2D square plate (i.e. domain)

dt = 0.0001; #defines Delta t measured in seconds (s)
#note: the finer the mesh the smaller the value of dt should be 
#functionally this means how much time has elapsed following the completion of a single iteration 
# in this can iteration is a single cycle of the for loop that containing the governing equation for an unsteady heat diffusion equation 
alpha= dt/(h*h) #defines the term (Gamma Delta t/ h^2) where Gamma (or the diffusion coefficient) is equal to 1 



#steps 3 and 4 
y = np.zeros((npoints,npoints)) #sets up a two 2d array that contain all zeros 
#this also serves as a way to intialize the interior points of the domain and set the boundary condition for the top row, leftmost column, and the rightmost column
y_new = np.zeros((npoints,npoints)) #sets up a two 2d array that contain all zeros 

y[0,:] = 1 #sets the boundary condition of the bottom row for y at a temperature of 1 C
y_new[0,:] = 1 #sets the boundary condition of the bottom row for y_new at a temperature of 1 C

y_transient = [],[],[]; #sets a variable that stores three different arrays 
#the purpose of this variable is to record the value of y at every timestep (i.e. after every iteration)



#steps 5 and 6 
error_mag = 1; #initialize the variable error_mag (just to make sure its value is greater than error_req)
error_req = 1e-6; # establishes the threshold value the error must fall below to break the following loop  
iterations = 0; # initializes the variable iterations as having a value of 0 
#the variable iterations is there to determine the number of cycles the code will have to do to reach a stable solution
error_track = []; #tracks how the error changes with time 
# [] indicate the error_track is an array of flexible size (right now it has a size of 0)

while (error_mag > error_req): # keeps doing the calculations until the error is less than the threshold value 
    
    # the range ensures that the for loop focuses on the points interior to the domain and does not modify the boundary conditions
    for index1 in range(1,npoints-1): #this for loop acts on points interior to the domain along the x-axis
        for index2 in range(1,npoints-1): #this for loop acts on points interior to the domain along the y-axis
            #uses the governing equation to update y_new for points interior to the domain 
            y_new[index1,index2] = y[index1,index2] + alpha*(y[(index1-1),index2] + y[(index1+1),index2] +y[index1,(index2-1)] + y[index1,(index2+1)] - (4*y[index1,index2]))
    
    for index1 in range(0,npoints): #this for loop acts on all points in the domain along the x-axis
        for index2 in range(0,npoints): #this for loop acts on all points in the domain along the y-axis
            y_transient[1].append(y_new[index1, index2]) #adds an aditional position to the second array y_transient and slots the value of the y_new[index1,index2] in that position 
    y_transient[0].append(iterations) #adds an aditional position to the first array y_transient and slots the value of the iterations in that position 
    
    # increments the variable iterations by 1 after the for loop sucessfully updates all the points interior to the domain 
    iterations = iterations + 1
    
    error_mag = 0; # resets the error to zero so that the error may be appropriately calculated
    
    # updates the error by summing the abosolute difference between the new and old values
    for index1 in range(1,npoints-1): #this for loop acts on points interior to the domain along the x-axis
        for index2 in range(1,npoints-1): #this for loop acts on points interior to the domain along the y-axis
            error_mag = error_mag + (abs(y[index1,index2]-y_new[index1,index2]))
            
    # records the error in an array 
    error_track.append(error_mag) #adds an aditional position to the array error_mag and slots the value of the error in that position 
            
    #the following loop is for the benifit of the coder to verify how the values of the code are proceeding 
    # the following if statement activates the chunk of code underneath if the iterations are a divisor of 1000
    #In other words if the iterations reaches a value that is a multiple of 1000 complete the following chunk of code
    if math.remainder(iterations,1000) == 0 : 
        print(iterations) # print the number of iterations
        print(error_mag) # print the value of the error 

    # updates the value of y by copying the values of y_new         
    for index1 in range(1,npoints-1): #this for loop completes the command below for points interior to the domain
        for index2 in range(1,npoints-1): #this for loop completes the command below for points interior to the domain
            y[index1,index2] = y_new[index1,index2]  #this transfers the value stored in y_new for particular set of indices to the same set of indices in y
    
print("the number of cycles this problem required to converge to a solution is " +str(iterations)+ " ")
time1 = iterations*dt;
print("the time this problem required to converge to a solution is " +str(time1)+ "s")



# Step 7 (Although not officially noted, this step is required to be able to present the final solution of the problem)

a = np.reshape(y_transient[1],(iterations,npoints,npoints)) #reshapes y_transient into an object that contains x npoints by npoints arrays where x is the number of iterations

iteratorArray = np.zeros((npoints,npoints)) #this creates an array of npoints where every term in the array is initialized as 0 

#this for acts on all vertices marking the grid over the 2D plate (i.e. domain)
for index3 in range(0,npoints):
   iteratorArray[index3,:] = index3 #this sets every value of each of the terms in array equal to its index (i.e. [0,1,2,..,n])
   for index4 in range(0,npoints):
       iteratorArray[:,index4] = index4 #this sets every value of each of the terms in array equal to its index (i.e. [0,1,2,..,n]) 

x1_dom = iteratorArray[index3,:] * h # multiply all elements of the array by the length of the edge of a cell (i.e. Delta x)
y1_dom = iteratorArray[index4,:] * h # multiply all elements of the array by the length of the edge of a cell (i.e. Delta y)
#The end result is that each position in the vector x_dom corresponds to a position along a row of the 2D square plate (i.e. domain) 
#The end result is that each position in the vector y_dom corresponds to a position along a column of the 2D square plate (i.e. domain) 
#in other words the result is a discretized version of points along the 2D square plate (i.e. domain)

X, Y = np.meshgrid(x1_dom, y1_dom)

lst = np.zeros(iterations) #sets up a two 2d array that contain all zeros

#sets the value found in array equal to its index where th index corresponds to a particular number of iterations
for index3 in range(0,iterations):
    lst[index3] = index3  

time = dt * (lst) # multiplies all elements in lst by the time it takes to complete a single iteration (dt). 
#The end result is that each position in the array 'time' corresponds to the time it takes to complete a specific number of iterations 
#said number of iterations is equal to the index of the array 'time' (for example t[9] is equal to the amount of time it takes to complete 9 iterations) 

#this graph is to show how quickly the error falls below the threshold value 
#which is a proxy for how quickly the unsteady heat diffusion reaches a steady state solution 
n = plt.scatter(time,error_track, color="red",s=1);
plt.xlabel("Time (s)")
plt.ylabel("error (no units)")
plt.title("how error changes with respect to time")
plt.show(n)


#the following code takes a particular timestep (i.e. iteration) and relates it to particular array (i.e. y_[index1,index2])
#what this means is that the following code can determine the temperature distribution of the 2D square plate for any timestep 
#determines what the temperature distribution is for a given timestep 

timestep_selected = iterations-1; #takes the timestep to be the last element in lst array 
y_timestep = a[timestep_selected,:,:] #selects the temperature distribution for the selected timstep 
#the timesteps points to an index attached to a particular npoints by npoints array 
#the a pulls out the 1 x npoints x npoints matrix and assigns it to y_timestep
y_timestep = np.reshape(y_timestep,(npoints,npoints)) #y_timestep is reshaped in a npoints x npoints array 
m = plt.contourf(X, Y, y_timestep,20,cmap = 'plasma'); #plots the countour plot for that particular timestep 
cbar = plt.colorbar(m)
plt.xlabel("Length of the plate along the x-axis (m)")
plt.ylabel("Length of the plate along the y-axis (m)")
plt.title("The temperature distribution at " +str(timestep_selected*dt)+ "s")



