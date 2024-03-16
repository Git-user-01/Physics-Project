import numpy as np
import matplotlib.pyplot as plt

#This code will find the final temperature distribution of a 2D square plate undergoing an unsteady convection diffusion equation  
#the boundary conditions are such that the bottom row and the leftmost column of the plate is held steady at 1 C where C = Celsius 
#the rightmost column and the top row of the plate are both held at 0 C 

#step 1
dom_length = 1; #defines both the length and width of the 2D square plate (i.e. domain) as being 1 meter



#step 2
# here the mesh of a 2d square plate (i.e. domain) will be defined
n_points = 34; #defines both the number of columns and rows for a 2d square plate
h1 = dom_length/(n_points-1); #defines the length of the edge of each cell in the mesh
#the value of h corresponds to delta x and delta y of the 2D square plate (i.e. domain)
x = np.linspace(0,dom_length,n_points) #x coordinate of each point 
y = np.linspace(0,dom_length,n_points) #y coordinate of each point 

dt = 0.0001; #defines Delta t measured in seconds (s)
#functionally this means how much time has elapsed following the completion of a single iteration 
# in this case an iteration is a single cycle of the 'for loop' containing the governing equation for an unsteady heat diffusion equation 
rho = 1; # the density of the fluid is 1 
u = 0; # sets the horizontal speed of the flow of the liquid 
v = 0; # sets the vertical speed of the flow of the liquid 
gamma = 1; # set the value of the diffusion coefficient
P = (rho*u*h1)/gamma; #sets the Peclet number 



#steps 3 and 4 
T = np.zeros((n_points,n_points)) #sets up a two 2d array that contain all zeros 
#this also serves as a way to intialize the interior points of the domain and set the boundary condition for the top row and the rightmost column
T_new = np.zeros((n_points,n_points)) #sets up a two 2d array that contain all zeros 

T[0,:] = 1; #sets the boundary condition of the bottom row for y at a temperature of 1 C
T[:,0] = 1; #sets the boundary condition of the leftmost column for y at a temperature of 1 C
T_new[0,:] = 1; #sets the boundary condition of the bottom row for y_new at a temperature of 1 C
T_new[:,0] = 1; #sets the boundary condition of the leftmost column for y at a temperature of 1 C



#steps 5 and 6 
error = 1;  #initialize the variable error_mag (just to make sure its value is greater than error_req)
error_req = 1e-6; # establishes the threshold value the error must fall below to break the following loop  
error_track = []; #tracks how the error changes with time 
iterations = 0; # initializes the variable 'iterations' as having a value of 0 
#the variable 'iterations' is there to determine the number of cycles the code will have to do to reach a stable solution

while error > error_req: # keeps doing the calculations until the error is less than the threshold value 
    
    # the range ensures that the for loop focuses on the points interior to the domain and does not modify the boundary conditions
    for i in range(1,n_points-1): #this for loop acts on points interior to the domain along the x-axis
        for j in range(1,n_points-1): #this for loop acts on points interior to the domain along the y-axis
            #uses the governing equation to update y_new for points interior to the domain 
            a_E = gamma - (rho*u*h1)/2; 
            a_W = gamma + (rho*u*h1)/2; 
            a_N = gamma - (rho*v*h1)/2; 
            a_S = gamma + (rho*v*h1)/2; 
            a_P =  (rho*u*h1)/2 - (rho*u*h1)/2 + (rho*v*h1)/2 - (rho*v*h1)/2 + gamma + gamma + gamma + gamma; 
            T_new[i,j] = ((a_E*T[i+1,j]) + (a_W*T[i-1,j]) + (a_N*T[i,j+1]) + (a_S*T[i,j-1]))/a_P; 
    
    # increments the variable iterations by 1 after the for loop sucessfully updates all the points interior to the domain 
    iterations = iterations + 1;
    
    error = 0; # resets the error to zero so that the error may be appropriately calculated
    
    # updates the error by summing the abosolute difference between the new and old values
    for i in range(1,n_points-1): #this for loop acts on points interior to the domain along the x-axis
        for j in range(1,n_points-1): #this for loop acts on points interior to the domain along the y-axis
            error = error + abs(T[i,j] - T_new[i,j]); 
    
    # records the error in an array 
    error_track.append(error) #adds an aditional position to the array error_track and slots the value of the error in that position 
    
    # updates the value of y by copying the values of y_new 
    for i in range(1,n_points-1): #this for loop completes the command below for points interior to the domain
        for j in range(1,n_points-1): #this for loop completes the command below for points interior to the domain
            T[i,j] = T_new[i,j] #this transfers the value stored in y_new for particular set of indices to the same set of indices in y

print("the number of cycles this problem required to converge to a solution is " +str(iterations)+ " ")
time1 = iterations*dt;
print("the time this problem required to converge to a solution is " +str(time1)+ "s")



# Step 7 (Although not officially noted, this step is required to be able to present the final solution of the problem)

iteratorArray1 = np.zeros((n_points,n_points)) #this creates an array of npoints where every term in the array is initialized as 0 

#this for acts on all vertices marking the grid over the 2D plate (i.e. domain)
for index3 in range(0,n_points):
   iteratorArray1[index3,:] = index3 #this sets every value of each of the terms in array equal to its index (i.e. [0,1,2,..,n])
   for index4 in range(0,n_points):
       iteratorArray1[:,index4] = index4 #this sets every value of each of the terms in array equal to its index (i.e. [0,1,2,..,n]) 

x1_dom = iteratorArray1[index3,:] * h1 # multiply all elements of the array by the length of the edge of a cell (i.e. Delta x)
y1_dom = iteratorArray1[index4,:] * h1 # multiply all elements of the array by the length of the edge of a cell (i.e. Delta y)
#The end result is that each position in the vector x_dom corresponds to a position along a row of the 2D square plate (i.e. domain) 
#The end result is that each position in the vector y_dom corresponds to a position along a column of the 2D square plate (i.e. domain) 
#in other words the result is a discretized version of points along the 2D square plate (i.e. domain)

X, Y = np.meshgrid(x1_dom, y1_dom) #creates the mesh to be used 

lst = np.zeros(iterations) #sets up an array that contain all zeros

#sets the value found in array equal to its index where the index corresponds to a particular number of iterations
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



# step 8 (this is not mentioned in the document but answers the second question posed in the poster)

#the following set of arrays are data about how the 'time to converge' changes based on the value of the diffusion coefficient 
gamma1 = [0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
times2 = [0.1013,0.1145,0.2156,0.2636,0.2861,0.2976,0.3040,0.3079,0.3104,0.312,0.3131] #time to converge for varying values of the diffusion coefficient when u=v=1
times3 = [0.0416,0.0479,0.1145,0.1732,0.2156,0.2443,0.2636,0.2768,0.2861,0.2921,0.2976] #time to converge for varying values of the diffusion coefficient when u=v=2
times4 = [0.0239,0.0277,0.0697,0.1145,0.1554,0.1891,0.2156,0.236,0.2516,0.2636,0.273] #time to converge for varying values of the diffusion coefficient when u=v=3

#the following set of arrays are a truncated form of the data above 
gamma11 = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
times21 = [0.1145,0.2156,0.2636,0.2861,0.2976,0.3040,0.3079,0.3104,0.312,0.3131] #time to converge for varying values of the diffusion coefficient when u=v=1
times31 = [0.0479,0.1145,0.1732,0.2156,0.2443,0.2636,0.2768,0.2861,0.2921,0.2976] #time to converge for varying values of the diffusion coefficient when u=v=2
times41 = [0.0277,0.0697,0.1145,0.1554,0.1891,0.2156,0.236,0.2516,0.2636,0.273] #time to converge for varying values of the diffusion coefficient when u=v=3

# the following lines of code draw a graph that shows the relationship between the diffusion coefficient and the 'time to converge'
w = plt.scatter(gamma11,times21,color="red", label = 'u = v = 1 m/s');
f = plt.scatter(gamma11,times31,color="blue", label = 'u = v = 2 m/s');
g = plt.scatter(gamma11,times41,color="purple", label = 'u = v = 3 m/s');
plt.xlabel("diffusion coefficient (m^2/s)")
plt.ylabel("Time to converge (s)")
plt.title("how 'time to converge' changes with respect to "+str('\u0393')+" ")
plt.legend(loc='lower right')
plt.show()

#plots the countour plot for the stable solution 
m = plt.contourf(X, Y, T,20,cmap = 'plasma'); 
cbar = plt.colorbar(m)
plt.xlabel("Length of the plate along the x-axis (m)")
plt.ylabel("Length of the plate along the y-axis (m)")
plt.title("The temperature distribution where gamma is " +str(time1)+" s")



    

