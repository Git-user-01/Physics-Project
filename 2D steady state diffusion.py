import numpy as np
import matplotlib.pyplot as plt

#This code will find the final temperature distribution of a 2D square plate undergoing steady state heat diffusion 
#the boundary conditions are such that the bottom row of the plate is held steady at 1 C where C = Celsius 
#The leftmost and rightmost column of the plate are both held at 0 C 
#and the top row of the plate is similarly held at 0 C 

#step 1
dom_size = 1; #defines both the length and width of the 2D square plate (i.e. domain) as being 1 meter



#step 2 
# here the mesh of a 2d square plate (i.e. domain) will be defined 
npoints = 34; #defines both the number of columns and rows for a 2d square plate 
h = dom_size/(npoints-1); #defines the length of the edge of each cell in the mesh 
#the value of h corresponds to delta x and delta y of the 2D square plate (i.e. domain)



#step 3 and 4 
y = np.zeros((npoints,npoints)) #sets up a two 2d array that contain all zeros 
#this also serves as a way to intialize the interior points of the domain and set the boundary condition for the top row, leftmost column, and the rightmost column
y_new = np.zeros((npoints,npoints)) #sets up a two 2d array that contain all zeros

y[0,:] = 1 #sets the boundary condition of the bottom row for y at a temperature of 1 C
y[:,0] = 1 
y_new[0,:] = 1 #sets the boundary condition of the bottom row for y_new at a temperature of 1 C
y_new[:,0] = 1 

#step 5 and 6 
error_mag = 1; #initialize the variable error_mag (just to make sure its value is greater than error_req)
error_req = 1e-6; # establishes the threshold value the error must fall below to break the following loop  
iterations = 0; # initializes the variable iterations as having a value of 0 
#the variable iterations is there to determine the number of cycles the code will have to do to reach a stable solution

while (error_mag > error_req): # keeps doing the calculations until the error is less than the threshold value 
    
    # the range ensures that the for loop focuses on the points interior to the domain and does not modify the boundary conditions
    for index1 in range(1,npoints-1): #this for loop acts on points interior to the domain along the x-axis
        for index2 in range(1,npoints-1): #this for loop acts on points interior to the domain along the y-axis
            y_new[index1,index2] = 0.25*(y[index1-1,index2] + y[index1+1,index2] +y[index1,index2-1] + y[index1,index2+1])
            #uses the governing equation to update y_new for points interior to the domain 
    
    # increments the variable iterations by 1 after the for loop sucessfully updates all the points interior to the domain 
    iterations = iterations + 1
        
    error_mag = 0; # resets the error to zero so that the error may be appropriately calculated
    
    # updates the error by summing the abosolute difference between the new and old values
    for index1 in range(1,npoints-1): #this for loop acts on points interior to the domain along the x-axis
        for index2 in range(1,npoints-1): #this for loop acts on points interior to the domain along the y-axis
            error_mag = error_mag + (abs(y[index1,index2]-y_new[index1,index2]))
        
    # updates the value of y by copying the values of y_new 
    for index1 in range(1,npoints-1): #this for loop completes the command below for points interior to the domain
        for index2 in range(1,npoints-1): #this for loop completes the command below for points interior to the domain
            y[index1,index2] = y_new[index1,index2] # this transfers the value stored in y_new for particular set of indices to the same set of indices in y

print("the number of cycles this problem required to converge to a solution is " +str(iterations)+ " ")



# Step 7 (Although not officially noted, this step is required to be able to present the final solution of the problem)

iteratorArray1 = np.zeros((npoints,npoints)) #this creates an array of npoints where every term in the array is initialized as 0 

#this for acts on all vertices marking the grid over the 2D plate (i.e. domain)
for index3 in range(0,npoints):
   iteratorArray1[index3,:] = index3 #this sets every value of each of the terms in array equal to its index (i.e. [0,1,2,..,n]) 
   for index4 in range(0,npoints):
       iteratorArray1[:,index4] = index4 #this sets every value of each of the terms in array equal to its index (i.e. [0,1,2,..,n]) 

x_dom = iteratorArray1[index3,:] * h # multiply all elements of the array by the length of the edge of a cell (i.e. Delta x)
y_dom = iteratorArray1[index4,:] * h # multiply all elements of the array by the length of the edge of a cell (i.e. Delta y)
#The end result is that each position in the vector x_dom corresponds to a position along a row of the 2D square plate (i.e. domain) 
#The end result is that each position in the vector y_dom corresponds to a position along a column of the 2D square plate (i.e. domain) 
#in other words the result is a discretized version of points along the 2D square plate (i.e. domain)

X, Y = np.meshgrid(x_dom, y_dom)

plt.contourf(X, Y, y,20,cmap = 'plasma'); #contour plot of how temperature varies with location 
plt.xlabel("Length of the plate along the x-axis (m)")
plt.ylabel("Length of the plate along the y-axis (m)")
plt.title("The temperature distribution at "+str(iterations*0.0001)+"s")
plt.colorbar()

print(iterations*0.0001)







