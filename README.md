# ConvexHull3D
This project consists of a parallel implementation of a divide-and-conquer algorithm for finding the Convex Hull in 3D. 

## Algorithm
The implementation is based on Chan’s [minimalist 3-d convex hull algorithm](http://tmc.web.engr.illinois.edu/ch3d/ch3d.pdf), an O(n logn) time topdown divide-and-conquer algorithm for computing convex hulls in three dimensions . This algorithm works by treating the three-dimensional problem as a two-dimensional kinetic problem. This means that the 3-d set S of n points is seen as a moving image, or movie, of 2-d projected points as the points are rotated at an angle **α** about the x-axis. The main idea behind this algorithm’s correctness is that affine transformations preserve points, straight-lines and planes, as well as convexity. Hence, points p<sub>i</sub> whose projected points p'<sub>i</sub> belong to the projected 2-d hulls, also belong to the 3-d convex hull.

As points are rotated, the kinetic 2-d hull changes by a point p'<sub>j</sub> being inserted between end-points p'<sub>i</sub>p'<sub>k</sub> of an edge, or by a point p'<sub>j</sub> being deleted, destroying two edges (p'<sub>i</sub>, p'<sub>j</sub>), (p'<sub>j</sub>, p'<sub>k</sub>) and creating a new one (p'<sub>i</sub>, p'<sub>k</sub>). For both of these events, a 3-d hull triangle -or facet - (p<sub>i</sub>p<sub>j</sub>p<sub>k</sub>) is identified. A point can be inserted or deleted at most once. Thus, by keeping track of the insertion and deletion events of the kinetic 2-d movie, each identified as a triplet of points, the 3-d convex hull can be constructed in linear time.

The divide-and-conquer approach can be applied by sorting the points *P* by increasing x-coordinates, and recursively solving the kinetic 2-d problem for the left and right hulls L and R, consisting of the points p<sub>1</sub>,...,p<sub>⌊n/2⌋</sub> and p<sub>⌊n/2⌋+1</sub>,...,p<sub>n</sub> respectively. Merging the two 3-d convex hulls L and R, whose kinetic 2-d event lists are now known, is done by keeping track of the lower (or upper) tangent joining the projected hulls. As the points are rotated, the merged hull changes by individual changes of the hulls L and R or by the lower (upper) tangent no longer supporting the neighbors of the end-points from below (above).

## Implementation
The algorithm was implemented by using 3 arrays as the main data structures: an array of points P, and two integer arrays A and B, representing event lists. The basic data structure for the algorithm is the Point, which is implemented as a struct. The struct attributes are the x, y and z coordinates, as well as two integers prev and next which represent if the point forms a 3-d hull facet with other points identified by their index in the local points array P. This emulates a doubly linked list without the hassle of reconstructing pointers each time data is shared among processors.

The algorithm’s main data structure consists of the event lists. An event in the theoretical scheme is represented as a triplet of points forming a 3-d hull triangle. However, to avoid data redundancy and costly send-receive operations, in this implementation an event is represented by a simple integer. This integer corresponds to a an entry in the local points array P. As the algorithm progresses, each point next and prev ”pointers” are updated to represent edges. As such, an event consisting of a triplet of points (p<sub>i</sub>, p<sub>j</sub>, p<sub>k</sub>) can be simplified as only considering the middle element p<sub>j</sub>, identified by its index j.

Array A always contains the final list of events representing a 3-d hull. On the other hand, array B is used as an auxiliary data structure holding the event lists of two hulls to be merged, each occupying one of the halves of B. For the actual Hull algorithm, these arrays form part of the arguments, and get modified in-place. To preserve the above conditions on arrays A and B, the two event arrays are swapped at each recursive call.

As can be implied from the above description, the points array P consists of a list of Point structs. It should be noted that the later is not a native data type which can simply be sent and received with MPI. An MPI data type had to be compiled for the purposes of sending and receiving points as usually done with atomic types such as INT or FLOAT.

Note that the implementation relies not only on event lists, but on point lists as well. At each communication step, a child process sends its local points array PLocal along with its event lists ALocal and BLocal, as the 3-d hull representation depends on both. For keeping the indices consistent in the parent process when receiving data from a child process, a mapping routine is executed before the sending takes place, so that each reference to an index is displaced an amount equal to the child’s number of points.

## Structure
The implementation is self contained, in that all the functionality is contained in the C++ file. The program expects the name of a points file as its input, with an additional flag indicating if the points are to be sorted or not (if they are already sorted when fed into the program). To run from the command line, the command is: 

```
mpirun -np 8 ./divide_and_conquer_MPI normal_origin_4096.in -sort
```

The input and output files get generated according to the following project file structure:

Insert image...

The input files consist of coordinates x, y and z in the format:

```
<number of points n>
{p_1.x, p_1.y, p_1.z}
{p_2.x, p_2.y, p_2.z}
...
{p_n.x, p_n.y, p_n.z}
```

The output are the events, or facets, of the 3-d hull represented as triplets of integers which make reference to the positions in a list formed by the sorted input points. An example output file looks like the following:

```
{2700, 2753, 2973}
{2753, 2973, 3349}
{1084, 2700, 2753}
...
{1084, 2669, 2700}
```

## Simulations
Additionally to the C++ code, a Mathematica notebook is included. This was used to visually asses the correctness of the implementation by plotting each triangle output from the main algorithm. The notebook was also used to simulate the minimalist algorithm in action, by rotating a set of 3D points and plotting triangle facets on the fly. This is particularly useful for exposition purposes. 
