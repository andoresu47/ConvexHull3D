/* 
 * MPI parallel implementation of a O(n log n) divide-and-conquer algorithm to find the convex hull in 3D, 
 * based on the kinetic-2d approach proposed in the article: "A Minimalist’s Implementation 
 * of the 3-d Divide-and-Conquer Convex Hull Algorithm", by Timothy M. Chan. 
 * 
 * input: coordinates file "points.in"
 * 		n
 *		{x_1, y_2, z_3}
 *		...
 *		{x_n, y_n, z_n}
 * 
 * output: indices of facets
 * 		{i_1, j_1, k_1}
 * 		{i_2, j_2, k_2}
 * 		...
 * 
 * This code assumes the following of the input points:
 * 	- No 3 points are colinear.
 * 	- No 4 points lie on the same plane. 
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>  
#include <stddef.h>


const double INF = 1e30;
const int NIL = -1;

/**
 * Point data structure, with emulated pointers to next and previous points representing edges.
 */ 
typedef struct Point{ 
	double x, y, z; 
	int prev, next;
} Point;

/*
 * Function to remove the file extension of a file, thus returning only the filename. 
 */
char *remove(char* mystr) {
    char *retstr;
    char *lastdot;
    if (mystr == NULL)
         return NULL;
    if ((retstr = (char*)malloc(strlen(mystr) + 1)) == NULL)
        return NULL;
    strcpy (retstr, mystr);
    lastdot = strrchr (retstr, '.');
    if (lastdot != NULL)
        *lastdot = '\0';
    return retstr;
}

/*
 * Function for determining the inequality/equality relationship between two Points based on their X coordinate.  
 */
int Comparator(const void * a, const void * b){
	Point point1 = *((Point*)a);
    Point point2 = *((Point*)b);

    if (point1.x > point2.x){
        return 1;
    }
    else if (point1.x == point2.x){
        return 0;
    }
    else{ 
        return -1;
	}
}

/*
* Function to determine if any of a triplet of point pointers is null. 
*/
inline bool hasnil(int p, int q, int r){
	return (p == NIL || q == NIL || r == NIL);
}

/*
* Function to determine if points p-q-r form a clockwise/left turn, 
* or a counterclockwise/right turn. Returns < 0 and > 0 respectively. 
*/
inline double turn(Point *p, Point *q, Point *r){  
	return (q->x - p->x) * (r->y - p->y) - 
			(r->x - p->x) * (q->y - p->y);
}

/*
* Function to determine the time when p-q-r will change between a left/right turn (ie. become colinear).
*/
inline double time(int pIndex, int qIndex, int rIndex, Point *head){ 
	if (pIndex == NIL || qIndex == NIL || rIndex == NIL){
		return INF;
	}
	Point *p = (head + pIndex);
	Point *q = (head + qIndex);
	Point *r = (head + rIndex);
	return ((q->x - p->x) * (r->z - p->z) - 
			(r->x - p->x) * (q->z - p->z)) / turn(p, q, r);
}

/*
 * Function to insert delete points from the linked-list-like point structure.
 */
void act(int pointIndex, Point *head){
	Point *prevP, *nextP;
	
	prevP = head + (head + pointIndex)->prev;
	nextP = head + (head + pointIndex)->next;
	
	if ((head + pointIndex)->prev == NIL){
		if ((head + pointIndex)->next != NIL){
			nextP->prev = pointIndex;
		}
	}
	else{
		// Insert point
		if (prevP->next != pointIndex){
			prevP->next = pointIndex;
			if ((head + pointIndex)->next != NIL){
				nextP->prev = pointIndex;
			}
		}
		// Delete point
		else{
			prevP->next = (head + pointIndex)->next;
			nextP->prev = (head + pointIndex)->prev;
		}
	}
}

/*
 * Function to re-index a point array to match its parent indices when combined. 
 */
void reindexPoints(Point *head, int size, int offset){
	int i, p, n;
	for(i = 0; i < size; i++){
		p = (head + i)->prev;
		n = (head + i)->next;
		if(p != NIL){
			(head + i)->prev = p + offset;
		}
		if(n != NIL){
			(head + i)->next = n + offset;
		}
	}
}

/*
 * Function to re-index an event array to match its reindexed point array. 
 */
void reindexEvents(int *A, int offset){
	int i;
	for(i = 0; A[i] != NIL; i++){
		A[i] = A[i] + offset;
	}
}

/*
 * Function to get indices of array given two pointers: head and displaced.
 */
int getIndex(Point *head, Point *displaced){
	return (displaced - head);
}

/*
 * Function to compute the maximum depth of a processor in the tree.
 */
int maxDepth(int treeHeight, int rank){
	int h, i, t;
	for(h = 0, i = treeHeight; h < treeHeight; h++, i--){
		t = rank / pow(2, h);
		if((t % 2) == 1){
			break;
		}
	}
	return i;
}

/* 
 * Function to swap pointers to integer arrays.
 */
void swapArrays(int **pa, int **pb)
{
    int *pc; 
    pc  = *pa;
    *pa = *pb;
    *pb = pc;
}

/*
 * Merging step for the divide-and-conquer algorithm to find the convex hull in 3D. 
 */
void merge(bool bottom, Point *list, int n, int *A, int *B, Point *head) {	
	// u: end of left hull based on x coordinate
	// v: beginning of right hull based on x coordinate
	// mid: beginning of right hull based on x coordinate
	int u, v, mid;  
	double t[6], oldt, newt;  
	int i, j, k, l, minl, prevL, nextL, prevR, nextR;
	
	// End of list for u
	u = getIndex(head, list + n/2 -1);
	mid = v = getIndex(head, list + n/2);
  
	// Find initial bridge
	// Start in the middle and move the bridge vertices out to each side.
	for ( ; ; ){  
		// Lower hull: turn must be negative
		if(bottom){
			if (!hasnil(u, v, (head + v)->next) && turn((head + u), (head + v), (head + ((head + v)->next))) < 0){
				v = (head + v)->next;
			}
			else if (!hasnil((head + u)->prev, u, v) && turn((head + ((head + u)->prev)), (head + u), (head + v)) < 0){
				u = (head + u)->prev;  
			}
			else{
				break;
			}
		}
		// Upper hull: turn must be positive
		else{
			if (!hasnil(u, v, (head + v)->next) && turn((head + u), (head + v), (head + ((head + v)->next))) > 0){
				v = (head + v)->next;
			}
			else if (!hasnil((head + u)->prev, u, v) && turn((head + ((head + u)->prev)), (head + u), (head + v)) > 0){
				u = (head + u)->prev;  
			}
			else{
				break;
			}
		}
	}

	// Merge by tracking bridge uv over time
	// Progress through time in an infinite loop until no more insertion/deletion events occur
	for (i = k = 0, j = n/2*2, oldt = -INF; ; oldt = newt) {  
		// To avoid buffer underflow
		if (B[i] == NIL){
			prevL = NIL;
			nextL = NIL;
		}
		else{
			prevL = (head + B[i])->prev;
			nextL = (head + B[i])->next;
		}
		
		// To avoid buffer underflow
		if (B[j] == NIL){
			prevR = NIL;
			nextR = NIL;
		}
		else{
			prevR = (head + B[j])->prev;
			nextR = (head + B[j])->next;
		}
		
		// Compute time values
		t[0] = time(prevL, B[i], nextL, head);  
		t[1] = time(prevR, B[j], nextR, head);    
		t[2] = time(u, (head + u)->next, v, head);  
		t[3] = time((head + u)->prev, u, v, head);
		t[4] = time(u, (head + v)->prev, v, head); 
		t[5] = time(u, v, (head + v)->next, head);
		
		for (newt = INF, l = 0; l < 6; l++){
			if (t[l] > oldt && t[l] < newt){ 
				minl = l; 
				newt = t[l]; 
			}
		}
		
		if (newt == INF){
			// No events found, the whole hull is merged
			break;
		}
		
		switch (minl) {
			// Left side
			case 0:  
				if ((head + B[i])->x < (head + u)->x){
					A[k++] = B[i]; 
				}
				act(B[i++], head);
				break;
			// Right side
			case 1:  
				if ((head + B[j])->x > (head + v)->x){
					A[k++] = B[j];  
				}
				act(B[j++], head);
				break;
			// Bridge
			case 2:  
				u = (head + u)->next;
				A[k++] = u;  
				break;
			// Bridge
			case 3:  
				A[k++] = u;
				u = (head + u)->prev; 
				break;
			// Bridge
			case 4:  
				v = (head + v)->prev;
				A[k++] = v;  
				break;
			// Bridge
			case 5:  
				A[k++] = v;  
				v = (head + v)->next;
				break;
		}
	}
	
	// Mark the end of the merged hull
	A[k] = NIL;
	
	/*
	* Go back in time to update pointers.
	* The vertices that are on the merged 2d hull at t=inf are not connected to their neighbours anymore, when the event
	* where they were added/removed to the 2d hull occured, so we have to go back and reconnect them (to those neighbours) 
	* so that each A[k] can be used as a face, along with A[k].prev and A[k].next
	*/

	// Connect bridge endpoints
	(head + u)->next = v;  
	(head + v)->prev = u;  
	// Loop from end of A[k] to start
	for (k--; k >= 0; k--) 
		 // On the left or right hull
		if ((head + A[k])->x <= (head + u)->x || (head + A[k])->x >= (head + v)->x) {
			act(A[k], head);
			if (A[k] == u) 
				u = (head + u)->prev; 
			else if (A[k] == v) 
				v = (head + v)->next;
		}
		// Inside the bridge (so it was a bridge endpoint)
		else { 
			// Put the point between the current two bridge endpoints (u and v)
			(head + u)->next = A[k]; 
			(head + A[k])->prev = u; 
			(head + v)->prev = A[k]; 
			(head + A[k])->next = v;
			
			// Make it a bridge endpoint
			// Left side
			if ((head + A[k])->x < (head + mid)->x){
				u = A[k]; 
			}
			// Right side
			else{
				v = A[k];
			}
		}
}

/*
 * Main divide-and-conquer algorithm to find the convex hull in 3D 
 */
void hull(bool bottom, Point *list, int n, int *A, int *B, Point *head) {  
	// Base case is a single point
	if (n == 1) { 
		// Remove the point from the list
		list->next = NIL; 
		list->prev = NIL;
		// Return no faces for the hull
		A[0] = NIL;
		return; 
	}
		
	// Recurse on left and right sides, swapping the A and B event arrays
	hull(bottom, list, n/2, B, A, head);  							// build left side
	hull(bottom, list + n/2, n-n/2, B+n, A+n, head);		        // build right side

	merge(bottom, list, n, A, B, head);
}

/*
 * Print hull faces.
 */
void printHull(int *A, Point *head)
{
	printf("\nFaces on the hull:\n");
	for (int i = 0; A[i] != NIL; i++) { 
		printf("Face %d: {%d %d %d}\n", i, (head + A[i])->prev, A[i], (head + A[i])->next);
		act(A[i], head);
	}
}

/*
 * Function that implements the parallel divide-and-conquer algorithm for finding the 3D convex hull based on kinetic 2D hulls.
 */
void divideAndConquer3DHull(int height, int rank, int n, int localNumPoints, int localArraysize, 
							Point *PLocalL, Point *PLocalU, int *ALocalL, int *ALocalU, int *BLocalL, int *BLocalU, 
							MPI_Datatype mpi_point_type, MPI_Comm comm){
	// Local variables							
	int parent, rightChild, rightChildSize, myHeight;
	
	// Compute the local lower and upper 3D Convex Hulls
	hull(true, PLocalL, localNumPoints, ALocalL, BLocalL, PLocalL);
	hull(false, PLocalU, localNumPoints, ALocalU, BLocalU, PLocalU);
	
	// Merge local 3D convex hulls recursively
	myHeight = 0;
	while (myHeight < height) { 	// not yet at top
        parent = (rank & (~(1 << myHeight)));
        
        if (parent == rank) { 		// left child
		    rightChild = (rank | (1 << myHeight));
			rightChildSize = n / pow(2, maxDepth(height, rightChild));
                        
			swapArrays(&ALocalL, &BLocalL);
			swapArrays(&ALocalU, &BLocalU);
			
  		    // Receive arrays from right child
  		    MPI_Recv(&(PLocalL[rightChildSize]), rightChildSize, mpi_point_type, rightChild, 0, comm, MPI_STATUS_IGNORE);
			MPI_Recv(&(BLocalL[rightChildSize * 2]), rightChildSize * 2, MPI_INT, rightChild, 0, comm, MPI_STATUS_IGNORE);
			MPI_Recv(&(PLocalU[rightChildSize]), rightChildSize, mpi_point_type, rightChild, 1, comm, MPI_STATUS_IGNORE);
			MPI_Recv(&(BLocalU[rightChildSize * 2]), rightChildSize * 2, MPI_INT, rightChild, 1, comm, MPI_STATUS_IGNORE);
			
  		    // Merge convex hulls 
  		    merge(true, PLocalL, rightChildSize * 2, ALocalL, BLocalL, PLocalL);
			merge(false, PLocalU, rightChildSize * 2, ALocalU, BLocalU, PLocalU);
			
            myHeight++;

        } else { // right child
			// Re-index the arrays to match parent indexing
			reindexPoints(PLocalL, localArraysize, localArraysize);
			reindexEvents(ALocalL, localArraysize);
			reindexPoints(PLocalU, localArraysize, localArraysize);
			reindexEvents(ALocalU, localArraysize);
			
			// Send local point and event arrays to parent
            MPI_Send(PLocalL, localArraysize, mpi_point_type, parent, 0, comm);
			MPI_Send(ALocalL, localArraysize * 2, MPI_INT, parent, 0, comm);
			MPI_Send(PLocalU, localArraysize, mpi_point_type, parent, 1, comm);
			MPI_Send(ALocalU, localArraysize * 2, MPI_INT, parent, 1, comm);
			
            if(myHeight != 0){
				// Free memory
				free(PLocalL);  
				free(ALocalL);
				free(BLocalL);
				free(PLocalU);  
				free(ALocalU);
				free(BLocalU);
			}
            myHeight = height;
			return;
        }
    }
}
 

int main(int argc, char **argv){
	// Local variables
	int n, i, localArraysize, localNumPoints, height;
	Point *P, *PLocalL, *PLocalU; 
	int *ALocalL, *ALocalU, *BLocalL, *BLocalU;
	double startTime, localTime, totalTime;
	bool flag;
	
	// For printing results in testing phase
	int turn = 0;
	
	// Initialize MPI
	int rank, numProcs;
	MPI_Status status;
	
	// Create a type of struct Point
	MPI_Datatype mpi_point_type;
    MPI_Datatype types[5] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT};
	int blocklengths[5] = {1, 1, 1, 1, 1};
	MPI_Aint offsets[5] = {offsetof(Point, x), 
						   offsetof(Point, y), 
						   offsetof(Point, z), 
						   offsetof(Point, prev), 
						   offsetof(Point, next)};
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);	
	
	MPI_Type_create_struct(5, blocklengths, offsets, types, &mpi_point_type);
	MPI_Type_commit(&mpi_point_type);
	
	// Compute total height of tree
    height = log2(numProcs);
	
	// If process 0, allocate memory for global points array and fill with values
	if (rank == 0){		
	
		if(argc > 1){
			char path[128] = "../input/raw/";
			strcat(path, argv[1]);
            
            //printf("Reading input from %s\n", path);
		
			// Read input points from file
			FILE * infile;
			infile = fopen (path, "r");
	
			fscanf(infile, "%d\n", &n);
	
			P = (Point*) malloc (n * sizeof (Point));
	
			for(i = 0; i < n; i++){
				fscanf(infile, "{%lf, %lf, %lf}\n", &P[i].x, &P[i].y, &P[i].z);
				P[i].prev = P[i].next = NIL;
			}
		
			fclose(infile);
		
			// If input points were not already sorted, they are sorted and written to a new file. 
			if(argc > 2){
				char path2[128] = "../input/sorted/";
				strcat(path2, argv[1]);
                
                //printf("Writing sorted points to %s\n", path2);
				
				if(strcmp(argv[2], "-sort") == 0){
					qsort (P, n, sizeof(Point), Comparator);
			
					FILE * auxfile;
					auxfile = fopen (path2, "w");
					fprintf(auxfile, "%d\n", n);
					for(i = 0; i < n; i++){
						fprintf(auxfile, "{%lf, %lf, %lf}\n", P[i].x, P[i].y, P[i].z);
					}
					fclose(auxfile);
				}
			}
		}
		else{
			// If no input file is passed the program aborts
			printf("No input file provided\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
	}
	
	// Every processor needs to have the total number of points n
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
				
	// Allocate memory for local point array based on maximum depth of corresponding processor	
	localArraysize = n / pow(2, maxDepth(height, rank));
	PLocalL = (Point *) malloc(localArraysize * sizeof(Point));
	PLocalU = (Point *) malloc(localArraysize * sizeof(Point));
    
    // Scatter to fill with equally-sized value chunks
	localNumPoints = n / numProcs;
	MPI_Scatter(P, localNumPoints, mpi_point_type, PLocalL, localNumPoints, mpi_point_type, 0, MPI_COMM_WORLD);
	
	// Free memory
	if(rank == 0){
		free(P);
	}
	
	// Copy lower hull point array into the upper hull array
	memcpy(PLocalU, PLocalL, localArraysize * sizeof(Point));
		
	// Allocate memory for local main and temporal event arrays, for both lower and upper hulls
	ALocalL = (int *) malloc (2 * localArraysize * sizeof (int));
	ALocalU = (int *) malloc (2 * localArraysize * sizeof (int));
	BLocalL = (int *) malloc (2 * localArraysize * sizeof (int));
	BLocalU = (int *) malloc (2 * localArraysize * sizeof (int));
		
	// Find 3D Convex Hull
	startTime = MPI_Wtime();				// Start timing
	divideAndConquer3DHull(height, rank, n, localNumPoints, localArraysize, 
							PLocalL, PLocalU, ALocalL, ALocalU, BLocalL, BLocalU, 
							mpi_point_type, MPI_COMM_WORLD);
	localTime = MPI_Wtime() - startTime;	// End timing
	//printf("Process %d took %f seconds \n", rank, localTime);
	    
	// Compute the parallel execution time
    MPI_Reduce(&localTime, &totalTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	
	// As process 0 contains the final merged 3D convex hull, 
	// it writes the final result to file. 
	if (rank == 0){
		int *CLocalL, *CLocalU;
        char path3[128] = "../output/";
        strcat(path3, remove(argv[1]));
        strcat(path3, "_MPI.out");
        
		//printf("Execution time: %f seconds \n", totalTime);
		printf("%f", totalTime);
        
        //printf("Writing output faces to %s\n", path3);
		
		FILE * outfile;
		outfile = fopen (path3, "w");
		
		// To know which array holds the final result as a consequence of the A-B swapping.
		if(height % 2 == 0){
			CLocalL = ALocalL;
			CLocalU = ALocalU;
		}
		else{
			CLocalL = BLocalL;
			CLocalU = BLocalU;
		}
		
		for (int i = 0; CLocalL[i] != NIL; i++) { 
			fprintf(outfile, "{%d, %d, %d}\n", (PLocalL + CLocalL[i])->prev, CLocalL[i], (PLocalL + CLocalL[i])->next);
			act(CLocalL[i], PLocalL);
		}
		for (int i = 0; CLocalU[i] != NIL; i++) { 
			fprintf(outfile, "{%d, %d, %d}\n", (PLocalU + CLocalU[i])->prev, CLocalU[i], (PLocalU + CLocalU[i])->next);
			act(CLocalU[i], PLocalU);
		}
		fclose(outfile);
		
		// Free memory
		free(PLocalL);  
		free(ALocalL);
		free(BLocalL);
		free(PLocalU);  
		free(ALocalU);
		free(BLocalU);
	}
	
	MPI_Type_free(&mpi_point_type);
	MPI_Finalize();
	exit(0);
}
