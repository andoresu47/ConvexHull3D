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
 
#define MIN(a,b) ((a) < (b) ? (a) : (b))

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>  
#include <stddef.h>


const double INF = 1e30;
const int NIL = -1;

// Point data structure
typedef struct Point{ 
	double x, y, z; 
	int prev, next;
} Point;

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
	// Insert point
	if ((head + (head + pointIndex)->prev)->next != pointIndex){
		(head + (head + pointIndex)->prev)->next = pointIndex;
		(head + (head + pointIndex)->next)->prev = pointIndex;
	}
	// Delete point
	else{
		(head + (head + pointIndex)->prev)->next = (head + pointIndex)->next;
		(head + (head + pointIndex)->next)->prev = (head + pointIndex)->prev;
	}
}

/*
 * Function to map local to global indices.
 */
int getGlobalIndex(int numPts, int numProc, int rank, int idx){
	int L = numPts / numProc;
    int R = numPts % numProc;
    int mu_prefix = rank * L + MIN(rank, R);
    int mu;
	
	return idx + mu;
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
	int u, v, mid;  
	double t[6], oldt, newt;  
	int i, j, k, l, minl;
	
	// End of list for u
	u = getIndex(head, list + n/2 -1);
	mid = v = getIndex(head, list + n/2);
  
	// Find initial bridge
	// Start in the middle and move the bridge vertices out to each side.
	for ( ; ; ){  
		int p1, q1, r1, 
			p2, q2, r2;
			
		p1 = u; q1 = v; r1 = (head + v)->next;
		p2 = (head + u)->prev; q2 = u; r2 = v;
		
		// Lower hull
		if(bottom){
			if (!hasnil(p1, q1, r1) && turn((head + p1), (head + q1), (head + r1)) < 0){
				v = (head + v)->next;
			}
			else if (!hasnil(p2, q2, r2) && turn((head + p2), (head + q2), (head + r2)) < 0){
				u = (head + u)->prev;  
			}
			else{
				break;
			}
		}
		// Upper hull
		else{
			if (!hasnil(p1, q1, r1) && turn((head + p1), (head + q1), (head + r1)) > 0){
				v = (head + v)->next;
			}
			else if (!hasnil(p2, q2, r2) && turn((head + p2), (head + q2), (head + r2)) > 0){
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
		t[0] = time((head + B[i])->prev, B[i], (head + B[i])->next, head);  
		t[1] = time((head + B[j])->prev, B[j], (head + B[j])->next, head);    
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
	hull(bottom, list + n/2, n-n/2, B+n/2*2, A+n/2*2, head);		// build right side

	merge(bottom, list, n, A, B, head);
}

/*
 * Print faces on the hull
 */
void printHull(int *A, Point *head)
{
	printf("\nFaces on the hull:\n");
	for (int i = 0; A[i] != NIL; i++) { 
		printf("Face %d: {%d %d %d}\n", i, (head + A[i])->prev, A[i], (head + A[i])->next);
		act(A[i], head);
	}
}


int main(int argc, char **argv){
	/* Local variables */
	int n, i, localArraysize, localNumPoints, height;
	int parent, rightChild, rightChildSize, myHeight;	
	Point *P, *PLocal; 
	int *ALocal, *BLocal;
	
	// For printing results in testing phase
	int turn = 0;
	
	/* Initialize MPI */
	int rank, numProcs;
	MPI_Status status;
	
	/* Create a type of struct Point */
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
	
	/* Compute total height of tree */
    height = log2(numProcs);
	
	/* If process 0, allocate memory for global points array and fill with values */
	if (rank == 0){
		// Read input points from file
		FILE * infile;
		infile = fopen ("points.in", "r");
	
		fscanf(infile, "%d\n", &n);
	
		P = (Point*) malloc (n * sizeof (Point));
	
		for(i = 0; i < n; i++){
			fscanf(infile, "{%lf, %lf, %lf}\n", &P[i].x, &P[i].y, &P[i].z);
			P[i].prev = P[i].next = NIL;
		}
		
		fclose(infile);
	
		// If input points were not already sorted, they are sorted and written to a new file. 
		if(argc > 1){
			if(strcmp(argv[1], "-sort") == 0){
				qsort (P, n, sizeof(Point), Comparator);
			
				FILE * auxfile;
				auxfile = fopen ("points_srtd.in", "w");
				fprintf(auxfile, "%d\n", n);
				for(i = 0; i < n; i++){
					fprintf(auxfile, "{%lf, %lf, %lf}\n", P[i].x, P[i].y, P[i].z);
				}
				fclose(auxfile);
			}
		}
	}
	
	/* Every processor needs to have the total number of points n */
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
				
	/* Allocate memory for local point array based on maximum depth of corresponding processor */	
	localArraysize = n / pow(2, maxDepth(height, rank));
	PLocal = (Point *) malloc(localArraysize * sizeof(Point));
	
	/* Scatter to fill with equally-sized value chunks */
	localNumPoints = n / numProcs;
	MPI_Scatter(P, localNumPoints, mpi_point_type, PLocal, localNumPoints, mpi_point_type, 0, MPI_COMM_WORLD);
		
	/* Allocate memory for local main and temporal event arrays */
	ALocal = (int *) malloc (2 * localArraysize * sizeof (int));
	BLocal = (int *) malloc (2 * localArraysize * sizeof (int));
	
	/* Compute the local 3D Convex Hull */
	hull(true, PLocal, localNumPoints, ALocal, BLocal, PLocal);
	
	/* Merge local 3D convex hulls recursively */
	myHeight = 0;
	while (myHeight < height) { 	// not yet at top
        parent = (rank & (~(1 << myHeight)));

        if (parent == rank) { 		// left child
		    rightChild = (rank | (1 << myHeight));
			rightChildSize = n / pow(2, maxDepth(height, rightChild));
			
			swapArrays(&ALocal, &BLocal);
			
  		    /* Receive arrays from right child */
  		    MPI_Recv(&(PLocal[rightChildSize]), rightChildSize, mpi_point_type, rightChild, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(&(BLocal[rightChildSize * 2]), rightChildSize * 2, MPI_INT, rightChild, 0, MPI_COMM_WORLD, &status);

  		    // Merge convex hulls 
  		    merge(true, PLocal, rightChildSize * 2, ALocal, BLocal, PLocal);
			
            myHeight++;

        } else { // right child
			// Re-index the arrays to match parent indexing
			reindexPoints(PLocal, localArraysize, localArraysize);
			reindexEvents(ALocal, localArraysize);
			
			// Send local point and event arrays to parent
            MPI_Send(PLocal, localArraysize, mpi_point_type, parent, 0, MPI_COMM_WORLD);
			MPI_Send(ALocal, localArraysize * 2, MPI_INT, parent, 0, MPI_COMM_WORLD);
            if(myHeight != 0){
				// Free memory
				free(PLocal);  
				free(ALocal);
				free(BLocal);
			}
            myHeight = height;
        }
    }
	
	// Process 0 contains the final merged 3D convex hull. 
	// Process 0 writes the final result to file. 
	if (rank == 0){
		FILE * outfile;
		outfile = fopen ("output_DC_MPI.out", "w");
		
		for (int i = 0; ALocal[i] != NIL; i++) { 
			fprintf(outfile, "{%d, %d, %d}\n", (PLocal + ALocal[i])->prev, ALocal[i], (PLocal + ALocal[i])->next);
			act(ALocal[i], PLocal);
		}
		fclose(outfile);
	}
	
	MPI_Type_free(&mpi_point_type);
	MPI_Finalize();
	exit(0);
}
