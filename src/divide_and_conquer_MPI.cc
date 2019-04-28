/* 
 * Parallel implementation of a O(n log n) divide-and-conquer algorithm to find the convex hull in 3D, 
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
#include <stddef.h>


const double INF = 1e30;
const int NIL = -1;

// Point data structure
struct Point{ 
	double x, y, z; 
	int prev, next;
};

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


int main(int argc, char **argv){
	/* Local variables */
	int n, i;  
	
	/* Initialize MPI */
	int rank, size;
	MPI_Status status;
	
	/* For creating a type of struct Point */
	MPI_Datatype mpi_point_type;
    MPI_Datatype types[5] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT};
	int blocklengths[5] = {1, 1, 1, 1, 1};
	MPI_Aint offsets[5] = {offsetof(Point, x), 
						   offsetof(Point, y), 
						   offsetof(Point, z), 
						   offsetof(Point, prev), 
						   offsetof(Point, next)};
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);	
	
	MPI_Type_create_struct(5, blocklengths, offsets, types, &mpi_point_type);
	MPI_Type_commit(&mpi_point_type);
	
	if (rank == 0){
		// Read input points from file
		FILE * infile;
		infile = fopen ("points.in", "r");
	
		fscanf(infile, "%d\n", &n);
		printf("%d\n", n);
	
		Point *P = (Point*) malloc (n * sizeof (Point));
	
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

		/* Print points 
		for(i = 0; i < n; i++){
			printf("%lf %lf %lf\n", P[i].x, P[i].y, P[i].z);
		}
		*/
		
		MPI_Send(&(P[0]), 1, mpi_point_type, 1, 0, MPI_COMM_WORLD);
		
		printf("Rank %d: sent structure point: {%lf, %lf, %lf}, [%d, %d]\n", rank, P[0].x, P[0].y, P[0].z, P[0].prev, P[0].next);
	}
	if(rank == 1){
		Point recv;
		MPI_Recv(&recv, 1, mpi_point_type, 0, 0, MPI_COMM_WORLD, &status);
		printf("Rank %d: Received structure point: {%lf, %lf, %lf}, [%d, %d]\n", rank, recv.x, recv.y, recv.z, recv.prev, recv.next);
	}
	
	/* Every processor needs to have the total number of points n */
	//MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
			
	/* Compute local indices for data distribution */
	//I = (int) ((n + size + rank - 1)/ size); 
	
	/* For mapping local to global indices */
	// int L = n / size;
	// int R = n % size;
	// int mu_prefix = rank * L + MIN(rank, R);
	// int mu;
	
	// PLocal = (Point *) malloc(I * sizeof(Point));
	
	//// MPI_Scatterv(P, I, displs, MPI_, recvbuf, recvcount, recvtype, 0, MPI_COMM_WORLD);
	
	MPI_Type_free(&mpi_point_type);
	MPI_Finalize();
	exit(0);
}
