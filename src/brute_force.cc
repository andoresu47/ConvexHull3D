/* Brute-force algorithm for 3D convex hull. 
 * O(n^4) time complexity.
 */ 

#include<stdlib.h>
#include <string>
#include <fstream>
#include <iostream>
//#include <mpi.h>

using namespace std;

struct Point { double x, y, z; };

inline double faceOrient(Point p, Point q, Point r, Point s){
	double uvi, uvj, uvk;
	
	uvi = (q.y-p.y) * (r.z-p.z) - (q.z-p.z) * (r.y-p.y);
    uvj = (q.z-p.z) * (r.x-p.x) - (q.x-p.x) * (r.z-p.z);
    uvk = (q.x-p.x) * (r.y-p.y) - (q.y-p.y) * (r.x-p.x);
	
	return (s.x - p.x) * uvi + (s.y - p.y) * uvj + (s.z - p.z) * uvk;
}

inline string toString(Point p){
	std::string x = "{" + std::to_string(p.x) + ", " + std::to_string(p.y) + ", " + std::to_string(p.z) + "}";
	return x;
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

main() {

  int n, i, j, k, l;
  std::cin >> n;
  
  ofstream outfile ("output.out");

  Point *P = new Point[n]; 
  Point centroid = Point();
  for (i = 0; i < n; i++) { 
	std::cin >> P[i].x; 
	std::cin >> P[i].y; 
	std::cin >> P[i].z; 
	
	centroid.x += P[i].x;
	centroid.y += P[i].y;
	centroid.z += P[i].z;
  }
  
  centroid.x = centroid.x / n;
  centroid.y = centroid.y / n;
  centroid.z = centroid.z / n;
  
  //for (i = 0; i < n; i++) { std::cout << P[i].x << " " << P[i].y << " " << P[i].z << "\n"; }
  //std::cout << centroid.x << " " << centroid.y << " " << centroid.z << "\n";
  
  qsort (P, n, sizeof(Point), Comparator);
  
  //for (i = 0; i < n; i++) { std::cout << P[i].x << " " << P[i].y << " " << P[i].z << "\n"; }
    
  for (i = 0; i < n; i++){
    for (j = 0; j < n; j++){
      for (k = 0; k < n; k++) {  
	  if (i != j && j != k && k != i){
        if (faceOrient(P[i],P[j],P[k],centroid) >= 0){
			continue;
		}
        for (l = 0; l < n; l++){  // check facet ijk
          if (faceOrient(P[i],P[j],P[k],P[l]) > 0){			  
			  break;
		  }
		}
        if (l == n){
			//std::cout << i << " " << j << " " << k << "\n";
			outfile << toString(P[i]) << "\n";
			outfile << toString(P[j]) << "\n";
			outfile << toString(P[k]) << "\n";
		}
	  }
      }
	}
  }

  delete P;
  outfile.close();
}
