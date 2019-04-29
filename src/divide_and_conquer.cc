/* 
 * Implementation of a O(n log n) divide-and-conquer algorithm to find the convex hull in 3D, 
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


#include <iostream>
#include <string>
#include <fstream>

using namespace std;

struct Point {
	double x, y, z;
	Point *prev, *next;
	
	/*
	* If the present point is still in the linked list, it gets dropped out (not pointed at anymore),
	* otherwise, it is inserted into the list, in between its own prev and next pointers. 
	*/
	void act() {  
		// insert
		if (prev->next != this){ 
			prev->next = this;
			next->prev = this;  
		}
		// delete
		else { 
			prev->next = next; 
			next->prev = prev; 
		}  
	}
};

const double INF = 1e99;
static Point nil = {INF, INF, INF, 0, 0};
Point *NIL = &nil;

/*
* Function to determine if any of a triplet of point pointers is null. 
*/
inline bool hasnil(Point *p, Point *q, Point *r){
	return (p == NIL || q == NIL || r == NIL);
}

/*
* Function to determine if points p-q-r form a clockwise/left turn, 
* or a counterclockwise/right turn. Returns < 0 and > 0 respectively. 
*/
inline double turn(Point *p, Point *q, Point *r) {  
	return (q->x - p->x) * (r->y - p->y) - 
			(r->x - p->x) * (q->y - p->y);
}

/*
* Function to determine the time when p-q-r will change between a left/right turn (ie. become colinear)
*/
inline double time(Point *p, Point *q, Point *r) { 
	if (p == NIL || q == NIL || r == NIL){
		return INF;
	}
	return ((q->x - p->x) * (r->z - p->z) - 
			(r->x - p->x) * (q->z - p->z)) / turn(p, q, r);
}

/*
* Function to provide a string representation of a point. 
*/
inline string toString(Point *p){
	std::string x = "{" + std::to_string(p->x) + ", " + std::to_string(p->y) + ", " + std::to_string(p->z) + "}";
	return x;
}

// Mergesort
Point *sort(Point P[], int n) {  

	Point *a, *b, *c, head;

	if (n == 1){
		P[0].next = NIL; 
		return P;
	}
	a = sort(P, n/2);
	b = sort(P+n/2, n-n/2);
	c = &head;
	do
		if (a->x < b->x){
			c = c->next = a; 
			a = a->next; 
		}
		else { 
			c = c->next = b; 
			b = b->next; 
		}
	while (c != NIL);
	return head.next;
}

// Main algorithm
void hull(bool bottom, Point *list, int n, Point **A, Point **B) {  

	Point *u, *v, *mid;  
	double t[6], oldt, newt;  
	int i, j, k, l, minl;

	// Base case is a single point
	if (n == 1) { 
		// Remove the point from the list
		list->next = NIL; 
		list->prev = NIL;
		// Return no faces for the hull
		A[0] = NIL;
		return; 
	}
	
	// Find the middle of the current problem set
	for (u = list, i = 0; i < n/2-1; u = u->next, i++) ;
	mid = v = u->next;
	
	// Recurse on each side of the middle, swapping the A and B arrays
	hull(bottom, list, n/2, B, A);  				// build left side
	hull(bottom, mid, n-n/2, B+n/2*2, A+n/2*2);		// build right side

	// ##### Merging step #####
  
	// Find initial bridge
	// Start in the middle and move the bridge vertices out to each side.
	for ( ; ; ){  
		// Lower hull
		if(bottom){
			if (!hasnil(u, v, v->next) && turn(u, v, v->next) < 0) v = v->next;
			else if (!hasnil(u->prev, u, v) && turn(u->prev, u, v) < 0) u = u->prev;  
			else break;
		}
		// Upper hull
		else{
			if (!hasnil(u, v, v->next) && turn(u, v, v->next) > 0) v = v->next;
			else if (!hasnil(u->prev, u, v) && turn(u->prev, u, v) > 0) u = u->prev;  
			else break;
		}
	}

	// Merge by tracking bridge uv over time
	for (i = k = 0, j = n/2*2, oldt = -INF; ; oldt = newt) {  
		t[0] = time(B[i]->prev, B[i], B[i]->next);  
		t[1] = time(B[j]->prev, B[j], B[j]->next);    
		t[2] = time(u, u->next, v);  
		t[3] = time(u->prev, u, v);
		t[4] = time(u, v->prev, v); 
		t[5] = time(u, v, v->next);
		
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
				if (B[i]->x < u->x){
					A[k++] = B[i]; 
				}
				B[i++]->act();  
				break;
			// Right side
			case 1:  
				if (B[j]->x > v->x){
					A[k++] = B[j];  
				}
				B[j++]->act();  
				break;
			// Bridge
			case 2:  
				u = u->next;
				A[k++] = u;  
				break;
			// Bridge
			case 3:  
				A[k++] = u;
				u = u->prev; 
				break;
			// Bridge
			case 4:  
				v = v->prev;
				A[k++] = v;  
				break;
			// Bridge
			case 5:  
				A[k++] = v;  
				v = v->next;
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
	u->next = v;  
	v->prev = u;  
	// Loop from end of A[k] to start
	for (k--; k >= 0; k--) 
		 // On the left or right hull
		if (A[k]->x <= u->x || A[k]->x >= v->x) {
			A[k]->act();
			if (A[k] == u) 
				u = u->prev; 
			else if (A[k] == v) 
				v = v->next;
		}
		// Inside the bridge (so it was a bridge endpoint)
		else { 
			// Put the point between the current two bridge endpoints (u and v)
			u->next = A[k]; 
			A[k]->prev = u; 
			v->prev = A[k]; 
			A[k]->next = v;
			
			// Make it a bridge endpoint
			// Left side
			if (A[k]->x < mid->x){
				u = A[k]; 
			}
			// Right side
			else{
				v = A[k];
			}
		}
}

main() {

	int n, i;  
	char left_bracket, comma_1, comma_2, right_bracket;
  
	// Input
	ifstream infile ("points.in");
  
	infile >> n;
  
	Point *P = new Point[n];  
  
	for (i = 0; i < n; i++){
		infile >> left_bracket >> P[i].x >> comma_1 >> P[i].y >> comma_2 >> P[i].z >> right_bracket;
	}
  
	infile.close();  
  
	ofstream outfile ("output_DC.out");

	Point *list;
	Point **A = new Point *[2*n], 
			**B = new Point *[2*n];
			
	// Compute the lower hull 
	list = sort(P, n);
	hull(true, list, n, A, B);
	  
	// Output for lower hull facets
	for (i = 0; A[i] != NIL; A[i++]->act()){
		outfile << "{" << A[i]->prev-P << ", " << A[i]-P << ", " << A[i]->next-P << "}\n";
	}
	
	// Compute the upper hull
	list = sort(P, n);
	hull(false, list, n, A, B);
	
	// Output for upper hull facets
	for (i = 0; A[i] != NIL; A[i++]->act()){ 
		outfile << "{" << A[i]->prev-P << ", " << A[i]-P << ", " << A[i]->next-P << "}\n";
	}
  
	delete A;  
	delete B;  
	delete P;
  
	outfile.close();
}
