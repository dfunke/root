/*
 * delaunayTriangulation.cxx
 *
 *  Created on: Jun 30, 2014
 *      Author: dfunke
 *
 *  This test creates a TGraph2D, fills it with 3 points and then performs
 *  the Delaunay triangulation of them.
 *
 *  Because of bug ROOT-XXX the data structures were not properly initialized and no triangle was found
 */

#include "TGraph2D.h"
#include "TGraphDelaunay.h"

int main(int argc, char **argv) {

	const int EXP = 2;

	TGraph2D graph;
	graph.SetPoint(0, 0.1, 0.2, 0);
	graph.SetPoint(1, 0.5, 0.4, 0);
	graph.SetPoint(2, 0.3, 0.3, 0);
	graph.SetPoint(3, 0.6, 0.1, 0);


	TGraphDelaunay delaunay(&graph);

	delaunay.FindAllTriangles();

	if(delaunay.GetNdt() == EXP)
		return 0;
	else {
		printf("Expected: %i\t Gotten: %i\n", EXP, delaunay.GetNdt());
		return 4;
	}

}

