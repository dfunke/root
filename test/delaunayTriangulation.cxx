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

void printDelaunay(const TGraphDelaunay & gd){

	auto graph = gd.GetGraph2D();

	for(const auto & face : gd){
		printf("[%u](%f,%f) - [%u](%f,%f) - [%u](%f,%f)\n",
				face.vertex(0)->info(), graph->GetX()[face.vertex(0)->info()], graph->GetY()[face.vertex(0)->info()],
				face.vertex(1)->info(), graph->GetX()[face.vertex(1)->info()], graph->GetY()[face.vertex(1)->info()],
				face.vertex(2)->info(), graph->GetX()[face.vertex(2)->info()], graph->GetY()[face.vertex(2)->info()]);
	}

}

int main() {

	const int EXP = 2;
	const bool VERBOSE = true;

	TGraph2D graph;
	graph.SetPoint(0, 0.1, 0.2, 0);
	graph.SetPoint(1, 0.5, 0.4, 0);
	graph.SetPoint(2, 0.3, 0.2, 0);
	graph.SetPoint(3, 0.6, 0.1, 0);


	TGraphDelaunay delaunay(&graph);

	delaunay.FindAllTriangles();

	if(delaunay.GetNdt() == EXP){
		if(VERBOSE) printDelaunay(delaunay);

	return 0;
	} else {
		printf("Expected: %i\t Gotten: %i\n", EXP, delaunay.GetNdt());
		if(VERBOSE) printDelaunay(delaunay);

		return 4;
	}

}

