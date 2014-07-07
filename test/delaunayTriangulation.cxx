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
#include "TGraphDelaunay2D.h"

#include "delaunayTriangulation_bug.h"

void printDelaunay(const TGraphDelaunay2D & gd){

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
	const bool VERBOSE = false;

	TGraph2D * graph = getGraph();

	auto h = graph->GetHistogram("");

	TGraphDelaunay2D delaunay(graph);

	for (int i = 0; i < 100; i++) {
		Double_t pt = 50 + 0.001;
		Double_t eta = 1. + 0.01 * i + 0.001;
		Double_t res = graph->Interpolate(pt, eta);
		Double_t res2 = delaunay.Eval(pt, eta);
		std::cout << eta << " " << res << "  " << res2 << std::endl;
	}

	if(delaunay.GetNdt() == EXP){
		if(VERBOSE) printDelaunay(delaunay);

	return 0;
	} else {
		printf("Expected: %i\t Gotten: %i\n", EXP, delaunay.GetNdt());
		if(VERBOSE) printDelaunay(delaunay);

		return 4;
	}

}

