// @(#)root/hist:$Id: TGraphDelaunay2D.cxx,v 1.00
// Author: Olivier Couet, Luke Jones (Royal Holloway, University of London)

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TMath.h"
#include "TGraph2D.h"
#include "TGraphDelaunay2D.h"

#include <thread>

ClassImp(TGraphDelaunay2D)


//______________________________________________________________________________
//
// TGraphDelaunay2D generates a Delaunay triangulation of a TGraph2D. This
// triangulation code derives from an implementation done by Luke Jones
// (Royal Holloway, University of London) in April 2002 in the PAW context.
//
// This software cannot be guaranteed to work under all circumstances. They
// were originally written to work with a few hundred points in an XY space
// with similar X and Y ranges.
//
// Definition of Delaunay triangulation (After B. Delaunay):
// For a set S of points in the Euclidean plane, the unique triangulation DT(S)
// of S such that no point in S is inside the circumcircle of any triangle in
// DT(S). DT(S) is the dual of the Voronoi diagram of S. If n is the number of
// points in S, the Voronoi diagram of S is the partitioning of the plane
// containing S points into n convex polygons such that each polygon contains
// exactly one point and every point in a given polygon is closer to its
// central point than to any other. A Voronoi diagram is sometimes also known
// as a Dirichlet tessellation.
//Begin_Html
/*
<img src="gif/dtvd.gif">
<br>
<a href="http://www.cs.cornell.edu/Info/People/chew/Delaunay.html">This applet</a>
gives a nice practical view of Delaunay triangulation and Voronoi diagram.
*/
//End_Html


//______________________________________________________________________________
TGraphDelaunay2D::TGraphDelaunay2D()
            : TNamed("TGraphDelaunay2D","TGraphDelaunay2D")
{
   // TGraphDelaunay2D default constructor

   fGraph2D      = 0;
   fX            = 0;
   fY            = 0;
   fZ            = 0;
   fNpoints      = 0;
   fOffsetX      = 0;
   fOffsetY      = 0;
   fScaleFactorX = 0;
   fScaleFactorY = 0;
   fZout         = 0.;
   fNdt          = 0;
   fInit         = kFALSE;
   fXNmin        = 0.;
   fXNmax        = 0.;
   fYNmin        = 0.;
   fYNmax        = 0.;
}


//______________________________________________________________________________
TGraphDelaunay2D::TGraphDelaunay2D(TGraph2D *g)
            : TNamed("TGraphDelaunay2D","TGraphDelaunay2D")
{
   // TGraphDelaunay2D normal constructor

   fGraph2D      = g;
   fX            = fGraph2D->GetX();
   fY            = fGraph2D->GetY();
   fZ            = fGraph2D->GetZ();
   fNpoints      = fGraph2D->GetN();
   fOffsetX      = 0;
   fOffsetY      = 0;
   fScaleFactorX = 0;
   fScaleFactorY = 0;
   fZout         = 0.;
   fNdt          = 0;
   fInit         = kFALSE;
   fXNmin        = 0.;
   fXNmax        = 0.;
   fYNmin        = 0.;
   fYNmax        = 0.;
}

//______________________________________________________________________________
Double_t TGraphDelaunay2D::Eval(Double_t x, Double_t y)
{
   // Return the z value corresponding to the (x,y) point in fGraph2D

   // Initialise the Delaunay algorithm if needed.
   // CreateTrianglesDataStructure computes fXoffset, fYoffset,
   // fXScaleFactor and fYScaleFactor;
   // needed in this function.
   FindAllTriangles();

   // Find the z value corresponding to the point (x,y).
   Double_t xx, yy;
   xx = linear_transform(x, fOffsetX, fScaleFactorX); //xx = xTransformer(x);
   yy = linear_transform(y, fOffsetY, fScaleFactorY); //yy = yTransformer(y);
   Double_t zz = _interpolateNormalized(xx, yy);

   // Wrong zeros may appear when points sit on a regular grid.
   // The following line try to avoid this problem.
   if (zz==0) zz = _interpolateNormalized(xx+0.0001, yy);

   return zz;
}

//______________________________________________________________________________
void TGraphDelaunay2D::FindAllTriangles()
{

   if (fInit) return; else fInit = kTRUE;

   // Function used internally only. It creates the data structures needed to
   // compute the Delaunay triangles.

   // Offset fX and fY so they average zero, and scale so the average
   // of the X and Y ranges is one. The normalized version of fX and fY used
   // in Interpolate.
   Double_t xmax = fGraph2D->GetXmax();
   Double_t ymax = fGraph2D->GetYmax();
   Double_t xmin = fGraph2D->GetXmin();
   Double_t ymin = fGraph2D->GetYmin();

   fOffsetX      = -(xmax+xmin)/2.;
   fOffsetY      = -(ymax+ymin)/2.;

   fScaleFactorX = 1./(xmax-xmin);
   fScaleFactorY = 1./(ymax-ymin);

   //xTransformer = std::bind(linear_transform, std::placeholders::_1, Xoffset, XScaleFactor);
   //yTransformer = std::bind(linear_transform, std::placeholders::_1, Yoffset, YScaleFactor);

   fXNmax        = linear_transform(xmax, fOffsetX, fScaleFactorX); //xTransformer(xmax);
   fXNmin        = linear_transform(xmin, fOffsetX, fScaleFactorX); //xTransformer(xmin);

   fYNmax        = linear_transform(ymax, fOffsetY, fScaleFactorY); //yTransformer(ymax);
   fYNmin        = linear_transform(ymin, fOffsetY, fScaleFactorY); //yTransformer(ymin);

   _normalizePoints(); // call backend specific point normalization

   _findTriangles(); // call backend specific triangle finding

   fNdt = fTriangles.size();

}

//______________________________________________________________________________
void TGraphDelaunay2D::SetMarginBinsContent(Double_t z)
{
   // Sets the histogram bin height for points lying outside the convex hull ie:
   // the bins in the margin.

   fZout = z;
}

//______________________________________________________________________________
// backend specific implementations
#ifdef HAS_CGAL

void TGraphDelaunay2D::_normalizePoints() {
	for (Int_t n = 0; n < fNpoints; n++) {
		//Point p(xTransformer(fX[n]), yTransformer(fY[n]));
		Point p(linear_transform(fX[n], fOffsetX, fScaleFactorX),
				linear_transform(fY[n], fOffsetY, fScaleFactorY));

		fNormalizedPoints.insert(std::make_pair(p, n));
	}
}

void TGraphDelaunay2D::_findTriangles() {
	fCGALdelaunay.insert(fNormalizedPoints.begin(), fNormalizedPoints.end());

	std::transform(fCGALdelaunay.finite_faces_begin(),
			fCGALdelaunay.finite_faces_end(), std::back_inserter(fTriangles),
			[] (const Delaunay::Face face) -> Triangle {

				Triangle tri;

				auto transform = [&] (const uint i) {
					tri.x[i] = face.vertex(i)->point().x();
					tri.y[i] = face.vertex(i)->point().y();
					tri.idx[i] = face.vertex(i)->info();
				};

				transform(0);
				transform(1);
				transform(2);

				return tri;

			});
}

Double_t TGraphDelaunay2D::_interpolateNormalized(Double_t xx, Double_t yy)
{
   // Finds the Delaunay triangle that the point (xi,yi) sits in (if any) and
   // calculate a z-value for it by linearly interpolating the z-values that
   // make up that triangle.

   // initialise the Delaunay algorithm if needed
    FindAllTriangles();

   //coordinate computation
	Point p(xx, yy);

	std::vector<std::pair<Point, Coord_type> > coords;
	auto nn = CGAL::natural_neighbor_coordinates_2(fCGALdelaunay, p,
			std::back_inserter(coords));

	//std::cout << std::this_thread::get_id() << ": Found " << coords.size() << " points" << std::endl;

	if(!nn.third) //neighbor finding was NOT successfull, return standard value
		return fZout;

	//printf("found neighbors %u\n", coords.size());

	Coord_type res = CGAL::linear_interpolation(coords.begin(), coords.end(),
			nn.second, Value_access(fNormalizedPoints, fZ));

	//std::cout << std::this_thread::get_id() << ": Result " << res << std::endl;

   return res;
}

#else

// fallback to triangle library
#include "triangle.h"

void TGraphDelaunay2D::_normalizePoints() {
	for (Int_t n = 0; n < fNpoints; n++) {
		fXN.push_back(linear_transform(fX[n], fOffsetX, fScaleFactorX));
		fYN.push_back(linear_transform(fY[n], fOffsetY, fScaleFactorY));
	}
}

void TGraphDelaunay2D::_findTriangles() {

	auto initStruct = [] (triangulateio & s) {
							  s.pointlist = nullptr;                                               /* In / out */
							  s.pointattributelist = nullptr;                                      /* In / out */
							  s.pointmarkerlist = nullptr;                                          /* In / out */
							  s.numberofpoints = 0;                                            /* In / out */
							  s.numberofpointattributes = 0;                                   /* In / out */

							  s.trianglelist = nullptr;                                             /* In / out */
							  s.triangleattributelist = nullptr;                                   /* In / out */
							  s.trianglearealist = nullptr;                                         /* In only */
							  s.neighborlist = nullptr;                                             /* Out only */
							  s.numberoftriangles = 0;                                         /* In / out */
							  s.numberofcorners = 0;                                           /* In / out */
							  s.numberoftriangleattributes = 0;                                /* In / out */

							  s.segmentlist = nullptr;                                              /* In / out */
							  s.segmentmarkerlist = nullptr;                                        /* In / out */
							  s.numberofsegments = 0;                                          /* In / out */

							  s.holelist = nullptr;                        /* In / pointer to array copied out */
							  s.numberofholes = 0;                                      /* In / copied out */

							  s.regionlist = nullptr;                      /* In / pointer to array copied out */
							  s.numberofregions = 0;                                    /* In / copied out */

							  s.edgelist = nullptr;                                                 /* Out only */
							  s.edgemarkerlist = nullptr;            /* Not used with Voronoi diagram; out only */
							  s.normlist = nullptr;                /* Used only with Voronoi diagram; out only */
							  s.numberofedges = 0;                                             /* Out only */
						};

	auto freeStruct = [] (triangulateio & s) {
							  free(s.pointlist);                                               /* In / out */
							  free(s.pointattributelist);                                      /* In / out */
							  free(s.pointmarkerlist);                                          /* In / out */

							  free(s.trianglelist);                                             /* In / out */
							  free(s.triangleattributelist);                                   /* In / out */
							  free(s.trianglearealist);                                         /* In only */
							  free(s.neighborlist);                                             /* Out only */

							  free(s.segmentlist);                                              /* In / out */
							  free(s.segmentmarkerlist);                                        /* In / out */

							  free(s.holelist);                        /* In / pointer to array copied out */

							  free(s.regionlist);                      /* In / pointer to array copied out */

							  free(s.edgelist);                                                 /* Out only */
							  free(s.edgemarkerlist);            /* Not used with Voronoi diagram; out only */
							  free(s.normlist);                /* Used only with Voronoi diagram; out only */
						};

	struct triangulateio in, out;
	initStruct(in); initStruct(out);

	/* Define input points. */

	in.numberofpoints = fNpoints;
	in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));

	for (uint i = 0; i < fNpoints; ++i) {
		in.pointlist[2 * i] = fXN[i];
		in.pointlist[2 * i + 1] = fYN[i];
	}

	triangulate("zQN", &in, &out, nullptr);

	for(uint t = 0; t < out.numberoftriangles; ++t){
		Triangle tri;

		auto transform = [&] (const uint v) {
			//each triangle as numberofcorners vertices ( = 3)
			tri.idx[v] = out.trianglelist[t*out.numberofcorners + v];

			//printf("triangle %u vertex %u: point %u/%i\n", t, v, tri.idx[v], out.numberofpoints);

			//pointlist is [x0 y0 x1 y1 ...]
			tri.x[v] = in.pointlist[tri.idx[v] * 2 + 0];

			//printf("\t x: %f\n", tri.x[v]);

			tri.y[v] = in.pointlist[tri.idx[v] * 2 + 1];

			//printf("\t y: %f\n", tri.y[v]);
		};

		transform(0);
		transform(1);
		transform(2);

		//see comment in header for CGAL fallback section
		tri.invDenom = fabs(1 / ( (tri.y[1] - tri.y[2])*(tri.x[0] - tri.x[2]) + (tri.x[2] - tri.x[1])*(tri.y[0] - tri.y[2]) ));

		fTriangles.push_back(tri);

	}

	//sort triangles by their centroid
	std::sort(fTriangles.begin(), fTriangles.end(),
			[] (const Triangle & t1, const Triangle & t2) -> bool{

		//compute centroid of both triangles
		double t1_cx = 0; double t1_cy = 0;
		double t2_cx = 0; double t2_cy = 0;

		for(uint i = 0; i < 3; ++i){
			t1_cx += t1.x[i];
			t1_cy += t1.y[i];

			t2_cx += t2.x[i];
			t2_cy += t2.y[i];
		}

		//we do not need the division, the ordering is unchanged
		//t1_cx /= 3; t1_cy /= 3;
		//t2_cx /= 3; t2_cy /= 3;

		return t1_cx < t2_cx || (t1_cx == t2_cx && t1_cy < t2_cy);

	});

	freeStruct(in); freeStruct(out);
}

Double_t TGraphDelaunay2D::_interpolateNormalized(Double_t xx, Double_t yy)
{
   // Finds the Delaunay triangle that the point (xi,yi) sits in (if any) and
   // calculate a z-value for it by linearly interpolating the z-values that
   // make up that triangle.

   // initialise the Delaunay algorithm if needed
    FindAllTriangles();

    //see comment in header for TriangleSupplement structure
    auto bayCoords = [&] (const uint t) -> std::tuple<double, double, double> {
    	double la = fabs(( (fTriangles[t].y[1] - fTriangles[t].y[2])*(xx - fTriangles[t].x[2])
    					 + (fTriangles[t].x[2] - fTriangles[t].x[1])*(yy - fTriangles[t].y[2]) ) * fTriangles[t].invDenom);
    	double lb = fabs(( (fTriangles[t].y[2] - fTriangles[t].y[0])*(xx - fTriangles[t].x[2])
    			         + (fTriangles[t].x[0] - fTriangles[t].x[2])*(yy - fTriangles[t].y[2]) ) * fTriangles[t].invDenom);

    	return std::make_tuple(la, lb, (1 - la - lb));
    };

    auto inTriangle = [] (const std::tuple<double, double, double> & coords) -> bool {
    	return std::get<0>(coords) >= 0 && std::get<1>(coords) >= 0 && std::get<2>(coords) >= 0;
    };

    uint t = fNdt / 2;
    uint step = t;

    auto coords = bayCoords(t);

    while(!inTriangle(coords) && step > 0){
    	step /= 2;

    	double t1_cx = 0; double t1_cy = 0;
		for(uint i = 0; i < 3; ++i){
			t1_cx += fTriangles[t].x[i];
			t1_cy += fTriangles[t].y[i];
		}

		t1_cx /= 3; t1_cy /= 3;

		if(xx < t1_cx || (xx == t1_cx && yy < t1_cy))
			t -= step;
		else
			t += step;
    }

    //we found the triangle -> interpolate using the barycentric interpolation
    if(inTriangle(coords))
    	return std::get<0>(coords) * fZ[fTriangles[t].idx[0]]
    		 + std::get<1>(coords) * fZ[fTriangles[t].idx[1]]
    		 + std::get<2>(coords) * fZ[fTriangles[t].idx[2]];

    printf("Could not find a triangle for point (%f,%f)\n", xx, yy);

    //no triangle found return standard value
   return fZout;
}

#endif
