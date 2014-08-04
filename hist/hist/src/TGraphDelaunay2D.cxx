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
   Double_t zz = InterpolateNormalized(xx, yy);

   // Wrong zeros may appear when points sit on a regular grid.
   // The following line try to avoid this problem.
   if (zz==0) zz = InterpolateNormalized(xx+0.0001, yy);

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

   for (Int_t n=0; n<fNpoints; n++) {
	   //Point p(xTransformer(fX[n]), yTransformer(fY[n]));
	   Point p(linear_transform(fX[n], fOffsetX, fScaleFactorX),
			   linear_transform(fY[n], fOffsetY, fScaleFactorY));

	   fNormalizedPoints.insert(std::make_pair(p, n));
   }

   fCGALdelaunay.insert(fNormalizedPoints.begin(), fNormalizedPoints.end());

   fNdt = fCGALdelaunay.number_of_faces();

   std::transform(fCGALdelaunay.finite_faces_begin(), fCGALdelaunay.finite_faces_end(), std::back_inserter(fTriangles),
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

//______________________________________________________________________________
Double_t TGraphDelaunay2D::InterpolateNormalized(Double_t xx, Double_t yy)
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

	std::cout << std::this_thread::get_id() << ": Found " << coords.size() << " points" << std::endl;

	if(!nn.third) //neighbor finding was NOT successfull, return standard value
		return fZout;

	//printf("found neighbors %u\n", coords.size());

	Coord_type res = CGAL::linear_interpolation(coords.begin(), coords.end(),
			nn.second, Value_access(fNormalizedPoints, fZ));

	std::cout << std::this_thread::get_id() << ": Result " << res << std::endl;

   return res;
}

//______________________________________________________________________________
void TGraphDelaunay2D::SetMarginBinsContent(Double_t z)
{
   // Sets the histogram bin height for points lying outside the convex hull ie:
   // the bins in the margin.

   fZout = z;
}
