// @(#)root/hist:$Id: TGraphDelaunay2D.h,v 1.00
// Author: Olivier Couet, Luke Jones (Royal Holloway, University of London)

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGraphDelaunay2D
#define ROOT_TGraphDelaunay2D

//for testing purposes HAS_CGAL can be [un]defined here
//#define HAS_CGAL

//for testing purposes THREAD_SAFE can [un]defined here
//#define THREAD_SAFE

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGraphDelaunay2D                                                     //
//                                                                      //
// This class uses the Delaunay triangles technique to interpolate and  //
// render the data set.                                                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TNamed
#include "TNamed.h"
#endif

#include <map>
#include <vector>
#include <set>
#include <functional>

#ifdef HAS_CGAL
	/* CGAL uses the name PTR as member name in its Handle class
	 * but its a macro defined in mmalloc.h of ROOT
	 * Safe it, disable it and then re-enable it later on*/
	#pragma push_macro("PTR")
	#undef PTR

   #include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
   #include <CGAL/Delaunay_triangulation_2.h>
   #include <CGAL/Triangulation_vertex_base_with_info_2.h>
   #include <CGAL/Interpolation_traits_2.h>
   #include <CGAL/natural_neighbor_coordinates_2.h>
   #include <CGAL/interpolation_functions.h>

	#pragma pop_macro("PTR")
#else
	// fallback to triangle library
	#include "triangle.h"
#endif

#ifdef THREAD_SAFE
	#include<atomic> //atomic operations for thread safety
#endif

class TGraph2D;
class TView;

class TGraphDelaunay2D : public TNamed {

public:

	struct Triangle {
		Double_t x[3];
		Double_t y[3];
		UInt_t idx[3];

		#ifndef HAS_CGAL
		Double_t invDenom; //see comment below in CGAL fall back section
		#endif
	};

	typedef std::vector<Triangle> Triangles;

private:
   TGraphDelaunay2D(const TGraphDelaunay2D&); // Not implemented
   TGraphDelaunay2D& operator=(const TGraphDelaunay2D&); // Not implemented

protected:

   //typedef std::function<Double_t(Double_t)>                  Transformer;

   Int_t       fNdt;         //!Number of Delaunay triangles found
   Int_t       fNpoints;     //!Number of data points in fGraph2D

   Double_t   *fX;           //!Pointer to fGraph2D->fX
   Double_t   *fY;           //!Pointer to fGraph2D->fY
   Double_t   *fZ;           //!Pointer to fGraph2D->fZ

   Double_t    fXNmin;       //!Minimum value of fXN
   Double_t    fXNmax;       //!Maximum value of fXN
   Double_t    fYNmin;       //!Minimum value of fYN
   Double_t    fYNmax;       //!Maximum value of fYN

   //Transformer xTransformer; //!transform x values to mapped space
   //Transformer yTransformer; //!transform y values to mapped space

   Double_t    fOffsetX;      //!Normalization offset X
   Double_t    fOffsetY;      //!Normalization offset Y

   Double_t    fScaleFactorX; //!Normalization factor X
   Double_t    fScaleFactorY; //!Normalization factor Y

   Double_t    fZout;        //!Histogram bin height for points lying outside the convex hull

#ifdef THREAD_SAFE

   enum class Initialization : char {UNINITIALIZED, INITIALIZING, INITIALIZED};
   std::atomic<Initialization> fInit; //!Indicate initialization state

#else
   Bool_t      fInit;        //!True if FindAllTriangels() has been performed
#endif

   TGraph2D   *fGraph2D;     //!2D graph containing the user data

   Triangles   fTriangles;   //!Triangles of Triangulation

#ifdef HAS_CGAL

   //Functor class for accessing the function values/gradients
   	template< class PointWithInfoMap, typename ValueType >
   	struct Data_access : public std::unary_function< typename PointWithInfoMap::key_type,
   				 std::pair<ValueType, bool> >
   	{

   	  Data_access(const PointWithInfoMap& points, const ValueType * values)
   			  : _points(points), _values(values){};

   	  std::pair< ValueType, bool>
   	  operator()(const typename PointWithInfoMap::key_type& p) const {
   		typename PointWithInfoMap::const_iterator mit = _points.find(p);
   		if(mit!= _points.end())
   		  return std::make_pair(_values[mit->second], true);
   		return std::make_pair(ValueType(), false);
   	  };

   	  const PointWithInfoMap& _points;
   	  const ValueType * _values;
   	};

   	typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
   	typedef CGAL::Triangulation_vertex_base_with_info_2<uint, K> Vb;
   	typedef CGAL::Triangulation_data_structure_2<Vb>             Tds;
   	typedef CGAL::Delaunay_triangulation_2<K, Tds>               Delaunay;
   	typedef CGAL::Interpolation_traits_2<K>                      Traits;
   	typedef K::FT                                                Coord_type;
   	typedef K::Point_2                                           Point;
   	typedef std::map<Point, Vb::Info, K::Less_xy_2>              PointWithInfoMap;
   	typedef Data_access< PointWithInfoMap, Double_t >            Value_access;

   Delaunay fCGALdelaunay; //! CGAL delaunay triangulation object
   PointWithInfoMap fNormalizedPoints; //! Normalized function values

#else // HAS_CGAL
   //fallback to triangle library

   /* Using barycentric coordinates for inTriangle test and interpolation
    *
    * Given triangle ABC and point P, P can be expressed by
    *
    * P.x = la * A.x + lb * B.x + lc * C.x
    * P.y = la * A.y + lb * B.y + lc * C.y
    *
    * with lc = 1 - la - lb
    *
    * P.x = la * A.x + lb * B.x + (1-la-lb) * C.x
    * P.y = la * A.y + lb * B.y + (1-la-lb) * C.y
    *
    * Rearranging yields
    *
    * la * (A.x - C.x) + lb * (B.x - C.x) = P.x - C.x
    * la * (A.y - C.y) + lb * (B.y - C.y) = P.y - C.y
    *
    * Thus
    *
    * la = ( (B.y - C.y)*(P.x - C.x) + (C.x - B.x)*(P.y - C.y) ) / ( (B.y - C.y)*(A.x - C.x) + (C.x - B.x)*(A.y - C.y) )
    * lb = ( (C.y - A.y)*(P.x - C.x) + (A.x - C.x)*(P.y - C.y) ) / ( (B.y - C.y)*(A.x - C.x) + (C.x - B.x)*(A.y - C.y) )
    * lc = 1 - la - lb
    *
    * We save the inverse denominator to speedup computation
    *
    * invDenom = 1 / ( (B.y - C.y)*(A.x - C.x) + (C.x - B.x)*(A.y - C.y) )
    *
    * P is in triangle (including edges if
    *
    * 0 <= [la, lb, lc] <= 1
    *
    * The interpolation of P.z is
    *
    * P.z = la * A.z + lb * B.z + lc * C.z
    *
    */

   std::vector<Double_t> fXN; //! normalized X
   std::vector<Double_t> fYN; //! normalized Y

   /* To speed up localisation of points a grid is layed over normalized space
    *
    * A reference to triangle ABC is added to _all_ grid cells that include ABC's bounding box
    */

   static const UInt_t fNCells = 25; //! number of cells to divide the normalized space
   Double_t fXCellStep; //! inverse denominator to calculate X cell = fNCells / (fXNmax - fXNmin)
   Double_t fYCellStep; //! inverse denominator to calculate X cell = fNCells / (fYNmax - fYNmin)
   std::set<UInt_t> fCells[(fNCells+1)*(fNCells+1)]; //! grid cells with containing triangles

   inline uint cell(uint x, uint y) const {
	   return x*(fNCells+1) + y;
   }

   inline int cellX(double x) const {
	   return (x - fXNmin) * fXCellStep;
   }

   inline int cellY(double y) const {
	   return (y - fYNmin) * fYCellStep;
   }

#endif //HAS_CGAL

public:

   TGraphDelaunay2D();
   TGraphDelaunay2D(TGraph2D *g);

   Double_t  Eval(Double_t x, Double_t y);
   void      FindAllTriangles();

   TGraph2D *GetGraph2D() const {return fGraph2D;}
   Double_t  GetMarginBinsContent() const {return fZout;}
   Int_t     GetNdt() const {return fNdt;}
   Double_t  GetXNmin() const {return fXNmin;}
   Double_t  GetXNmax() const {return fXNmax;}
   Double_t  GetYNmin() const {return fYNmin;}
   Double_t  GetYNmax() const {return fYNmax;}

   void      SetMarginBinsContent(Double_t z=0.);

   Triangles::const_iterator begin() const { return fTriangles.begin(); }
   Triangles::const_iterator end()  const { return fTriangles.end(); }

   ClassDef(TGraphDelaunay2D,1)  // Delaunay triangulation

private:

   inline double linear_transform(double x, double offset, double factor){
	   return (x+offset)*factor;
   }

   void _normalizePoints();
   void _findTriangles();
   Double_t  _interpolateNormalized(Double_t x, Double_t y);

};

#endif
