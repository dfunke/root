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
#include <functional>

class TGraph2D;
class TView;

//CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Interpolation_traits_2.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/interpolation_functions.h>

class TGraphDelaunay2D : public TNamed {

public:

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

	//typedef std::function<Double_t(Double_t)>                  Transformer;

private:
   TGraphDelaunay2D(const TGraphDelaunay2D&); // Not implemented
   TGraphDelaunay2D& operator=(const TGraphDelaunay2D&); // Not implemented

protected:

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

   Bool_t      fInit;        //!True if CreateTrianglesDataStructure() and FindHull() have been performed

   TGraph2D   *fGraph2D;     //!2D graph containing the user data

   Delaunay fCGALdelaunay; //! CGAL delaunay triangulation object
   PointWithInfoMap fNormalizedPoints; //! Normalized function values

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

   Delaunay::Finite_faces_iterator begin() const { return fCGALdelaunay.finite_faces_begin(); }
   Delaunay::Finite_faces_iterator end()  const { return fCGALdelaunay.finite_faces_end(); }

   ClassDef(TGraphDelaunay2D,1)  // Delaunay triangulation

private:

   inline double linear_transform(double x, double offset, double factor){
	   return (x+offset)*factor;
   }

   Double_t  InterpolateNormalized(Double_t x, Double_t y);

};

#endif
