// @(#)root/sql:$Id$
// Author: Sergey Linev  20/11/2005

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//________________________________________________________________________
//
// TSQLClassInfo class containes info about tables specific to one class and
// version. It provides names of table for that class. For each version of
// class not more than two tables can exists. Normal table has typically
// name like TH1_ver4 and additional table has name like TH1_raw4
// List of this objects are kept by TSQLFile class
//
//________________________________________________________________________

#include "TSQLClassInfo.h"

#include "TObjArray.h"


ClassImp(TSQLClassColumnInfo)

////////////////////////////////////////////////////////////////////////////////
/// default constructor

TSQLClassColumnInfo::TSQLClassColumnInfo() :
   TObject(),
   fName(),
   fSQLName(),
   fSQLType()
{
}

////////////////////////////////////////////////////////////////////////////////
/// normal constructor

TSQLClassColumnInfo::TSQLClassColumnInfo(const char* name,
                                         const char* sqlname,
                                         const char* sqltype) :
   TObject(),
   fName(name),
   fSQLName(sqlname),
   fSQLType(sqltype)
{
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TSQLClassColumnInfo::~TSQLClassColumnInfo()
{
}


ClassImp(TSQLClassInfo)

////////////////////////////////////////////////////////////////////////////////
/// default constructor

TSQLClassInfo::TSQLClassInfo() :
   TObject(),
   fClassName(),
   fClassVersion(0),
   fClassId(0),
   fClassTable(),
   fRawTable(),
   fColumns(0),
   fRawtableExist(kFALSE)
{
}

////////////////////////////////////////////////////////////////////////////////
/// normal constructor of TSQLClassInfo class
/// Sets names of tables, which are used for that version of class

TSQLClassInfo::TSQLClassInfo(Long64_t classid,
                             const char* classname,
                             Int_t version) :
   TObject(),
   fClassName(classname),
   fClassVersion(version),
   fClassId(classid),
   fClassTable(),
   fRawTable(),
   fColumns(0),
   fRawtableExist(kFALSE)
{
   fClassTable.Form("%s_ver%d", classname, version);
   fRawTable.Form("%s_raw%d", classname, version);
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TSQLClassInfo::~TSQLClassInfo()
{
   if (fColumns!=0) {
      fColumns->Delete();
      delete fColumns;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// assigns new list of columns

void TSQLClassInfo::SetColumns(TObjArray* columns)
{
   if (fColumns!=0) {
      fColumns->Delete();
      delete fColumns;
   }
   fColumns = columns;
}

////////////////////////////////////////////////////////////////////////////////
/// set current status of class tables

void TSQLClassInfo::SetTableStatus(TObjArray* columns, Bool_t israwtable)
{
   SetColumns(columns);
   fRawtableExist = israwtable;
}

////////////////////////////////////////////////////////////////////////////////
/// Search for column of that name
/// Can search either for full column name (sqlname = kFALSE, default)
/// or for name, used as column name (sqlname = kTRUE)
/// Return index of column in list (-1 if not found)

Int_t TSQLClassInfo::FindColumn(const char* name, Bool_t sqlname)
{
   if ((name==0) || (fColumns==0)) return -1;

   TIter next(fColumns);

   TSQLClassColumnInfo* col = 0;

   Int_t indx = 0;

   while ((col = (TSQLClassColumnInfo*) next()) != 0) {
      const char* colname = sqlname ? col->GetSQLName() : col->GetName();
      if (strcmp(colname, name)==0) return indx;
      indx++;
   }

   return -1;
}
