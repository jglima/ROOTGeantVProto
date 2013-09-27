#ifndef __CINT__
#include "Vc/vector.h"
#include "TPointStruct.h"
typedef Vc::double_v vd; // short for vector double
#endif

#include "TGeoBBox.h"



void TGeoBBox::DistFromOutsideSOA_v(StructOfCoord const & point, StructOfCoord const & dir, Int_t iact, 
			Double_t const * step, Double_t *safe, Double_t * distance , Int_t np) const 
{
  TGeoBBox::DistFromOutsideS_v(point, dir, this->TGeoBBox::GetDX(), this->TGeoBBox::GetDY(), this->TGeoBBox::GetDZ(), this->TGeoBBox::GetOrigin(),  step, distance, np );
}





// // this is the actual TGeoBBox::DistFromOutsideSOA_Vc doing the computation with possibility of early return
void TGeoBBox::DistFromOutsideS_v(StructOfCoord const & point, StructOfCoord const & dir, 
                          double dx, double dy, double dz, const double *origin, const double * stepmax, double * distance, int np)
{
#ifndef __CINT__
   vd zero=0.;
   vd vcorigin[3]={origin[0],origin[1],origin[2]};
   int tailsize =np % Vc::double_v::Size; 
   for(volatile unsigned int i = 0; i < np-tailsize; i+= Vc::double_v::Size) 
    {
      // this is a lot of memory movement ( maybe a cast would be better )?
      vd x(&point.x[i]);
      vd y(&point.y[i]);
      vd z(&point.z[i]);
      vd dirx(&dir.x[i]);
      vd diry(&dir.y[i]);
      vd dirz(&dir.z[i]);
      vd s(&stepmax[i]);
      vd dist(0.);
      TGeoBBox::DistFromOutsideSOA_Vc( x, y, z, dirx, diry, dirz, dx, dy, dz, vcorigin, s, dist );
      dist.store(&distance[i]);
   }
   // do the tail part for the moment, we just call the old static version
   for(unsigned int i = 0; i < tailsize; ++i)
     {
       double p[3]={point.x[np-tailsize+i], point.y[np-tailsize+i], point.z[np-tailsize+i]};
       double d[3]={dir.x[np-tailsize+i], dir.y[np-tailsize+i], dir.z[np-tailsize+i]};

       //
       distance[np-tailsize+i]=TGeoBBox::DistFromOutside(p, d, dx, dy, dz, origin, stepmax[np-tailsize+i] );
     }
#endif
}


// SOA version____________________________________________________________________________
void TGeoBBox::DistFromInsideS_v(const StructOfCoord & __restrict__ point,const StructOfCoord & __restrict__ dir, 
				      Double_t dx, Double_t dy, Double_t dz, const Double_t *__restrict__ origin, const Double_t *__restrict__ stepmax, Double_t *__restrict__ distance, int npoints )
{
#ifndef __INTEL_COMPILER 
  const double * x = (const double *) __builtin_assume_aligned (point.x, 32); 
  const double * y = (const double *) __builtin_assume_aligned (point.y, 32); 
  const double * z = (const double *) __builtin_assume_aligned (point.z, 32); 
  const double * dirx = (const double *) __builtin_assume_aligned (dir.x, 32); 
  const double * diry = (const double *) __builtin_assume_aligned (dir.y, 32); 
  const double * dirz = (const double *) __builtin_assume_aligned (dir.z, 32); 
  double * dist = (double *) __builtin_assume_aligned (distance, 32); 
#else
  const double * x = (const double *) point.x; 
  const double * y = (const double *) point.y; 
  const double * z = (const double *) point.z; 
  const double * dirx = (const double *) dir.x; 
  const double * diry = (const double *) dir.y; 
  const double * dirz = (const double *) dir.z; 
#endif

  // this and the previous should be the same; here I have done manual inlining
#pragma vector aligned
  for(size_t k=0; k<npoints; ++k) //@EXPECTVEC
    {
      Double_t s,smin,saf[6];
      Double_t newpt[3];

      newpt[0] = x[k] - origin[0];
      saf[0] = dx + newpt[0];
      saf[1] = dx - newpt[0];
      newpt[1] = y[k] - origin[1];
      saf[2] = dy + newpt[1];
      saf[3] = dy - newpt[1];
      newpt[2] = z[k] - origin[2];
      saf[4] = dz + newpt[2];
      saf[5] = dz - newpt[2];
      
      smin=TGeoShape::Big();
      double sx, sy, sz;
      double tiny=1e-20;
      sx = (dirx[k]>0)? (saf[1]/(dirx[k]+tiny)):(-saf[0]/(dirx[k]-tiny));
      sy = (diry[k]>0)? (saf[3]/(diry[k]+tiny)):(-saf[2]/(diry[k]-tiny));
      sz = (dirz[k]>0)? (saf[5]/(dirz[k]+tiny)):(-saf[4]/(dirz[k]-tiny));
      //      sx = (saf[1]/(std::fabs(dirx[k])+tiny));
      //      sy = (saf[3]/(std::fabs(diry[k])+tiny));
      //      sz = (saf[5]/(std::fabs(dirz[k])+tiny));

      smin = sx;
      smin = (sy < smin)? sy : smin;
      smin = (sz < smin)? sz : smin;
      distance[k] = (smin < 0)? 0 : smin;
    }
}

