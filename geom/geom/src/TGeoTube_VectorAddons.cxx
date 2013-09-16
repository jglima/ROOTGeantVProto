
#ifndef __CINT__
#include "Vc/vector.h"
#include "PointStruct.h"
#endif

#include "TGeoTube.h"

//#ifdef VEC_EXTENSIONS // REPEATED: WHERE CAN ALL THIS GO

typedef Vc::double_v vd; // short for vector double
typedef Vc::double_m vdm; // short for double mask 
typedef Vc::int_v vi; // short for vector integer

static vd tol_v = 1.E-10;

struct VecUtil{
static Vc::double_m IsSameWithinTolerance(Vc::double_v const & a, Vc::double_v const & b )
{
  Vc::double_m c = Vc::abs(a-b) < tol_v;
  return c;
}
}; 


//_____________________________________________________________________________ 
inline                    
void DistToTube_Vc(Vc::double_v const & r2_v, Vc::double_v const & n2_v, Vc::double_v const & rdotn_v, Double_t radius ,Vc::double_v & b_v, Vc::double_v & delta_v) 
{
  vd t1_v = 1./n2_v;
  vd radius_v(radius);
  vd t3_v = r2_v-(radius_v*radius_v);
  b_v = t1_v*rdotn_v;
  vd c_v = t1_v*t3_v;
  delta_v = b_v*b_v-c_v;
  vdm delta_m = delta_v > Vc::Zero;
  delta_v(delta_m) = Vc::sqrt(delta_v);
  delta_v(!delta_m) = -1;
}


Vc::double_v DistFromInsideS_Vc(Vc::double_v const & x_v, Vc::double_v const & y_v, Vc::double_v const & z_v, Vc::double_v const & dirx_v, Vc::double_v const & diry_v, Vc::double_v const & dirz_v, Double_t rmin, Double_t rmax, Double_t dz)
{
  // Compute distance from inside point to surface of the tube (static)                                                                                                                                  
  // Boundary safe algorithm.                                                                                                                                                                            
  // compute distance to surface                                                                                                                                                                         
  // Do Z                                                                                                                                                                                                
  vd const dz_v(dz);

  vd s_v(1.E30); // vector with the results                                                                                                                                                              

  vd dirzsign_v(1.);
  dirzsign_v(dirz_v < 0) = -1.;

 vd sz_v = (dz_v * dirzsign_v - z_v)/dirz_v;
  s_v(sz_v < 0) = 0.;
  vdm done_m = sz_v < 0;

  if ( done_m.isFull() ) return s_v;

  // Distance to R                                                                                                                                                                                       
  vd n2_v = dirx_v*dirx_v + diry_v*diry_v;
  s_v(Vc::abs(n2_v) < tol_v) = sz_v;
  done_m |= Vc::abs(n2_v) < tol_v;

  vd r2_v = x_v*x_v + y_v*y_v;
  vd rdotn_v = x_v*dirx_v + y_v*diry_v;
  vdm rdotn_m = rdotn_v > 0.;
  vd const rmin2_v = rmin*rmin;
  vd const rmax2_v = rmax*rmax;

  vd b_v(0.);
  vd d_v(0.);
  vd sr_v(1.E30);

  // distance to the inner cylinder                                                                                                                                                                      
  if (rmin > 0.)
    {
      vdm r2_m = r2_v <= (rmin2_v + tol_v);
      s_v(!done_m && r2_m && !rdotn_m) = 0.;
      done_m |= r2_m && !rdotn_m;

      DistToTube_Vc(r2_v, n2_v, rdotn_v, rmin, b_v, d_v);
      sr_v = -b_v - d_v;
      vdm sr_m = (sr_v > 0.) && (d_v > 0.) && !rdotn_m && !r2_m;
      s_v(!done_m && sr_m) = Vc::min(sz_v, sr_v);
      done_m |= sr_m;
    }
  if ( done_m.isFull() ) return s_v;

  //distance to the outer cylinder                                                                                                                                                                       

  vdm r2_m = r2_v >= (rmax2_v - tol_v);
  s_v(!done_m && r2_m && rdotn_m) = 0.;
  done_m |= r2_m && rdotn_m;

  DistToTube_Vc(r2_v, n2_v, rdotn_v, rmax, b_v, d_v);
  sr_v = -b_v + d_v;
  vdm sr_m = (sr_v > 0.) && (d_v > 0.);
  s_v(!done_m && sr_m) = Vc::min(sz_v, sr_v);
  done_m |= sr_m;

  s_v(!done_m) = 0.;
  return s_v;
}


// //_____________________________________________________________________________
// inline
// Vc::double_v DistFromInsideS_Vc(Vc::double_v const & x_v, Vc::double_v const & y_v, Vc::double_v const & z_v, Vc::double_v const & dirx_v, Vc::double_v const & diry_v, Vc::double_v const & dirz_v, Double_t rmin, Double_t rmax, Double_t dz) 
// {
//   // Compute distance from inside point to surface of the tube (static)
//   // Boundary safe algorithm.
//   // compute distance to surface
//   // Do Z

//   /* early check if really inside */
//   vd sz_v(1.E30); // was TGeoShape::Big(), sz_v keeps distance to bottom and top plate 
//   vd const dz_v(dz);

//   vd dirzsign_v(1.);

//   static const vd tiny_v(1.E-20);

//   dirzsign_v(dirz_v < 0)=-1.;
//   // do distance to top or bottom
//   sz_v = (dz_v * dirzsign_v - z_v)/(dirz_v + tiny_v);
  
//   // if in reallity outside ( might skip this )
//   // if (sz_v <= 0 ) return 0.0;
//   sz_v( sz_v < 0 )=0.;
//   // todo
//   //  done_m = sz_v < 0;
//   // later you say done_m |= ...

//   // do distance to R
//   vd nsq_v=dirx_v*dirx_v + diry_v*diry_v;

//   // this has to be checked 
//   //  if (Vc::abs(nsq_v)<TGeoShape::Tolerance()) return sz_v;

//   vd rsq_v = x_v*x_v+y_v*y_v;
//   vd rdotn_v = x_v*dirx_v+ y_v*diry_v;
//   vd b_v,d_v;
//   vd sr_v(1.E30);

//   // inner cylinder only if needed
//   if( rmin>0 ) 
//     {
//       // Protection in case point is actually outside the tube

//       /* switch off for the moment */
//       //  if (rsq <= rmin*rmin+TGeoShape::Tolerance()) {
//       //	if (rdotn<0) return 0.0;
//       // } 

//       // else {
//       DistToTube_Vc(rsq_v,nsq_v,rdotn_v,rmin,b_v,d_v);
//       sr_v( (rdotn_v < 0.) & (d_v > 0.) ) = -b_v - d_v;

//       vdm allhitinner_m = sr_v>0.;

//       // update done here
//       if (allhitinner_m.isFull()) return Vc::min(sz_v,sr_v);
//     }
 
//   // outer cylinder:

//   // check disabled for the moment
//   //  if (rsq >= rmax*rmax-TGeoShape::Tolerance()) {
//   //  if (rdotn>=0) return 0.0;
//   // }

//   DistToTube_Vc(rsq_v,nsq_v,rdotn_v,rmax,b_v,d_v);
//   sr_v(d_v>0) = -b_v + d_v;
//   return Vc::min(sz_v,sr_v);
//   // }
//   // return 0.;
// }
//*/

void TGeoTube::DistFromInsideS_v(StructOfCoord const & pointi, StructOfCoord const & diri, Double_t rmin, Double_t rmax, Double_t dz, Double_t * distance, Int_t np)
{
  int tailsize = np % Vc::double_v::Size; 
  
  for(volatile unsigned int i = 0; i < np-tailsize; i+= Vc::double_v::Size)
    {
      vd x_v(&pointi.x[i]); // will copy a certain number of x's into vc vector x;
      vd y_v(&pointi.y[i]);   
      vd z_v(&pointi.z[i]);
      vd dirx_v(&diri.x[i]);
      vd diry_v(&diri.y[i]);
      vd dirz_v(&diri.z[i]);
      
      vd dist_v(0.);
      
      dist_v = DistFromInsideS_Vc(x_v, y_v, z_v, dirx_v, diry_v, dirz_v, rmin, rmax, dz);
      
      // write the results for each particle
      dist_v.store(&distance[i]);
    }
  
  // do the tail part for the moment, we just call the old static version
  for(unsigned int i = 0; i < tailsize; ++i)
    {
      Double_t point[3] = {pointi.x[np-tailsize+i], pointi.y[np-tailsize+i], pointi.z[np-tailsize+i]};
      double dir[3] = {diri.x[np-tailsize+i], diri.y[np-tailsize+i], diri.z[np-tailsize+i]};
      distance[np-tailsize+i] = TGeoTube::DistFromInsideS(point, dir, rmin, rmax, dz);
    }
}



// implicit kernel function_____________________________________________________________________________                                    
inline                                                         
Vc::double_v DistFromOutsideS_Vc(Vc::double_v const & x_v, Vc::double_v const & y_v, Vc::double_v const & z_v, Vc::double_v const & dirx_v, Vc::double_v const & diry_v, Vc::double_v const & dirz_v, Double_t rmin, Double_t rmax, Double_t dz) 
{
  // Static method to compute distance from outside point to a tube with given parameters
  
  vdm done_m(true); // which elements of the vector are ready to be returned
  vd s_v(1.E30);
  
  // check Z planes
  
  vd dz_v(dz);
  vd zi_v = dz_v - Vc::abs(z_v);
  vdm inz_m = zi_v >= Vc::Zero;
  
  done_m = !inz_m && (z_v*dirz_v >= Vc::Zero); // particle outside the z-range and moving away 
  if( done_m.isFull() ) return s_v;
  
  s_v(!done_m) = -zi_v/Vc::abs(dirz_v);
  vd xi_v = x_v + s_v*dirx_v;
  vd yi_v = y_v + s_v*diry_v;
  vd ri2_v = xi_v*xi_v + yi_v*yi_v;
  vd rmin2_v = rmin*rmin;
  vd rmax2_v = rmax*rmax;
  
  done_m |= !inz_m && (rmin2_v <= ri2_v) && (ri2_v <= rmax2_v);
  if( done_m.isFull() ) return s_v;
  
  // TO DO: CHECK IF IT'S INSIDE (BOUNDARY WITHIN MACHINE PRECISION, TOO SPECIFIC)

  vd r2_v = x_v*x_v + y_v*y_v;
  vd n2_v = dirx_v*dirx_v + diry_v*diry_v;
  vd rdotn_v = x_v*dirx_v + y_v*diry_v;
  vd b_v(0.);
  vd d_v(0.);
  vdm inrmax_m = (r2_v - rmax2_v) <= tol_v; 
  // vdm inrmin_m = (rmin2_v - r2_v) <= tol_v;

  // Check outer cylinder (only r>rmax has to be considered)
  vdm still_m = Vc::abs(n2_v) < tol_v;
  s_v(!done_m && still_m) = 1.E30;
  done_m |= still_m;
  if ( done_m.isFull() ) return s_v;
  
  DistToTube_Vc(r2_v, n2_v, rdotn_v, rmax, b_v, d_v);

  s_v(!done_m) = -b_v - d_v;
  zi_v = z_v + s_v*dirz_v;
  done_m |= !inrmax_m && (d_v > 0) && (s_v > 0) && (Vc::abs(zi_v) <= dz_v);
  if ( done_m.isFull() ) return s_v;
  
  // check inner cylinder
  if(rmin > 0)
    {
      DistToTube_Vc(r2_v, n2_v, rdotn_v, rmin, b_v, d_v);
      s_v(!done_m) = -b_v + d_v;
      zi_v = z_v + s_v*dirz_v;
      done_m |= (d_v > 0) && (s_v > 0) && (Vc::abs(zi_v) <= dz_v);
    }

  s_v(!done_m) = 1.E30;
  return s_v;
  
}


void TGeoTube::DistFromOutsideS_v(StructOfCoord const & pointi, StructOfCoord const & diri, Double_t rmin, Double_t rmax, Double_t dz, Double_t * distance, Int_t np)
{
  int tailsize = np % Vc::double_v::Size; 
  
  for(volatile unsigned int i = 0; i < np-tailsize; i+= Vc::double_v::Size)
    {
      vd x_v(&pointi.x[i]); // will copy a certain number of x's into vc vector x;
      vd y_v(&pointi.y[i]);   
      vd z_v(&pointi.z[i]);
      vd dirx_v(&diri.x[i]);
      vd diry_v(&diri.y[i]);
      vd dirz_v(&diri.z[i]);
      
      vd dist_v(0.);
      
      dist_v = DistFromOutsideS_Vc(x_v, y_v, z_v, dirx_v, diry_v, dirz_v, rmin, rmax, dz);
      
      // write the results for each particle
      dist_v.store(&distance[i]);
    }
  
  // do the tail part for the moment, we just call the old static version
  for(unsigned int i = 0; i < tailsize; ++i)
    {
      Double_t point[3] = {pointi.x[np-tailsize+i], pointi.y[np-tailsize+i], pointi.z[np-tailsize+i]};
      double dir[3] = {diri.x[np-tailsize+i], diri.y[np-tailsize+i], diri.z[np-tailsize+i]};
      distance[np-tailsize+i] = TGeoTube::DistFromOutsideS(point, dir, rmin, rmax, dz);
    }
}


