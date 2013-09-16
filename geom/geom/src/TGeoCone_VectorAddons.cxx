// implementation copied from Geant-V geometry development branch

#ifndef __CINT__
#include "Vc/vector.h"
#include "PointStruct.h"
#endif

#include "TGeoCone.h"

//#ifdef VEC_EXTENSIONS // REPEATED: WHERE CAN ALL THIS GO

typedef Vc::double_v vd; // short for vector double
typedef Vc::double_m vdm; // short for double mask 
typedef Vc::int_v vi; // short for vector integer

static vd tol_v = 1.E-10;

inline                      
void DistToCone_v(Vc::double_v const & x_v, Vc::double_v const & y_v, Vc::double_v const & z_v, Vc::double_v const & dirx_v, Vc::double_v const & diry_v, Vc::double_v const & dirz_v, Double_t dz, Double_t rmin, Double_t rmax, Vc::double_v & b_v, Vc::double_v & delta_v)
{
  delta_v = -1.;
  if (dz<0) return;
  
  vdm done_m(false);
  
  vd ro0_v = 0.5 * (rmin + rmax);
  vd tz_v = 0.5 * (rmax - rmin)/dz;
  vd r2_v = x_v * x_v + y_v * y_v;
  vd rc_v = ro0_v + z_v * tz_v;

  vd a_v = dirx_v * dirx_v + diry_v * diry_v - tz_v * tz_v * dirz_v * dirz_v;
  b_v = x_v * dirx_v + y_v * diry_v - tz_v * rc_v * dirz_v;
  vd c_v = r2_v - rc_v * rc_v;
  
  vdm a_m = Vc::abs(a_v) < tol_v;
  done_m = a_m && (Vc::abs(b_v) < tol_v);
  if ( done_m.isFull() ) return;
  
  b_v(!done_m && a_m) = 0.5 * c_v/b_v;
  delta_v(!done_m && a_m) = Vc::Zero;
  done_m |= a_m;
  if ( done_m.isFull() ) return;

  a_v = 1./a_v;
  b_v(!done_m) *= a_v;
  c_v *= a_v;
  delta_v(!done_m) = b_v * b_v - c_v;
  delta_v(!done_m && (delta_v > Vc::Zero)) = Vc::sqrt(delta_v);
  done_m |= delta_v > Vc::Zero;
  delta_v(!done_m) = -1.;
  return;
}


//_____________________________________________________________________________                                                              inline                      
Vc::double_v DistFromOutsideS_Vc(Vc::double_v const & x_v, Vc::double_v const & y_v, Vc::double_v const & z_v, Vc::double_v const & dirx_v, Vc::double_v const & diry_v, Vc::double_v const & dirz_v, Double_t dz, Double_t rminl, Double_t rmaxl, Double_t rminh, Double_t rmaxh)
{
  vd s_v(1.E30);

  if (dz<=0) return s_v;
  
  // Distance to the z-planes
  vdm done_m(false);
  vd dz_v(dz);
  vdm inz_m = Vc::abs(z_v) <= dz_v;

  done_m = !inz_m && ((z_v * dirz_v) >=0); // Point outside the z-range and moving away
  if( done_m.isFull() ) return s_v;

  s_v(!done_m) = (Vc::abs(z_v) - dz)/Vc::abs(dirz_v);
  vd xp_v = x_v + s_v * dirx_v;
  vd yp_v = y_v + s_v * diry_v;
  vd r2_v = xp_v * xp_v + yp_v * yp_v;

  vd rminl2_v = rminl * rminl;
  vd rmaxl2_v = rmaxl * rmaxl;
  vd rminh2_v = rminh * rminh;
  vd rmaxh2_v = rmaxh * rmaxh;

  vdm reach_m = (z_v <= -dz_v) && (r2_v >= rminl2_v) && (r2_v <= rmaxl2_v);
  reach_m |= (z_v >= dz_v) && (r2_v >= rminh2_v) && (r2_v <= rmaxh2_v);

  done_m |= reach_m;
  if ( done_m.isFull() ) return s_v;

  // Determine in which region in the r-range the particle is
  r2_v = x_v * x_v + y_v * y_v;

  Double_t dzinv = 1./dz;
  Double_t ro1 = 0.5*(rminl+rminh);
  Bool_t hasrmin = (ro1>0)?kTRUE:kFALSE;
  vdm inrmin_m(true);

  if (hasrmin)  {
    vd tg1_v = 0.5*(rminh-rminl)*dzinv;
    vd rin_v(ro1);
    rin_v += tg1_v*z_v;
    inrmin_m = (rin_v <= Vc::Zero) || (r2_v >= rin_v*(rin_v - tol_v));
  } else {
    vd tg1_v(0.);
    vd rin_v(0.);
  }

  vd ro2_v = 0.5 * (rmaxl + rmaxh);
  vd tg2_v = 0.5 * (rmaxh - rmaxl) * dzinv;
  vd rout_v = tg2_v * z_v + ro2_v;
  vdm inrmax_m = (rout_v > Vc::Zero) && (r2_v<rout_v*(rout_v+tol_v));

  vdm in_m = inz_m && inrmin_m && inrmax_m;
  
  // TO DO: CASE INSIDE CONE (BOUNDARY WITHIN MACHINE PRECISION)
  
  s_v(!done_m) = 1.E30;

  // Compute distance to inner cone
  vd b_v(0.);
  vd delta_v(0.);
  vd d_v(1.E30);
  vd zp_v(0.);
  vdm zp_m(false);

   // hasrmin part
  if(hasrmin)
    {
      DistToCone_v(x_v, y_v, z_v, dirx_v, diry_v, dirz_v, dz, rminl, rminh, b_v, delta_v);
      done_m |= (delta_v < 0) && !inrmin_m;

      d_v = -b_v + delta_v;
      d_v((d_v < 0) && !inrmin_m) = -b_v - delta_v;
      zp_v = z_v + d_v * dirz_v;
      zp_m = (d_v > 0) && (delta_v > 0) && (Vc::abs(zp_v) <= dz_v);

      s_v(!done_m && zp_m) = d_v;
      done_m |= zp_m && !inrmin_m;
    }
  
  done_m |= inrmax_m;
  if (done_m.isFull() ) return s_v;

  // We can cross outer cone, both solutions possible
  // compute distance to outer cone

  DistToCone_v(x_v, y_v, z_v, dirx_v, diry_v, dirz_v, dz, rmaxl, rmaxh, b_v, delta_v);
  done_m |= delta_v < 0;
  if( done_m.isFull() ) return s_v;
  
  d_v = -b_v - delta_v;
  vdm d_m = (d_v>0) && (d_v<s_v);
  zp_v = z_v + d_v * dirz_v;
  zp_m = (Vc::abs(zp_v) <= dz_v) && d_m;
  
  s_v(!done_m && zp_m) = d_v;
  done_m |= zp_m;
  
  d_v = -b_v + delta_v;
  done_m |= (d_v <= 0) || (d_v > s_v);
  if( done_m.isFull() ) return s_v;

  zp_v = z_v + d_v * dirz_v;
  zp_m = Vc::abs(zp_v) <= dz_v;
  
  s_v(!done_m && zp_m) = d_v;

  return s_v;
}


void TGeoCone::DistFromOutsideS_v(StructOfCoord const & pointi, StructOfCoord const & diri, Double_t dz, Double_t rminl, Double_t rmaxl, Double_t rminh, Double_t rmaxh, Double_t * distance, Int_t np)
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
      
      dist_v = DistFromOutsideS_Vc(x_v, y_v, z_v, dirx_v, diry_v, dirz_v, dz, rminl, rmaxl, rminh, rmaxh);
      
      // write the results for each particle
      dist_v.store(&distance[i]);
    }
  
  // do the tail part for the moment, we just call the old static version
  for(unsigned int i = 0; i < tailsize; ++i)
    {
      Double_t point[3] = {pointi.x[np-tailsize+i], pointi.y[np-tailsize+i], pointi.z[np-tailsize+i]};
      double dir[3] = {diri.x[np-tailsize+i], diri.y[np-tailsize+i], diri.z[np-tailsize+i]};
      distance[np-tailsize+i] = TGeoCone::DistFromOutsideS(point, dir, dz, rminl, rmaxl, rminh, rmaxh);
    }
}

