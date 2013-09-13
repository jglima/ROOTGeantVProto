#ifndef __CINT__
#include "Vc/vector.h"
#endif

typedef Vc::double_v vd; // short for vector double

void TGeoBBox::DistFromOutside_v(StructOfCoord const & point, StructOfCoord const & dir, Int_t iact, 
			Double_t const * step, Double_t *safe, Double_t * distance , Int_t np) const 
{
  TGeoBBox_v::DistFromOutsideS_v(point, dir, this->TGeoBBox::GetDX(), this->TGeoBBox::GetDY(), this->TGeoBBox::GetDZ(), this->TGeoBBox::GetOrigin(),  step, distance, np );
}


 // this is the actual kernel doing the computation with possibility of early return
void kernel( Vc::double_v const & x , Vc::double_v const & y, Vc::double_v const & z, Vc::double_v const & dirx, Vc::double_v const & diry, Vc::double_v const & dirz, double dx, double dy, double dz, Vc::double_v const origin[3], Vc::double_v const & stepmax, Vc::double_v & distance )
 {
#ifndef __CINT__
   Vc::double_v in=Vc::double_v(1.);
   Vc::double_v saf[3];
   Vc::double_v newpt[3];
   Vc::double_v tiny(1e-20);
   Vc::double_v big(1e30);
   Vc::double_v faraway=Vc::double_v(0.); // initializing all components to zero
   Vc::double_v par[3]={dx,dy,dz}; // very convenient

   newpt[0] = x-origin[0];
   saf[0] = Vc::abs(newpt[0])-par[0]; // can we profit here from abs function in array form?
   newpt[1] = y-origin[1];
   saf[1] = Vc::abs(newpt[1])-par[1];
   newpt[2] = z-origin[2];
   saf[2] = Vc::abs(newpt[2])-par[2];
   faraway(saf[0]>=stepmax | saf[1]>=stepmax | saf[2]>=stepmax)=1; 
   in(saf[0]<0. & saf[1]<0. & saf[2]<0.)=0;
   distance=big;

   if( faraway > Vc::double_v(0.) )
     {
       return; // return big
     }
  
   // proceed to analysis of hits
   Vc::double_v snxt[3];
   Vc::double_v hit0=Vc::double_v(0.);
   snxt[0] = saf[0]/Vc::abs(dirx+tiny); // distance to y-z face
   Vc::double_v coord1=newpt[1]+snxt[0]*diry; // calculate new y and z coordinate
   Vc::double_v coord2=newpt[2]+snxt[0]*dirz;
   hit0( saf[0] > 0 & newpt[0]*dirx < 0 & ( Vc::abs(coord1)<= par[1] & Vc::abs(coord2)<= par[2] ) ) = 1; // if out and right direction
   
   Vc::double_v hit1=Vc::double_v(0.);
   snxt[1] = saf[1]/Vc::abs(diry+tiny); // distance to x-z face
   coord1=newpt[0]+snxt[1]*dirx; // calculate new x and z coordinate
   coord2=newpt[2]+snxt[1]*dirz;
   hit1( saf[1] > 0 & newpt[1]*diry < 0 & ( Vc::abs(coord1)<= par[0] & Vc::abs(coord2)<= par[2] ) ) = 1; // if out and right direction

   Vc::double_v hit2=Vc::double_v(0.);
   snxt[2] = saf[2]/Vc::abs(dirz+tiny); // distance to x-y face
   coord1=newpt[0]+snxt[2]*dirx; // calculate new x and y coordinate
   coord2=newpt[1]+snxt[2]*diry;
   hit2( saf[2] > 0 & newpt[2]*dirz < 0 & ( Vc::abs(coord1)<= par[0] & Vc::abs(coord2)<= par[1] ) ) = 1; // if out and right direction

   distance( hit0>0 | hit1>0 | hit2>0 )=(hit0*snxt[0] + hit1*snxt[1] + hit2*snxt[2]);
   distance=in*distance;
   return;
#endif
}


// // this is the actual kernel doing the computation with possibility of early return
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
      DistFromOutsideS_kernel( x, y, z, dirx, diry, dirz, dx, dy, dz, vcorigin, s, dist );
      dist.store(&distance[i]);
   }
   // do the tail part for the moment, we just call the old static version
   for(unsigned int i = 0; i < tailsize; ++i)
     {
       double p[3]={point.x[np-tailsize+i], point.y[np-tailsize+i], point.z[np-tailsize+i]};
       double d[3]={dir.x[np-tailsize+i], dir.y[np-tailsize+i], dir.z[np-tailsize+i]};
       distance[np-tailsize+i]=TGeoBBox_v::DistFromOutsideS(p, d, dx, dy, dz, origin, stepmax[np-tailsize+i] );
     }
#endif
}


