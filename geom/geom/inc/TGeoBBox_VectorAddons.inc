// Automatically generated codefile containing addons for vector-interfaces ( as used by the Geant vector prototoype ) 

   virtual void Contains_l( Double_t const *point, Bool_t * isin , Int_t vecsize ) {
         for( int k=0;k < vecsize; ++k){
               isin[k]= TGeoBBox::Contains( (Double_t *) &point[3*k] );
         }
   }
   virtual void ContainsSOA_l( StructOfCoord  const & pointi, Bool_t * isin , Int_t vecsize ) {
         for( int k=0;k < vecsize; ++k){
               double point[3]; point[0]=pointi.x[k];point[1]=pointi.y[k];point[2]=pointi.z[k];
               isin[k]= TGeoBBox::Contains( point );
         }
   }
   virtual void Safety_l( Double_t const *point, Bool_t inside, Double_t * safe , Int_t vecsize ) {
         for( int k=0;k < vecsize; ++k){
               safe[k]= TGeoBBox::Safety( (Double_t *) &point[3*k], inside );
         }
   }
   virtual void SafetySOA_l( StructOfCoord  const & pointi, Bool_t inside, Double_t * safe , Int_t vecsize ) {
         for( int k=0;k < vecsize; ++k){
               double point[3]; point[0]=pointi.x[k];point[1]=pointi.y[k];point[2]=pointi.z[k];
               safe[k]= TGeoBBox::Safety( point, inside );
         }
   }
   virtual void DistFromInside_l( Double_t const *point, Double_t const *dir, Int_t /*iact*/ , Double_t const * step, Double_t * /*safe*/ , Double_t * dist, Int_t vecsize ) {
         for( int k=0;k < vecsize; ++k){
               dist[k]= TGeoBBox::DistFromInside( (Double_t *) &point[3*k], (Double_t *) &dir[3*k], 3, step[k] , 0 );
         }
   }
   virtual void DistFromInsideSOA_l( StructOfCoord const &  pointi, StructOfCoord const &  diri, Int_t /*iact*/ , Double_t const * step, Double_t * /*safe*/ , Double_t * dist, Int_t vecsize ) {
         for( int k=0;k < vecsize; ++k){
               double point[3]; point[0]=pointi.x[k];point[1]=pointi.y[k];point[2]=pointi.z[k];
               double dir[3]; dir[0]=diri.x[k];dir[1]=diri.y[k];dir[2]=diri.z[k];
               dist[k]= TGeoBBox::DistFromInside( point, dir, 3, step[k] , 0 );
         }
   }
   virtual void DistFromOutside_l( Double_t const *point, Double_t const *dir, Int_t /*iact*/, Double_t const * step, Double_t * /*safe*/ , Double_t * dist, Int_t vecsize ) {
         for( int k=0;k < vecsize; ++k){
               dist[k]= TGeoBBox::DistFromOutside( (Double_t *) &point[3*k], (Double_t *) &dir[3*k], 3, step[k] , 0 );
         }
   }
   virtual void DistFromOutsideSOA_l( StructOfCoord  const & pointi, StructOfCoord  const & diri, Int_t /*iact*/, Double_t const * step, Double_t * /*safe*/ , Double_t * dist, Int_t vecsize ) {
         for( int k=0;k < vecsize; ++k){
               double point[3]; point[0]=pointi.x[k];point[1]=pointi.y[k];point[2]=pointi.z[k];
               double dir[3]; dir[0]=diri.x[k];dir[1]=diri.y[k];dir[2]=diri.z[k];
               dist[k]= TGeoBBox::DistFromOutside( point, dir, 3, step[k] , 0 );
         }
   }

   virtual void DistFromOutsideSOA_v( StructOfCoord  const & pointi, StructOfCoord  const & diri, Int_t /*iact*/, Double_t const * step, Double_t * /*safe*/ , Double_t * dist, Int_t vecsize ) const;

   static void DistFromOutsideS_v(StructOfCoord const & point, StructOfCoord const & dir, 
                          double dx, double dy, double dz, const double *origin, const double * stepmax, double * distance, int np);


   virtual void DistFromInsideSOA_v( StructOfCoord const &  pointi, StructOfCoord const &  diri, Int_t /*iact*/ , Double_t const * step, Double_t * /*safe*/ , Double_t * dist, Int_t vecsize ) const {
     TGeoBBox::DistFromInsideS_v(pointi, diri, this->TGeoBBox::GetDX(), this->TGeoBBox::GetDY(), this->TGeoBBox::GetDZ(), this->TGeoBBox::GetOrigin(),  step, dist, vecsize );
   }


   static void DistFromInsideS_v(StructOfCoord const & point, StructOfCoord const & dir, 
                          double dx, double dy, double dz, const double *origin, const double * stepmax, double * distance, int np);

 // we are also exposing the elementary vector function
 // this is the actual kernel doing the computation with possibility of early return

//  static void DistFromOutsideSOA_Vc( Vc::double_v const & x , Vc::double_v const & y, Vc::double_v const & z, Vc::double_v const & dirx, Vc::double_v const & diry, Vc::double_v //const & dirz, double dx, double dy, double dz, Vc::double_v const origin[3], Vc::double_v const & stepmax, Vc::double_v & distance ) 

#ifndef __CINT__
static void DistFromOutsideSOA_Vc( Vc::double_v const & x , Vc::double_v const & y, Vc::double_v const & z, Vc::double_v const & dirx, Vc::double_v const & diry, Vc::double_v const & dirz, double dx, double dy, double dz, Vc::double_v const origin[3], Vc::double_v const & stepmax, Vc::double_v & distance ) 
{
   Vc::double_v in(1.);
   Vc::double_v saf[3];
   Vc::double_v newpt[3];
   Vc::double_v tiny(1e-20);
   Vc::double_v big(1e30);
   Vc::double_v faraway(0.); // initializing all components to zero
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
 }
#endif



