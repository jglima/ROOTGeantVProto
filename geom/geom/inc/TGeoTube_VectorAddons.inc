// Automatically generated codefile containing addons for vector-interfaces ( as used by the Geant vector prototoype ) 

   virtual void Contains_l( Double_t const *point, Bool_t * isin , Int_t vecsize ) {
         for( int k=0;k < vecsize; ++k){
               isin[k]= TGeoTube::Contains( (Double_t *) &point[3*k] );
         }
   }
   virtual void ContainsSOA_l( StructOfCoord  const & pointi, Bool_t * isin , Int_t vecsize ) {
         for( int k=0;k < vecsize; ++k){
               double point[3]; point[0]=pointi.x[k];point[1]=pointi.y[k];point[2]=pointi.z[k];
               isin[k]= TGeoTube::Contains( point );
         }
   }
   virtual void Safety_l( Double_t const *point, Bool_t inside, Double_t * safe , Int_t vecsize ) {
         for( int k=0;k < vecsize; ++k){
               safe[k]= TGeoTube::Safety( (Double_t *) &point[3*k], inside );
         }
   }
   virtual void SafetySOA_l( StructOfCoord  const & pointi, Bool_t inside, Double_t * safe , Int_t vecsize ) {
         for( int k=0;k < vecsize; ++k){
               double point[3]; point[0]=pointi.x[k];point[1]=pointi.y[k];point[2]=pointi.z[k];
               safe[k]= TGeoTube::Safety( point, inside );
         }
   }
   virtual void DistFromInside_l( Double_t const *point, Double_t const *dir, Int_t /*iact*/ , Double_t const * step, Double_t * /*safe*/ , Double_t * dist, Int_t vecsize ) {
         for( int k=0;k < vecsize; ++k){
               dist[k]= TGeoTube::DistFromInside( (Double_t *) &point[3*k], (Double_t *) &dir[3*k], 3, step[k] , 0 );
         }
   }
   virtual void DistFromInsideSOA_l( StructOfCoord const &  pointi, StructOfCoord const &  diri, Int_t /*iact*/ , Double_t const * step, Double_t * /*safe*/ , Double_t * dist, Int_t vecsize ) {
         for( int k=0;k < vecsize; ++k){
               double point[3]; point[0]=pointi.x[k];point[1]=pointi.y[k];point[2]=pointi.z[k];
               double dir[3]; dir[0]=diri.x[k];dir[1]=diri.y[k];dir[2]=diri.z[k];
               dist[k]= TGeoTube::DistFromInside( point, dir, 3, step[k] , 0 );
         }
   }
   virtual void DistFromOutside_l( Double_t const *point, Double_t const *dir, Int_t /*iact*/, Double_t const * step, Double_t * /*safe*/ , Double_t * dist, Int_t vecsize ) {
         for( int k=0;k < vecsize; ++k){
               dist[k]= TGeoTube::DistFromOutside( (Double_t *) &point[3*k], (Double_t *) &dir[3*k], 3, step[k] , 0 );
         }
   }
   virtual void DistFromOutsideSOA_l( StructOfCoord  const & pointi, StructOfCoord  const & diri, Int_t /*iact*/, Double_t const * step, Double_t * /*safe*/ , Double_t * dist, Int_t vecsize ) {
         for( int k=0;k < vecsize; ++k){
               double point[3]; point[0]=pointi.x[k];point[1]=pointi.y[k];point[2]=pointi.z[k];
               double dir[3]; dir[0]=diri.x[k];dir[1]=diri.y[k];dir[2]=diri.z[k];
               dist[k]= TGeoTube::DistFromOutside( point, dir, 3, step[k] , 0 );
         }
   }


// REAL VECTOR INTERFACES: THIS IMPLEMENTATION SHOULD BE GRADUALLY MOVED SOMEWHERE ELSE 

   virtual void DistFromInsideSOA_v( StructOfCoord const &  pointi, StructOfCoord const &  diri, Int_t /*iact*/ , Double_t const * step, Double_t * /*safe*/ , Double_t * dist, Int_t vecsize ) const {
         for( int k=0;k < vecsize; ++k){
               double point[3]; point[0]=pointi.x[k];point[1]=pointi.y[k];point[2]=pointi.z[k];
               double dir[3]; dir[0]=diri.x[k];dir[1]=diri.y[k];dir[2]=diri.z[k];
               dist[k]= TGeoTube::DistFromInside( point, dir, 3, step[k] , 0 );
         }
   }


// REAL VECTOR INTERFACES: THIS IMPLEMENTATION SHOULD BE GRADUALLY MOVED SOMEWHERE ELSE 

   virtual void DistFromOutsideSOA_v( StructOfCoord  const & pointi, StructOfCoord  const & diri, Int_t /*iact*/, Double_t const * step, Double_t * /*safe*/ , Double_t * dist, Int_t vecsize ) const {
         for( int k=0;k < vecsize; ++k){
               double point[3]; point[0]=pointi.x[k];point[1]=pointi.y[k];point[2]=pointi.z[k];
               double dir[3]; dir[0]=diri.x[k];dir[1]=diri.y[k];dir[2]=diri.z[k];
               dist[k]= TGeoTube::DistFromOutside( point, dir, 3, step[k] , 0 );
         }
   }
