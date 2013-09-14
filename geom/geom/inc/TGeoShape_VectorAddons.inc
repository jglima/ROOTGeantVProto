// Automatically generated codefile containing addons for vector-interfaces ( as used by the Geant vector prototoype ) 


   virtual void Contains_l( Double_t const *, Bool_t * , Int_t  ) {  }
   virtual void ContainsSOA_l( StructOfCoord  const & , Bool_t *  , Int_t  ) {   }
   virtual void Safety_l( Double_t const *, Bool_t , Double_t *  , Int_t  ) {   }
   virtual void SafetySOA_l( StructOfCoord const &  , Bool_t , Double_t *  , Int_t  ) {   }
   virtual void DistFromInside_l( Double_t const *, Double_t const *, Int_t , Double_t const * , Double_t * , Double_t * , Int_t  ) {   } 
   virtual void DistFromInsideSOA_l( StructOfCoord const &  , StructOfCoord  const & , Int_t , Double_t const * , Double_t * , Double_t * , Int_t  ) {   }
   virtual void DistFromOutside_l( Double_t const *, Double_t const *, Int_t , Double_t const * , Double_t * , Double_t * , Int_t  ) {   }
   virtual void DistFromOutsideSOA_l( StructOfCoord const &  , StructOfCoord const & , Int_t , Double_t const * , Double_t * , Double_t * , Int_t  ) {   }


// REAL VECTOR INTERFACES: THIS IMPLEMENTATION SHOULD BE GRADUALLY MOVED SOMEWHERE ELSE 

   virtual void DistFromInsideSOA_v( StructOfCoord const &  pointi, StructOfCoord const &  diri, Int_t /*iact*/ , Double_t const * step, Double_t * /*safe*/ , Double_t * dist, Int_t vecsize ) =0;
 

// REAL VECTOR INTERFACES: THIS IMPLEMENTATION SHOULD BE GRADUALLY MOVED SOMEWHERE ELSE 

   virtual void DistFromOutsideSOA_v( StructOfCoord  const & pointi, StructOfCoord  const & diri, Int_t /*iact*/, Double_t const * step, Double_t * /*safe*/ , Double_t * /* dist */, Int_t vecsize ) =0;
 
