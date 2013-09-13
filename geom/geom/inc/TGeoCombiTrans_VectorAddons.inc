virtual void         MasterToLocal_l(Double_t const *master, Double_t *local, unsigned int np) const {
for( unsigned int k=0; k< np; ++k)
  {
 this->TGeoCombiTrans::MasterToLocal(&master[3*k], &local[3*k]);
  }
}

virtual void         MasterToLocalVect_l(Double_t const *master, Double_t *local, unsigned int np) const {
  for( unsigned int k=0; k< np; ++k)
   {
	 this->TGeoCombiTrans::MasterToLocalVect(&master[3*k], &local[3*k]);
   }
}

virtual void         MasterToLocalBomb_l(Double_t const *master, Double_t *local, unsigned int np) const {
  for( unsigned int k=0; k< np; ++k)
    {
	 this->TGeoCombiTrans::MasterToLocalBomb(&master[3*k], &local[3*k]);
    }
}

virtual void         LocalToMaster_l(Double_t const *master, Double_t *local, unsigned int np) const {
  for( unsigned int k=0; k< np; ++k)
    {
	 this->TGeoCombiTrans::LocalToMaster(&master[3*k], &local[3*k]);
    }
}

virtual void         LocalToMasterVect_l(Double_t const *master, Double_t *local, unsigned int np) const {
  for( unsigned int k=0; k< np; ++k)
    {
	 this->TGeoCombiTrans::LocalToMasterVect(&master[3*k], &local[3*k]);
    }
}

virtual void         LocalToMasterBomb_l(Double_t const *master, Double_t *local, unsigned int np) const {
  for( unsigned int k=0; k< np; ++k)
    {
	 this->TGeoCombiTrans::LocalToMasterBomb(&master[3*k], &local[3*k]);
    }
}

//virtual void MasterToLocal_v(StructOfCoord const &master, StructOfCoord &local, Int_t np ) const;

//virtual void LocalToMaster_v(StructOfCoord const &local, StructOfCoord &master, Int_t np ) const;

//virtual void MasterToLocalVect_v(StructOfCoord const &master, StructOfCoord &local, Int_t np ) const;

//virtual void LocalToMasterVect_v(StructOfCoord const &local, StructOfCoord &master, Int_t np ) const;

//virtual void MasterToLocalCombined_v(StructOfCoord const &master, StructOfCoord &local, StructOfCoord const &mastervec, StructOfCoord &localvec, Int_t np ) const;

//virtual void LocalToMasterCombined_v(StructOfCoord const &local, StructOfCoord &master, StructOfCoord const &localvec, StructOfCoord &mastervec, Int_t np ) const;

