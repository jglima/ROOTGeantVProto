virtual void         MasterToLocal_l(Double_t const *master, Double_t *local, unsigned int np) const {
for( unsigned int k=0; k< np; ++k)
  {
 this->TGeoIdendity::MasterToLocal(&master[3*k], &local[3*k]);
  }
}
virtual void         MasterToLocalVect_l(Double_t const *master, Double_t *local, unsigned int np) const {
  for( unsigned int k=0; k< np; ++k)
   {
	 this->TGeoIdendity::MasterToLocalVect(&master[3*k], &local[3*k]);
   }
}
virtual void         MasterToLocalBomb_l(Double_t const *master, Double_t *local, unsigned int np) const {
  for( unsigned int k=0; k< np; ++k)
    {
	 this->TGeoIdendity::MasterToLocalBomb(&master[3*k], &local[3*k]);
    }
}
virtual void         LocalToMaster_l(Double_t const *master, Double_t *local, unsigned int np) const {
  for( unsigned int k=0; k< np; ++k)
    {
	 this->TGeoIdendity::LocalToMaster(&master[3*k], &local[3*k]);
    }
}
virtual void         LocalToMasterVect_l(Double_t const *master, Double_t *local, unsigned int np) const {
  for( unsigned int k=0; k< np; ++k)
    {
	 this->TGeoIdendity::LocalToMasterVect(&master[3*k], &local[3*k]);
    }
}
virtual void         LocalToMasterBomb_l(Double_t const *master, Double_t *local, unsigned int np) const {
  for( unsigned int k=0; k< np; ++k)
    {
	 this->TGeoIdendity::LocalToMasterBomb(&master[3*k], &local[3*k]);
    }
}

virtual void MasterToLocal_v(StructOfCoord const &master, StructOfCoord &local, unsigned int np ) const

virtual void LocalToMaster_v(StructOfCoord const &local, StructOfCoord &master, unsigned int np ) const

virtual void MasterToLocalVect_v(StructOfCoord const &master, StructOfCoord &local, unsigned int np ) const

virtual void LocalToMasterVect_v(StructOfCoord const &local, StructOfCoord &master, unsigned int np ) const

virtual void MasterToLocalCombined_v(StructOfCoord const &master, StructOfCoord &local, unsigned int np ) const

virtual void LocalToMasterCombined_v(StructOfCoord const &local, StructOfCoord &master, unsigned int np ) const

