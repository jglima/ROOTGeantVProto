ClassNames = [ "TGeoMatrix",
               "TGeoTranslation",                
	       "TGeoRotation", 
               "TGeoScale", 
               "TGeoCombiTrans",
	       "TGeoIdentity"
	      ]


def writeHeader( matrixtype, filehandle ):
    filehandle.write("virtual void         MasterToLocal_l(Double_t const *master, Double_t *local, unsigned int np) const {\n"+
     "for( unsigned int k=0; k< np; ++k)\n"+
     "  {\n"+
	" this->"+matrixtype+"::MasterToLocal(&master[3*k], &local[3*k]);\n"+
     "  }\n"+
   "}\n\n"+
   "virtual void         MasterToLocalVect_l(Double_t const *master, Double_t *local, unsigned int np) const {\n"+
   "  for( unsigned int k=0; k< np; ++k)\n"+
   "   {\n"
   "	 this->"+matrixtype+"::MasterToLocalVect(&master[3*k], &local[3*k]);\n"+
   "   }\n"+
   "}\n\n"+
   "virtual void         MasterToLocalBomb_l(Double_t const *master, Double_t *local, unsigned int np) const {\n"+
   "  for( unsigned int k=0; k< np; ++k)\n"+
   "    {\n"+
   "	 this->"+matrixtype+"::MasterToLocalBomb(&master[3*k], &local[3*k]);\n"+
   "    }\n"+
   "}\n\n"+

   "virtual void         LocalToMaster_l(Double_t const *master, Double_t *local, unsigned int np) const {\n"+
   "  for( unsigned int k=0; k< np; ++k)\n"+
   "    {\n"+
   "	 this->"+matrixtype+"::LocalToMaster(&master[3*k], &local[3*k]);\n"+
   "    }\n"+
   "}\n\n"
   "virtual void         LocalToMasterVect_l(Double_t const *master, Double_t *local, unsigned int np) const {\n"+
   "  for( unsigned int k=0; k< np; ++k)\n"+
   "    {\n"+
   "	 this->"+matrixtype+"::LocalToMasterVect(&master[3*k], &local[3*k]);\n"+
   "    }\n"+
   "}\n\n"+
   "virtual void         LocalToMasterBomb_l(Double_t const *master, Double_t *local, unsigned int np) const {\n"+
   "  for( unsigned int k=0; k< np; ++k)\n"+
   "    {\n"+
   "	 this->"+matrixtype+"::LocalToMasterBomb(&master[3*k], &local[3*k]);\n"+
   "    }\n"+
   "}\n\n" + 

#   // define some real vectorized interfaces
   "virtual void MasterToLocal_v(StructOfCoord const &master, StructOfCoord &local, Int_t np ) const;\n\n"+
   "virtual void LocalToMaster_v(StructOfCoord const &local, StructOfCoord &master, Int_t np ) const;\n\n"+
   "virtual void MasterToLocalVect_v(StructOfCoord const &master, StructOfCoord &local, Int_t np ) const;\n\n"+
   "virtual void LocalToMasterVect_v(StructOfCoord const &local, StructOfCoord &master, Int_t np ) const;\n\n"+
   "virtual void MasterToLocalCombined_v(StructOfCoord const &master, StructOfCoord &local, StructOfCoord const &mastervec, StructOfCoord &localvec, Int_t np ) const;\n\n"+
   "virtual void LocalToMasterCombined_v(StructOfCoord const &local, StructOfCoord &master, StructOfCoord const &localvec, StructOfCoord &mastervec, Int_t np ) const;\n\n")


def main():
    for shape in ClassNames:
        shapefilename=shape+"_VectorAddons.inc"
        filehandle = open(shapefilename, 'w')
        writeHeader(shape, filehandle)
        filehandle.close()

if __name__ == "__main__":
    main()
