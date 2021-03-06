      subroutine strip_cext(analfile,namleni,namleno,lenext)
      character*(*) analfile
      character*4 extinp
      character*6 extinp6
      namleno=namleni
c     Remove coordinate extension
      extinp=analfile(namleno-3:namleno)
      if (extinp .eq. '.pdb' .or. extinp .eq. '.PDB' .or.
     -    extinp .eq. '.crd' .or. extinp .eq. '.CRD' .or.
     -    extinp .eq. '.dat' .or. extinp .eq. '.DAT' .or.
     -    extinp .eq. '.mae' .or. extinp .eq. '.flt' .or.
     -    extinp .eq. '.slt') then
        lenext=4
      else if (extinp .eq. 'mol2') then
        lenext=5
      else
        extinp6=analfile(namleno-5:namleno)
        if (extinp6 .eq. '.pdbqs' .or. extinp6 .eq. '.PDBQS' .or.
     -      extinp6 .eq. '.pdbqt' .or. extinp6 .eq. '.PDBQT')
     -    lenext=6
      end if
      namleno=namleno-lenext
      return
      end
