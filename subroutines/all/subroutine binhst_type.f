      subroutine binhst_type(ihist,c,ibox,istuner,ieof,iout,
     -  maxat)
      dimension c(3,maxat)
c     Determine if MMC binary trajectory file has box and/or tuning info
      real*8 etoto
      dimension test5(5),test4(4),test3(3)
      character*4 chhd
c     print *,'BINHST_TYPE maxat=',maxat
c     Check for mistaken Charrm dcd
      read (ihist,end=888,err=888) chhd
      if (chhd .eq. 'CORD' .or. chhd .eq. 'VELD') then
        print *,'Trajectory is a Charmm DCD file'
        stop
      end if
      rewind ihist
      ieof=0
      istuner=0
      ibox=0
      call binhst_read(ihist,c,etoto,nmc,nwatr,ieof,0,maxat)
      read (ihist,end=100,err=100) test5
c     New header was found
      rewind ihist
      return
c     Tuning and/or box record found
100   rewind ihist
      call binhst_read(ihist,c,etoto,nmc,nwatr,ieof,0,maxat)
      read (ihist,end=200,err=200) test4
c     Tuning record was found - no box info
      istuner=1
      rewind ihist
      return
200   rewind ihist
c     There must be box information
      ibox=1
      call binhst_read(ihist,c,etoto,nmc,nwatr,ieof,0,maxat)
      read (ihist,end=888,err=888) test3
      read (ihist,end=300,err=300) test5
c     No tuning info
      rewind ihist
      return
300   rewind ihist
c     There is be tuning iformation too information
      istuner=1
888   ieof=1
      if (iout .gt. 0)
     -   write(iout,*) 'First box info record is in error'
      return
      end
