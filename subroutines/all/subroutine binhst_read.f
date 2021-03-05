      subroutine binhst_read(ihist,c,etoto,nmc,nwatr,ieof,iout,
     -  maxat)
      dimension c(3,maxat)
      real*8 etoto
c     Read header and coordinate record from a binary history file
      real*8 uusfacd,tesi(4)
c     print *,'BINHST_READ maxat=',maxat
      ieof=0
      read (ihist,end=888,err=888) nwatr,natomsr,uusfacd,nmc,nidmc,
     -  niaccp,ndaccp,ia0,ia1,etoto,tesi,cplpar
      nmolecr=nwatr+1
      nsvp=0
      if (nwatr .gt. 0) nsvp=(natomsr-(ia1-ia0+1))/nwatr
      if (nwatr .gt. 0) then
        read (ihist,end=999,err=999) ((c(k,ii),k=1,3),ii=ia0,ia1),
     -    (((c(k,ia1+(i-1)*nsvp+j),k=1,3),j=1,nsvp),i=1,nwatr)
      else
        read (ihist,end=999,err=999) ((c(k,ii),k=1,3),ii=ia0,ia1)
      end if
      return
888   ieof=1
      if (iout .gt. 0) write(iout,*) 'File header is truncated'
      return
999   ieof=2
      if (iout .gt. 0)
     -   write(iout,*) 'First coordinate record is in error'
      return
c2000  format(' nwat,nats=',2i11,' nmc=',i10,' ia0,ia1=',2i6,
c     -  ' nsv=',i2,/,
c     -  ' etot=',e13.5,' tesi=',4e12.5,/,' cplpar=',e12.5)
      end
