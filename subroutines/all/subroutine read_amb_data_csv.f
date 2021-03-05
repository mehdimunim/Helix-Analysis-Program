      subroutine read_amb_data_csv(data,err,ncol,ix1,ix2,resnames,
     -  nres,ipairs,ierr_mean,inpt,mxcol,mxrsd,mxrec)
      dimension  data(mxcol,mxrec),err(mxcol,mxrec),ix1(mxrec),
     -  ix2(mxrec)
      character*8 resnames(mxrsd)
      character*3 rn1,rn2
      character*1000 line
c     Read data from Amber energy analysis tables
c     print *,'READ_AMB_DATA_CSV ipairs,ncol=',ipairs,ncol
      len=1000
      nres=0
      lc=len
      do while (lc .gt. 1)
        call blankout(line,1,len)
        read (inpt,1000,end=999) line
        call lastchar(line,lc,len)
        if (lc .gt. 1) then
          nres=nres+1
c         write (6,*) 'NRES=',nres,' LINE=',line(1:20)
          if (nres .gt. mxrec) then
            write (6,1001) mxrec
            stop
          end if
          ic=1
          rn1=line(ic:ic+2)
          call findnextchar(',',line,ic,len)
          read (line(ic-4:ic-1),*) ix1(nres)
          ic=ic+1
          ic1=ic
          if (ipairs .eq. 0) then
            call findnextchar(',',line,ic,len)
            call blankout(resnames(ix1(nres)),1,8)
            resnames(ix1(nres))(1:3)=rn1
          else
            rn2=line(ic:ic+2)
            call findnextchar(',',line,ic,len)
            read (line(ic-4:ic-1),*) ix2(nres)
            call blankout(resnames(ix1(nres)),1,8)
            call blankout(resnames(ix2(nres)),1,8)
            resnames(ix1(nres))(1:3)=rn1
            resnames(ix2(nres))(1:3)=rn2
c           write (77,8977) nres,rn1,rn2,ix1(nres),ix2(nres),
c    -       resnames(ix1(nres)),resnames(ix2(nres))
c8977        format(' nres=',i4,' rn1,2=',a,'|',a,'| ix1,2=',2i5,
c    -         ' rn1,2=',a,'|',a,'|')
          end if
          ic=ic+1
          do nc=1,ncol
            call getnextcsv(line,i,data(nc,nres),2,ic,len)
            call getnextcsv(line,i,sd,2,ic,len)
            call getnextcsv(line,i,err_mean,2,ic,len)
            if (ierr_mean .eq. 1) err(nc,nres)=sd
            if (ierr_mean .eq. 2) err(nc,nres)=err_mean
          end do
c         write (77,8877) ncol,nc,nres,(data(i,nres),i=1,ncol)
c8877     format(' ncol,nc,nres=',i2,i3,i5,' data=',10f10.3)
        end if
      end do
999   return
      stop
1000  format(a)
1001  format(' ERROR: Number of residue-residue records exceeds ',
     -  'the limit (',i9,')',/,8x,'Redimension with larger MAXREC (and',
     -  ' possibly, larger MAXPHI')
      end
