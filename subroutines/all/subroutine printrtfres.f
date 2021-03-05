      subroutine printrtfres(if,il,nneig,ineig,resnam,q,atnam,potnam,
     -  nglist,iqok,ipotok,iout,outfile,namleno,maxng,maxat)
      dimension nneig(maxat),ineig(maxng,maxat),q(maxat),nglist(maxng)
      character*4 resnam,atnam(maxat),potnam(maxat)
      character*(*) outfile
      if (iqok .eq. 0) write (iout,1000) 'Charges were unavailable'
      if (ipotok .eq. 0) write (iout,1000) 'FP labels were unavailable'
c-----Print residue name, charge
      qsum=0.0
      qabssum=0.0
      if (iqok .eq. 1) then
        do ia=if,il
          qsum=qsum+q(ia)
          qabssum=qabssum+abs(q(ia))
        end do
      end if
      if (qabssum .eq. 0.0) write (iout,1000) 'All charges are zero ??'
c-----Print atoms
      write (iout,1003) resnam,qsum
      qsum=0.0
      do ia=if,il
        write (iout,1002) atnam(ia),potnam(ia),q(ia)
        qsum=qsum+q(ia)
        if (isinteger(qsum,4,qoffset) .eq. 1 .and. qabssum .gt. 0.001)
     -     write (iout,1004)
      end do
c-----Print bonds
      do ia=if,il
        nn=nneig(ia)
        if (nn .gt. 0) then
          nnl=0
          do ja=1,nn
            in=ineig(ja,ia)
            if (in .gt. ia) then
              nnl=nnl+1
              nglist(nnl)=in
            end if
          end do
          if (nnl .gt. 0) write (iout,1001)
     -      (atnam(ia),atnam(nglist(j)),j=1,nnl)
        end if
      end do
      write (6,2000) resnam,outfile(1:namleno),il-if+1,qsum
      return
1000  format('! ',a)
1001  format('BOND  ',12(a4,1x))
1002  format('ATOM ',a4,1x,a4,f10.5)
1003  format('RESI ',a4,f10.5)
1004  format('GROUP')
2000  format(' RTF file for residue ',a,' was printed a to file ',a,/,
     -  ' Number of atoms=',i4,' Charge sum=',f6.2)
      end
