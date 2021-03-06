      subroutine ramachandran(c,n,index,line,nconfig,pi,nresfound,
     -  iresno,ixres,ir1,ir2,ic1,iw0,maxng,mxrec)
      dimension c(3,n),index(n),iresno(n),ixres(n)
      character* 132 line(mxrec)
      parameter (MAXFRAMES=50000,MAXCOPY=600)
      common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
     -  increment,increment2,res(2,MAXFRAMES,MAXCOPY),
     -  x0(MAXCOPY),y0(MAXCOPY),nxselres,ixselres(MAXCOPY)
      parameter (MAXREC=200000,MAXRSD=70000,MAXNEIG=70,MAXPHI=400)
      parameter (IFILL1=MAXPHI*MAXPHI*MAXPHI-(7+2*MAXNEIG)*MAXREC)
      parameter (IFILL5=(MAXNEIG+6)*MAXREC+IFILL1-44*MAXRSD)
      common /nnwork/ ineig(MAXNEIG,MAXREC),nneig(MAXREC),
     -  nprossacc(6,6,MAXRSD),issprossacc(5,MAXRSD),
     -  ixypross(2,MAXRSD),isspross(MAXRSD),fill(IFILL5)
      character*1 prosscode(5)
      character*4 atnam
      data prosscode /'H','E','T','P','C'/
c     C - N - CA - C - N
c     i1  i2  i3   i4  i5
c         <--ires-->
c     print *,'RAMA n,maxng,mxrec=',n,maxng,mxrec
      call trajlimtest(ixres(n),MAXFRAMES)
      nresfound=0
      do ires=1,ixres(n)
        res(1,ires,maxpres)=999.9
        res(2,ires,maxpres)=999.9
      end do
      do ia=1,n
        atnam=line(index(ia))(ic1:ic1+3)
        call leftadjust4(atnam,atnam)
        if (atnam .eq. 'CA  ') then
c         Alpha carbon found
          i3=ia
          call ca_to_bb(i3,iresno,nneig,ineig,index,line,ic1,
     -      i1,i2,i4,i5,ires,iprotein,maxng,mxrec)
          if (iprotein .gt. 0) then
c           Full residue backbone found
c           write (iw0,*)'i1-5=',i1,i2,i3,i4,i5
c           write (iw0,*)'ires i1-5=',iresno(i1),iresno(i2),
c    -        iresno(i3),iresno(i4),iresno(i5)
            phirad=dihangl(c,i1,i2,i3,i4,0,mxrec)
            psirad=dihangl(c,i2,i3,i4,i5,0,mxrec)
            phi=phirad*(180.0/pi)
            psi=psirad*(180.0/pi)
            if (nframe .eq. 0) then
              write (iw0,1000) phi,psi,ires,line(index(i3))(ir1:ir2)
            else
              write (iw0,1000) phi,psi,ires,line(index(i3))(ir1:ir2),
     -          ' ',nframe
              if (nxselres .gt. 0) then
c               Save dial plot info
                do ix=1,nxselres
                  if (ires .eq. ixselres(ix)) then
                    call trajlimtest(nframe,MAXFRAMES)
                    res(1,nframe,2*ix-1)=cos(phirad)
                    res(2,nframe,2*ix-1)=sin(phirad)
                    res(1,nframe,2*ix)=cos(psirad)
                    res(2,nframe,2*ix)=sin(psirad)
                  end if
                end do
              end if
            end if
            nresfound=nresfound+1
            call trajlimtest(nresfound,MAXFRAMES)
            res(1,nresfound,maxpres)=phi
            res(2,nresfound,maxpres)=psi
c            write (iw0,2000) (c(k,i1),k=1,3),(c(k,i2),k=1,3),
c     -       (c(k,i3),k=1,3),(c(k,i4),k=1,3),(c(k,i5),k=1,3)
c2000         format(' c1=',3f8.3,' c2=',f8.3,' c3=',f8.3,' c4=',3f8.3,
c     -         ' c5=',f8.3)
c           Save and accumulate PROSS indices information
            ixypross(1,nresfound)=ixpross(phi)
            ixypross(2,nresfound)=ixpross(psi)
            nprossacc(ixypross(1,nresfound),
     -                ixypross(2,nresfound),nresfound)=
     -        nprossacc(ixypross(1,nresfound),
     -                  ixypross(2,nresfound),nresfound)+1
          end if
        end if
      end do
c     Accumulate PROSS indices information
      call zeroiti(isspross,0,nresfound)
c     Look for Alpha
      ir=1
      do while (ir .le. nresfound)
        ir0=ir
        ixy=ixypross(1,ir)+6*(ixypross(2,ir)-1)
        if (ixy .eq. 8 .or. ixy .eq. 14) then
          irr=ir
          do while (ixy .eq. 8 .or. ixy .eq. 14)
            irr=irr+1
            ixy=0
            if (irr .lt. nresfound)
     -        ixy=ixypross(1,irr)+6*(ixypross(2,irr)-1)
          end do
          if (irr-ir .ge. 4) then
            do i=ir,irr
              isspross(i)=1
            end do
          end if
          ir=irr+1
        end if
        if (ir .eq. ir0) ir=ir+1
      end do
c     Look for Beta
      ir=1
      do while (ir .le. nresfound)
        ir0=ir
        if (isspross(ir) .eq. 0) then
          ixy=ixypross(1,ir)+6*(ixypross(2,ir)-1)
          if (ixy .eq. 13 .or. ixy .eq. 24 .or. ixy .eq. 29 .or.
     -        ixy .eq. 30 .or. ixy .eq. 28 .or. ixy .eq. 36) then
            irr=ir
            do while (ixy .eq. 13 .or. ixy .eq. 24 .or. ixy .eq. 29 .or.
     -        ixy .eq. 30 .or. ixy .eq. 28 .or. ixy .eq. 36)
              irr=irr+1
              ixy=0
              if (irr .lt. nresfound)
     -          ixy=ixypross(1,irr)+6*(ixypross(2,irr)-1)
            end do
            if (irr-ir .ge. 2) then
              do i=ir,irr
                isspross(i)=2
              end do
            end if
            ir=irr+1
          end if
        end if
        if (ir .eq. ir0) ir=ir+1
      end do
c     Look for turn
      ir=1
      do while (ir .le. nresfound-1)
        ir0=ir
        if (isspross(ir) .eq. 0 .and. isspross(ir+1) .eq. 0) then
          ixy1=ixypross(1,ir)+6*(ixypross(2,ir)-1)
          ixy2=ixypross(1,ir+1)+6*(ixypross(2,ir+1)-1)
          if ((ixy1 .eq.  8 .and. (ixy2 .eq. 8 .or. ixy2 .eq. 14 .or.
     -        ixy2 .eq. 13)) .or. (ixy1 .eq.  14 .and. (ixy2 .eq. 8 .or.
     -        ixy2 .eq. 14 .or. ixy2 .eq. 13)) .or. (ixy1 .eq. 13 .and.
     -        (ixy2 .eq.  8 .or. ixy2 .eq. 14 .or. ixy2 .eq.  13)) .or.
     -        (ixy1 .eq. 30 .and. (ixy2 .eq. 20 .or. ixy2 .eq. 16 .or.
     -        ixy2 .eq. 17)) .or. (ixy1 .eq. 24 .and. (ixy2 .eq. 20 .or.
     -        ixy2 .eq. 16 .or. ixy2 .eq. 17)) .or. (ixy1 .eq. 20 .and.
     -        (ixy2 .eq. 20 .or. ixy2 .eq. 16 .or. ixy2 .eq.  17)) .or.
     -        (ixy1 .eq. 16 .and. (ixy2 .eq. 20 .or. ixy2 .eq. 16 .or.
     -        ixy2 .eq. 17)) .or. (ixy1 .eq. 17 .and. (ixy2 .eq. 20 .or.
     -        ixy2 .eq. 16 .or. ixy2 .eq. 17)) .or. (ixy1 .eq. 32 .and.
     -        (ixy2 .eq. 8  .or. ixy2 .eq. 14 .or. ixy2 .eq.  13)) .or.
     -        (ixy1 .eq.  4 .and. (ixy2 .eq. 8 .or. ixy2 .eq. 14 .or.
     -        ixy2 .eq. 13))) then
            isspross(ir)=3
            isspross(ir+1)=3
            ir=ir+2
          end if
        end if
        if (ir .eq. ir0) ir=ir+1
      end do
c     Look for Pi
      ir=1
      do while (ir .le. nresfound)
        if (isspross(ir) .eq. 0) then
          ixy=ixypross(1,ir)+6*(ixypross(2,ir)-1)
          if (ixy .eq. 30 .or. ixy .eq. 24) isspross(ir)=4
        end if
        ir=ir+1
      end do
c     Rest is Coil
      do ir=1,nresfound
        if (isspross(ir) .eq. 0) isspross(ir)=5
      end do
      do ir=1,nresfound
        if (isspross(ir) .lt. 1 .or. isspross(ir) .gt. 5)
     -    print *,'ir=',ir,' isspross=',isspross(ir)
        issprossacc(isspross(ir),ir)=issprossacc(isspross(ir),ir)+1
      end do
      write (iw0,1002) (prosscode(isspross(ir)),ir=1,nresfound)
      if (nresfound .eq. 0) then
        print *,'ERROR: no protein backbone found'
        stop
      else if (nconfig .le. 1) then
        write (6,1001) nresfound
      end if
      return
1000  format(2f10.3,' (phi and psi angles) residue#=',i6,' (',a,')',
     -  a,'Frame=',i5)
1001  format(' Number of protein residues found=',i5)
1002  format(' PROSS SS codes: ',50a1)
c1003  format(' WARNING: dial-plot residue #',i3,', residue',i6,' is ',
c     -  'not a protein',/,10x,'No dials will be plotted for it')
      end
