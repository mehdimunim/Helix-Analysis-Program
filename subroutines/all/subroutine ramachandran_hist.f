      subroutine ramachandran_hist(nres,resnames,nrescol,iw,
     -  idistw,isummw)
      character*(*) resnames(nres)
      parameter (MAXREC=200000,MAXRSD=70000,MAXNEIG=70,MAXPHI=400)
      parameter (IFILL1=MAXPHI*MAXPHI*MAXPHI-(7+2*MAXNEIG)*MAXREC)
      parameter (IFILL5=(MAXNEIG+6)*MAXREC+IFILL1-44*MAXRSD)
      common /nnwork/ ineig(MAXNEIG,MAXREC),nneig(MAXREC),
     -  nprossacc(6,6,MAXRSD),issprossacc(5,MAXRSD),
     -  ixypross(2,MAXRSD),isspross(MAXRSD),fill(IFILL5)
      character*1 prosscode(5)
      character*1 prossixy(36)
      dimension iang(6)
      data prosscode /'H','E','T','P','C'/
      data prossixy /'H','N','T','r','l','B', 'I','O','U','q','k','C',
     -               'J','P','V','p','j','D', 'K','Q','W','o','i','E',
     -               'L','R','X','n','h','F', 'G','M','S','m','g','A'/
c     print *,'RAMA hist nres=',nres
      if (idistw .gt. 0) then
        write (iw,*)
        do i=1,6
          iang(i)=-120+(i-1)*60
        end do
        do i=1,6
          do j=1,6
            ixy=i+6*(j-1)
            write (iw,1001) iang(i),iang(j),prossixy(ixy)
          end do
        end do
        write (iw,*)
        do ir=1,nres
          write (iw,1002) ir,resnames(ir)(1:nrescol),(iang(j),
     -      (nprossacc(i,j,ir),i=1,6),j=6,1,-1)
        end do
      end if
      if (isummw .gt. 0) then
        do ir=1,nres
          write (iw,1000) ir,resnames(ir)(1:nrescol),
     -      (prosscode(i),issprossacc(i,ir),i=1,5)
        end do
      end if
      return
1000  format(' Res',i5,' (',a,') PROSS distr:',5(1x,a1,':',i6))
1001  format(' PROSS SS codes and angles:',(' phi=',i4,' psi=',i4,
     -  ' code=',a1))
1002  format(' PROSS angular grid distribution for residue ',i5,
     -  ' (',a,'):',6(/,' psi=',i4,':',6i9))
      end
