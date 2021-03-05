      subroutine findbestcorrep(iout,ifcl,ilcl,index,mx2d)
      dimension index(mx2d)
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-
     -  (4*MAXBONDS+MAX2D+2*MAX2D*MAX2D))
      common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
     -  ihbpair(2,MAXBONDS),a1(MAXBONDS),a2(MAXBONDS),fill(IFILL2)
      ravmax=100000.0
      ixavmax=0
      rminmax=100000.0
      ixminmax=0
      nmem=ilcl-ifcl+1
      rxymax=0.0
      ixmax=0.
      iymax=0.
      do ixx=ifcl,ilcl
        ix=index(ixx)
        r2av=0.0
        r2max=0.0
        do iyy=ifcl,ilcl
          if (iyy .ne. ixx) then
            iy=index(iyy)
            r2=rmsd2d(ix,iy)
            r2av=r2av+r2
            if (r2 .gt. r2max) r2max=r2
            if (r2 .gt. rxymax) then
              rxymax=r2
              ixmax=ix
              iymax=iy
            end if
          end if
        end do
        if (r2max .lt. rminmax) then
          rminmax=r2max
          ixminmax=ix
        end if
        r2av=r2av/nmem
        if (r2av .lt. ravmax) then
          ravmax=r2av
          ixavmax=ix
        end if
      end do
      write (iout,1001) 'max',ixminmax,'max',rminmax
      write (iout,1001) 'avg',ixavmax,'avg',ravmax
      write (iout,1002) rxymax,ixmax,iymax
      write (iout,1003) ixminmax,ixavmax,rmsd2d(ixavmax,ixminmax)
      write (iout,1004)
      return
1001  format(' Cluster center based on the lowest ',a,' distance: ',i5,
     -  1x,a,' distance=',f8.3)
1002  format(' Cluster diameter=',f8.2,' (distance between ',i6,
     -  ' and ',i6,')')
1003  format(' Distance between the two center estimates',
     -  ' (',i5,',',i5,') =',f8.2)
1004  format(' NOTE: numbers above refer to the unfiltered bond number')
      end
