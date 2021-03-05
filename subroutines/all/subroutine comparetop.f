      subroutine comparetop(c,n,nneig,ineig,iatnum,innlist,nslt,nslv,
     -  cell,ncell,ioppbc,maxng,mxat)
      dimension c(3,mxat),nneig(mxat),ineig(maxng,mxat),iatnum(mxat),
     -  cell(3,ncell)
      common /connatdat/ ramax(99),ramax2(99),hlimfac,ianfg(99),
     -  namfcg(100),nrmw
c     Check if list of neighbors actually represent bonded pairs
c     print *,'COMPARETOP n,nslt,nslv,ncell=',n,nslt,nslv,ncell
      if (innlist .lt. 0) return
      if (innlist .eq. 0) then
        print *,'PROGRAM ERROR: no neighbour list in topology check'
        return
      end if
      discrepmax=0.0
      nfail=0
      do ia=1,n
        do jja=1,nneig(ia)
          ja=ineig(jja,ia)
          if (ja .gt. ia) then
            r2=dist2(c(1,ia),c(1,ja))
            rlm=amax1(ramax2(iatnum(ia)),ramax2(iatnum(ja)))*1.2
            if (iatnum(ia) .eq. 1 .or. iatnum(ja) .eq. 1)
     -         rlm=rlm*hlimfac
            if (r2 .gt. rlm) then
              nfail=nfail+1
              if (nfail .le. 25) write (6,1000) ia,ja,sqrt(r2)
              if (r2 .gt. discrepmax) discrepmax=r2
            end if
          end if
        end do
      end do
      if (nfail .gt. 0) write (6,1001) nfail,sqrt(discrepmax)
      percout=0.0
      if (ioppbc .gt. 0 .and. n-nslt .gt. 100) then
c       See how many solvents are outside the cell
        nout=0
        numslv=(n-nslt)/nslv
        do iw=1,numslv
          call genimdist(c(1,nslt+(iw-1)*nslv+1),cell,1,ncell,img,d2)
          if (img .gt. 1) nout=nout+1
        end do
        if (nout .gt. 0) write (6,1002) nout
        percout=float(100*nout)/float(numslv)
        if (percout .gt. 2.0)
     -    print *,'The PBC cell may have been incorrectly specified'
      end if
      if (nfail .gt. 0 .or. percout .gt. 2.0) then
        call askstop(1)
      else
        print *,'First frame topology and solvent PBC checks passed'
      end if
      return
1000  format(' Atoms',i6,' and ',i6,' are bonded in the input ',
     -  'structure but ',/,6x,'they are ',f6.2,' A apart in the ',
     -  'first frame of the trajectory')
1001  format(' There were ',i6,' discrepancies',/,
     -  ' The largest discrepancy=',f8.2,' A',/,
     -  ' Large discrepancies usually arise when the structure file',/,
     -  ' does not correspond to the trajectory.')
1002  format(i6,' solvents appear to be outside the periodic cell')
      end
