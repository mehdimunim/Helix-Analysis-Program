      subroutine writeconn(ioutpdb,ifres,ilres,line,index,inamcol1,
     -  inamcol2,c,n,iresbondcorr,nrrbond,nframe,iout,mxbonds,maxrsd,
     -  maxrec)
      dimension ifres(maxrsd),ilres(maxrsd),index(maxrec),c(3,maxrec),
     -  nrrbond(mxbonds)
      character* 132 line(maxrec)
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-
     -  (6*MAXBONDS+MAX2D+2*MAX2D*MAX2D))
      common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
     -  ihbtores(MAXBONDS),nusepair(MAXBONDS),nhb_atot(MAXBONDS),
     -  nhb_rtot(MAXBONDS),ibond(MAXBONDS),a1(MAXBONDS),fill(IFILL2)
      common /bondpairs/ ihbpair(2,MAXBONDS),ihb_pair_res(3,MAXBONDS)
      dimension cHe(3)
      lname=inamcol2-inamcol1+1
      if (lname .lt. 1) then
        print *,'PROGRAM ERROR inamcol1,2=',inamcol1,inamcol2
        stop
      end if
      nnoCA=0
      nrr=1
      nrw=0
      corrmax=0.0
      corrmin=1.0
      iconncorr=0
      iconnacorr=0
      if (iresbondcorr .eq. 1) then
        call askyn('Do you want to connect correlated bonds',39,1,-1,
     -    iconncorr,00,0)
        if (iconncorr .eq. 1) then
          call getreal(
     -      'Maximum (correlation) distance measure',38,1.0,corrmax,1,0)
          write (iout,2003) ' ','<',corrmax
          write (ioutpdb,2003) 'REMARK ','<',corrmax
        else
          call askyn('Do you want to connect anti-correlated bonds',44,
     -      1,-1,iconnacorr,00,0)
          if (iconnacorr .eq. 1) then
            call getreal(
     -        'Minimum (correlation) distance measure',38,1.0,corrmin,
     -        1,0)
            write (iout,2003) ' ','>',corrmin
            write (ioutpdb,2003) 'REMARK ','>',corrmin
          end if
        end if
      end if
      if (iconncorr+iconnacorr .gt. 0) then
        write (ioutpdb,1000) 'TER'
        do while (ihb_pair_res(3,nrr) .gt. 0 .and. nrr .le. MAXBONDS)
          irr1=ihb_pair_res(1,nrr)
          call findCA(line,index,ifres(irr1),ilres(irr1),inamcol1,
     -      inamcol2,iCA1,nnoCA,maxrec)
          irr2=ihb_pair_res(2,nrr)
          call findCA(line,index,ifres(irr2),ilres(irr2),inamcol1,
     -      inamcol2,iCA2,nnoCA,maxrec)
          do k=1,3
            cHe(k)=(c(k,iCA1)+c(k,iCA2))/2.0
          end do
c          write (6,9782) nrr,irr1,irr2,iCA1,iCA2,nusepair(nrr)
c9782      format(i4,' irr1,2=',2i5,' iCA1,2=',2i6,' nusepair=',i2)
          nrw=nrw+1
          write (ioutpdb,1001) n+nrw,nrw,cHe,
     -      float(nrrbond(nrw))/float(nframe)
          nrr=nrr+1
        end do
        write (ioutpdb,1000) 'TER'
c       Connect bonds where the correlation distance measure is > corrmax
c       or < corrmax
        do ir1=1,nrw
          do ir2=ir1+1,nrw
            if (rmsd2d(ir1,ir2) .gt. corrmin .or.
     -        rmsd2d(ir1,ir2) .lt. corrmax)
     -          write (ioutpdb,2001) n+ir1,n+ir2
          end do
        end do
        if (nnoCA .gt. 0) write (6,2000) nnoCA
        if (nnoCA .gt. 0) write (iout,2000) nnoCA
      else
        do while (ihb_pair_res(3,nrr) .gt. 0 .and. nrr .le. MAXBONDS)
          nrr=nrr+1
        end do
      end if
      write (ioutpdb,1000) 'REMARK CA - CA bonds'
      do irr=1,nrr
        irr1=ihb_pair_res(1,irr)
        call findCA(line,index,ifres(irr1),ilres(irr1),inamcol1,
     -    inamcol2,iCA1,nnoCA,maxrec)
        irr2=ihb_pair_res(2,irr)
        call findCA(line,index,ifres(irr2),ilres(irr2),inamcol1,
     -    inamcol2,iCA2,nnoCA,maxrec)
        write (ioutpdb,2001) iCA1,iCA2
        write (iout,2004) irr1,irr2,iCA1,iCA2
      end do
      write (iout,2002) nrr
      return
1000  format(a)
1001  format('ATOM',i7,' He   HE  ','X',i4,4x,3f8.2,1x,f5.3,' 0.000',
     -  11x,'He')
2000  format(' NOTE:',i4,' residues had no CA atom - first atom in ',
     -  'the residue is used')
2001  format('CONECT',2i5)
2002  format(' Number of CA-CA bonds created=',i4)
2003  format(a,' Bonds created between CA-CA pairs with correlation ',
     -  'measure ',a,f6.2)
2004  format(' CA - CA bonds between residues',i6,' and',i6,' (atom ',
     -  'numbers',2i7)
      end
