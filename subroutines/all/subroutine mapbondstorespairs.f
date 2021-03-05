      subroutine mapbondstorespairs(nbfound,nbfoundorig,nbresfound,
     -  nbresorig,nframe,ixres,indexres,bondname,lbondname,line,index,
     -  irescol1,irescol2,nresslt,percmind,percmaxd,minresdistd,
     -  maxresdistd,percmin,percmax,minresdist,maxresdist,ifres,iresno,
     -  isegno,nochange,indexbond,ifa,ila,itemp1,itemp2,itemp3,iout,
     -  mxbonds,maxrsd,maxrec)
      dimension index(maxrec),ixres(maxrec),indexbond(maxrsd),
     -  indexres(mxbonds),ifres(maxrec),iresno(maxrec),isegno(maxrec),
     -  ifa(mxbonds),ila(mxbonds),itemp1(maxrec),itemp2(maxrec),
     -  itemp3(maxrec)
      character* 132 line(maxrec)
      character*(*) bondname
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-
     -  (6*MAXBONDS+MAX2D+2*MAX2D*MAX2D))
      common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
     -  ihbtores(MAXBONDS),nusepair(MAXBONDS),nhb_atot(MAXBONDS),
     -  nhb_rtot(MAXBONDS),ibond(MAXBONDS),a1(MAXBONDS),fill(IFILL2)
      common /bondpairs/ ihbpair(2,MAXBONDS),ihb_pair_res(3,MAXBONDS)
c     print *,'MAPBONDSTORESPAIRS NBFOUND,NBFOUNDORIG,NBRESORIG=',
c    -  nbfound,nbfoundorig,nbresorig
      call indexit(indexbond,1,nbfoundorig,0)
      write (6,1005)
      write (iout,1007)
      write (iout,1022) nbresorig,bondname(1:lbondname)
      call askyn(
     -  'Do yo want to filter the residue-aggregated bonds',49,1,1,
     -  ifiltrr,83,0)
      if (ifiltrr .eq. 1) then
        call getfiltlims(percmind,percmaxd,minresdistd,maxresdistd,
     -    percmin,percmax,minresdist,maxresdist,nresslt,nochange,
     -    bondname,lbondname,' res-res ',9,iout)
        write (iout,1020) 'filtered out'
        ndel=0
        do irr=1,nbresorig
          ir1=ihb_pair_res(1,irr)
          ir2=ihb_pair_res(2,irr)
          perc=100.0*float(ihb_pair_res(3,irr))/float(nframe)
          irdiff=iabs(ir1-ir2)
          if (perc .lt. percmin .or. perc .gt. percmax .or.
     -      irdiff .lt. minresdist .or. irdiff .gt. maxresdist) then
            ndel=ndel+1
            ia1=ifres(ir1)
            ia2=ifres(ir2)
            write (iout,1002) irr,iresno(ia1),
     -        line(index(ia1))(irescol1:irescol2),isegno(ia1),
     -        iresno(ia2),line(index(ia2))(irescol1:irescol2),
     -        isegno(ia2),perc
          else
            ihb_pair_res(1,irr-ndel)=ihb_pair_res(1,irr)
            ihb_pair_res(2,irr-ndel)=ihb_pair_res(2,irr)
            ihb_pair_res(3,irr-ndel)=ihb_pair_res(3,irr)
          end if
        end do
        nbresfound=nbresorig-ndel
        write (6,1011) nbresfound,nbresorig
        write (iout,1011) nbresorig,nbresfound
        do irr=1,ndel
          ihb_pair_res(3,nbresfound+irr)=0
        end do
        ndel=0
        do ib=1,nbfoundorig
          ihb=indexbond(ib)
          ir1=ixres(ihbpair(1,ihb))
          ir2=ixres(ihbpair(2,ihb))
          ifound=0
          do irr=1,nbresfound
            if ((ir1 .eq. ihb_pair_res(1,irr) .and.
     -           ir2 .eq. ihb_pair_res(2,irr)) .or.
     -          (ir2 .eq. ihb_pair_res(1,irr) .and.
     -           ir1 .eq. ihb_pair_res(2,irr))) then
              ifound=1
              go to 100
            end if
          end do
100       if (ifound .eq. 0) then
            ndel=ndel+1
          else
            indexbond(ib-ndel)=indexbond(ib)
          end if
        end do
        nbfound=nbfoundorig-ndel
        write (6,1023) bondname(1:lbondname),nbfoundorig,nbfound
        write (iout,1023) bondname(1:lbondname),nbfoundorig,nbfound
      else
        nbfound=nbfoundorig
      end if
      nbuse=nbfound
      call zeroiti(ihbtores,0,nbfoundorig)
      call zeroiti(nusepair,0,nbresorig)
      call indexit(itemp1,1,nbresorig,0)
      do ibb=1,nbuse
        ib=indexbond(ibb)
        ir1=ixres(ihbpair(1,ib))
        ir2=ixres(ihbpair(2,ib))
        do irr=1,nbresfound
          if (ihb_pair_res(1,irr) .eq. ir1) then
            if (ihb_pair_res(2,irr) .eq. ir2) then
              ihbtores(ib)=irr
              nusepair(irr)=nusepair(irr)+1
              go to 200
            end if
          end if
        end do
        if (nochange .eq. 1) then
          write (6,1000) (ihbpair(k,ib),k=1,2),ir1,ir2
          write (iout,1000) (ihbpair(k,ib),k=1,2),ir1,ir2
        end if
200     continue
      end do
      do irr=1,nbresfound
        itemp2(irr)=ihb_pair_res(1,irr)
      end do
      call indexit(indexres,1,nbresfound,0)
      call mrgsrti(6,indexres,itemp2,nbresfound,ifa,ila,itemp1,itemp3,
     -  mxbonds)
      write (iout,1020)
     -  'kept - sequence numbers also refer to the bond track #'
      do ipp=1,nbresfound
        ip=indexres(ipp)
        ir1=ihb_pair_res(1,ip)
        ir2=ihb_pair_res(2,ip)
        ia1=ifres(ir1)
        ia2=ifres(ir2)
        perc=100.0*float(ihb_pair_res(3,ip))/float(nframe)
        write (iout,1002) ipp,iresno(ia1),
     -    line(index(ia1))(irescol1:irescol2),isegno(ia1),iresno(ia2),
     -    line(index(ia2))(irescol1:irescol2),isegno(ia2),perc,' ',
     -    nusepair(ipp)
      end do
      write (iout,*)
c      write (iout,8711) (ihbtores(i),i=1,nbuse)
c8711  format(' IHBTORES:',/,(15i5))
c      write (6,8712) (irrix(i),i=1,nbresfound)
c8712  format(' IRRIX:',/,(15i5))
      return
1000  format(' PROGRAM ERROR: bond between atoms',i6,i7,' (residue ',
     -  'indices',i5,i6,')'/,' is not found in the residue pair list')
1002  format(i6,' Bonded % of residues',i6,' (',a,') C/S',i2,' and',i6,
     -  ' (',a,') C/S',i2,'=',f5.1,a,' Nbonds=',i3)
1005  format(/,' Generating residue-contracted bond information')
1007  format(/,' === Residue-contracted bond information:')
1011  format(' Number of residue-residue pairs filtered down from',i6,
     -  ' to',i5)
1020  format(/,' List of residue-residue pairs ',a)
1022  format(' There are ',i4,' residue-residue pairs with ',a,' bonds')
1023  format(' After res-res filtering the number of ',a,' bonds',/,
     -  ' changed from',i6,' to',i5)
      end
