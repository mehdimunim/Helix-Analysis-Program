      subroutine mpxblist(n,i012,impx,rmpx,itemp,itempres,icc1,icc2,
     -  inc1,inc2,irc1,irc2,irnc1,irnc2,line,index,ixres,nhbdist,
     -  rhbdist,nbfound,nbresfound,nbonds,c,rmax,ifail,iout,mxbonds,
     -  mxrec)
      dimension index(n),ixres(n),i012(n),impx(n),rmpx(n),itemp(n),
     -  itempres(n),nhbdist(mxbonds),rhbdist(mxbonds),c(3,n)
      character*132 line(mxrec)
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-
     -  (6*MAXBONDS+MAX2D+2*MAX2D*MAX2D))
      common /nnwork/ rmsd2d(MAX2D,MAX2D),nng(MAX2D),ing(MAX2D,MAX2D),
     -  ihbtores(MAXBONDS),nusepair(MAXBONDS),nhb_atot(MAXBONDS),
     -  nhb_rtot(MAXBONDS),a1(MAXBONDS),a2(MAXBONDS),fill(IFILL2)
      common /bondpairs/ ihbpair(2,MAXBONDS),ihb_pair_res(3,MAXBONDS)
      parameter (MAXFRAMES=50000,MAXCOPY=600,MAXITEMS=2*MAXCOPY-2)
      common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
     -  increment,increment2,ires(MAXITEMS,MAXFRAMES),
     -  scres(2,MAXFRAMES),x0(MAXCOPY),y0(MAXCOPY),nxselres,
     -  ixselres(MAXCOPY)
c     print *,'MPXBLIST n,iout=',n,iout
      nbitmax=30*MAXITEMS
      if (nframe .eq. 0) write (iout,1000)
      if (nframe .gt. 0) write (iout,1000) ' ',nframe
      call zeroiti(itemp,0,n)
      call zeroiti(itempres,0,n)
      ifail=0
      nbonds=0
      do ia=1,n
        if (i012(ia) .eq. 1) then
          ja=impx(ia)
          if (ja .gt. 0) then
            if (impx(ja) .eq. ia) then
              dist_ia_ja=sqrt(dist2(c(1,ia),c(1,ja)))
              iskip=0
              if (rmax .lt. 9999.0) then
                if (dist_ia_ja .gt. rmax) iskip=1
              end if
              if (iskip .eq. 0) then
                call readint(line(index(ia)),irnc1,irnc2,irnia,2,1,ire)
                call readint(line(index(ja)),irnc1,irnc2,irnja,2,1,ire)
                write (iout,1001) line(index(ia))(icc1:icc2),ia,
     -            line(index(ia))(inc1:inc2),irnia,
     -            line(index(ia))(irc1:irc2),line(index(ja))(icc1:icc2),
     -            ja,line(index(ja))(inc1:inc2),irnja,
     -            line(index(ja))(irc1:irc2),sqrt(rmpx(ia))
                nbonds=nbonds+1
                if (nframe .gt. 0) then
c                 Save bond if not found yet
                  ihb=1
                  do while (ihb .le. nbfound .and.
     -              (ia .ne. ihbpair(1,ihb) .or. ja .ne.ihbpair(2,ihb)))
                    ihb=ihb+1
                  end do
                  if (ihb .gt. nbfound) then
                    nbfound=nbfound+1
                    if (ihb .lt. MAXBONDS) then
                      ihbpair(1,nbfound)=ia
                      ihbpair(2,nbfound)=ja
                    else
                      if (nframe .gt. 0) write (6,1003) nframe
                      if (nframe .gt. 0) write (iout,1003) nframe
                      write (6,1002) mxbonds
                      write (iout,1002) mxbonds
                      call askyn(
     -                  'Do you want to continue without tracking',40,0,
     -                   1,istopscan,112,0)
                      ifail=2*istopscan-1
                    end if
                  else
                    nhbdist(ihb)=nhbdist(ihb)+1
                  end if
                  rhbdist(ihb)=rhbdist(ihb)+dist_ia_ja
                  itemp(ihb)=1
                  ir1=ixres(ia)
                  ir2=ixres(ja)
                  ihb=1
                  do while (ihb .le. nbresfound .and.
     -              (ir1 .ne. ihb_pair_res(1,ihb) .or.
     -               ir2 .ne. ihb_pair_res(2,ihb)))
                    ihb=ihb+1
                  end do
                  if (ihb .gt. nbresfound) then
                    nbresfound=nbresfound+1
                    ihb_pair_res(1,nbresfound)=ir1
                    ihb_pair_res(2,nbresfound)=ir2
                    ihb_pair_res(3,nbresfound)=1
                    ihb_pair_res(3,nbresfound+1)=0
                    itempres(nbresfound)=1
                  else if (itempres(ihb) .eq. 0) then
                    ihb_pair_res(3,ihb)=ihb_pair_res(3,ihb)+1
                    itempres(ihb)=1
                  end if
                end if
              end if
            end if
          end if
        end if
      end do
      if (nbfound .le. nbitmax)
     -  call savebitc(ires(1,nframe),itemp,nbfound,30,MAXITEMS)
      return
1000  format(' List of mutually proximal heavy atom pair(s) found',a,
     -  ' at frame #',i6)
1001  format(1x,a,i5,1x,a4,1x,i4,',',a4,' - ',a,1x,i5,1x,a4,',',i4,1x,
     -  a4,' r=',f6.1,' A')
1002  format(' ERROR: maximum number of mutually proximal pairs to ',
     -  'store (',i5,') is exceeded',/,8x,'Reduce the participating ',/,
     -  8x,'atoms or increase the size of the array res in the common',
     -  ' block /analres/')
1003  format(' Frame number=',i5)
      end
