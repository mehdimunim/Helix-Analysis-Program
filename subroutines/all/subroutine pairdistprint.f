      subroutine pairdistprint(nframe,npairs,listpairdist,iclusterdist,
     -  iclustermem,ifstclst1,ifstclst2,ilstclst2,pairdistsum,
     -  pairdistsum2,pairdistwsum,npairdist,pairdistminmax,pairgrid,
     -  rmaxpair,line,index,inamcol1,inamcol2,irescol1,irescol2,
     -  iresncol1,iresncol2,inpcrdtyp,ioins,iout,nslt,maxddbin,
     -  maxddistr,maxcdlist,maxrec)
      dimension listpairdist(2,maxddistr),npairdist(maxddbin,maxddistr),
     -  iclustermem(maxcdlist),ifstclst1(maxddistr),
     -  ifstclst2(maxddistr),ilstclst2(maxddistr),
     -  pairdistsum(maxddistr),pairdistsum2(maxddistr),
     -  pairdistwsum(2,maxddistr),pairdistminmax(2,maxddistr),
     -  index(nslt)
      character*3 star3
      character*8 an1,an2,rn1,rn2,irn1,irn2
      character* 132 line(maxrec)
      character*80 pline
      dimension dist(100)
      data star3 /' * '/,nnamcol /0/,nrescol /0/,nresncol /0/,
     -  an1 /'        '/,an2 /'        '/,rn1 /'        '/,
     -  rn2 /'        '/,irn1 /'        '/,irn2 /'        '/
c     For now, only maxddbin=20 works
      if (nframe .gt. 1) write (iout,1003) rmaxpair,pairgrid
      nnamcol=inamcol2-inamcol1+1
      nrescol=irescol2-irescol1+1
      nresncol=iresncol2-iresncol1+1
c     print *,'PAIRDISTPRINT nrescol,iresncol1,iresncol2=',
c    -  nrescol,iresncol1,iresncol2
      do ip=1,npairs
        write (iout,*)
        davg=pairdistsum(ip)/nframe
c       dwavg=pairdistwsum(1,ip)/pairdistwsum(2,ip)
        dwavg=(pairdistwsum(2,ip)/nframe)**(-1.0/6.0)
        sd=sqrt(abs(pairdistsum2(ip)/nframe-davg**2))
        if (iclusterdist .eq. 0) then
          i1=listpairdist(1,ip)
          i2=listpairdist(2,ip)
        else
          pline(1:10)=' Cluster1:'
          ncol=10
          ip11=ifstclst1(ip)
          ip21=ifstclst2(ip)
          ip12=ip21-1
          ip22=ilstclst2(ip)
          if (inpcrdtyp .le. ioins) then
            nnamcol=inamcol2-inamcol1+1
            nrescol=irescol2-irescol1+1
            ncolinc=nnamcol+nrescol+7
            do i=ip11,ip12
              write (pline(ncol+1:ncol+5),1009) iclustermem(i)
              ncol=ncol+6
              pline(ncol:ncol)=' '
              pline(ncol+1:ncol+nnamcol)=
     -          line(index(iclustermem(i)))(inamcol1:inamcol2)
              ncol=ncol+nnamcol+1
              pline(ncol:ncol)=' '
              pline(ncol+1:ncol+nrescol)=
     -          line(index(iclustermem(i)))(irescol1:irescol2)
              ncol=ncol+nrescol+1
              pline(ncol:ncol)=' '
            end do
          else
            ncolinc=5
            do i=ip11,ip12
              write (pline(ncol+1:ncol+ncolinc),1009)iclustermem(i)
              ncol=ncol+ncolinc
            end do
          end if
          if (ncol+(ip22-ip21+1)*ncolinc+10 .gt. 80) then
            write (iout,1007) pline(1:ncol)
            ncol=0
          end if
          pline(ncol+1:ncol+10)=' Cluster2:'
          ncol=ncol+10
          if (inpcrdtyp .le. ioins) then
            do i=ip21,ip22
              write (pline(ncol+1:ncol+5),1009) iclustermem(i)
              ncol=ncol+6
              pline(ncol:ncol)=' '
              pline(ncol+1:ncol+nnamcol)=
     -          line(index(iclustermem(i)))(inamcol1:inamcol2)
              ncol=ncol+nnamcol+1
              pline(ncol:ncol)=' '
              pline(ncol+1:ncol+nrescol)=
     -          line(index(iclustermem(i)))(irescol1:irescol2)
              ncol=ncol+nrescol+1
              pline(ncol:ncol)=' '
            end do
          else
            do i=ip21,ip22
              write (pline(ncol+1:ncol+ncolinc),1009)iclustermem(i)
              ncol=ncol+ncolinc
            end do
          end if
          if (ncol .gt. 70) then
            write (iout,1007) pline(1:ncol)
            write (iout,1010) davg
          else
            write (pline(ncol+1:ncol+10),1010) davg
            write (iout,1007) pline(1:ncol+10)
          end if
        end if
        if (nframe .gt. 1) write (iout,1004) ip,davg,sd,
     -    (pairdistminmax(k,ip),k=1,2),dwavg
        if (inpcrdtyp .le. ioins) then
          an1(1:nnamcol)=line(index(i1))(inamcol1:inamcol2)
          an2(1:nnamcol)=line(index(i2))(inamcol1:inamcol2)
          rn1(1:nrescol)=line(index(i1))(irescol1:irescol2)
          rn2(1:nrescol)=line(index(i2))(irescol1:irescol2)
          irn1(1:nresncol)=line(index(i1))(iresncol1:iresncol2)
          irn2(1:nresncol)=line(index(i2))(iresncol1:iresncol2)
          write (iout,1000) i1,an1(1:nnamcol),rn1(1:nrescol),
     -      irn1(1:nresncol),i2,an2(1:nrescol),rn2,irn2(1:nresncol)
        end if
        if (nframe .gt. 1) then
          write (iout,1002) (npairdist(k,ip),k=1,maxddbin)
          write (iout,*)
          pdmax=0.0
          do ic=1,maxddbin
            dist(ic)=float(npairdist(ic,ip))/float(nframe)
            if (dist(ic) .gt. pdmax) pdmax=dist(ic)
          end do
          maxperc=100
          distfact=10.0
          if (pdmax .lt. 0.3) then
            maxperc=30
            distfact=0.3
          else if (pdmax .lt. 0.5) then
            maxperc=50
            distfact=0.5
          else if (pdmax .lt. 0.8) then
            maxperc=80
            distfact=0.8
          end if
          write (iout,1006)
          if (maxddbin .le. 20) then
            lbin=3
          else if (maxddbin .le. 30) then
            lbin=2
          else if (maxddbin .le. 60) then
            lbin=1
          else
            print *,'ERROR: too many bins to plot - distributions ',
     -        'will be only printed'
            lbin=0
          end if
c         if (maxddbin .le. 60) then
          if (maxddbin .eq. 20) then
            lastcol=12+lbin*maxddbin
            pline(lastcol:lastcol)='|'
            do il=10,1,-1
              call blankout(pline,1,lastcol-1)
              pline(11:11)='|'
              if (il .eq. 10) write (pline(5:9),1008) maxperc
              do ic=1,maxddbin
                if (dist(ic) .gt. distfact*float(il-1)/10.0)
     -            pline(12+(ic-1)*lbin:11+ic*lbin)=star3
              end do
              write (iout,1007) pline(1:72)
            end do
            write (iout,1006)
            write (iout,1005) ((rmaxpair/5.0)*i,i=1,5)
          else
            print *,'Sorry, only MAXDDBIN=20 works for now'
          end if
        else
          write (iout,1011) ip,davg,dwavg
        end if
      end do
      return
1000  format(' Atom',i7,' (',a,1x,a,1x,a,') - atom',i7,' (',a,1x,a,1x,a,
     -  ')')
1002  format(' Distribution :',10i6,/,(15x,10i6))
1003  format(' Distance distributions are limited to the 0 -',f5.1,
     -  ' A range; bin size=',f5.3,' A')
1004  format(' Pair',i4,' <d>=',f5.2,' SD=',f5.2,' Min=',f5.2,
     -  ' Max=',f5.2,' NMR(<d>)=',f5.2)
1005  format(11x,5f12.2)
1006  format(10x,'|',5(3('-+-'),'-|-'),'|')
1007  format(a)
1008  format(i3,' %')
1009  format(i5)
1010  format(' <d>=',f5.2)
1011  format(' Pair',i4,' <d>=',f5.2,' NMR(<d>)=',f5.2,' A')
      end
