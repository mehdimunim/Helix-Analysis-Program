      subroutine clusterdistcalc(c,nslt,npairs,iclustermem,
     -  ifstclst1,ifstclst2,ilstclst2,pairdistsum,pairdistsum2,
     -  pairdistwsum,npairdist,pairdistminmax,pairgrid,iout,
     -  maxddbin,maxpairs,maxclustermem)
      dimension c(3,nslt),iclustermem(maxclustermem),
     -  ifstclst1(maxpairs),ifstclst2(maxpairs),ilstclst2(maxpairs),
     -  npairdist(maxddbin,maxpairs),pairdistsum(maxpairs),
     -  pairdistsum2(maxpairs),pairdistwsum(2,maxpairs),
     -  pairdistminmax(2,maxpairs)
      data d2 /0.0/
      do ip=1,npairs
        d2max=100000.0
        do i1=ifstclst1(ip),ifstclst2(ip)-1
          do i2=ifstclst2(ip),ilstclst2(ip)
            d2=dist2(c(1,iclustermem(i1)),c(1,iclustermem(i2)))
            if (d2 .lt. d2max) d2max=d2
          end do
        end do
        d=sqrt(d2)
        if (d .lt. pairdistminmax(1,ip)) pairdistminmax(1,ip)=d
        if (d .gt. pairdistminmax(2,ip)) pairdistminmax(2,ip)=d
        pairdistsum(ip)=pairdistsum(ip)+d
        pairdistsum2(ip)=pairdistsum2(ip)+d2
        pairdistwsum(1,ip)=pairdistwsum(1,ip)+d/d2**3
        pairdistwsum(2,ip)=pairdistwsum(2,ip)+1.0/d2**3
        id=d/pairgrid+1.0
        if (id .gt. maxddbin) id=maxddbin
        npairdist(id,ip)=npairdist(id,ip)+1
        if (iout .gt. 0) write (iout,1000) ip,d
      end do
      return
1000  format(' Pair ',i3,' distance=',f5.2)
      end
