      subroutine pairdistcalc(c,nslt,npairs,listpairdist,pairdistsum,
     -  pairdistsum2,pairdistwsum,npairdist,pairdistminmax,pairgrid,
     -  iout,maxddbin,maxddistr)
      dimension c(3,nslt),listpairdist(2,maxddistr),
     -  npairdist(maxddbin,maxddistr),pairdistsum(maxddistr),
     -  pairdistsum2(maxddistr),pairdistwsum(2,maxddistr),
     -  pairdistminmax(2,maxddistr)
      do ip=1,npairs
        d2=dist2(c(1,listpairdist(1,ip)),c(1,listpairdist(2,ip)))
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
        if (iout .gt. 0) write (iout,1000) (listpairdist(k,ip),k=1,2),d
      end do
      return
1000  format(' Atom',i7,' - atom',i7,' distance=',f5.2)
      end
