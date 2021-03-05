      subroutine checkunphys(c,nslt,n,naslv,islvw,iatnum,ifchrg,isegno,
     -  idelseg,indexo,nneig,nneiga,nhbneig,ineig,nhneig,nnneig,ncneig,
     -  nsneig,npneig,ixres,line,irescol1,irescol2,inamcol1,inamcol2,
     -  index,nconfig,innlist,molresflag,ioppbc,cell,ncell,edge,ixyzhex,
     -  molsltlim,nmolslt,hblimfac,angmin,ctfac,bondminfac,maxdist,iles,
     -  iwchk,isltonly,indices,nbox,rlim,innlistread,nframe,radtodeg,
     -  nerr,maxrepconf,maxng,maxbox,maxrsd,maxrec)
      dimension nneig(n),ineig(maxng,n),iatnum(n),ifchrg(n),c(3,n),
     -  idelseg(1000),indexo(n),nhbneig(n),nneiga(n),nhneig(n),
     -  nnneig(n),ncneig(n),nsneig(n),npneig(n),ixres(n),
     -  molresflag(maxrsd),index(n),indices(maxbox,maxrec),nbox(maxrec),
     -  isegno(n),cell(3,ncell),edge(3),ixyzhex(3),molsltlim(3,nmolslt),
     -  rlim(maxng)
      character* 132 line(maxrec)
      call bondcheck(iwchk,1,nslt,iatnum,nneig,ineig,maxng,c,maxdist,
     -  line,irescol1,irescol2,inamcol1,inamcol2,index,rlim,nbonderr,
     -  maxrec)
      if (isltonly .eq. 0) nlim=n
      if (isltonly .eq. 1) nlim=nslt
      if (iles .eq. 1) then
c       Mark segments to be disregarded
        nmark=0
        call trnsfi(indexo,isegno,nslt)
        do ia=1,nslt
          if (idelseg(isegno(ia)) .eq. 1) nmark=nmark+1
          if (idelseg(isegno(ia)) .eq. 1) isegno(ia)=-1
        end do
      end if
      if (innlistread .eq. 0)
     -  call nnlist(nslt,islvw,naslv,nlim,iatnum,ifchrg,c,nneig,nneiga,
     -    nhbneig,ineig,nhneig,nnneig,ncneig,nsneig,npneig,line,
     -    irescol1,irescol2,inamcol1,inamcol2,index,nconfig,innlist,
     -    molresflag,hblimfac,angmin,0,indices,nbox,isegno,ixres,
     -    maxrepconf,1,nframe,radtodeg,0,maxbox,maxng,maxrsd,maxrec)
      call contactcheck(iwchk,nslt,n,iatnum,c,line,irescol1,irescol2,
     -  inamcol1,inamcol2,index,ctfac,bondminfac,isltonly,
     -  naslv,nneig,ineig,isegno,ioppbc,cell,ncell,edge,ixyzhex,
     -  molsltlim,nmolslt,ncontacterr,maxng,maxrec)
        if (iles .eq. 1) call trnsfi(isegno,indexo,nslt)
      nerr=nbonderr+ncontacterr
      return
      end
