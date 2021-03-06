      subroutine nnlist(nslt,islvw,nslv,n,iatnum,ifchrg,c,nneig,nneiga,
     -  nhbneig,ineig,nhneig,nnneig,ncneig,nsneig,npneig,line,irescol1,
     -  irescol2,inamcol1,inamcol2,index,nconfig,innlist,molresflag,
     -  hblimfac,angmin,ihbondcalc,indices,nbox,isegno,ixres,maxrepconf,
     -  nowarn,nframe,radtodeg,LEVTEST,maxbox,maxng,maxrsd,maxrec)
      dimension nneig(n),ineig(maxng,n),iatnum(n),c(3,n),nhbneig(n),
     -  nneiga(n),nhneig(n),nnneig(n),ncneig(n),nsneig(n),npneig(n),
     -  molresflag(maxrsd),index(n),indices(maxbox,maxrec),nbox(maxrec),
     -  isegno(n),ixres(n),ifchrg(n)
      character* 132 line(maxrec)
      character*4 atnami
      character*8 resnami
c     Set up neighbour list
      character*4 namfcg
      common /connatdat/ ramax(99),ramax2(99),hlimfac,ianfg(99),
     -  namfcg(100),nrmw
      common /clonedat/ nclone,iaclnf(1000),iaclnl(1000),ncopcln(1000)
c     print *,'NNLIST n,nslt,islvw,ihbondcalc=',n,nslt,islvw,ihbondcalc
      if (innlist .gt. 0) then
        return
      else
        innlist=1
      end if
      do i=1,99
        ramax2(i)=ramax(i)**2
      end do
      call nninit(nneig,nhbneig,ineig,nhneig,nnneig,ncneig,
     -  nsneig,npneig,1,n,1,maxng)
c     Generate connectivity for solvents
      numsolv=(n-nslt)/nslv
c     print *,'NNLIST n,nslt,nslv,numsolv=',n,nslt,nslv,numsolv
      do is=1,numsolv
        call nnlist0o(nslt+(is-1)*nslv+1,nslt+is*nslv,iatnum,c,nneig,
     -    nhbneig,ineig,nhneig,nnneig,ncneig,nsneig,npneig,line,
     -    irescol1,irescol2,inamcol1,inamcol2,index,maxng,hblimfac,
     -    maxrec)
      end do
      ntestclone=0
      nhb=nslt
      if (ihbondcalc .eq. 1) nhb=n
      if (nclone .eq. 0) then
        call nnlist0(1,nhb,nslt,islvw,iatnum,ifchrg,c,nneig,
     -    nhbneig,ineig,nhneig,nnneig,ncneig,nsneig,npneig,line,
     -    irescol1,irescol2,inamcol1,inamcol2,index,maxng,
     -    hblimfac,angmin,ihbondcalc,indices,nbox,ixres,isegno,ifail,
     -    nframe,radtodeg,maxbox,maxrec,LEVTEST)
      else
        call nnlist0(1,iaclnf(1)-1,nslt,islvw,iatnum,ifchrg,c,nneig,
     -    nhbneig,ineig,nhneig,nnneig,ncneig,nsneig,npneig,line,
     -    irescol1,irescol2,inamcol1,inamcol2,index,maxng,
     -    hblimfac,angmin,ihbondcalc,indices,nbox,ixres,isegno,ifail,
     -    nframe,radtodeg,maxbox,maxrec,LEVTEST)
        ifirst=iaclnf(1)
        ntestclone=ifirst-1
        do ic=1,nclone
c         Create list separately for the clones
          ilast=ifirst+iaclnl(ic)-iaclnf(ic)
          call nnlist0(ifirst,ilast,nslt,islvw,iatnum,ifchrg,c,nneig,
     -      nhbneig,ineig,nhneig,nnneig,ncneig,nsneig,npneig,line,
     -      irescol1,irescol2,inamcol1,inamcol2,index,maxng,hblimfac,
     -      angmin,ihbondcalc,indices,nbox,ixres,isegno,ifail,nframe,
     -      radtodeg,maxbox,maxrec,LEVTEST)
c        Copy cloned nn info
         do id=2,ncopcln(ic)
           do ia=ifirst,ilast
             ianew=ia+(id-1)*(ilast-ifirst+1)
             nneig(ianew)=nneig(ia)
             nneiga(ianew)=nneiga(ia)
             nhbneig(ianew)=nhbneig(ia)
             nhneig(ianew)=nhneig(ia)
             nnneig(ianew)=nnneig(ia)
             nsneig(ianew)=nsneig(ia)
             npneig(ianew)=npneig(ia)
             call trnsfi(ineig(1,ianew),ineig(1,ia),maxng)
           end do
         end do
         ifirst=ifirst+ncopcln(ic)*(ilast-ifirst+1)
         ilast=ifirst-1
        end do
        if (ilast .lt. nslt)
     -    call nnlist0(ilast+1,nslt,nslt,islvw,iatnum,ifchrg,c,nneig,
     -      nhbneig,ineig,nhneig,nnneig,ncneig,nsneig,npneig,line,
     -      irescol1,irescol2,inamcol1,inamcol2,index,maxng,hblimfac,
     -      angmin,ihbondcalc,indices,nbox,ixres,isegno,ifail,nframe,
     -      radtodeg,maxbox,maxrec,LEVTEST)
        nclone=0
      end if
      if (LEVTEST .gt. 0) then
c       Repeat for test
        ntest=nslt
        if (ntestclone .gt. 0) ntest=ntestclone
        print *,'NN comparison test  for ',ntest,' atoms'
        do i=1,ntest
          write (77,7811) i,(ineig(j,i),j=1,maxng)
        end do
        call nninit(nneig,nhbneig,ineig,nhneig,nnneig,ncneig,
     -    nsneig,npneig,1,ntest,1,maxng)
        call nnlist0o(1,ntest,iatnum,c,nneig,nhbneig,
     -    ineig,nhneig,nnneig,ncneig,nsneig,npneig,line,irescol1,
     -    irescol2,inamcol1,inamcol2,index,maxng,hblimfac,maxrec)
        do i=1,ntest
          write (78,7811) i,(ineig(j,i),j=1,maxng)
        end do
      end if
c     Save the pseudo-atom neighbours for the real atoms
      do i=1,n
        nneiga(i)=nneig(i)
      end do
c     Functional group search will not see the pseudo-atom neighbours
      do i=1,n
        if (iatnum(i) .ge. 88 .and. iatnum(i) .le. 90 .and.
     -      nneig(i) .gt. 0) then
c         Pseudo atom found
          nn=nneig(i)
          do in=1,nn
            ia=ineig(in,i)
            if (ia .lt. 88 .or. ia .gt. 90) then
c             Add i to the neighbour list of ia
              nneig(ia)=nneig(ia)+1
              ineig(nneig(ia),ia)=i
            end if
          end do
        end if
      end do
      nlim=n
      if (nslv .eq. 1) nlim=nslt
c     Check for unconnected atoms not already made a molecule
      if (nowarn .eq. 0) then
        nslvmsg=0
        do ia=1,nlim
          if (nneig(ia) .eq. 0 .and. molresflag(ixres(ia)) .lt. 2) then
            imol1=0
            if (ia .eq. 1) then
              if (isegno(ia) .ne. isegno(ia+1)) imol1=1
            else if (ia .eq. nlim) then
              if (isegno(ia) .ne. isegno(ia-1)) imol1=1
            else
              if (isegno(ia) .ne. isegno(ia-1) .and.
     -            isegno(ia) .ne. isegno(ia+1)) imol1=1
            end if
            if (imol1 .eq. 0 .and. ramax2(iatnum(ia)) .gt. 0.0) then
              if (nslvmsg .le. 25) then
                if (inamcol2 .ge. inamcol1) then
                  resnami='     '
                  resnami=line(index(ia))(irescol1:irescol2)
                  atnami=line(index(ia))(inamcol1:inamcol2)
                  if (resnami(1:3) .ne. 'HOH' .and.
     -                nconfig .le. maxrepconf)
     -              write (6,2003) ia,atnami,resnami,iatnum(ia)
                else
                  if (nconfig .le. maxrepconf)
     -              write (6,2002) ia,iatnum(ia)
                end if
                if (nslvmsg .eq. 25)
     -            print *,'Solvent related warnings are turned off'
              else
              end if
              if (ia .gt. nslt)nslvmsg=nslvmsg+1
            end if
          end if
        end do
        if (nslvmsg .gt. 25) write (6,2001) nslvmsg
      end if
      return
2001  format(' WARNING: A total of ',i7,' solvent atoms had no bonds')
2002  format(' WARNING: atom ',i6,' (atomic #=',i2,
     -  ') has no bonds at all')
2003  format(' WARNING: atom ',i6,' (',a4,1x,a6,', atomic #=',i2,
     -  ') has no bonds at all')
7811  format(i6,' in=',(20i6))
      end
