      subroutine princax(c,co,aw,awo,n,indxdel,evecs0,evals0,iprint,
     -  inputref,irefcall,radtodeg,LEVTEST)
      dimension c(3,n),co(3,n),aw(n),awo(n),indxdel(n),evecs0(3,3),
     -  evals0(3)
      parameter (MAXFRAMES=50000,MAXCOPY=600)
      common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
     -  increment,increment2,res(2,MAXFRAMES,MAXCOPY),
     -  x0(MAXCOPY),y0(MAXCOPY),nxselres,ixselres(MAXCOPY)
      common /pcadat/ evecsprev(3,3),angprev(3)
      real*8 ddd
      character*1 signlab(3)
      dimension com(3),index(3),ifa(3),ila(3),itemp(3),temp(3),
     -  overlap(3,3),evecs(3,3),evals(3),ang(3),evecstemp(3,3),
     -  tenssave(3,3)
      real*8 tensinert(3,3),diag(3),offdiag(3)
      data signlab /3*' '/
      if (iprint .eq. 0)  then
        iout=6
      else
        iout=iprint
      end if
c     write (iout,*)'PRINCAX inputref,irefcall=',inputref,irefcall
      write (iout,*)
      nfinal=n
      call trnsfr(co,c,3*n)
      call trnsfr(awo,aw,n)
      call extract(co,indxdel,3,n,nfinal)
      call extract(awo,indxdel,1,n,nfinal)
      call cofms(co,com,nfinal,aw)
      call zeroitd(tensinert,9)
      awsum=0.0
      do ia=1,nfinal
        do i=2,3
          do j=1,i-1
            tensinert(i,j)=tensinert(i,j)+
     -        aw(ia)*(co(i,ia)-com(i))*(co(j,ia)-com(j))
          end do
        end do
        tensinert(1,1)=tensinert(1,1)+aw(ia)*
     -    ((co(2,ia)-com(2))**2+(co(3,ia)-com(3))**2)
        tensinert(2,2)=tensinert(2,2)+aw(ia)*
     -    ((co(1,ia)-com(1))**2+(co(3,ia)-com(3))**2)
        tensinert(3,3)=tensinert(3,3)+aw(ia)*
     -    ((co(1,ia)-com(1))**2+(co(2,ia)-com(2))**2)
      end do
      do i=2,3
        do j=1,i-1
           tensinert(j,i)=tensinert(i,j)
        end do
      end do
      do i=1,3
        do j=1,3
           tenssave(j,i)=tensinert(i,j)
        end do
      end do
c     Find the eigenvectors a and eigenvalues mu of tensinert
      call dtred2(tensinert,3,3,diag,offdiag)
      call dtqli(diag,offdiag,3,3,tensinert,ierr)
      if (ierr .gt. 0) then
        write (6,2004)
        write (iout,2004)
        return
      end if
      do k=1,3
        evals(k)=diag(k)
      end do
      call indexit(index,1,3,0)
      call mrgsrt(6,index,evals,3,ifa,ila,itemp,temp,3)
      do k=1,3
        do l=1,3
          evecs(k,l)=tensinert(l,index(k))
        end do
      end do
      if (LEVTEST .gt. 1) then
c       Test the eigenvectors after sort
        do i=1,3
          do j=1,3
            s=0.0
            do k=1,3
              s=s+tenssave(i,k)*evecs(j,k)
            end do
            evecstemp(i,j)=s/evals(j)
          end do
        end do
c       The columns of evecstemp should be also the eigenvectors
        if (LEVTEST .gt. 2) then
          write (iout,7777) 'evecstest',((evecstemp(i,k),i=1,3),k=1,3)
          write (iout,7777) 'evecs',((evecs(k,i),i=1,3),k=1,3)
        end if
        do i=1,3
          rmin=100.0
          rmax=-100.0
          do j=1,3
            if (evecstemp(j,i) .ne. 0.0) then
              r=evecs(i,j)/evecstemp(j,i)
              if (rmin .gt. r) rmin=r
              if (rmax .lt. r) rmax=r
            end if
          end do
          if (abs(rmax-rmin) .gt. 0.01)
     -      write (iout,*) 'eval failure=',rmax-rmin
        end do
      end if
c     The rows of the matrix evecs are the eigenvectors
c      write (6,1000) diag,evals,tensinert,evecs
c1000  format(' diag=',3f15.5,/,' evals=',3f15.5,/,
c     -  ' tensinert=',/,3(3f10.5,/),' evecs=',/,3(3f10.5,/))
      if (irefcall .eq. 1) then
        call trnsfr(evals0,evals,3)
        call trnsfr(evecs0,evecs,9)
        write (iout,2000) (i,evals(i),(evecs(i,k),k=1,3),' ',i=1,3)
        write (6,2000) (i,evals(i),(evecs(i,k),k=1,3),' ',i=1,3)
        call trnsfr(evecsprev,evecs0,9)
        return
      end if
      if (nframe .le. 1) then
        if (inputref .eq. 0) then
          call trnsfr(evals0,evals,3)
          call trnsfr(evecs0,evecs,9)
          call zeroit(ang,3)
        end if
      end if
      call indexit(index,1,3,0)
      if (nframe .gt. 1 .or. inputref .eq. 1) then
        if (LEVTEST .gt. 0) then
c         Now check overlap with previous orientation
          call overlapcheck(evecsprev,evecs,overlap,index,nneg)
          write (iout,7777) 'overlap with previous frame',
     -      ((overlap(k,i),i=1,3),k=1,3)
          write (iout,*) 'index=',index
          if (nneg .gt. 0) write (iout,2003) (overlap(index(i),i),i=1,3)
          call overlapcheck(evecs0,evecs,overlap,index,nneg)
          write (iout,7777) 'overlap with reference frame',
     -      ((overlap(k,i),i=1,3),k=1,3)
          write (iout,*) 'index=',index
          if (nneg .gt. 0) write (iout,2003) (overlap(index(i),i),i=1,3)
        end if
        call overlapcheck(evecs0,evecs,overlap,index,nneg)
        do i=1,3
          ang(i)=radtodeg*
     -      dacoscheck(ddd,abs(overlap(index(i),i)),0,6,'PRINCAX')
        end do
      end if
      if (nframe .gt. 0) then
        call trajlimtest(nframe,MAXFRAMES)
        res(1,nframe,1)=evecs(1,1)
        res(2,nframe,1)=evecs(1,2)
        res(1,nframe,2)=evecs(1,3)
        res(2,nframe,2)=evecs(2,1)
        res(1,nframe,3)=evecs(2,2)
        res(2,nframe,3)=evecs(2,3)
        res(1,nframe,4)=evecs(3,1)
        res(2,nframe,4)=evecs(3,2)
        res(1,nframe,5)=evecs(3,3)
        res(1,nframe,6)=ang(1)
        res(2,nframe,6)=ang(2)
        res(1,nframe,7)=ang(3)
      end if
      if (nframe .le. 1 .and. inputref .eq. 0)
     -   call indexit(index,1,3,0)
      write (iout,2000)
     -  (i,evals(index(i)),(evecs(i,k),k=1,3),signlab(i),i=1,3)
      if (nframe .gt. 0) write (iout,2001) ang
      if (irefcall .eq. 0) call trnsfr(evecsprev,evecs,9)
      return
2000  format(' Evalue ',i1,'=',f15.4,' Principal axes=',3f9.5,1x,a)
2001  format(' Angles between initial and current principal axes=',
     -  3f8.2,' deg')
2003  format(' Negative overlap found:',3f6.1)
2004  format(' Calculation aborted due to diagonalization failure')
7777  format(1x,a,':',/,(3f10.5))
      end
