      subroutine compare_rmsf(resnames,nrescol,siglev,bfac,plotfile,
     -  lplotfile,irmin,irmax,mxrsd)
      dimension bfac(mxrsd)
      character*8 resnames(mxrsd)
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      parameter (MAXRSD=70000)
      parameter (IFILL2=MAXPHI*MAXPHI*MAXPHI-4*MAXRSD)
      common /nnwork/ rmsf1(MAXRSD),rmsf2(MAXRSD),sd1(MAXRSD),
     -  sd2(MAXRSD),fill(IFILL2)
      common /colorinfo/ ncolcode,maxcolcode
      character*(*) plotfile
      common /t_test/ t_test_CI_1(5),t_test_CI_2(5),t_test_table(5,30)
      character*200 inpfile1,inpfile2
c     print *,'COMPARE_RRDIST maxcolcode,mxrsd=',maxcolcode,mxrsd
      call read_rmsf(rmsf1,sd1,41,'first',5,irmin,irmax,
     -  inpfile1,linpfile1,MAXRSD)
      call read_rmsf(rmsf2,sd2,42,'second',6,irmin2,irmax2,
     -  inpfile2,linpfile2,MAXRSD)
      if (irmin .ne. irmin2 .or. irmax .ne. irmax2) then
        print *,'Residue ranges differ - can not run the comparison'
        return
      end if
      lplotfile=0
      call openfile(42,0,'RMSF difference',15,'new',
     -  plotfile,lplotfile,notfnd,0,1,1,0,0)
      call zeroit(bfac,irmax)
c     Both SD values are computed over a sample of 10 block averages
      n1=10
      n2=10
      do ir=irmin,irmax
        rmsf_diff=rmsf1(ir)-rmsf2(ir)
        s12=sd1(ir)**2/n1+sd2(ir)**2/n2
        t=abs(rmsf_diff)/sqrt(s12)
        dof=s12**2/((sd1(ir)**2/n1)**2/(n1-1)+(sd2(ir)**2/n2)**2/(n2-1))
        idof=dof
        idof=min0(30,max0(1,idof))
        isig=1
        sig=1.0
        do while (t .gt. t_test_table(isig,idof) .and. isig .lt. 5)
          sig=t_test_CI_2(isig)
          isig=isig+1
        end do
        if (isig .eq. 5) then
          if (t .gt. t_test_table(isig,idof)) then
            sig=t_test_CI_2(isig)
          else
            isig=isig-1
          end if
        end if
        write (42,1001) ir,resnames(ir)(1:nrescol),rmsf1(ir),sd1(ir),
     -    rmsf2(ir),sd2(ir),rmsf_diff,t,dof,t_test_CI_2(isig)
        if (siglev .gt. 0.0 .and. sig .gt. siglev) rmsf_diff=0.0
        bfac(ir)=rmsf_diff
      end do
      return
1001  format(i6,1x,a,' F1=',f5.1,' sd=',f5.1,' F2=',f5.1,
     -  ' sd2=',f5.1,' D=',f5.1,' t=',f5.1,' df=',f4.1,' p<',f5.3)
      end
