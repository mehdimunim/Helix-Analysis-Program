      subroutine dialps(iout,labdial,llabdial,title,ltitle,remark,
     -  lremark,nrdiv,ndials,ndprow,nframesaved,pi,ipspage,npspageinc,
     -  nfravgd,idincr,noopen,noclose,iconndial,ioutpr,mappdf,ipdfgrd)
      character*(*) labdial(ndials)
      character*80 title
      character*(*) remark
      dimension llabdial(ndials)
      parameter (MAXFRAMES=50000,MAXCOPY=600)
      common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
     -  increment,increment2,res(2,MAXFRAMES,MAXCOPY),
     -  x0(MAXCOPY),y0(MAXCOPY),nxselres,ixselres(MAXCOPY)
      character*200 trajnam,trajnam2
      common /trajname/ trajnam,trajnam2,ltrajnam,ltrajnam2,ifirsttraj,
     -  ifirsttraj2,ilasttraj,ilasttraj2,incrementtraj,incrementtraj2
c     print *,'DIALPS idincr=',idincr,' ndials=',ndials
      if (ndials .eq. 0) then
        print *,'ERROR: no torsions were selected for dial plots'
        return
      end if
      print *,'Drawing ',ndials,' dials for ',nframesaved,' frames'
      nlt=min0(80,ltrajnam)
c     xm=500.0
      xm=550.0
      ym=650.0
      nddone=0
      if (ndprow .eq. 0) ndprow=xm/sqrt(xm*ym/float(ndials))+1
      edge=0.95*xm/ndprow
      maxrows=(ym-50.0)/(edge*1.2)
      maxrows=max0(1,maxrows)
      npages=(ndials-1)/(ndprow*maxrows)+1+npspageinc
      if (nfravgd .eq. 0)
     -  call getint('Number of snapshots to average in the dials',43,
     -    1,1,nframe,nfravgd,28)
      if (noopen .eq. 0) then
        call openps(iout,xm,ym,title(1:ltitle),ltitle,remark,lremark,
     -    trajnam,nlt,trajnam,0,npages,ipspage)
      else
        call plothead(iout,xm,ym,title,ltitle,remark,lremark,trajnam,
     -    nlt,trajnam,0)
      end if
      ym=ym-100
      if (nfravgd .gt. 1) then
        write (iout,1000) xm*0.04,ym*0.95+100,nfravgd
        ymfac=0.85
      else
        ymfac=0.88
      end if
      npage=0
      do while (nddone .lt. ndials)
        nddo=min0(ndials-nddone,ndprow*maxrows)
        call partwindow(xm,ym,ymfac,nddo,edge,x0,y0,nddone,ndprow,nrows,
     -    MAXCOPY)
        ipspage=ipspage+1
        write (iout,1003) ipspage
        ifirst=1
        ilast=0
        do id=nddone+1,nddone+nddo
          if (id .eq. nddone+nddo .or. mod(id,ndprow) .eq. 0) ilast=1
          call drawdial(iout,edge,nframesaved,nfravgd,id,labdial(id),
     -      llabdial(id),nrdiv,ndprow,ifirst,ilast,idincr,
     -      iconndial,ioutpr,mappdf,ipdfgrd,pi)
          ifirst=0
        end do
        nddone=nddone+nddo
        write (iout,*) 'showpage'
        if (nddone .lt. ndials) write (iout,*) '25 25 translate'
      end do
      if (noclose .eq. 0) close (iout)
      return
1000  format(2f6.1,' m',/,
     -  '(Plots were averaged over',i4,' frames) show')
1003  format('%%Page: 1 ',i4)
      end
