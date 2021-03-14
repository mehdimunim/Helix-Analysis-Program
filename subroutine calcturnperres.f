      subroutine calcturnperres(turnperres,nres,incrot,perpvec,axisdir,
     -  anglechangeref,irefang,pi,MAXHX)
      real*8 perpvec(3,MAXHX),axisdir(3)
      dimension anglechangeref(MAXHX)
      nresu=nres-2*incrot
      a11=(2*nresu**3+3*nresu**2+nresu)/6
      a21=(nresu*(nresu+1))/2
      a12=a21
      a22=nresu
      b1=0.0
      b2=0.0
      turna=0.0
      nturn=0
      turnchangeav=0.0
      do ir=2+incrot,nres-incrot
        call angcomp(perpvec(1,ir-1),axisdir,perpvec(1,ir),turnchange)
        if (turnchange .lt. 0.0) turnchange=turnchange+2.0*pi
        if (irefang .eq. 1) then
c         Compare with turn angle in the reference conformation
          if (abs(turnchange-anglechangeref(ir)) .gt.
     -        abs((turnchange-2.0*pi)-anglechangeref(ir)))
     -      turnchange=turnchange-2.0*pi
cd77          write (77,*) 'ir,turnchange,ref=',
cd77     -      ir,turnchange*180.0/pi,anglechangeref(ir)*180.0/pi
        end if
        turnchangeav=turnchangeav+turnchange
        turna=turna+turnchange
cd77        write (77,*) 'ir,turna,turnchange=',ir,
cd77     -    turna*180.0/pi,turnchange*180.0/pi
        b1=b1+(ir-incrot)*turna
        b2=b2+turna
      end do
      turnperres=(a22*b1-a12*b2)/(a22*a11-a12*a21)
      turnchangeav=turnchangeav/float(nres-2*incrot-1)
cd77      write (77,*) 'turnchangeav=',turnchangeav,turnchangeav*180.0/pi
      return
      end
