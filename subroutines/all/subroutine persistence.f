      subroutine persistence(nhbdist,nhbpers,itf,itl,maxlenon,
     -  maxlenoff,iframeunit,framefac,nframe,iuselaston,iout)
      character*6 frunit
      common /frameunit/ lfrunit(4),frunit(4)
c     write (iout,*) 'PERSISTENCE ITF,ITL=',itf,itl,
c    -  ' IUSELASTON=',iuselaston
      if (iuselaston .eq. 1) lentrack=itl-itf+1
      if (iuselaston .eq. 0) lentrack=nframe-itf+1
      rlenb=float(nhbdist)/float(nhbpers+1)
      rlenoff=float(lentrack-nhbdist)/float(nhbpers)
      percbreak=100.0*float(lentrack-nhbdist)/float(lentrack)
      write (iout,1003) nhbpers,percbreak,lentrack,nhbdist,itf,itl
      write (iout,1005) rlenb,rlenoff,maxlenon,maxlenoff,frunit(1)
      if (iframeunit .eq. 2) write (iout,1007) rlenb*framefac,
     -  rlenoff*framefac,maxlenon*framefac,maxlenoff*framefac,
     -  frunit(iframeunit)(1:lfrunit(iframeunit))
      if (iframeunit .gt. 2) write (iout,1006) rlenb*framefac,
     -  rlenoff*framefac,maxlenon*framefac,maxlenoff*framefac,
     -  frunit(iframeunit)(1:lfrunit(iframeunit))
      return
1003  format('   Nbreak=',i5,' %off=',f5.1,' Ltrack=',i5,' Non=',i5,
     -  ' First on=',i5,' Last on=',i5)
1005  format(' <Lon>=',f9.1,' <Loff>=',f9.1,' Lon(max)=',i9,
     -  ' Loff(max)=',i9,1x,a)
1006  format(' <Lon>=',f9.4,' <Loff>=',f9.4,' Lon(max)=',f9.4,
     -  ' Loff(max)=',f9.4,1x,a)
1007  format(' <Lon>=',f9.2,' <Loff>=',f9.2,' Lon(max)=',f9.2,
     -  ' Loff(max)=',f9.2,1x,a)
      end
