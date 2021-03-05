      subroutine writeout_h(dir,ip,fp,rms,message,iprint)
      real*8 dir(3),ip(3),fp(3),rms
      integer iprint,kprint
      character*60 message

c     don't waste time
      kprint=mod(iprint,10)
      if (kprint .eq. 0) return

 10   format(a6,1x,3(f12.6,1x))
 20   format('     RMS deviation: ',f11.6)
      if (message(1:5) .ne. '     ') write(*,*) '  '//message
      write(*,10) '     D',dir
      write(*,10) '     I',ip
      write(*,10) '     F',fp
      write(*,20) rms
      end
