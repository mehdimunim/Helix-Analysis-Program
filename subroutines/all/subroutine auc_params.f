      subroutine auc_params(iauctype,lastframeinp,loffmin,nreusemax,
     -  nseg_scr,iaucw,nframe,iframeunit,framefac,label,llabel,
     -  iuselaston,iout,mxframes)
c     Establish auc caculation details
      character*(*) label
      character*6 frunit
      common /frameunit/ lfrunit(4),frunit(4)
      character*1 ans
      dimension lauc_type(5)
      character*60 auc_type(5)
      data auc_type /
     -  'AUC data generation is ended at the last ON state',
     -  'AUC data generation is ended at the last frame of the track',
     -  'Track is padded with additional zeros',
     -  'Track is padded by reusing the last stretch of the track',
     -  'Track is padded by random 0/1 string, p(1)=fraction on'/,
     -  lauc_type /49,59,37,56,54/
      call quiz(ans,iauctype,'o',' ',0,
     -  'autocorrelation calculation mode',32,0,5,6,000)
      write (iout,1013) auc_type(iauctype)(1:lauc_type(iauctype))
      lastframeinp=0
      if (iauctype .gt. 2) then
        call getint('Frame number to pad the tracks to',33,mxframes,1,
     -    mxframes,lastframeinp,000)
        write (iout,1014) lastframeinp
      end if
      if (iauctype .gt. 3) then
        write (6,1015)
        call getint(
     -    'Length of the off stretch to shift to zero padding',50,
     -    mxframes,1,mxframes,loffmin,000)
        if (loffmin .lt. mxframes) write (iout,1016) loffmin
        if (iauctype .eq. 4) then
          call getint('Maximum number of frames to reuse',33,nframe,
     -      1,nframe,nreusemax,000)
          if (nreusemax .lt. nframe) write (iout,1017) nreusemax
          call getint(
     -      'Number of segments to scramble before reusing a track',
     -      53,0,1,100,nseg_scr,000)
          if (nseg_scr .gt. 0) write (iout,1012) nseg_scr
        end if
      end if
      call askyn(
     -  'Do you want to write the full autocorrelation function file',
     -  55,1,-1,iaucw,0,0)
      write (iout,1011) label(1:llabel)
      if (iuselaston .eq. 1) write (iout,1000)
      if (iuselaston .eq. 0) write (iout,1001)
      write (iout,1002)
      if (iframeunit .eq. 1) then
        write (iout,1010) nframe
        if (iaucw .eq. 1) write (iout,1004) max0(1,nframe/20)
      else
        write (iout,1008) nframe*framefac,
     -    frunit(iframeunit)(1:lfrunit(iframeunit))
        if (iaucw .eq. 1) write (iout,1009) 10*framefac,
     -    frunit(iframeunit)(1:lfrunit(iframeunit)),
     -    max0(1,nframe/20)*framefac,
     -    frunit(iframeunit)(1:lfrunit(iframeunit))
      end if
      return
1000  format(' Ltrack: #of frames between first on to the last on')
1001  format(' Ltrack: #of frames between first on to the last frame')
1002  format(' Non: # frames the bond was on; Nbreak: Number of breaks',
     -  /,' <Lon>,<Loff>: Average length of a contiguous on and ',
     -  'off track, resp.',/,' Lon(max),Loff(max): longest contiguous ',
     -  'on and off track, resp.',/,' %break: ',
     -  '100*(Ltrack-Non)/Ltrack')
1004  format(' AUC(1-10): Autocorrelation for the first 10 frames',
     -  /,' AUC(incr): Autocorrelation at ',i6,' frame intervals')
1008  format(' Length of the trajectory analyzed=',f10.2,1x,a)
1009  format(' AUC(1-10): Bond autocorrelation for the first',f8.2,1x,a,
     -  /,' AUC(incr): Bond autocorrelation at',f8.2,1x,a,' intervals')
1010  format(' Number of frames analyzed=',i5)
1011  format(/,' === ',a,' autocorrelation results')
1012  format(' Reused tracks will be scrambled in ',i4,' segments')
1013  format(1x,a)
1014  format(' Tracks will be extended until frame number ',i6)
1015  format(' If the trajectory ended well before the last frame then',
     -  ' 0 padding will be used')
1016  format(' If the end of the trajectory contains an all 0 stretch ',
     -  'longer than',i7,' than',/,' the padding will switch to 0s')
1017  format(' Reused tracks will be limited to the last',i6,' frames ',
     -  'of each track')
      end
