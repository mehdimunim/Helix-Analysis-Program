      subroutine filenamenum(inpfile,namleni,numfile,namlenn,n,inout)
      character*200 inpfile,numfile
      character*1 separatorchar
      common /filenuminfo/ iaskunderscore,separatorchar
c     inout=-1: input file; inout=+1: output file; inout=+2: unpacking output f
      iunderscore=0
c     print *,'FNAMNUM inout=',inout,' iaskunderscore=',iaskunderscore,
c    -  ' separatorchar=',separatorchar
c     Find last '.'
      do ic=namleni,1,-1
        if (inpfile(ic:ic) .eq. '.') then
          nl1=ic
          go to 100
        end if
      end do
      nl1=namleni+1
100   nl0=nl1-1
      if (nl0 .gt. 1) then
c       Skip '.1.' if present
        if (inpfile(nl0-1:nl0+1) .eq. '.1.' .or.
     -      inpfile(nl0-1:nl0+1) .eq. '_1.') nl1=nl1-2
        if (inout .eq. 2 .and. iaskunderscore .eq. 1) then
c         Unpacking traj or configs - always ask separator character
          print *,'Filenumber separator is . (e.g., x.2.pdb)'
          call askyn('Use _ as the filenumber separator instead',
     -      41,1,-1,iunderscore,0,0)
          iaskunderscore=-1
        end if
        if (separatorchar .eq. '.' .and. iaskunderscore .eq. 1) then
          if (inpfile(nl0-1:nl0-1) .eq. '_') then
            separatorchar='_'
            iunderscore=1
          else if (inpfile(nl0-1:nl0-1) .ne. '.') then
            if (inout .eq. -1) then
              call askyn(
     -          'Are the file numbers separated by _ (e.g., x_2.pdb)',
     -          51,1,-1,iunderscore,0,0)
            else
              print *,'Filenumber separator is . (e.g., x_2.pdb)'
              call askyn('Use _ as the filenumber separator instead',
     -          41,1,-1,iunderscore,0,0)
              iaskunderscore=-1
            end if
          end if
        end if
      end if
      if (iabs(iaskunderscore) .eq. 1) then
        if (iunderscore .eq. 1) separatorchar='_'
        if (inout .gt. 0)
     -    write (6,2000) separatorchar,inpfile(nl1:namleni)
        iaskunderscore=0
      end if
c     Insert numeral
      numfile(1:nl0)=inpfile(1:nl0)
      numfile(nl1:nl1)=separatorchar
      nl1=nl1+1
      call writeint(numfile,nl1,n,lenw)
      nl1=nl1-1
c     Add part of the filename after the '.' (if any)
      if (namleni .gt. nl0) then
        numfile(nl1+1:nl1+namleni-nl0)=inpfile(nl0+1:namleni)
        nl1=nl1+namleni-nl0
      end if
      namlenn=nl1
      return
2000  format(' Output file names should be of the form x',a1,'*',a)
      end
