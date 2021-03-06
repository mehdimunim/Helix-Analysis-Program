      subroutine getfiltlims(percmind,percmaxd,minresdistd,
     -  maxresdistd,percmin,percmax,minresdist,maxresdist,
     -  nresslt,nochange,label,llabel,bondtyp,lbondtyp,iout)
      character*(*) label,bondtyp
c     print *,'GETFILTLIMS label=',label,' bondtyp=',bondtyp
      call getreal('MINimum percentage presence',27,percmind,percmin,1,
     -  7)
      call getreal('MAXimum percentage presence',27,percmaxd,percmax,1,
     -  7)
      if (percmin .gt. 0.0 .or. percmax .lt. 100.0) then
        nochange=0
        write (iout,1000) bondtyp(1:lbondtyp),label(1:llabel),percmin,
     -    percmax
      end if
      call getint('MINimum residue-residue sequence distance',41,
     -  0,minresdistd,0,minresdist,8)
      call getint('MAXimum residue-residue sequence distance',41,
     -  maxresdistd,1,nresslt-1,maxresdist,8)
      if (minresdist .gt. 0 .or. maxresdist .lt. nresslt-1) then
        nochange=0
        write (iout,1001) bondtyp(1:lbondtyp),label(1:llabel),
     -    minresdist,maxresdist
      end if
1000  format(' Plotting and correlation calculation will exclude ',a,a,
     -  ' bonds',/,5x,'that are present in less than ',f6.1,'% of the ',
     -  'time or',/,5x,'more than ',f6.1,'% of the time')
1001  format(' Plotting and correlation calculation will exclude ',a,a,
     -  ' bonds',/,5x,'whose residue-residue distance is less than ',i4,
     -  ' or greater than ',i5)
      end
