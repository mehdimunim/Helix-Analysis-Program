      subroutine datprt(iout,version,iprtver,mark,lmark,
     -  hostname,lhostname,iheadnode,ihostonly)
c#    MMC routine 252 lstmod: 04/24/05
c*****Prints the date and the time
      use iflport
      character*8 version
      character*(*) mark
      character*100 hostname
      character*12 today
      common /today_date/ ltoday,today
c      : Intel Fortran code
c     C@ AB : Absoft Fortran code
c     C@ G7 : Gnu G77 Fortran code
c     C@ UG : SGI IRIX Fortran code
c     C@ HP : Hewlett-Packard Fortran code
c     C@ AX : IBM AIX code
c     C@ UX : Generic Unix code
C@AB      dimension idayspm(12)
c     INTEL compiler: on Linux  only, requires -Vaxlib compilation option
      integer idtvalue(8),lmonths(12)
C@AB      character*8 date_dat
C@AB      character*10 time_dat
C@AB      character*5 zone_dat
      character*8 date_dat
      character*10 time_dat
      character*5 zone_dat
C@G7      character*3 mon
C@G7      character*24 adate
C@G7      external fdate
C@G7      external etime
C@UG      character*24 fdate
C@UG      external fdate
C@AX      character*24 fdate_
C@AX      external fdate_
C@HP      character*24 fdate
C@HP      external fdate
      character*5 months(12)
      data months/'Jan.','Feb.','March','April','May','June',
     -   'July','Aug.','Sep.','Oct.','Nov.','Dec.'/
      data lmonths /4,4,5,5,3,7*4/
C@AB      data idayspm /31,28,31,30,31,30,31,31,30,31,30,31/
c     Find out host name
      lhostname=0
      iheadnode=0
      call blankout(hostname,1,100)
C@AB      call getenv('HOST',VALUE=hostname)
C@UG      call getenv('HOST',hostname)
C@G7      call getenv('HOST',hostname)
      istat=hostnam(hostname)
      call lastchar(hostname,lhostname,100)
      if (lhostname .gt. 0) then
        if (hostname(1:lhostname) .eq. 'minerva2' .or.
     -      hostname(1:lhostname) .eq. 'login1') iheadnode=1
      else
        iheadnode=-1
      end if
      if (ihostonly .eq. 1) return
      if (iprtver .eq. 1) write (iout,1000) mark(1:lmark),version
      call zeroiti(idtvalue,0,8)
      ltoday=0
      ief=0
      ief=1
      iab=0
C@AB      iab=1
c     Print time and date
C@UX      call system('date')
C@UG      write (iout,1008) mark(1:lmark),fdate()
C@AX      write (iout,1008) mark(1:lmark),fdate_()
C@G7      call fdate(adate)
C@G7      read (adate,1012) mon,iday,ihour,imin,isec,iyear
C@G7      write (iout,1018) mark(1:lmark),mon,iday,iyear
C@G7      write (today,1002) mon,iday,iyear
C@G7      ltoday=12
C@G7      write (iout,1011) mark(1:lmark),ihour,imin,isec
      if (ief+iab .gt. 0) then
        call date_and_time(date_dat,time_dat,zone_dat,idtvalue)
C@AB        call date_and_time(date_dat,time_dat,zone_dat,idtvalue)
        write (iout,1009) mark(1:lmark),
     -    months(idtvalue(2))(1:lmonths(idtvalue(2))),idtvalue(3),
     -    idtvalue(1)
        write (iout,1011) mark(1:lmark),idtvalue(5),idtvalue(6),
     -    idtvalue(7)
        write (today,1001) idtvalue(2),idtvalue(3),idtvalue(1)
        ltoday=10
        do ic=1,ltoday
          if (today(ic:ic).eq. ' ') today(ic:ic)='0'
        end do
      end if
      return
1000  format(a,'SIMULAID Version: ',a)
1001  format(i2,'/',i2,'/',i4)
C@G71002  format(a3,' ',i2,', ',i4)
C@UG1008  format(a,'Date: ',a)
C@AX1008  format(a,'Date: ',a)
1009  format(a,'Date: ',a,i3,', ',i4)
1011  format(a,'Time:',i5,' hours,',i3,' minutes,',i3,' seconds')
C@G71012  format(4x,a3,4(1x,i2),1x,i4)
C@G71018  format(a,'Date: ',a3,1x,i2,', ',i4)
      end
