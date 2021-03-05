      subroutine setdatapath
c     Find out if any of the stored datapaths has a pdb_nam.dat
      character*100 datapath,datapaths
      common /environment/ npaths,ldatapath,ldatapaths(5),
     -  datapath,datapaths(5)
      ldatapath=0
c     Try stored paths
      do ip=1,npaths
        call trypath(datapaths(ip),ldatapaths(ip),datapath,ldatapath)
        if (ldatapath .gt. 0) then
          write (6,1000) datapath(1:ldatapath)
          return
        end if
      end do
      if (ldatapath .eq. 0) then
        call blankout(datapaths(5),1,100)
        call getenv('PWD',datapaths(5))
        call lastchar(datapaths(5),ldatapathtry,100)
        call trypath(datapaths(5),ldatapathtry,datapath,ldatapath)
        if (ldatapath .gt. 0) then
          write (6,1000) datapath(1:ldatapath)
          return
        end if
      end if
      write (6,1001)
100   call getname(datapaths(5),ldatapathtry,
     -  'Path to directory of Simulaid-provided conversion files',55,
     -  100,' ',0,1,000,0)
      if (ldatapathtry .gt. 0) then
        call trypath(datapaths(5),ldatapathtry,datapath,ldatapath)
        if (ldatapath .gt. 0) then
          write (6,1000) datapath(1:ldatapath)
          return
        else
          print *,'No conversion files were found in',
     -      datapath(1:ldatapath)
          call askyn('Do you want to try a different path',35,1,1,itry,
     -      0,0)
          if (itry .gt. 0) go to 100
        end if
      end if
      return
1000  format(' Conversion files found in directory ',a)
1001  format(' NOTE: no valid datapath is found to Simulaid-provided',
     -  ' conversion files',/,' - the conversion files are assumed to ',
     -  'be in the current directory',/,' You can edit the 4th entry ',
     -  'of the character array datapaths and the ',/,'correponding ',
     -  'array entry of ldatapath or enter it here (or just hit ENTER ',
     -  'if the coversion files are in the current directory or you ',
     -  'will not need them)')
      end
