      subroutine read_rmsf(rmsf,sd,inpt,lab,llab,irmin,irmax,
     -  inpfile,linpfile,mxres)
      dimension rmsf(mxres),sd(mxres)
      character*(*) lab
      character*200 inpfile,prompt,line
      prompt(1:llab)=lab(1:llab)
      prompt(llab+1:llab+5)=' RMSF'
      lprompt=llab+5
      linpfile=0
      call openfile(inpt,0,prompt,lprompt,'old',inpfile,linpfile,notfnd,
     -  0,1,1,0,0)
      istart=1
      irmin=0
      ir=0
      call zeroit(rmsf,mxres)
      call zeroit(sd,mxres)
      do while (.true.)
        read (inpt,1000,end=999) line
        if (line(7:10) .eq. 'RMSF') then
          read (line(1:6),*,end=888,err=888) ir
          read (line(11:16),*,end=888,err=888) rmsf_ir
          read (line(21:26),*,end=888,err=888) sd_ir
          if (irmin .eq. 0) irmin=ir
          if (ir .gt. 0 .and. ir .le. mxres) then
            rmsf(ir)=rmsf_ir
            sd(ir)=sd_ir
          else
            write (6,1002) ir,1,mxres
          end if
        end if
      end do
999   close (inpt)
      irmax=ir
      write (6,1003) irmax
      return
888   print *,'ERROR reading record'
      print *,line(1:66)
      close (inpt)
      write (6,1001)
      return
c    1 RMSF 16.96 SD=  1.99 decide=Uncorrelated
c123456789012345678901234567890
cRes #     1(ILE ) - Res #     1(ILE ): <d>=  0.0000 A sd=  0.0000
1000  format(a)
1001  format(' Residue range found: [',i6,',',i6,']')
1002  format(' ERROR: residue # read (',i6,') is outside the range [',
     -  i6,',',i6,']')
1003  format(' Number of RMSF records found=',i6)
      end
