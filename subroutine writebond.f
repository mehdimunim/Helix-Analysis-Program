      subroutine writebond(makb,makc,iout_bond,iout_conn,nneig,ineig,
     -  nslt,maxng)
      dimension nneig(nslt),ineig(maxng,nslt)
      character*80 linewr
c     print *,'WRITEBOND iout_bond,iout_conn=',iout_bond,iout_conn,
c    -  ' nslt=',nslt,' maxng=',maxng
      if (makb .gt. 0) then
        nb2=0
        do ia=1,nslt
          nb2=nb2+nneig(ia)
c         write (6,8693) ia,(ineig(j,ia),j=1,nneig(ia))
c8693     format(i4,' nn:',10i4)
        end do
        if (mod(nb2,2) .ne. 0) then
          print *,'PROGRAM ERROR: nneig sum is not even:',nb2
          stop
        end if
        nb=nb2/2
        write (iout_bond,*) 'MAKB ',nb,' ~'
        linewr(73:74)=' ~'
        llw=74
        nw=0
        iw=0
        do ia=1,nslt
          do jaa=1,nneig(ia)
            ja=ineig(jaa,ia)
            if (ja .gt. ia) then
              nw=nw+1
              iw=mod(nw-1,4)*18
              write (linewr(iw+1:iw+18),1109) ia,ja
              if (nw .eq. nb)call blankout(linewr,iw+19,80)
              if (mod(nw,4) .eq. 0 .or. nw .eq. nb) then
                write (iout_bond,1000) linewr(1:llw)
              end if
            end if 
          end do
        end do
      end if 
      if (makc .gt. 0) then
        do ia=1,nslt
          do jaa=1,nneig(ia)
            ja=ineig(jaa,ia)
            if (ja .gt. ia) write (iout_conn,1110) ia,ja
          end do
        end do
      end if 
      return
1000  format(a)
1109  format(2i9)
1110  format('CONECT',4i5)
      end
