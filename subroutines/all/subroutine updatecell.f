      subroutine updatecell(inptrajtyp,edge)
      dimension edge(3)
      real*8 xtlabc,xtlabc0
      common /boxdat/ xtlabc(6),xtlabc0(6),box(3),box0(3),edge_gen(3,3),
     -  cell0(3,27),cell(3,27),cellalt(3,27),
     -  ncell,ioppbc,noboxinfoar,noboxinfow,noboxrep,
     -  istuner,iboxtypfound,ixcrd(3),ixang(3),ixyzhex(3),
     -  ixyzhextraj(3),isizewarn
      character*1 xyz
      common /axislab/ xyz(3)
      real*8 cellfac
      dimension cellfac(3),ixyzhextrajcharmm(3)
      data ixyzhextrajcharmm /3,2,1/
      if (ioppbc .eq. 5) then
c       For skewed hexagons, proportionality does not hold
        if (ixyzhextraj(1) .eq. 0) then
          call trnsfi(ixyzhextraj,ixyzhextrajcharmm,3)
          write (6,1000) xyz(ixyzhextraj(1)),'hexagon prism axis'
          write (6,1000) xyz(ixyzhextraj(2)),
     -      'axis crossing a hexagon vertex'
          print *,'If this is wrong, modify the subroutine updatecell'
        end if
        if (inptrajtyp .eq. 1) then
          do k=1,3
            edge(k)=xtlabc(ixcrd(ixyzhextraj(k)))
          end do
        else
          do k=1,3
            edge(k)=box(ixyzhextraj(k))
          end do
        end if
        call crorgn(edge,edge_gen,ioppbc,3,ncell,cell,cellalt,
     -    ixyzhex,rins,rc)
      else
        if (inptrajtyp .eq. 1) then
          do k=1,3
            cellfac(k)=xtlabc(ixcrd(k))/xtlabc0(ixcrd(k))
          end do
        else
          do k=1,3
            cellfac(k)=box(k)/box0(k)
          end do
        end if
        do ic=1,ncell
          do k=1,3
            cell(k,ic)=cell0(k,ic)*cellfac(k)
          end do
        end do
      end if
c      write (6,7671) xtlabc0,xtlabc,
c     -  ((cell0(k,ii),k=1,3),ii=ncell-1,ncell),
c     -  ((cell(k,ii),k=1,3),ii=ncell-1,ncell)
c7671   format(' xtlabc0=',6f12.8,/,' xtlabc=',6f12.8,/,
c     -  ' cell0i',6f12.8,' cell=',6f12.8)
       return
1000  format(' NOTE: trajectory assumes ',a,' to be the ',a)
      end
