      subroutine ca_to_bb(icaa,iresno,nneig,ineig,index,line,
     -  ic1,icna,ina,ica,inca,ires,iprotein,maxng,mxrec)
      character* 132 line(mxrec)
      dimension iresno(mxrec),index(mxrec),nneig(mxrec),
     -  ineig(maxng,mxrec)
      character*4 atnam
c     Find the atom indices for C(n-1),N,CA,C,N(n+1)
      iprotein=0
      ires=iresno(icaa)
      ina=0
      ica=0
c     print *,'ic1,maxng,mxrec=',ic1,maxng,mxrec
c     print *,'icaa,ires=',icaa,ires
      do j=1,nneig(icaa)
        ja=ineig(j,icaa)
        if (iresno(ja) .ge. ires) then
          atnam=line(index(ja))(ic1:ic1+3)
          call leftadjust4(atnam,atnam)
          if (atnam .eq. 'N   ') ina=ja
          if (atnam .eq. 'C   ') ica=ja
        end if
      end do
c     print *,' ina,ica=',ina,ica
      if (ina*ica .gt. 0) then
c       N-CA-C found
        icna=0
        inca=0
        do j=1,nneig(ina)
          ja=ineig(j,ina)
          if (iresno(ja) .eq. ires-1) then
            atnam=line(index(ja))(ic1:ic1+3)
            call leftadjust4(atnam,atnam)
            if (atnam .eq. 'C   ') icna=ja
          end if
        end do
c       print *,' icna   =',icna
        do j=1,nneig(ica)
          ja=ineig(j,ica)
          if (iresno(ja) .gt. ires) then
            atnam=line(index(ja))(ic1:ic1+3)
            call leftadjust4(atnam,atnam)
            if (atnam .eq. 'N   ') inca=ja
          end if
        end do
c       print *,' inca   =',inca
        if (icna*inca .gt. 0) iprotein=1
      end if
      return
      end
