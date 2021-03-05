      subroutine bondord(iatnum,mmtype,nattot,nneig,ineig,nhneig,ibnd,
     -  maxneig,c,index,nntemp,nvalx,inc1,inc2,irc1,irc2,line,nconfig,
     -  maxrepconf,maxrec)
c     Establish the bond orders for Macromodel format (as best as possible)
      dimension c(3,nattot),iatnum(nattot),mmtype(nattot),nneig(nattot),
     -  ineig(maxneig,nattot),nhneig(nattot),ibnd(maxneig,nattot),
     -  index(nattot),nntemp(nattot),nvalx(nattot)
      character* 132 line(maxrec)
      character*2 iatnm2
      common /atnam/ aw(99),nval(99),nvalmax(99),mmofan(99),
     -  mmatno(64),iatnm2(99)
      character*4 dbonds(3,100),dbres(100),thisres,thisresl,
     -  atnaml,atnam1l
      character*8 resnam
      dimension ifres(100),ilres(100)
      real*8 cosa
      data jmin /0/
c     For residue dbres(1,i) there is a double bond between atoms
c     dbres(2,i) and dbres(3,i)
c     All names leftadjusted
      data dbonds /'TYR ','CG  ','CD1 ','TYR ','CE1 ','CZ  ',
     -  'TYR ','CE2 ','CD2 ','PHE ','CG  ','CD1 ',
     -  'PHE ','CE1 ','CZ  ','PHE ','CE2 ','CD2 ',
     -  'HSD ','CG  ','CD2 ',
     -  'TRP ','CG  ','CD1 ','TRP ','CD2 ','CE2 ',
     -  'TRP ','CZ2 ','CH2 ','TRP ','CZ3 ','CE3 ',
     -  'ADE ','C4  ','C5  ','GUA ','C4  ','C5  ',
     -  'CYT ','C5  ','C6  ','THY ','C5  ','C6  ',
     -  'URA ','C5  ','C6  ',
     -  252*'    '/,resnam /'     '/
c     Scan dbond list to find a list of special residues
      thisres='****'
      ndres=0
      do i=1,100
        if (thisres .ne. dbonds(1,i)) then
c         New residue in the double bond list
          thisres=dbonds(1,i)
          ndres=ndres+1
          dbres(ndres)=thisres
          ifres(ndres)=i
          if (ndres .gt. 1) ilres(ndres-1)=i-1
        end if
      end do
      ndres=ndres-1
      if (nconfig .eq. 1) then
        print  *,'Built in double bond list:'
        do i=1,ndres
          write (6,1101) i,dbres(i),
     -      (dbonds(2,j),dbonds(3,j),j=ifres(i),ilres(i))
        end do
      end if
c     Find hybridization
      ncgeneric=0
      do i=1,nattot
        mmtype(i)=0
c        write (77,9857) i,iatnum(i),(ineig(j,i),j=1,nneig(i))
c9857    format(i5,' atnum=',i3,' ing=',20i5)
        if (iatnum(i) .eq. 6 .or. iatnum(i) .eq. 7) then
          if (nneig(i) .eq. 4 .or.
     -        (iatnum(i) .eq. 7 .and. nneig(i) .eq. 3)) then
            mmtype(i)=3
          else if (nneig(i) .gt. 1) then
            angav=0
            do in=1,nneig(i)
c             Calculate angle betwen in-th and (in+1)th neighbor
              d1s=0.0
              d2s=0.0
              d12s=0.0
              in1=in+1
              if (in1 .gt. nneig(i)) in1=1
              do k=1,3
                dx1=c(k,ineig(in,i))-c(k,i)
                dx2=c(k,ineig(in1,i))-c(k,i)
                d1s=d1s+dx1*dx1
                d2s=d2s+dx2*dx2
                d12s=d12s+dx1*dx2
              end do
              cosa=dble(d12s)/sqrt(d1s*d2s)
              ang=(180.0/3.141592)*dacoscheck(cosa,ccc,1,6,'BONDORD')
              angav=angav+ang
            end do
            angav=angav/nneig(i)
            if (angav .gt. 150.0) then
              mmtype(i)=1
            else if (angav .gt. 117.0) then
              mmtype(i)=2
            else
              mmtype(i)=3
            end if
c            write (199,7171) i,line(index(i))(inc1:inc2),angav,mmtype(i)
c7171        format(i5,2x,a4,' angle=',f8.3,' mmtype(i)=',i2)
          else if (nneig(i) .gt. 0) then
            if (iatnum(ineig(1,i)) .eq. 16) then
              mmtype(i)=6
            else
              mmtype(i)=14
              ncgeneric=ncgeneric+1
            end if
          end if
        end if
      end do
c     See end-node carbons and nitrogens
      do i=1,nattot
        if (nneig(i) .eq. 1 .and.
     -    (iatnum(i) .eq. 6 .or. iatnum(i) .eq. 7)) then
          if (mmtype(ineig(1,i)) .eq. 1 .or.
     -        mmtype(ineig(1,i)) .eq. 2) then
            if (iatnum(i) .eq. 6) mmtype(i)=2
            if (iatnum(i) .eq. 7) mmtype(i)=1
          end if
          if (mmtype(ineig(1,i)) .eq. 3) mmtype(i)=3
        end if
      end do
c     Find bond orders
      do i=1,nattot
        if (nval(iatnum(i)) .le. 0)
     -    write (6,2005) i,iatnum(i),nval(iatnum(i))
c       nvalx is the excess valence
        nvalx(i)=nval(iatnum(i))-nneig(i)
        if (nvalx(i) .lt. 0) then
          if (nvalx(i)+(nvalmax(iatnum(i))-nval(iatnum(i))) .lt. 0 .and.
     -        nconfig .lt. maxrepconf) write (6,2001) i,line(index(i))
     -        (inc1:inc2),nneig(i),nval(iatnum(i)),iatnum(i)
          nvalx(i)=0
        else if (iatnum(i) .eq. 6 .and. nneig(i) .gt. 1 .and.
     -    nneig(i) .lt. 4 .and. nhneig(i) .eq. 0) then
c         For carbons than might be united atoms, deduce united atoms
c         Do end node atoms last
          if (nneig(i) .eq. 3) then
            if (mmtype(i) .eq. 3) nvalx(i)=nvalx(i)-1
          else if (nneig(i) .eq. 2) then
            if (mmtype(i) .eq. 2) nvalx(i)=nvalx(i)-1
            if (mmtype(i) .eq. 3) nvalx(i)=nvalx(i)-2
          end if
        end if
      end do
c     Look at end-node carbons
      do i=1,nattot
        if (iatnum(i) .eq. 6 .and. nneig(i) .eq. 1) then
          if (mmtype(ineig(1,i)) .eq. 1 .or.
     -        mmtype(ineig(1,i)) .eq. 2) nvalx(i)=nvalx(i)-2
          if (mmtype(ineig(1,i)) .eq. 3) nvalx(i)=nvalx(i)-3
        end if
      end do
c     Initialize bond orders to one
      do i=1,nattot
        nntemp(i)=nneig(i)
        do j=1,nneig(i)
          ibnd(j,i)=1
        end do
      end do
c     Scan the atom list to add the explicitly given double bonds
      thisres='****'
      do i=1,nattot
        resnam='        '
        resnam=line(index(i))(irc1:irc2)
        imod=0
        if (resnam(1:4) .ne. thisres) then
c         New residue
          thisres=resnam(1:4)
          call leftadjust4(thisres,thisresl)
          if (thisresl(1:2)  .eq. 'HS') thisresl='HSD '
c         See if there is double bond in the list for this residue
          do ir=1,ndres
            if (thisresl .eq. dbres(ir)) then
              imod=ir
              go to 1100
            end if
          end do
        end if
1100    if (imod .gt. 0) then
          call leftadjust4(line(index(i))(inc1:inc2),atnaml)
          do irl=ifres(imod),ilres(imod)
c           Check if atom (i) is in the dblist
            if (atnaml .eq. dbonds(2,irl)) then
c             Check the bonds of atom i against the list
              do j=1,nneig(i)
                in1=ineig(j,i)
                call leftadjust4(line(index(in1))(inc1:inc2),atnam1l)
                if (atnam1l .eq. dbonds(3,irl)) then
c                 Match found
                  ibnd(j,i)=ibnd(j,i)+1
                  nvalx(i)=nvalx(i)-1
c                 print *,'Added i,in1=',i,in1
                  do jj=1,nneig(in1)
                    if (ineig(jj,in1) .eq. i) then
                      ibnd(jj,in1)=ibnd(jj,in1)+1
                      nvalx(in1)=nvalx(in1)-1
                      go to 1200
                    end if
                  end do
                  print *,'Cant add reverse bond i,in1=',i,in1
                end if
              end do
            end if
          end do
        end if
1200    continue
      end do
      if (nconfig .eq. 1) write (6,2002)
c     Take care of atoms with one neighbours first, in an iterative fashion
      naddtot=0
      nvalerr=0
300   nadd=0
      do i=1,nattot
        if (nntemp(i) .eq. 1) then
          ii=i
200       in1=ineig(1,ii)
          nn1=nntemp(in1)
          do j=1,nn1
            if (ineig(j,in1) .eq. ii) jn1=j
          end do
          if (nvalx(ii) .gt. 0 .and. nvalx(in1) .gt. 0) then
c           Multiple terminal bond found, add to bond orders
            ndel=min0(nvalx(ii),nvalx(in1))
            nadd=nadd+ndel
            nvalx(ii)=nvalx(ii)-ndel
            nvalx(in1)=nvalx(in1)-ndel
            ibnd(1,ii)=ibnd(1,ii)+ndel
            ibnd(jn1,in1)=ibnd(jn1,in1)+ndel
c           Drop atom from further consideration
            if (nvalx(ii) .gt. 0) then
              nvalx(ii)=-nvalx(ii)
              nvalerr=nvalerr+1
            end if
            call swapng(ineig,ibnd,in1,jn1,nn1,nattot,maxneig)
            nntemp(ii)=0
            nntemp(in1)=nntemp(in1)-1
            if (nntemp(in1) .eq. 1) then
c             Immediate neighbour became an end node - repeat
              ii=in1
              go to 200
            end if
          end if
        end if
      end do
      naddtot=naddtot+nadd
      if (nadd .gt. 0) go to 300
      if (naddtot .gt. 0 .and. nconfig .le. maxrepconf)
     -   print *,naddtot,' tree bond orders added'
c     At this point only atoms in loops are left
      nadd=0
      do 310 i=1,nattot
        if (nvalx(i) .gt. 0) then
c         Distribute additional bond orders
          nn=nneig(i)
c         Find the neighbor with the lowest bond order
          ibomin=1000
          do j=1,nn
            if (nvalx(ineig(j,i)) .gt. 0 .and. ibomin .gt.
     -          ibnd(j,i)) then
              ibomin=ibnd(j,i)
              jmin=j
            end if
          end do
          if (ibomin .lt. 1000) then
            nadd=nadd+1
            ibnd(jmin,i)=ibnd(jmin,i)+1
            jn=ineig(jmin,i)
            nnj=nneig(jn)
            do k=1,nnj
              if (ineig(k,jn) .eq. i) ibnd(k,jn)=ibnd(k,jn)+1
            end do
            nvalx(i)=nvalx(i)-1
            nvalx(jn)=nvalx(jn)-1
            if (nvalx(i) .eq. 0) go to 310
          end if
        end if
310   continue
      if (nadd .gt. 0 .and. nconfig .le. maxrepconf)
     -  print *,nadd,' loop bond orders added'
c     Now see atoms with nvalmax > nval (like N)
      nadd=0
      do i=1,nattot
        if (nvalmax(iatnum(i)) .gt. nval(iatnum(i))) then
          if (nneig(i) .lt. nvalmax(iatnum(i))) then
            ibosum=0
            do in=1,nneig(i)
              ibosum=ibosum+ibnd(in,i)
            end do
            if (nvalx(i) .ge. 0) nvalx(i)=nvalmax(iatnum(i))-ibosum
            if (ibosum .lt. nvalmax(iatnum(i))) then
c             Check if any neighbor has unsatisfied valence
              do in=1,nneig(i)
                i1=ineig(in,i)
                if (nvalx(i1) .gt. 0) then
c                 Increase bond order of i-i1)
                  ibnd(in,i)=ibnd(in,i)+1
                  do in1=1,nneig(i1)
                    if (ineig(in1,i1) .eq. i)
     -                ibnd(in1,i1)=ibnd(in1,i1)+1
                  end do
                  nvalx(i1)=nvalx(i1)-1
                  if (nvalx(i) .gt. 0) nvalx(i)=nvalx(i)-1
                  nadd=nadd+1
                end if
              end do
            end if
          end if
        end if
      end do
      if (nadd .gt. 0) print *,nadd,' N+ type bond orders added'
      nswap=0
400   nvxe=0
      do i=1,nattot
        if (nvalx(i) .gt. nvalmax(iatnum(i))-nval(iatnum(i))) then
          nvxe=nvxe+1
        end if
      end do
      if (nvxe .gt. 0) then
c       See if loop double bonds can be swapped around
        nadd=0
        do i=1,nattot
c          write (77,9571) i,(ineig(j,i),j=1,nneig(i))
c9571      format(i4,' ing=',20i5)
          if (nvalx(i) .gt. 0 .and. iatnum(i) .ne. 7) then
c           Scan neighbor for atoms with multiple bond
            do in=1,nntemp(i)
              i1=ineig(in,i)
c             write (77,*) 'i,in,i1=',i,in,i1
              do in1=1,nntemp(i1)
                i2=ineig(in1,i1)
c                write (77,*) 'i,in1,i2=',i,in1,i2
                if (ibnd(in1,i1) .gt. 1 .and. i2 .ne. i) then
c                 Multiple bond was found in the neighborhood: i1-i2
c                 Check if i2 has a neighbor with unsatisfied valence
                  do in2=1,nntemp(i2)
                    i3=ineig(in2,i2)
c                   write (77,*) 'i,in2,i3=',i,in2,i3
                    if (i3 .ne. i1 .and. nvalx(i3) .gt. 0) then
c                     Decrement bond order of i1-i2; increment i-i1 and i2-i3
c                     First move the neig info to their new place
                      call swapng(ineig,ibnd,i2,in2,nntemp(i2),nattot,
     -                  maxneig)
                      do in3=1,nntemp(i3)
                        if (ineig(in3,i3) .eq. i2) call
     -                    swapng(ineig,ibnd,i3,in3,nntemp(i3),nattot,
     -                      maxneig)
                      end do
c                       write (77,*) 'i,i3,nntemp(i3)=',i,i3,nntemp(i3)
                      ibnd(nntemp(i3),i3)=ibnd(nntemp(i3),i3)+1
                      nvalx(i3)=nvalx(i3)-1
                      nntemp(i3)=nntemp(i3)-1
                      call swapng(ineig,ibnd,i,in1,nntemp(i),nattot,
     -                  maxneig)
                      do in21=1,nntemp(i1)
                        if (ineig(in21,i1) .eq. i) call
     -                    swapng(ineig,ibnd,i,in21,nntemp(i1),nattot,
     -                  maxneig)
                      end do
                      ibnd(nntemp(i),i)=ibnd(nntemp(i),i)+1
                      nvalx(i)=nvalx(i)-1
                      nntemp(i)=nntemp(i)-1
                      nadd=nadd+1
                      go to 410
                    end if
                  end do
                end if
              end do
            end do
          end if
410       continue
        end do
        nswap=nswap+nadd
        if (nadd .gt. 0) go to 400
      end if
      if (nswap .gt. 0) print *,nswap,' bond orders added by swapping'
      nvxe=0
      do i=1,nattot
        if (nvalx(i) .gt. nvalmax(iatnum(i))-nval(iatnum(i))) then
          nvxe=nvxe+1
        end if
      end do
      if (nconfig .lt. maxrepconf) write (6,2000) nvxe
c      do i=1,nattot
c        write (100,7711) i,iatnum(i),nneig(i),nvalx(i),
c     -      (ineig(j,i),ibnd(j,i),j=1,nneig(i))
c7711    format(i4,' atno=',i3,' nn=',i2,' nvalx=',i3,
c     -    ' ineig,iibnd=',6(i4,i2,2x))
c      end do
c     Finally, set atomtypes
      do i=1,nattot
        if (iatnum(i) .eq. 6) then
          nunsatur=mmtype(i)+1-nneig(i)
          if (nunsatur .gt. 0) then
c           Unsaturated carbon --> united atom
            if (mmtype(i) .eq. 3) mmtype(i)=3+nunsatur
            if (mmtype(i) .eq. 2) mmtype(i)=6+nunsatur
            if (mmtype(i) .eq. 1) mmtype(i)=9+nunsatur
          end if
c         Saturated carbon's atomtype is just the hybridization number
        else if (iatnum(i) .eq. 7) then
          if (mmtype(i) .eq. 3) then
            if (nneig(i) .eq. 4) then
              mmtype(i)=32
            else if (nneig(i) .eq. 3) then
              mmtype(i)=26
            else if (nneig(i) .eq. 2) then
              mmtype(i)=27
            else if (nneig(i) .eq. 1) then
              mmtype(i)=28
            end if
          else if (mmtype(i) .eq. 2) then
            if (nneig(i) .eq. 3) then
              mmtype(i)=31
            else if (nneig(i) .eq. 2) then
              mmtype(i)=25
            else if (nneig(i) .eq. 1) then
              mmtype(i)=30
            end if
          else if (mmtype(i) .eq. 1) then
            mmtype(i)=24
          end if
        else if (iatnum(i) .eq. 8) then
          if (nneig(i) .eq. 1 .and. ibnd(1,i) .eq. 2) then
            mmtype(i)=15
          else
            mmtype(i)=mmofan(iatnum(i))
          end if
        else if (iatnum(i) .eq. 1) then
          if (nneig(i) .lt. 1) then
            write (6,2003) i
            mmtype(i)=mmofan(iatnum(i))
          else
            if (iatnum(ineig(1,i)) .eq. 8) then
              mmtype(i)=42
            else if (iatnum(ineig(1,i)) .eq. 7) then
              mmtype(i)=43
            else if (iatnum(ineig(1,i)) .eq. 6 .or.
     -          iatnum(ineig(1,i)) .eq. 15 .or.
     -          iatnum(ineig(1,i)) .eq. 16) then
              mmtype(i)=41
            else
              mmtype(i)=mmofan(iatnum(i))
            end if
          end if
        else
          mmtype(i)=mmofan(iatnum(i))
        end if
      end do
      do i=1,nattot
        if (mmtype(i) .eq.  0) write (6,2004) i,iatnum(i),nneig(i)
      end do
      if (ncgeneric .gt. 0)
     -  print *,'WARNING:',ncgeneric,' carbon types were set to generic'
      return
c1201  format(1x,i3,6(1x,i5,1x,i1))
1101  format(i4,1x,a4,3x,10(' (',a4,'- ',a4,')'))
2000  format(' Number of atoms with unsatisfied valence=',i6)
2001  format(' WARNING: no of neighbours of atom',i5,'(',a4,')=',i2,
     -  ' exceeds valence (',i1,') ',' atno=',i4)
2002  format(' Bond orders will be determined by a simple algorithm',/,
     -  '   - possibly some multiple bonds will not be recognized',/,
     -  '   - for example, united atoms may not recognized as such')
2003  format(' WARNING: Hydrogen ',i6,' has no bond')
2004  format(' ERROR: Atom ',i5,' has no atomtype. ',/,
     -  8x,'Atomic number=',i3,' Number of neighbors=',i3)
2005  format(' ERROR: Atom ',i5,' atomic number ',i3,' has invalid ',
     -  'valence:',i3)
      end
