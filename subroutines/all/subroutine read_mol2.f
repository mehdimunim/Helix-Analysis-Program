      subroutine read_mol2(iomol2,line,nlines,natsmol2,iresno,iatnum,
     -  title,ltitle,nnmol2,inmol2,iseqncol1,iseqncol2,inamcol1,
     -  irescol1,iresncol1,iresncol2,iccol1,iccol2,ipotcol1,iqcol1,
     -  iqcol2,lenrec,cmol2,qmol2,index,iout,nerr,maxng,maxat)
      dimension iresno(maxat),iatnum(maxat),nnmol2(maxat),
     -  inmol2(maxng,maxat),cmol2(3,maxat),qmol2(maxat),index(maxat)
      character*(*) title
      character* 132 line(maxat)
      character*1 xyz
      common /axislab/ xyz(3)
c     character*2 bondord(maxng,maxat)
c     character*4 resnamprev
      character*80 lineinp
c     print *,'READMOL2 iomol2,maxng,maxat=',iomol2,maxng,maxat
c     print *,'READMOL2 inamcol1=',inamcol1,' irescol1=',irescol1,
c    -  ' iresncol1=',iresncol1
      nerr=0
      natsmol2=0
      nres=0
c     resnamprev='XYZU'
      call blankout(lineinp,1,80)
      nlines=0
      naterr=0
      nlines_mol=0
      do while (lineinp(1:17) .ne. '@<TRIPOS>MOLECULE')
        read (iomol2,1000,end=777) lineinp
        nlines=nlines+1
        line(nlines)=lineinp
      end do
      nlines_mol=nlines
      do while (lineinp(1:13) .ne. '@<TRIPOS>ATOM')
        call blankout(lineinp,1,80)
        read (iomol2,1000,end=777) lineinp
        nlines=nlines+1
        line(nlines)(1:80)=lineinp
      end do
c     print *,line(nlines_mol+2)(1:60)
777   if (nlines_mol .eq. 0) then
        print *,'ERROR: no MOLECULE record is found'
        stop
      end if
      call lastchar(line(nlines_mol+1),ltitle,80)
      title(1:ltitle)=line(nlines_mol+1)(1:ltitle)
      read (line(nlines_mol+2),*,err=771) nats,nbonds
c     print *,'nats,nbonds=',nats,nbonds
771   do while (lineinp(1:13) .ne. '@<TRIPOS>BOND')
        call blankout(lineinp,1,80)
        read (iomol2,1000,end=888) lineinp
        nlines=nlines+1
        line(nlines)(1:80)=lineinp
        if (lineinp(1:13) .ne. '@<TRIPOS>BOND') then
          natsmol2=natsmol2+1
          index(natsmol2)=nlines
          call blankout(line(nlines),1,80)
          if (natsmol2 .gt. maxat) then
            write (6,2001) maxat
            nerr=1
            return
          end if
          call blankout(line(nlines),1,80)
          ic=1
c         Atom id
          call nextchar(lineinp,ic,lenrec)
          ic1=ic
          call nextblank(lineinp,ic,lenrec)
          nsp=iseqncol2-iseqncol1+1-(ic-ic1)
          line(nlines)(iseqncol1+nsp:iseqncol2)=lineinp(ic1:ic-1)
c         Atom name
          call nextchar(lineinp,ic,lenrec)
          ic1=ic
          call nextblank(lineinp,ic,lenrec)
          ic2=ic-1
c         Check for blank within the atom name
          icc=ic
          call nextchar(lineinp,icc,lenrec)
          if (icc-ic1 .le. 4) then
c           Fill in the blank
            ic2=icc
            do ic=ic1,ic2
              if (lineinp(ic:ic) .eq. ' ') lineinp(ic:ic)='_'
            end do
          end if
          line(nlines)(inamcol1:inamcol1+ic2-ic1)=lineinp(ic1:ic2)
          iatnum(natsmol2)=ianum(lineinp(ic1:ic2),1,ic2-ic1+1)
c         Coordinates
          lencoord=(iccol2-iccol1+1)/3
          do k=1,3
            call nextchar(lineinp,ic,lenrec)
            if (ic .eq. lenrec) go to 555
            ic1=ic
            call nextblank(lineinp,ic,lenrec)
            if (ic-ic1 .lt. 3) then
              write (6,2002) xyz(k),lineinp(ic1:ic-1)
            end if
            read (lineinp(ic1:ic-1),*,err=555) cmol2(k,natsmol2)
            write (line(nlines)
     -        (iccol1+(k-1)*lencoord:iccol1+k*lencoord-1),1001)
     -        cmol2(k,natsmol2)
          end do
c         Atomtype
          call nextchar(lineinp,ic,lenrec)
          if (ic .eq. lenrec) go to 555
          ic1=ic
          call nextblank(lineinp,ic,lenrec)
          line(nlines)(ipotcol1:ipotcol1+ic-ic1-1)=lineinp(ic1:ic-1)
c         Substance id (residue number)
          call nextchar(lineinp,ic,lenrec)
          ic1=ic
          call nextblank(lineinp,ic,lenrec)
          read (lineinp(ic1:ic-1),*,err=555) nres
c         Residue (substructure) name
          call nextchar(lineinp,ic,lenrec)
          if (ic .eq. lenrec) go to 555
          ic1=ic
          call nextblank(lineinp,ic,lenrec)
          line(nlines)(irescol1:irescol1+ic-ic1-1)=lineinp(ic1:ic-1)
c         if (resnamprev(1:ic-ic1) .ne. lineinp(ic1:ic-1)) then
c           New residue
c           nres=nres+1
c           resnamprev(1:ic-ic1)=lineinp(ic1:ic-1)
c         end if
          write (line(nlines)(iresncol1:iresncol2),1003) nres
          iresno(natsmol2)=nres
c         Charge
          call nextchar(lineinp,ic,lenrec)
          if (ic .eq. lenrec) go to 555
          ic1=ic
          call nextblank(lineinp,ic,lenrec)
          read (lineinp(ic1:ic-1),*,err=555) qmol2(natsmol2)
          write (line(nlines)(iqcol1:iqcol2),1002) qmol2(natsmol2)
c         Status bit - not read for now
        end if
        go to 556
555     call lastchar(lineinp,iclast,80)
        write (6,2000) 'atom',natsmol2,lineinp(1:iclast)
        naterr=naterr+1
        if (naterr .eq. 25) then
          call askyn('There were 25 errors so far. Do you want to stop',
     -      48,1,+1,istop,000,0)
          if (istop .eq. 1) stop
        endif
556     continue
      end do
888   write (iout,*) 'Read ',natsmol2,' atoms'
      nerr=nerr+naterr
      if (nats .ne. natsmol2) then
        write (6,2003) natsmol2,nats
        call askstop(1)
      end if
      call zeroiti(nnmol2,0,natsmol2)
      nbmol2=0
      do ibn=1,nbonds
        call blankout(lineinp,1,80)
        read (iomol2,1000,end=999,err=666) lineinp
        nlines=nlines+1
        line(nlines)=lineinp
        read (lineinp,*,end=999) ib,ia,ja
c       Bond order label is ignored for now
c       ic=1
c       do ifnd=1,4
c         call nextchar(lineinp,ic,lenrec)
c         ic1=ic
c         call nextblank(lineinp,ic,lenrec)
c       end do
        if (nnmol2(ia) .lt. maxng) then
          nnmol2(ia)=nnmol2(ia)+1
          inmol2(nnmol2(ia),ia)=ja
c         bondord(nnmol2(ia),ia)(1:ic-ic1)=lineinp(ic1:ic-1)
        else
          print *,'ERROR: atom ',ia,' has more than ',maxng,
     -      ' neighbors'
          nerr=nerr+1
        end if
        if (nnmol2(ja) .lt. maxng) then
          nnmol2(ja)=nnmol2(ja)+1
          inmol2(nnmol2(ja),ja)=ia
c         bondord(nnmol2(ja),ja)(1:ic-ic1)=lineinp(ic1:ic-1)
        else
          print *,'ERROR: atom ',ja,' has more than ',maxng,
     -      ' neighbors'
          nerr=nerr+1
        end if
        nbmol2=nbmol2+1
      end do
999   write (iout,*) 'Read ',nbmol2,' bonds'
      return
666   call lastchar(lineinp,iclast,80)
      write (6,2000) 'bond',nbmol2,lineinp(1:iclast)
      nerr=nerr+1
      return
1000  format(a)
1001  format(f12.5)
1002  format(f7.4)
1003  format(i4)
2000  format(' ERROR: Invalid .mol2 ',a,' record, natsmol2=',i5,
     -  ' record:',/,1x,a)
2001  format(' ERROR: maximum number of atoms (',i9,') is exceeded',/,
     -  ' Recompile the program with parameter MAXREC set larger')
2002  format(' ERROR: ',a1,'-coordinate has too few digits:',a)
2003  format(' WARNING: # of atoms read (',i6,') differs from the # ',
     -  'specified (',i6,')')
      end
