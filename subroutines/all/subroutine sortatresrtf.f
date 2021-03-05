      subroutine sortatresrtf(line,n,nslt,index,indexn,indexo,isegno,
     -  itemp1,irescol1,irescol2,iresncol1,iresncol2,isegcol1,isegcol2,
     -  inamcol1,inamcol2,inpcrdtyp,nconfig,maxrec)
      character* 132 line(maxrec)
      dimension isegno(n),index(n),indexo(n),indexn(maxrec),itemp1(n)
c*****Sort atoms within a residue to follow an RTF order
      character*4 ires,reso,resn,atnam,sego,segn
      character*8 convdat
      character*200 sfilename
      common /savedat/ mxresdat,maxcondat,ifst(1000),ilst(1000),
     -  nres,iresgen,lsfilename,ires(1000),convdat(7,10000),sfilename
      character*5 crdext
      common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
     -  iommod,iommc,iommc4,iogro,iomol2,iomae,iocif,ioxxx,ioins,
     -  ionxyz,iosxyz,iosxyzrq,iograsp,iofull,lext(19),crdext(19)
c     print *,'sortatresrtf n=',n
      nresproc=0
      natproc=0
      reso='    '
      sego='    '
      ireso=0
      nresnotfound=0
      ifaslt=1
      call trnsfi(indexo,index,n)
      call zeroiti(indexn,0,n)
      lench=irescol2-irescol1+1
      lens=isegcol2-isegcol1+1
      numlen=iresncol2-iresncol1+1
      do ia=1,n
        resn=line(index(ia))(irescol1:irescol2)
        segn=line(index(ia))(isegcol1:isegcol2)
        if (numlen .eq. 5) then
          read (line(index(ia))(iresncol1:iresncol2),1000,err=999) iresn
        else if (numlen .eq. 4) then
          read (line(index(ia))(iresncol1:iresncol2),1001,err=999) iresn
        else
          print *,'Invalid resnum length=',numlen
          stop
        end if
        if (ia .eq. n .or. iresn .ne. ireso .or.
     -      segn(1:lens) .ne. sego(1:lens)) then
c         End of residue found
          ilaslt=ia-1
          if (ia .eq. n) ilaslt=n
          nats=ilaslt-ifaslt+1
          if (nresproc .gt. 0) then
c           Find the residue
            call findname(reso,ires,1,nres,ix,lench)
c           print *,'nresproc,reso,ix=',nresproc,reso,ix
            if (ix .gt. 0) then
c             Check the number of atoms
              ifadb=ifst(ix)
              iladb=ilst(ix)
              natsdb=iladb-ifadb+1
              if (nats .ne. natsdb) write (6,2000) reso,nats,natsdb
c             Sort the residue
c              write (77,7712) ifadb,iladb,
c      -         (convdat(2,i)(1:4),i=ifadb,iladb)
c7712          format('ifadb,iladb=',2i4,' convdat2=',/(10a5,/))
              nnotfound=0
              natsdbfound=0
              do iaa=ifaslt,ilaslt
                atnam=line(index(iaa))(inamcol1:inamcol2)
                call leftadjust4(atnam,atnam)
                if (idigit(atnam(1:1),1) .eq. 1)
     -            call regularpdb(atnam,atnam,-1)
c               write (77,*) 'iaa,resnam,atnam=',iaa,reso,'*',atnam,'*'
                do ira=ifadb,iladb
                  if (atnam .eq. convdat(2,ira)(1:4)) then
                    indexn(ira-ifadb+ifaslt)=iaa
                    natsdbfound=natsdbfound+1
                    go to 100
                  end if
                end do
                write (6,2001) iaa,
     -            line(index(iaa))(inamcol1:inamcol2),reso
                nnotfound=nnotfound+1
                itemp1(nnotfound)=iaa
100             continue
              end do
              if (natsdbfound .lt. natsdb)
     -          write (6,2003) reso,natsdbfound,natsdb
              if (nnotfound .gt. 0) then
c6577  format(1x,a,' ifaslt,ilaslt=',2i5,' indexn:',(/,10i5))
c                write (6,6577) 'Before compacting',ifaslt,ilaslt,
c     -            (indexn(iaa),iaa=ifaslt,ilaslt)
c               First compact the index list
                ndel=0
                do iaa=ifaslt,ifaslt+iladb-ifadb+1
                  if (indexn(iaa) .gt. 0) then
                    indexn(iaa-ndel)=indexn(iaa)
                    if (ndel .gt. 0) indexn(iaa)=0
                  else
                    ndel=ndel+1
                  end if
                end do
                nmoved=0
                do iaa=ifaslt,ilaslt
                  if (indexn(iaa) .eq. 0) then
                    nmoved=nmoved+1
                    indexn(ilaslt-nnotfound+nmoved)=itemp1(nmoved)
                  end if
                end do
              end if
            else if (n .le. nslt) then
              nresnotfound=nresnotfound+1
              call indexit(indexn,ifaslt,ilaslt,0)
            end if
          end if
          nresproc=nresproc+1
          sego=segn
          reso=resn
          ireso=iresn
          ifaslt=ilaslt+1
        end if
      end do
      if (nresnotfound .gt. 0) print *,
     -  'The number of residues not found in the database=',nresnotfound
c     Now, indexn contains the new order (from 1 to n) or zero for
c     database atoms that were found no match.
c6711  format(1x,a,/,(30i4))
c      write (6,6711) 'index',(index(i),i=1,n)
c      write (6,6711) 'indexn',(indexn(i),i=1,n)
c     Check for missed atoms
      ndel=0
      do ia=1,n
        if (indexn(ia) .eq. 0) then
          ndel=ndel+1
          print *,'PROGRAM ERROR: atom ',ia,' is misplaced'
        end if
      end do
      if (ndel .gt. 0) then
        print *,'PROGRAM ERROR:',ndel,' atoms were not taken care of'
        print *,'- aborting rearrangement'
        return
      end if
c      write (6,6711) 'Cond indexn',(indexn(i),i=1,n)
      do ia=1,n
        index(ia)=indexo(indexn(ia))
c       isegno should be the same within residues
        isegno(ia)=isegno(indexn(ia))
      end do
c      write (6,6711) 'Cond index',(index(i),i=1,n)
      if (ndel .gt. 0) print *,'Number of atom records deleted=',ndel
      if (inpcrdtyp .eq. iommod) then
c       Macromodel, rerrange the connectivity
        do ia=1,n
          do ib=1,6
            ic1=6+(ib-1)*8
            ic2=ic1+4
c           write (6,888) ia,ib,ic1,ic2,iold
c888        format(' ia,ib=',2i3,' ic1,2=',2i3,' iold=',i3)
            read (line(index(ia))(ic1:ic2),1000) iold
            if (iold .gt. 0)
     -        write (line(index(ia))(ic1:ic2),1000) indexn(iold)
          end do
        end do
      end if
      if (nconfig .eq. 1) then
        print *,'Number of residues processed=',nresproc-1
        if (nresnotfound .gt. 0) then
           print *,'WARNING: ',nresnotfound,
     -     ' solute residues were not found in the RTF info file'
          call askstop(0)
        end if
      end if
      return
999   write (6,2002) ia,index(ia),line(index(ia))(iresncol1:iresncol2)
      stop
1000  format(i5)
1001  format(i4)
2000  format(' WARNING: number of atoms in the input residue ',a4,/,
     -  ' differs from the number in the RTF info file:',i5,' vs',i5)
2001  format(' No match was found for atom ',i4,2x,a4,' residue ',a4,/,
     -  ' - it will be moved to the end of this residue')
2002  format(' ERROR: invalid residue number for atom no',i5,
     - ' line no',i5,':',a)
2003  format(' WARNING: not all atoms in the input residue ',a4,/,
     -  ' were matched to the RTF info file:',i5,' vs',i5)
      end
