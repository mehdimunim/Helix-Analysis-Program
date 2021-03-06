      subroutine extend_nnlist(nneig,ineig,n14neig,nslt,maxng,maxrec)
      dimension nneig(maxrec),n14neig(maxrec),ineig(maxng,maxrec)
      call getint('Minimum number of chemical bond to separate',43,
     -  3,1,4,nbondsep,118)
      do ia=1,nslt
        n14neig(ia)=nneig(ia)
        nnprevlev=0
        if (nbondsep .gt. 2) then
          nerr=0
          do lev=1,nbondsep-2
c           Include (lev+1)th neighbors to the excluded list
            len=n14neig(ia)
            do in=nnprevlev+1,len
              inn=ineig(in,ia)
              do in2=1,nneig(inn)
                inn2=ineig(in2,inn)
c               Check for duplicates
                ifound=0
                do iaa=1,len
                  if (inn2 .eq. ia .or. ineig(iaa,ia) .eq. inn2)
     -              ifound=1
                end do
                if (ifound .eq. 0) then
                  if (n14neig(ia) .lt. maxng) then
                    n14neig(ia)=n14neig(ia)+1
                    ineig(n14neig(ia),ia)=inn2
                  else
                    write (6,1000) ia,maxng
                    nerr=nerr+1
                  end if
                end if
              end do
            end do
            nnprevlev=len
          end do
        end if
      end do
      return
1000  format(' ERROR: number of 1-3 or 1-4 neighbors of atom ',i6,
     -  'exeeds limit (',i3,')',/,
     -  8x,'- recompile with increased MAXNEIG')
      end
