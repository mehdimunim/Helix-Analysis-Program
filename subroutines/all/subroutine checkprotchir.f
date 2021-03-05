      subroutine checkprotchir(c,n,line,index,ifres,nneig,ineig,
     -  inamcol1,inamcol2,irescol1,irescol2,iresncol1,iresncol2,
     -  nminus0,nplus0,nminus,nplus,ild,iverb,maxneig,maxrsd,maxrec)
      dimension c(3,maxrec),nneig(n),ineig(maxneig,n),index(maxrec),
     -  ifres(maxrec)
      character*1 ild(maxrsd)
      character* 132 line(maxrec)
      character*8 atnam
      lnam=inamcol2-inamcol1+1
      ncheck=0
      nminus=0
      nplus=0
      nplusminus=0
      nhgen=0
      iresno=0
      read (line(index(1))(iresncol1:iresncol2),*) iresnoprev
c     print *,'CHPROTC n,maxrec,iverb,maxneig=',n,maxrec,iverb,maxneig
      icafound=0
      do ia=1,n
        atnam(1:lnam)=line(index(ia))(inamcol1:inamcol2)
        if (lnam .eq. 4) atnam(5:8)='    '
        call leftadjustn(atnam,atnam,8)
        read (line(index(ia))(iresncol1:iresncol2),*) iresno
        if (atnam(1:4) .eq. 'CA  ') then
          icafound=1
          ic=0
          in=0
          icb=0
          iha=0
          nfound=0
          ica=ia
          ierr=0
          do iaa=1,nneig(ia)
            ing=ineig(iaa,ia)
            atnam(1:lnam)=line(index(ing))(inamcol1:inamcol2)
            if (lnam .eq. 4) atnam(5:8)='    '
            call leftadjustn(atnam,atnam,8)
            call lookforneig(atnam(1:4),'C   ',ing,ic,nfound,ica,ierr)
            call lookforneig(atnam(1:4),'N   ',ing,in,nfound,ica,ierr)
            call lookforneig(atnam(1:4),'CB  ',ing,icb,nfound,ica,ierr)
            call lookforneig(atnam(1:4),'HA  ',ing,iha,nfound,ica,ierr)
          end do
          if (ierr .eq. 0) then
            if (nfound .eq. 3 .and. iha .eq. 0) then
c             HA is missing only - generate it
              iha=maxrec
              nfound=nfound+1
              do k=1,3
                c(k,maxrec)=c(k,ica)+
     -            (3.0*c(k,ica)-c(k,icb)-c(k,in)-c(k,ic))
              end do
              nhgen=nhgen+1
            end if
            if (nfound .eq. 4) then
              ncheck=ncheck+1
              call checkchir(c,maxrec,ica,in,iha,icb,ic,isg)
              if (isg .lt. 0) then
                nminus=nminus+1
                ild(ncheck)='D'
              else if (isg .gt. 0)  then
                nplus=nplus+1
                ild(ncheck)='L'
              end if
            else
              nplusminus=nplusminus+1
              ncheck=ncheck+1
              ild(ncheck)='-'
            end if
          else
            nplusminus=nplusminus+1
            ncheck=ncheck+1
            ild(ncheck)='e'
          end if
        end if
        if (iresno .ne. iresnoprev) then
          if (icafound .eq. 0) then
c           Non AA residue
            nplusminus=nplusminus+1
            ncheck=ncheck+1
            ild(ncheck)='-'
          end if
          iresnoprev=iresno
          icafound=0
        end if
      end do
      if (nhgen .gt. 0)
     -  print *,'Number of HAs generated for chirality check=',nhgen
      if (nminus0 .ge. 0 .and. nplus0 .ge .0) then
        nchch=iabs(nminus-nminus0)+iabs(nplus-nplus0)
        if (nchch .gt. 0)
     -    print *,'PROGRAM ERROR:',nchch,' CA atom(s) changed chirality'
      else
        nchch=0
      end if
      if (iverb .gt. 0 .or. nchch .gt. 0) then
        if (nplus+nminus .gt. 0) then
          write (6,1000) nplus+nminus,nplus,nminus,nplusminus
          if (nplus*nminus .ne. 0) then
c           Both L and D was found
            call print1charlist(ild,ncheck,line,index,ifres,irescol1,
     -        irescol2,maxrec)
          end if
          if (nminus .gt. 0) write(6,1001)
        else
          print *,'No chiral CA atom was found'
        end if
      end if
      return
1000  format(i6,' chiral CAs were found, ',i5,' in L and ',i5,' in D ',
     -  'conformation',/,i6,' achiral CAs were found (glycine)')
1001  format(' NOTE: chiral peptide center in D conformation was found')
      end
