      subroutine contactcheck(iout,nslt,n,iatnum,c,line,irescol1,
     -  irescol2,inamcol1,inamcol2,index,ctlimfac,bondminfac,
     -  isltonly,naslv,nneig,ineig,isegno,ioppbc,cell,ncell,edge,
     -  ixyzhex,molsltlim,nmolslt,nerr,maxng,maxrec)
      character* 132 line(maxrec)
      dimension iatnum(n),c(3,n),index(n),nneig(n),ineig(maxng,n),
     -  isegno(n),cell(3,ncell),edge(3),ixyzhex(3),molsltlim(3,nmolslt)
      common /connatdat/ ramax(99),ramax2(99),hlimfac,ianfg(99),
     -  namfcg(100),nrmw
      character*2 iatnm2
      common /atnam/ aw(99),nval(99),nvalmax(99),mmofan(99),
     -  mmatno(64),iatnm2(99)
      dimension rij(3)
      character*4 pbclab
      character*4 atnami
      character*8 resnami
      dimension rpbc(3)
c     Check solute topology
      nbshort=0
      rminshort=1000.0
      nerr=0
      do i=1,nslt
c       Check bonds if they are too short
        if ((iatnum(i) .lt. 88 .or. iatnum(i) .gt. 90)
     -      .and. isegno(i) .ge. 0) then
          if (inamcol2 .ge. inamcol1) then
            resnami='     '
            resnami=line(index(i))(irescol1:irescol2)
            atnami=line(index(i))(inamcol1:inamcol2)
          end if
          do jj=1,nneig(i)
            j=ineig(jj,i)
            if (j .gt. i) then
              r2=dist2(c(1,i),c(1,j))
              call decidebondcut(iatnum(i),iatnum(j),rlim)
              rlim=rlim*bondminfac
              if (r2 .lt. rlim) then
c               Too short bond found.
                nerr=nerr+1
                r=sqrt(r2)
                if (inamcol2 .lt. inamcol1) then
                  write (iout,2002) i,iatnm2(iatnum(i)),j,
     -              iatnm2(iatnum(j)),r,sqrt(rlim)
                else
                  write (iout,2003) i,resnami,atnami,j,
     -              line(index(j))(irescol1:irescol2),
     -              line(index(j))(inamcol1:inamcol2),
     -              r,sqrt(rlim)
                end if
                nbshort=nbshort+1
                if (rminshort .gt. r) rminshort=r
              end if
            end if
          end do
        end if
      end do
      if (ioppbc .ge. 0) call distminsetup(edge,ioppbc)
      ncsltintra=0
      ncsltinter=0
      ncsltinterpbc=0
      rminsltintra=1000.0
      rminsltinter=1000.0
      rminsltinterpbc=1000.0
      do im=1,nmolslt
c       Check for close approach
        do i=molsltlim(1,im),molsltlim(2,im)
          if ((iatnum(i) .lt. 88 .or. iatnum(i) .gt. 90)
     -        .and. isegno(i) .ge. 0) then
            if (inamcol2 .ge. inamcol1) then
              resnami='     '
              resnami=line(index(i))(irescol1:irescol2)
              atnami=line(index(i))(inamcol1:inamcol2)
            end if
c           Check solute intramolecular NB distances
            do j=i+1,molsltlim(2,im)
              if ((iatnum(j) .lt. 88 .or. iatnum(j) .gt. 90)
     -            .and. isegno(j) .ge. 0) then
                ibonded=0
                do in=1,nneig(i)
                  if (j .eq. ineig(in,i)) ibonded=1
                end do
                if (ibonded .eq. 0) then
                  r2=dist2(c(1,i),c(1,j))
                  call decidebondcut(iatnum(i),iatnum(j),rlim)
                  rlimct=rlim*ctlimfac
                  if (r2 .le. rlimct) then
c                   Contact between atoms found.
                    nerr=nerr+1
                    r=sqrt(r2)
                    if (inamcol2 .lt. inamcol1) then
                      write (iout,2000) 'Slt NB ',i,
     -                  iatnm2(iatnum(i)),j,iatnm2(iatnum(j)),
     -                  r,sqrt(rlimct)
                    else
                      write (iout,2001) 'Slt NB ',i,resnami,atnami,j,
     -                  line(index(j))(irescol1:irescol2),
     -                  line(index(j))(inamcol1:inamcol2),
     -                  r,sqrt(rlimct)
                    end if
                    ncsltintra=ncsltintra+1
                    if (rminsltintra .gt. r) rminsltintra=r
                  end if
                end if
              end if
            end do
c           Check solute intermolecular distances
            do jm=im+1,nmolslt
              do j=molsltlim(1,jm),molsltlim(2,jm)
                if ((iatnum(j) .lt. 88 .or. iatnum(j) .gt. 90)
     -              .and. isegno(j) .ge. 0) then
                  pbclab='    '
                  if (ioppbc .ge. 0) then
c                   Check for images
                    call arrdiff(c(1,i),c(1,j),rpbc,3)
                    call distmincalc(ioppbc,cell,ncell,ixyzhex,edge,
     -                rpbc(1),rpbc(2),rpbc(3),icmin,r2)
                    if (icmin .gt. 1) pbclab=' PBC'
                  else
                    r2=dist2(c(1,i),c(1,j))
                  end if
                  call decidebondcut(iatnum(i),iatnum(j),rlim)
                  rlimct=rlim*ctlimfac
                  if (r2 .le. rlimct) then
c                   Contact between atoms found.
                    nerr=nerr+1
                    r=sqrt(r2)
                    if (inamcol2 .lt. inamcol1) then
                      write (iout,2000) 'Slt-Slt',i,
     -                  iatnm2(iatnum(i)),j,iatnm2(iatnum(j)),
     -                  r,sqrt(rlimct),pbclab
                    else
                      write (iout,2001) 'Slt-Slt',i,resnami,
     -                  atnami,j,line(index(j))(irescol1:irescol2),
     -                  line(index(j))(inamcol1:inamcol2),
     -                  r,sqrt(rlimct),pbclab
                    end if
                    if (pbclab .eq. ' PBC') then
                      ncsltinterpbc=ncsltinterpbc+1
                      if (rminsltinterpbc .gt. r) rminsltinterpbc=r
                    else
                      ncsltinter=ncsltinter+1
                      if (rminsltinter .gt. r) rminsltinter=r
                    end if
                    if (ioppbc .gt. 0) then
                      call arrdiff(c(1,i),c(1,j),rij,3)
                      call genimdist(rij,cell,1,ncell,icminn,rmin2)
                      if (abs(rmin2-r2) .gt. 0.1 .or.
     -                   icmin .eq. 1 .and. icminn .gt. 1 .or.
     -                   icmin .gt. 1 .and. icminn .eq. 1)
     -                  write (iout,2009) sqrt(rmin2),r2,icmin,icminn
                    end if
                  end if
                end if
              end do
            end do
          end if
        end do
      end do
      if (isltonly .eq. 1) go to 999
c     Check for solute-solvent contacts
      nsv=(n-nslt)/naslv
      ncsltslv=0
      ncsltslvpbc=0
      rminsltslv=1000.0
      rminsltslvpbc=1000.0
      do ia=1,nslt
        do jm=1,nsv
          pbclab='    '
          icmin=1
          jm0=nslt+(jm-1)*naslv
          call arrdiff(c(1,ia),c(1,jm0+1),rpbc,3)
          if (ioppbc .ge. 0) then
c           Check for images
            call distmincalc(ioppbc,cell,ncell,ixyzhex,edge,
     -        rpbc(1),rpbc(2),rpbc(3),icmin,r2)
            if (icmin .gt. 1) pbclab=' PBC'
          else
            r2=rpbc(1)**2+rpbc(2)**2+rpbc(3)**2
          end if
          isvclose=1
          if (r2 .lt. 25.0) then
c           General image search to get neighbor cell no (could be improved)
            if (pbclab .eq. ' PBC') call genimdist(rpbc,cell,
     -        1,ncell,icminw,rmin2)
            do ja=2,naslv
              call arrdiff(c(1,ia),c(1,jm0+ja),rpbc,3)
              if (pbclab .eq. ' PBC') then
                do k=1,3
                  rpbc(k)=rpbc(k)-cell(k,icminw)
                end do
              end if
              r2j=rpbc(1)**2+rpbc(2)**2+rpbc(3)**2
              if (r2j .lt. r2) then
                r2=r2j
                isvclose=ja
              end if
            end do
            ja=jm0+isvclose
            call decidebondcut(iatnum(ia),iatnum(ja),rlim)
            rlimct=rlim*ctlimfac
            if (r2 .lt. rlimct) then
              r=sqrt(r2)
              nerr=nerr+1
              if (inamcol2 .lt. inamcol1) then
                write (iout,2005) ia,iatnm2(ia),jm,ja,iatnm2(ja),
     -           r,sqrt(rlimct),pbclab
              else
                write (iout,2006) ia,line(index(ia))(irescol1:irescol2),
     -            line(index(ia))(inamcol1:inamcol2),ja,
     -            iatnm2(iatnum(ja)),jm,r,sqrt(rlimct),pbclab
              end if
              if (pbclab .eq. ' PBC') then
                ncsltslvpbc=ncsltslvpbc+1
                if (rminsltslvpbc .gt. r) rminsltslvpbc=r
              else
                ncsltslv=ncsltslv+1
                if (rminsltslv .gt. r) rminsltslv=r
              end if
              if (ioppbc .gt. 0) then
                call arrdiff(c(1,ia),c(1,ja),rij,3)
                call genimdist(rij,cell,1,ncell,icminn,rmin2)
                if (abs(rmin2-r2) .gt. 0.1 .or.
     -             icmin .eq. 1 .and. icminn .gt. 1 .or.
     -             icmin .gt. 1 .and. icminn .eq. 1)
     -            write (iout,2009) sqrt(rmin2),r2,icmin,icminn
              end if
            end if
          end if
        end do
      end do
c     Check for solvent-solvent contacts
      ncslvslv=0
      ncslvslvpbc=0
      rminslvslv=1000.0
      rminslvslvpbc=1000.0
      do im=1,nsv
        im0=nslt+(im-1)*naslv
        do jm=im+1,nsv
          jm0=nslt+(jm-1)*naslv
          call arrdiff(c(1,im0+1),c(1,jm0+1),rpbc,3)
          icmin=1
          pbclab='    '
          if (ioppbc .ge. 0) then
c           Check for images
            call distmincalc(ioppbc,cell,ncell,ixyzhex,edge,
     -        rpbc(1),rpbc(2),rpbc(3),icmin,r2)
            if (icmin .gt. 1) pbclab=' PBC'
          else
            r2=rpbc(1)**2+rpbc(2)**2+rpbc(3)**2
          end if
          isvclose=1
          jsvclose=1
          if (r2 .lt. 25.0) then
            if (pbclab .eq. ' PBC') call genimdist(rpbc,cell,
     -        1,ncell,icminw,rmin2)
            do ia=1,naslv
              do ja=1,naslv
                call arrdiff(c(1,im0+ia),c(1,jm0+ja),rpbc,3)
                if (pbclab .eq. ' PBC') then
                  do k=1,3
                    rpbc(k)=rpbc(k)-cell(k,icminw)
                  end do
                end if
                r2j=rpbc(1)**2+rpbc(2)**2+rpbc(3)**2
                if (r2j .lt. r2) then
                  r2=r2j
                  isvclose=ia
                  jsvclose=ja
                end if
              end do
              iac=im0+isvclose
              jac=jm0+jsvclose
              call decidebondcut(iatnum(iac),iatnum(jac),rlim)
              rlimct=rlim*ctlimfac
              r=sqrt(r2)
              if (r2 .lt. rlimct) then
                 nerr=nerr+1
                write (iout,2007) iac,iatnm2(iatnum(iac)),im,jac,
     -            iatnm2(iatnum(jac)),jm,r,sqrt(rlimct),pbclab
                if (pbclab .eq. ' PBC') then
                  ncslvslvpbc=ncslvslvpbc+1
                  if (rminslvslvpbc .gt. r) rminslvslvpbc=r
                else
                  ncslvslv=ncslvslv+1
                  if (rminslvslv .gt. r) rminslvslv=r
                end if
                if (ioppbc .gt. 0) then
                  call arrdiff(c(1,iac),c(1,jac),rij,3)
                  call genimdist(rij,cell,1,ncell,icminn,rmin2)
                  if (abs(rmin2-r2) .gt. 0.1 .or.
     -               icmin .eq. 1 .and. icminn .gt. 1 .or.
     -               icmin .gt. 1 .and. icminn .eq. 1)
     -              write (iout,2009) sqrt(rmin2),r2,icmin,icminn
                end if
              end if
            end do
          end if
        end do
      end do
999   if (rminshort .eq. 1000.0) rminshort=0.0
      if (nbshort .gt. 0) then
        write (iout,2008) 'short bonds',nbshort,rminshort
        write (6,2008) 'short bonds',nbshort,rminshort
      else
        write (iout,*) 'No short bonds found'
      end if
      if (ncsltintra .gt. 0) then
        write (iout,2008) 'solute intramolecular clashes',
     -    ncsltintra,rminsltintra
        write (6,2008) 'solute intramolecular clashes',
     -    ncsltintra,rminsltintra
        else
          write (iout,*) 'No solute intramolecular clashes found'
        end if
      if (nmolslt .gt. 1) then
        if (ncsltinter .gt. 0) then
          write (iout,2008) 'solute intermolecular non-PBC clashes',
     -      ncsltinter,rminsltinter
          write (6,2008) 'solute intermolecular non-PBC clashes',
     -      ncsltinter,rminsltinter
        else
          write (iout,*) 'No solute intramolecular clashes found'
        end if
        if (ioppbc .ge. 0)  then
          if (ncsltinterpbc .gt. 0) then
            write (iout,2008) 'solute intermolecular PBC clashes',
     -        ncsltinterpbc,rminsltinterpbc
            write (6,2008) 'solute intermolecular PBC clashes',
     -        ncsltinterpbc,rminsltinterpbc
          else
            write (iout,*) 'No solute intermolecular PBC clashes found'
          end if
        end if
      end if
      if (isltonly .eq. 0) then
        if (ncsltslv .gt. 0) then
          write (iout,2008) 'solute -solvent non-PBC clashes',
     -      ncsltslv,rminsltslv
          write (6,2008) 'solute-solvent non-PBC clashes',
     -      ncsltslv,rminsltslv
        else
          write (iout,*) 'No solute-solvent non-PBC clashes found'
        end if
        write (iout,2008) 'solvent-solvent non-PBC clashes',
     -    ncslvslv,rminslvslv
        if (ioppbc .ge. 0) then
          if (ncsltslvpbc .gt. 0) then
            write (iout,2008) 'solute -solvent PBC clashes',
     -        ncsltslvpbc,rminsltslvpbc
            write (6,2008) 'solute-solvent PBC clashes',
     -        ncsltslvpbc,rminsltslvpbc
          else
            write (iout,*) 'No solute-solvent PBC clashes found'
          end if
          if (ncslvslvpbc .gt. 0) then
            write (iout,2008) 'solvent-solvent PBC clashes',
     -        ncslvslvpbc,rminslvslvpbc
            write (6,2008) 'solvent-solvent PBC clashes',
     -        ncslvslvpbc,rminslvslvpbc
          else
            write (iout,*) 'No solvent-solvent PBC clashes found'
          end if
        end if
      end if
      return
2000  format(1x,a,' clash: d[',i5,' (',a2,') - ',i5,'(',a2,')]=',
     -  f4.2,' A (<',f4.2,')',a)
2001  format(1x,a,' clash: d[',i5,' (',a4,1x,a5,') - ',i5,
     -  ' (',(a4,1x,a4),')]=',f4.2,' A (<',f4.2,')',a)
2002  format(' Too short bond: d[',i5,' (',a2,') - ',i5,'(',a2,')]=',
     -  f6.2,' A (<',f4.2,')')
2003  format(' Too short bond: d[',i5,' (',a4,1x,a5,') - ',
     -  i5,'(',a4,1x,a5,')]=',f4.2,' A (<',f4.2,')')
2005  format(' Solvent clash: d(',i5,1x,a2,' - ',
     -  i7,1x,a2,' slv #',i6,')=',f4.2,' A (<',f4.2,')',a)
2006  format(' Solvent clash: d(',i5,1x,a4,1x,a4,' - ',
     -  i7,1x,a2,' slv #',i6,')=',f4.2,' A (<',f4.2,')',a)
2007  format(' Solvent clash: d(',i7,1x,a2,' #',i6,' - ',
     -  i7,1x,a2,' #',i6,')=',f4.2,' A (<',f4.2,')',a)
2008  format(' !!! Number of ',a,' found=',i7,' Min=',f4.2,' A')
2009  format(' PROGRAM ERROR: r(check)=',f9.4,' r(calc)=',f9.4,
     -  ' icmin,n=',2i3)
      end
