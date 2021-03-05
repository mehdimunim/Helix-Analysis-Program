      subroutine hblist(c,n,nslt,islvw,naslv,iatnum,ifchrg,isegno,nneig,
     -  nneiga,nhbneig,ineig,nhneig,nnneig,ncneig,nsneig,npneig,ixres,
     -  nconfig,hbf0,angm0,molresflag,hblimfac,angmin,iwhbl,ic1,ic2,ir1,
     -  ir2,irn1,irn2,blankline,line,index,indexn,ibnd,indexo,isolvent,
     -  iprint,icallnn,nosameseg,nframe,radtodeg,maxrepconf,
     -  maxng,maxbox,maxrsd,maxrec)
      dimension c(3,n),iatnum(n),ifchrg(n),isegno(n),nneig(n),
     -  ineig(maxng,n),nhbneig(n),nneiga(n),nhneig(n),nnneig(n),
     -  ncneig(n),nsneig(n),npneig(n),ixres(n),index(n),indexn(n),
     -  ibnd(maxbox,maxrec),molresflag(maxrsd),indexo(maxrec)
      character* 132 line(maxrec),linep(25),blankline
      character*1 r1,bbsc(4)
      common /connatdat/ ramax(99),ramax2(99),hlimfac,ianfg(99),
     -  namfcg(100),nrmw
      character*8 resnam
      character*8 atname,atname1,atnamea,atname1a
      dimension rhbl(25),nh12(4)
      data bbsc/'B','S','V','?'/
      data pi /3.2141692/
c     print *,'HBLIST n,nslt,ic1,ic2=',n,nslt,ic1,ic2
      innl=0
      nhbats=nslt
      if (isolvent .eq. 1) nhbats=n
c     print *,'NHB=',nhbats
c     print *,'hbf0,angm0=',hbf0,angm0
c     print *,'hblimfac,angmin=',hblimfac,angmin
c     print *,'isegno(1),(n),icallnn=',isegno(1),isegno(n),icallnn
      if (hbf0 .ne. hblimfac .or. angm0 .ne. angmin .or.
     -    isegno(1) .ne. isegno(n) .or. icallnn .eq. 1)
     -  call nnlist(nslt,islvw,naslv,nhbats,iatnum,ifchrg,c,nneig,
     -    nneiga,nhbneig,ineig,nhneig,nnneig,ncneig,nsneig,npneig,line,
     -    ir1,ir2,ic1,ic2,index,nconfig,innl,molresflag,hblimfac,angmin,
     -    1,ibnd,indexo,isegno,ixres,maxrepconf,1,nframe,radtodeg,
     -    0,maxbox,maxng,maxrsd,maxrec)
      if (iprint .eq. 0) return
      rDHmin=200.0
      rDHmax=0.0
      nhats=0
c     Determine if backbone,sidechain or other
      atname='        '
      atname1='        '
      resnam='        '
      do i=1,nhbats
        resnam(1:ir2-ir1+1)=line(index(i))(ir1:ir2)
        call changeprot(resnam,r1,2)
        if (i .gt. nslt) then
          indexn(i)=3
        else if (r1 .eq. '*') then
          indexn(i)=4
        else
          indexn(i)=2
          atname(1:ic2-ic1+1)=line(index(i))(ic1:ic1+3)
          call leftadjustn(atname,atnamea,8)
          if (ineig(1,i) .gt. 0) then
            atname1(1:ic2-ic1+1)=line(index(ineig(1,i)))(ic1:ic1+3)
            call leftadjustn(atname1,atname1a,8)
            if (atnamea(1:3) .eq. 'H   ' .or.
     -          atnamea(1:3) .eq. 'HN ') then
              if (atname1a(1:3) .eq. 'N  ') indexn(i)=1
            else if (atnamea(1:3) .eq. 'O  ') then
              if (atname1a(1:3) .eq. 'C  ') indexn(i)=1
            end if
          end if
c          write (77,8871) i,atnamea,atname1a,indexn(i)
c8871      format(i5,' atnamea=',a,' atname1a=',a,' indexn(i)=',i2)
        end if
      end do
      call checkhblist(n,ineig,nhbneig,maxng)
      write (iwhbl,*) 'List of hydrogen bond(s) found:'
      nhbonds=0
      do ia=1,nhbats
        if (iatnum(ia) .eq. 1) then
c         Atom ia is always the donor H
          ibb=indexn(ia)
c         ibb,jbb: 1,2,3: Backbone/sidechain/solvent/?
          nhats=nhats+1
          nhbia=0
c         ihb0 is the heavy atom of the donor H
          call get_heavyat(ia,nneig,ineig,ixres,nframe,ihb0,maxng,
     -      maxrec)
          do ja=1,nhbneig(ia)
            ihb=ineig(maxng+1-ja,ia)
            if (nosameseg .eq. 0 .or. isegno(ia) .ne. isegno(ihb)) then
              nhbia=nhbia+1
              jbb=indexn(ihb)
              if (ihb .gt. n)
     -          print *,'ERROR: illegal ihb (',ihb,') at ia=',ia
              call angdistw(c(1,ia),c(1,ihb),c(1,ihb0),rHB,rb,rab,ang)
              ang=ang*(180.0/pi)
              rhbl(ja)=rHB
              call readint(line(index(ia)),irn1,irn2,irna,2,1,irerr)
              call readint(line(index(ihb)),irn1,irn2,irnb,2,1,irerr)
              linep(ja)=blankline
              write (linep(ja),2040) ia,line(index(ia))(ic1:ic2),irna,
     -          line(index(ia))(ir1:ir2),ihb,line(index(ihb))(ic1:ic2),
     -          irnb,line(index(ihb))(ir1:ir2),rHB,ang,rab,
     -          bbsc(ibb),bbsc(jbb)
              if (rDHmin .gt. rHB) then
                rDHmin=rHB
                imn1=ia
                imn2=ineig(maxng+1-ja,ia)
              end if
              if (rDHmax .lt. rHB) then
                rDHmax=rHB
                imx1=ia
                imx2=ineig(maxng+1-ja,ia)
              end if
              nhbonds=nhbonds+1
            end if
          end do
          do ja=1,min0(25,nhbneig(ia))
            call writeline(iwhbl,linep(ja),1,79,0)
          end do
        else if (ifchrg(ia) .gt. 0) then
c         Cation - only H-bonds to water
          nhbng=nhbneig(ia)
          do ja=1,nhbneig(ia)
            ihb=ineig(maxng+1-ja,ia)
            if (ihb .le. nslt) print *,'PROGRAM ERROR: atom ',ia,
     -        ' is a cation but is ','H-bonded to a non-solvent'
            if (iatnum(ihb) .ne. 8) print *,'PROGRAM ERROR: atom ',ia,
     -        ' is H-bonded to a cation but is not an oxygen'
            nnh=0
            do in=1,nneig(ihb)
              jn=ineig(in,ihb)
              if (iatnum(jn) .eq. 1) then
                nnh=nnh+1
                nh12(nnh)=jn
              end if
            end do
            if (nnh .ne. 2) then
              if (nframe .eq. 0) write (6,2000) ihb,nnh
              if (nframe .gt. 0) write (6,2000) ihb,nnh,' ',nframe
            end if
c           write (77,*) 'ia,ihb,nh12=',ia,ihb,nh12
            call angdistw(c(1,ihb),c(1,ia),c(1,nh12(1)),rHB,rb,rab,ang)
            ang=ang*radtodeg
            if (ang .lt. angmin) then
              call angdistw(c(1,ihb),c(1,ia),c(1,nh12(2)),rHB,rb,rab,
     -          ang)
              ang=ang*radtodeg
            end if
            call readint(line(index(ia)),irn1,irn2,irna,2,1,irerr)
            call readint(line(index(ihb)),irn1,irn2,irnb,2,1,irerr)
            linep(ja)=blankline
            write (linep(ja),2041) ia,line(index(ia))(ic1:ic2),irna,
     -        line(index(ia))(ir1:ir2),ihb,line(index(ihb))(ic1:ic2),
     -        irnb,line(index(ihb))(ir1:ir2),rHB,ang
            if (rDHmin .gt. rHB) then
              rDHmin=rHB
              imn1=ia
              imn2=ineig(maxng+1-ja,ia)
            end if
            if (rDHmax .lt. rHB) then
              rDHmax=rHB
              imx1=ia
              imx2=ineig(maxng+1-ja,ia)
            end if
            nhbonds=nhbonds+1
          end do
          do ja=1,min0(25,nhbneig(ia))
            call writeline(iwhbl,linep(ja),1,79,0)
          end do
        end if
      end do
      if (nconfig .le. maxrepconf) then
        if (nhbonds .gt. 0) then
          write (6,2049) nconfig,nhbonds
          write (6,2043) 'Shortest',rDHmin,imn1,imn2
          write (6,2043) 'Longest ',rDHmax,imx1,imx2
          write (iwhbl,2049) nconfig,nhbonds
          write (iwhbl,2043) 'Shortest',rDHmin,imn1,imn2
          write (iwhbl,2043) 'Longest ',rDHmax,imx1,imx2
        else
          if (nhats .eq. 0) then
            print *,'WARNING: structure contains NO hydrogens'
          else
            print *,'No hydrogen bonds found - try with a larger factor'
          end if
        end if
      end if
      return
2000  format(' WARNING: water oxygen ',i6,' has',i3,' hydrogens bonded',
     -  ' to it',a,' Nframe=',i6)
2040  format(i5,1x,a4,i4,1x,a4,' -',i5,1x,a4,i4,1x,a4,
     -  ' rHA=',f4.2,' a(D-H..A)=',f6.1,' rDA=',f4.2,1x,a1,'-',a1)
2041  format(i5,1x,a4,i4,1x,a4,' -',i5,1x,a4,i4,1x,a4,
     -  ' rO+=',f4.2,' a(+..O-H)=',f6.1)
2043  format(1x,a8,' A...H bond=',f8.3,' between atoms',i5,' and',i5)
2049  format(' Frame ',i6,' Number of hydrogen bonds found=',i5)
      end
