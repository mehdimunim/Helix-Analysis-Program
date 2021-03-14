      subroutine printbondlist(iangpr,itorpr,n1,n2,nstta,cslt,nneig,
     -  ineig,line,inamcol1,inamcol2,irescol1,irescol2,iresncol1,
     -  iresncol2,index,cv,ichiral,maxneig,radtodeg,iout,maxrec)
c#    MMC routine  98/b lstmod: 06/21/05
c*****List bonds, angles, torsions as requested
      character* 132 line(maxrec)
      dimension cslt(3,nstta),nneig(nstta),ineig(maxneig,nstta),
     -  index(nstta),cv(nstta),ichiral(nstta)
      character*8 lab
      data lab /' CHIRAL '/
      dimension ang(28),iang(2,28)
      do ia=n1,n2
        if (nneig(ia) .ne. 1) then
          do in=1,nneig(ia)
            cv(in)=sqrt(dist2(cslt(1,ia),cslt(1,ineig(in,ia))))
            if (cv(in) .lt. 0.1) then
              write (iout,2265) ineig(in,ia)
            end if
          end do
          llab=1
          if (ichiral(ia) .gt. 0) llab=8
          write (iout,2148) ia,line(index(ia))(inamcol1:inamcol2),
     -      line(index(ia))(iresncol1:iresncol2),
     -      line(index(ia))(irescol1:irescol2),lab(1:llab),nneig(ia)
          write (iout,2149) (ineig(in,ia),
     -      line(index(ineig(in,ia)))(inamcol1:inamcol2),
     -      line(index(ineig(in,ia)))(irescol1:irescol2),
     -      line(index(ineig(in,ia)))(iresncol1:iresncol2),cv(in),
     -      in=1,nneig(ia))
          if (iangpr .eq. 1 .and. nneig(ia) .ge. 2) then
c           Generate angle list
            nnlim=min0(8,nneig(ia))
            nang=0
            do in1=1,nnlim
              do in2=in1+1,nnlim
                nang=nang+1
                ang(nang)=radtodeg*angleijk(cslt,nstta,ineig(in1,ia),ia,
     -            ineig(in2,ia),iout)
                iang(1,nang)=ineig(in1,ia)
                iang(2,nang)=ineig(in2,ia)
              end do
            end do
            write (iout,2142) (iang(1,i),ia,iang(2,i),ang(i),i=1,nang)
          end if
          if (itorpr .eq. 1) then
c           Generate torsion angle list
            do in=1,nneig(ia)
              inj=ineig(in,ia)
              if (inj .gt. in) then
                nnlim1=min0(8,nneig(ia))
                nnlim2=min0(8,nneig(inj))
                nang=0
                do in1=1,nnlim1
                  inn1=ineig(in1,ia)
                  if (inn1 .ne. inj) then
                    do in2=in1+1,nnlim2
                      inn2=ineig(in2,inj)
                      if (inn2 .ne. ia) then
                        nang=nang+1
c                       ang(nang)=dihangl(cslt(1,inn1),cslt(1,ia),
c    -                    cslt(1,inj),cslt(1,inn2),1,iout)
                        ang(nang)=dihangl(cslt,inn1,ia,inj,inn2,0,
     -                    maxrec)
                        iang(1,nang)=inn1
                        iang(2,nang)=inn2
                      end if
                    end do
                  end if
                end do
                write (iout,2182) (' ',iang(1,i),ia,inj,iang(2,i),
     -            radtodeg*ang(i),i=1,nang)
              end if
            end do
          end if
        end if
      end do
      return
2142  format(2('  Angle(',i5,'-',i5,'-',i5,')=',f8.3,' deg'))
2148  format(i6,' (',a,1x,a,1x,a,')',a,'-',i3,' neighbours:')
2149  format((7x,2(i5,' (',a,1x,a,1x,a,') r=',f6.4,' A')))
2182  format(2(a,' Torsion(',i5,'-',i5,'-',i5,'-',i5,')=',f7.2,' deg'))
2265  format(' ----- WARNING: unphysically short bond with atom',i7)
      end
