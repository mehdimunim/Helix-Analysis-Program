      subroutine stat14(c,n,nneig,ineig,index,line,nconfig,pi,iatnm2,
     -  ir1,ir2,ic1,ic2,inpcrdtyp,ioins,iw0,maxrepconf,maxng,maxrec)
      dimension nneig(n),ineig(maxng,n),c(3,n),index(n)
      character*2 iatnm2(99)
      character* 132 line(maxrec)
      rmin=1000.0
      rmax=-rmin
      do ia=1,n
        do j=1,nneig(ia)
          ja=ineig(j,ia)
          if (ia .lt. ja) then
c           Bond ia-ja found
            do ii=1,nneig(ia)
              iaa=ineig(ii,ia)
              do jj=1,nneig(ja)
                jaa=ineig(jj,ja)
                if (jaa .ne. ia .and. iaa .ne. ja
     -              .and. iaa .ne. jaa) then
c                 iaa-ia-ja-jaa is a 1-4  chain
                  r2=sqrt(dist2(c(1,iaa),c(1,jaa)))
                  if (r2 .gt. rmax) then
                    rmax=r2
                    iamax=iaa
                    jamax=jaa
                  else if (r2 .lt. rmin) then
                    rmin=r2
                    iamin=iaa
                    jamin=jaa
                  end if
                  tors=dihangl(c,iaa,ia,ja,jaa,1,maxrec)*(180.0/pi)
                  if (inpcrdtyp .le. ioins) then
                    write (iw0,2042) iaa,line(index(iaa))(ic1:ic2),
     -               line(index(iaa))(ir1:ir2),
     -               ia,line(index(ia))(ic1:ic2),
     -               line(index(ia))(ir1:ir2),
     -               ja,line(index(ja))(ic1:ic2),
     -               line(index(ja))(ir1:ir2),
     -               jaa,line(index(jaa))(ic1:ic2),
     -               line(index(jaa))(ir1:ir2),r2,tors
                  else
                    write (iw0,2039) iaa,iatnm2(iaa),ia,iatnm2(ia),
     -               ja,iatnm2(ja),jaa,iatnm2(jaa),r2,tors
                  end if
                end if
              end do
            end do
          end if
        end do
      end do
      if (nconfig .le. maxrepconf)
     -  write (6,2056) rmin,iamin,jamin,rmax,iamax,jamax
      return
2039  format(3(i5,' (',a2,') - '),i5,' (',a2,'):',/,
     -  ' 1-4 distance=',f8.3,' torsion angle=',f10.3)
2042  format(3(i5,' (',a4,1x,a4,') - '),i5,' (',a4,1x,a4,'):',/,
     -  '1-4 distance=',f5.2,' torsion angle=',f7.1)
2056  format(' Shortest 1-4 distance=',f6.2,' between atoms ',2i6,/,
     -  ' Longest 1-4 distance= ',f6.2,' between atoms ',2i6)
      end
