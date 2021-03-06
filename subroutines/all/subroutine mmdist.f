      subroutine mmdist(c,n,atw,iatnum,nmolslt,molsltlim,c2,dist,
     -  ignoreh,iw0,iw1)
      dimension c(3,n),atw(n),iatnum(n),molsltlim(3,n),c2(3,n),dist(n)
      real*8 c0(3),atwsum
c     print *,'MMDIST nmolslt,nmc,iw0,iw1=',nmolslt,nmc,iw0,iw1
      do im=1,nmolslt
        call zeroitd(c0,3)
        atwsum=0.d0
        do ia=molsltlim(1,im),molsltlim(2,im)
          atwsum=atwsum+atw(ia)
          do k=1,3
            c0(k)=c0(k)+atw(ia)*c(k,ia)
          end do
          do k=1,3
            c0(k)=c0(k)+atw(ia)*c(k,ia)
          end do
        end do
        do k=1,3
          c2(k,im)=c0(k)/atwsum
        end do
      end do
      do im=1,nmolslt
        do imm=1,nmolslt
          dist(imm)=sqrt(dist2(c2(1,im),c2(1,imm)))
        end do
        write (iw0,1000) im,(dist(imm),imm=1,nmolslt)
        do imm=1,nmolslt
          rdist2=999999.0
          do iam=molsltlim(1,im),molsltlim(2,im)
            if (ignoreh .eq. 0) then
              do iamm=molsltlim(1,imm),molsltlim(2,imm)
                d2=dist2(c(1,iam),c(1,iamm))
                if (d2 .lt. rdist2) rdist2=d2
              end do
            else
              if (iatnum(iam) .gt. 1) then
                do iamm=molsltlim(1,imm),molsltlim(2,imm)
                  if (iatnum(iamm) .gt. 1) then
                    d2=dist2(c(1,iam),c(1,iamm))
                    if (d2 .lt. rdist2) rdist2=d2
                  end if
                end do
              end if
            end if
          end do
          dist(imm)=sqrt(rdist2)
        end do
        write (iw1,1000) im,(dist(imm),imm=1,nmolslt)
      end do
      return
1000  format(' im=',i4,' D=',5f12.5,/,(11x,5f12.5))
      end
