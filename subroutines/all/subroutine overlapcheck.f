      subroutine overlapcheck(evecsprev,evecs,overlap,index,nneg)
      dimension evecsprev(3,3),evecs(3,3),overlap(3,3),index(3)
      dimension absoverlap(3,3)
      data imax1 /0/,jmax1 /0/,imax2 /0/,jmax2 /0/
      call zeroit(overlap,9)
      do i=1,3
        do j=1,3
          do k=1,3
            overlap(i,j)=overlap(i,j)+evecsprev(i,k)*evecs(j,k)
          end do
          absoverlap(i,j)=abs(overlap(i,j))
        end do
      end do
      call zeroiti(index,0,3)
      omax=0.0
      do i=1,3
        do j=1,3
          if (omax  .lt. absoverlap(i,j)) then
            omax=absoverlap(i,j)
            imax1=i
            jmax1=j
          end if
        end do
      end do
      index(imax1)=jmax1
      omax=0.0
      do i=1,3
        if (i .ne. imax1) then
          do j=1,3
            if (j .ne. jmax1) then
              if (omax  .lt. absoverlap(i,j)) then
                omax=absoverlap(i,j)
                imax2=i
                jmax2=j
              end if
            end if
          end do
        end if
      end do
      index(imax2)=jmax2
      do i=1,3
        if (index(i) .eq. 0) index(i)=6/(jmax1*jmax2)
      end do
      nneg=0
      do i=1,3
        if (overlap(index(i),i) .lt. 0.0) nneg=nneg+1
      end do
      return
      end
