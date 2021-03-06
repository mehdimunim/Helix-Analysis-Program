      subroutine checkonedir(c,ctemp,itemp,ifst,ilst,edge,uxyz,
     -  vertex,cmin,cmax,c0,vol,iaxis,ishift,n)
      dimension c(3,n),ctemp(3,n),itemp(n),edge(3),uxyz(3,3),
     -  vertex(3,8),cmin(3),cmax(3),c0(3)
      dimension shift(3)
      print *,'CHECKONEDIR iaxis=',iaxis,' ifst,ilst=',ifst,ilst
      ishift=0
      nats=ilst-ifst+1
      do ix=-1,1,2
c       Translate back and forth along the x axis
        if (ishift .eq. 0) then
c         Generate translation along the iaxis-th axis
          do k=1,3
            shift(k)=edge(k)*uxyz(k,iaxis)
          end do
c         print *,'SHIFT=',shift
          call shiftmol(c(1,ifst),nats,shift,ctemp,float(ix))
          call trnsfr(ctemp(1,nats+1),vertex,24)
          call extension(ctemp,itemp,0,1,nats+8,cmin,cmax,c0,1,0,volnew)
          print *,'VOL,VOLNEW=',vol,volnew,' ix,iaxis=',ix,iaxis
          if (volnew/vol .lt. 0.9) then
            vol=volnew
            call trnsfr(c(1,ifst),ctemp,3*nats)
            print *,'MOVED ifst,ilst=',ifst,ilst
            ishift=1
          end if
        end if
      end do
      return
      end
