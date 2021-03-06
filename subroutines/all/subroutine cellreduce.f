      subroutine cellreduce(it,c,ih,n,nnh,edge,edgen,edge_gen,ioppbc,
     -  npbc,cell,ncell,cellalt,ixyzhex,closmn,volnew,nwnew)
      dimension c(3,n),ih(n),edge(3),edgen(3),edge_gen(3,3),cell(3,27),
     -  cellalt(3,27),ixyzhex(3)
c     When it=0, find out if any of the cell sizes can be further reduced
c     without affecting the closest approach
c     When it=1, scale down the cell until the minimum approach is reached
      dimension edgenn(3),delxyz(3)
c     print *,'CELLRED  edge=',edge
c     print *,'CELLRED  closmn=',closmn
      if (closmn .eq. 0.0) then
        call crorgn(edge,edge_gen,ioppbc,3,ncell,cell,cellalt,
     -    ixyzhex,rinscr,rcirc)
        closmin=distminimg(c,ih,n,nnh,edge,ioppbc,cell,ncell,ixyzhex,
     -    inn,jnn)
      else
        closmin=closmn
      end if
      call trnsfr(edgen,edge,3)
      call trnsfr(edgenn,edge,3)
      call zeroit(delxyz,3)
      delsum=0.0
      icmax=npbc
      if (it .eq. 1) icmax=1
      do ic=1,icmax
        del=0.01
        if (icmax .gt. 1) then
          edgenn(ic)=edge(ic)-del
        else
          do k=1,3
            edgenn(k)=edge(k)*(1.0-del/edge(1))
          end do
        end if
        call crorgn(edgenn,edge_gen,ioppbc,3,ncell,cell,cellalt,
     -    ixyzhex,rinscr,rcirc)
        dnew=distminimg(c,ih,n,nnh,edgenn,ioppbc,cell,ncell,ixyzhex,
     -    inn,jnn)
c       print *,'Cellred ic,closmin,dnew=',ic,closmin,dnew
        if (dnew .ge. closmin-1.e-4) then
c         Successful reduction - see how far it can be taken
          edgen(ic)=edgenn(ic)
          del=0.5
          nred=0
          maxred=25
          do while (del .gt. 0.005 .and. nred .lt. maxred)
            nred=nred+1
            if (icmax .gt. 1) then
              edgenn(ic)=edgen(ic)-del
            else
              do k=1,3
                edgenn(k)=edgen(k)*(1.0-del/edge(1))
              end do
            end if
c           write (6,8788) edgenn
c8788        format(' edgenn=',3f10.5)
            call crorgn(edgenn,edge_gen,ioppbc,3,ncell,cell,cellalt,
     -        ixyzhex,rinscr,rcirc)
            dnew=distminimg(c,ih,n,nnh,edgenn,ioppbc,cell,ncell,ixyzhex,
     -        inn,jnn)
c             print *,'dnew,del=',dnew,del
            if (dnew .gt. closmin-1.e-4) then
              if (icmax .gt. 1) then
                edgen(ic)=edgenn(ic)
                delxyz(ic)=delxyz(ic)+del
              else
                call trnsfr(edgen,edgenn,3)
                delsum=delsum+del
              end if
              del=2.0*del
            else
              del=del/2.0
            end if
          end do
        else if (icmax .gt. 1) then
          edgenn(ic)=edgenn(ic)+del
        end if
      end do
      if (icmax .gt. 1) then
        write (6,1000) (delxyz(k),k=1,3)
      else
        delfac=(1.0-delsum/edge(1))
        write (6,1001) delfac
      end if
      call prtcell(ioppbc,edgen,edge_gen,0.0,volnew,nwnew,0)
      return
1000  format(5x,' Anisotropic cell contraction achieved by:',3f6.2,' A')
1001  format(5x,' Isotropic cell contraction by a factor of ',f6.4)
      end
