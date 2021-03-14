      subroutine clusterplot(ips,xtraj,value,ifirst,ilast,ncl,ixclst,
     -  nframe,xtrajlab,lxtrajlab,ipspage,noclose,mx2d)
      dimension xtraj(mx2d),value(mx2d),ifirst(mx2d),ilast(mx2d),
     -  ixclst(mx2d)
      character*(*) xtrajlab
c     Plot the history of cluster membership
c     print *,'CLUSTERPLOT ncl=',ncl,' nframe=',nframe,' ips=',ips
      if (nframe .gt. mx2d) then
        write (6,1000) nframe,mx2d
        return
      end if
      do ic=1,ncl
        do ia=ifirst(ic),ilast(ic)
          if (ixclst(ia) .lt. 1 .or. ixclst(ia) .gt. nframe) then
            print *,'CLUSTER ERROR ic,ia,ixclst(ia)=',ic,ia,ixclst(ia)
          else
            value(ixclst(ia))=ic
          end if
        end do
      end do
      iprt=0
      call plot2fun(ips,1,xtraj,value,value,nframe,0.0,0.0,0,0.0,
     -  1.0,ncl+1,0.0,0.0,0,'Cluster membership',18,' ',0,
     -  xtrajlab,lxtrajlab,'Cluster #',9,' ',0,1,iprt,
     -  6,1,1,0,0,0,0,ipspage,noclose,1,1)
      return
1000  format(' PROGRAM ERROR in CLUSTERPLOT: nframe (',i8,') > ',
     -  'max2d (',i5,')')
      end
