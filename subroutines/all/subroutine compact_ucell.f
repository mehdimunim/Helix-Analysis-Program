      subroutine compact_ucell(c,ctemp,itemp,molsltlim,n,nmolslt,
     -  vertex,cmin,cmax,c0,edge,uxyz,nshift)
      dimension c(3,n),ctemp(3,n),itemp(n),molsltlim(3,nmolslt),
     -  edge(3),uxyz(3,3),vertex(3,8),cmin(3),cmax(3),c0(3)
      print *,'COMPACT_UCELL n,nmolslt=',n,nmolslt,' VOL=',vol
      icompact=1
      nshift=0
      do while (icompact .gt. 0)
        icompact=0
        do im=1,nmolslt
          ifst=molsltlim(1,im)
          ilst=molsltlim(2,im)
          nats=ilst-ifst+1
c         print *,'IFST,ILST=',ifst,ilst
          call trnsfr(ctemp,c(1,ifst),3*nats)
          call trnsfr(ctemp(1,nats+1),vertex,24)
          call extension(ctemp,itemp,0,1,nats+8,cmin,cmax,c0,1,0,vol)
          call checkonedir(c,ctemp,itemp,ifst,ilst,edge,uxyz,
     -      vertex,cmin,cmax,c0,vol,1,ixshift,n)
          call checkonedir(c,ctemp,itemp,ifst,ilst,edge,uxyz,
     -      vertex,cmin,cmax,c0,vol,2,iyshift,n)
          call checkonedir(c,ctemp,itemp,ifst,ilst,edge,uxyz,
     -      vertex,cmin,cmax,c0,vol,3,izshift,n)
          ishift=ixshift+iyshift+izshift
          print *,'ixshift,iyshift,izshift=',ixshift,iyshift,izshift
          icompact=icompact+ishift
          nshift=nshift+ishift
        end do
      end do
c     print *,'COMP NSHIFT=',nshift
      return
      end
