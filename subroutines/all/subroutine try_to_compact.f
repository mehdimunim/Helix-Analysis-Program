      subroutine try_to_compact(c,ctemp,molsltlim,edge,nslt,nmolslt,
     -  icompact)
      dimension c(3,nslt),ctemp(3,nslt),molsltlim(3,nmolslt),edge(3)
      dimension shift(3)
      do im=1,nmolslt
        ifst=molsltlim(1,im)
        ilst=molsltlim(2,im)
        icomp=1
        do while (icomp .gt. 0)
          icomp=0
          call countout_rect(ifst,ilst,c,edge,nouttot,kmax,nslt)
          if (float(nouttot)/float(ilst-ifst+1) .gt. 0.1) then
c           Translate
            if (kmax .lt. 0) then
              kmax=-kmax
              fact_s=1.0
            else
              fact_s=-1.0
            end if
            call zeroit(shift,3)
            shift(kmax)=edge(kmax)
            call shiftmol(c(1,ifst),ilst-ifst+1,shift,ctemp(1,ifst),
     -        fact_s)
c           Check if ctemp is more compact
            call countout_rect(ifst,ilst,ctemp,edge,nouttotnew,kmax,
     -        nslt)
            if (nouttotnew .lt. nouttot) then
              call trnsfr(c(1,ifst),ctemp(1,ifst),3*(ilst-ifst+1))
              icompact=1
              icomp=1
            end if
          end if
        end do
      end do
      return
      end
