      subroutine volcalc(nrand,c,na,isegno,ian,nsegslt,molsltlim,
     -  icellno,rvdw,rvdwsolv2,rvdwsolv,rvdwsolv22,vmac,vshell,vint,
     -  vmacsd,vshellsd,vintsd,vfrstsh,vfrstshsd,rsolv,ixcell,index,
     -  ifirst,ilast,itemp1,itemp2,nconfig,iout,levout,maxrsd,maxrec)
c     Perform the Monte Carlo estimate of ligand and interfacial volumes
      dimension c(3,na),isegno(na),ian(na),rvdw(na),icellno(na),
     -  rvdwsolv2(na),rvdwsolv(na),rvdwsolv22(na),molsltlim(3,maxrsd),
     -  ixcell(na),index(maxrec),ifirst(maxrec),ilast(maxrec),
     -  itemp1(maxrec),itemp2(maxrec)
      parameter (MAXINTVOL=10)
      dimension corner(3),edge(3),rand(3),r(3),ri(3),rj(3),xyzmin(3),
     -  xyzmax(3),ecell(3),rvdwan(99),rijexmn(1000),
     -  rij2mnmx(1000),iamn(1000),ing(MAXINTVOL),
     -  nintij(MAXINTVOL,MAXINTVOL),nxyz(3),ixyz(3)
      character*80 line
      data rvdwan /1.2,4*0.0,2.0,1.5,1.4,1.35,5*0.0,1.9,1.85,1.8,
     -  17*0.0,1.95,64*0.0/
      rsolv2=2.0*rsolv
      rvdwmax=2.0
      padding=rvdwmax+rsolv
      spacing=rvdwmax+2.0*rsolv
      call cellpart(c,ian,itemp1,1,na,padding,spacing,corner,ecell,edge,
     -  xyzmin,xyzmax,vtot,nxyz,ixyz,icellno,itemp2,index,ifirst,ilast,
     -  ntotcell,levout,iout,maxrec)
c     Get a list of filled cells
      ncell=0
      do ic=1,ntotcell
        if (ifirst(ic) .ne. 0) then
          ncell=ncell+1
          ixcell(ncell)=ic
        end if
      end do
c     Mark neighbors of filled cells
      call zeroiti(itemp1,0,ntotcell)
      do ic=1,ncell
        icell=ixcell(ic)-1
        ixyz(1)=mod(icell,nxyz(1))
        icell=(icell-ixyz(1))/nxyz(1)
        ixyz(2)=mod(icell,nxyz(2))
        ixyz(3)=(icell-ixyz(2))/nxyz(2)
        do ix=max0(0,ixyz(1)-1),min0(nxyz(1)-1,ixyz(1)+1)
          do iy=max0(0,ixyz(2)-1),min0(nxyz(2)-1,ixyz(2)+1)
            do iz=max0(0,ixyz(3)-1),min0(nxyz(3)-1,ixyz(3)+1)
              icng=1+ix+nxyz(1)*iy+nxyz(1)*nxyz(2)*iz
              itemp1(icng)=1
            end do
          end do
        end do
      end do
      nrancell=0
      do ic=1,ntotcell
        if (itemp1(ic) .gt. 0) then
          nrancell=nrancell+1
          itemp1(nrancell)=ic
        end if
      end do
      if (nconfig .lt. 2) then
        write (iout,1003) ecell,ntotcell,ncell,nrancell,rsolv2
        write (6,1003) ecell,ntotcell,ncell,nrancell,rsolv2
        if (nsegslt .gt.  MAXINTVOL)  then
          write (iout,1005) nsegslt,MAXINTVOL
          write (6,1005) nsegslt,MAXINTVOL
        end if
      end if
      if (levout .gt. 0) write (iout,7000) nxyz,ntotcell,nrancell,ncell
      do ia=1,na
        rvdw(ia)=rvdwan(ian(ia))
        rvdwsolv(ia)=rvdw(ia)+rsolv
        rvdwsolv2(ia)=(rvdw(ia)+rsolv)**2
        rvdwsolv22(ia)=(rvdw(ia)+2.0*rsolv)**2
        if (levout .gt. 4)
     -    write (iout,7001) ia,ian(ia),rvdw(ia),isegno(ia)
      end do
      nmac=0
      nshell=0
      nfrstsh=0
      nint=0
      nmacprev=0
      nshellprev=0
      nfrstshprev=0
      nintprev=0
      nrandprev=0
      vmacsm=0.0
      vshellsm=0.0
      vfrstshsm=0.0
      vintsm=0.0
      vmacsm2=0.0
      vshellsm2=0.0
      vfrstshsm2=0.0
      vintsm2=0.0
      if (nsegslt .le. MAXINTVOL) call zeroiti(nintij,0,MAXINTVOL**2)
      nrand10=max0(1,nrand/10)
      do ir=1,nrand
        call randpx(3,rand)
        do ic=1,nrancell
          icell=itemp1(ic)-1
c         Shift random copy to cell icell and examine just its neighbor cells
c         Deconvolute icell
          ixyz(1)=mod(icell,nxyz(1))
          icell=(icell-ixyz(1))/nxyz(1)
          ixyz(2)=mod(icell,nxyz(2))
          ixyz(3)=(icell-ixyz(2))/nxyz(2)
          if (levout .gt. 3) write (iout,7003) ir,ic,icell,ixyz
          do k=1,3
            r(k)=corner(k)+(ixyz(k)+rand(k))*ecell(k)
          end do
          do im=1,nsegslt
            rijexmn(im)=100000.0
            rij2mnmx(im)=100000.0
          end do
          call zeroiti(iamn,0,nsegslt)
          do ix=max0(0,ixyz(1)-1),min0(nxyz(1)-1,ixyz(1)+1)
            do iy=max0(0,ixyz(2)-1),min0(nxyz(2)-1,ixyz(2)+1)
              do iz=max0(0,ixyz(3)-1),min0(nxyz(3)-1,ixyz(3)+1)
                icn=1+ix+nxyz(1)*iy+nxyz(1)*nxyz(2)*iz
                if (ifirst(icn) .gt. 0) then
c                 Cell is not empty
                  do iaa=ifirst(icn),ilast(icn)
                    ia=index(iaa)
                    rij2=dist2(r,c(1,ia))
                    if (levout .gt. 4)
     -                write (iout,7004) iaa,ia,ix,iy,iz,icn,rij2
                    im=isegno(ia)
                    if (rij2 .lt. rij2mnmx(im)) then
                      rijex=sqrt(rij2)-rvdw(ia)
                      if (rijex .lt. rijexmn(im)) then
                        rijexmn(im)=rijex
                        rij2mnmx(im)=(rvdwmax+rijex)**2
                        iamn(im)=ia
                      end if
                    end if
                  end do
                end if
              end do
            end do
          end do
          nng=0
          rijexmin=10000.0
          do im=1,nsegslt
            if (iamn(im) .gt. 0) then
              if (rijexmn(im) .lt. rijexmin) then
                rijexmin=rijexmn(im)
                iamin=iamn(im)
              end if
              rij2mnmx(im)=(rvdw(iamn(im))+rijexmn(im))**2
            end if
c           Check if potential interface neighbor
            if (iamn(im) .gt. 0) then
              if (rij2mnmx(im) .le. rvdwsolv22(iamn(im))) then
                nng=nng+1
                ing(nng)=im
              end if
            end if
          end do
          if (rijexmin .le. 0.0) then
c           Inside solute's vdW region
            nmac=nmac+1
          else if (rijexmin .le. rsolv) then
c           Within solvent-excluded shell
            nshell=nshell+1
            if (levout .eq. 2) then
              call blankout(line,1,80)
              line(1:4)='ATOM'
              line(22:22)='S'
              line(13:16)='He  '
              line(18:20)='SHL'
              write (line(7:11),2004) nshell
              write (line(23:26),2003) 1
              write (line(31:54),2002) r
              write (line(55:60),2001) 1.0
              write (line(61:66),2001) 0.0
c             write (77,2000) line(1:80)
            end if
          else if (rijexmin .le. rsolv2) then
c           Within first solvation shell
            nfrstsh=nfrstsh+1
          end if
          if (nng .gt. 1) then
c           Possible interface
            nrint=0
            do in=1,nng
              ia=iamn(ing(in))
              do jn=1,in-1
                ja=iamn(ing(jn))
                ijint=0
                if (rij2mnmx(in) .lt. rvdwsolv2(ia) .and.
     -              rij2mnmx(jn) .lt. rvdwsolv2(ja)) then
                  ijint=1
                else
                  rij2=dist2(c(1,ia),c(1,ja))
                  if (rij2 .le. (rvdwsolv(ia)+rvdwsolv(ja))**2) then
c                   Calculate the scalar product
                    rijsum=0.0
                    do k=1,3
                      ri(k)=c(k,ia)-r(k)
                      rj(k)=c(k,ja)-r(k)
                      rijsum=rijsum+ri(k)*rj(k)
                    end do
                    ac=rijsum/sqrt(abs(rij2mnmx(in)*rij2mnmx(jn)))
c                   Calculate the scalar threshold value
                    acmax=(rvdwsolv2(ia)+rvdwsolv2(ja)-rij2)/
     -                (2.0*rvdwsolv(ia)*rvdwsolv(ja))
                    if (ac .le. acmax) ijint=1
                  end if
                end if
                if (ijint .eq. 1) then
                  nrint=nrint+1
                  if (nsegslt .le. MAXINTVOL)
     -              nintij(in,jn)=nintij(in,jn)+1
                end if
              end do
            end do
            if (nrint .gt. 0) then
              nint=nint+1
              if (levout .eq. 3) then
                call blankout(line,1,80)
                line(1:4)='ATOM'
                line(22:22)='I'
                line(13:16)='He  '
                line(18:20)='INT'
                write (line(7:11),2004) nint
                write (line(23:26),2003) 1
                write (line(31:54),2002) r
                write (line(55:60),2001) 1.0
                write (line(61:66),2001) 0.0
                write (77,2000) line(1:80)
              end if
            end if
          end if
        end do
        if (mod(ir,nrand10) .eq. 0) then
          if (nconfig .eq. 1 .and. ir+nrand10 .le. nrand) then
            vshell=vtot*float(nshell)/float(ntotcell*ir)
            vfrstsh=vtot*float(nshell+nfrstsh)/float(ntotcell*ir)
            vmac=vtot*float(nmac)/float(ntotcell*ir)
            vint=vtot*float(nint)/float(ntotcell*ir)
            write (iout,1000) ir,vmac,vshell,vfrstsh,vint
          end if
          rndiff=float(ntotcell*(ir-nrandprev))
          vm10=vtot*float(nmac-nmacprev)/rndiff
          vmacsm=vmacsm+vm10
          vmacsm2=vmacsm2+vm10**2
          vs10=vtot*float(nshell-nshellprev)/rndiff
          vshellsm=vshellsm+vs10
          vshellsm2=vshellsm2+vs10**2
          vf10=vtot*float((nshell+nfrstsh)-(nshellprev+nfrstshprev))/
     -      rndiff
          vfrstshsm=vfrstshsm+vf10
          vfrstshsm2=vfrstshsm2+vf10**2
          vi10=vtot*float(nint-nintprev)/rndiff
          vintsm=vintsm+vi10
          vintsm2=vintsm2+vi10**2
          nmacprev=nmac
          nshellprev=nshell
          nfrstshprev=nfrstsh
          nintprev=nint
          nrandprev=ir
        end if
      end do
      vshell=vtot*float(nshell)/float(ntotcell*nrand)
      vfrstsh=vtot*float(nshell+nfrstsh)/float(ntotcell*nrand)
      vmac=vtot*float(nmac)/float(ntotcell*nrand)
      vint=vtot*float(nint)/float(ntotcell*nrand)
      vmacsd=sqrt(abs(vmacsm2/10.0-(vmacsm/10.0)**2)/9.0)
      vshellsd=sqrt(abs(vshellsm2/10.0-(vshellsm/10.0)**2)/9.0)
      vfrstshsd=sqrt(abs(vfrstshsm2/10.0-(vfrstshsm/10.0)**2)/9.0)
      vintsd=sqrt(abs(vintsm2/10.0-(vintsm/10.0)**2)/9.0)
      write (iout,1000) nrand,vmac,vshell,vfrstsh,vint,' ',
     -  vmacsd,vshellsd,vfrstshsd,vintsd
      if (nsegslt .le. MAXINTVOL .and. nsegslt .gt. 2) then
        write (iout,1001) (in,in=1,nsegslt)
        do in=1,nsegslt
          write (iout,1002) molsltlim(1,in),molsltlim(2,in),in,
     -      (vtot*float(nintij(in,jn))/float(ntotcell*nrand),jn=1,in)
        end do
      end if
      return
1000  format(' Nr=',i7,' Vslt=',f10.1,' Vxsh=',f9.1,' Vfsh=',f9.1,
     -  ' Vint=',f9.1,' A**2',a,/,
     -  11x,'   SD=',f10.1,'   SD=',f9.1,'   SD=',f9.1,'   SD=',f9.1)
1001  format(' Interface volume between solute molecules',/,
     -  3x,10i7)
1002  format(i6,' - ',i6,i3,10f7.1)
1003  format(' Copies of random points will be placed in ',
     -  'rectangles of dimension',/,f10.5,2(' * ',f10.5),' A**3',/,
     -  ' There are',i6,' such cells',i5,' of them contains solute',/,
     -  ' Translated copies will be placed into all occupied cells and',
     -  ' its neighbors,',/,' a total of',i6,' cells',
     -  ' Nr: result for Nr random numbers',/
     -  ' Vslt: estimate of solute volume (in A^3)',/,
     -  ' Vxsh: estimate of solvent-excludes shell volume (in A^3)',/,
     -  ' Vfsh: estimate of first solvation shell volume (in A^3)',/,
     -  '       shell thickness=',f5.2,' A',/,
     -  ' Vint: estimate of the interface volume volume (in A^3)',/)
1005  format(' NOTE: The number of solute molecules (',i5,') exceeds',
     -  i4,/,' - pairwise interface volumes will not be calculated')
2000  format(a)
2001  format(f6.3)
2002  format(3f8.3)
2003  format(i4)
2004  format(i5)
7000  format(' nyxz=',3i3,' ntotcell=',i4,' nrancell=',i4,' ncell=',i4)
7001  format(' ia,ian=',2i4,' rvdw(ia)=',f8.2,'isegno=',i4)
7003  format(' ir,ic,icell=',3i5,' ixyz=',3i3)
7004  format(' iaa,ia=',2i5,' ixyz=',3i3,' icn=',i5,' rij2=',e12.5)
      end
