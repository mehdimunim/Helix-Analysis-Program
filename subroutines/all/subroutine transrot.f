      subroutine transrot(c,ct,atw,atwtemp,rot,xyzmin,c0,edge_gen,
     -  molsltlim,nmolslt,nmolsltnoion,minresflag,n,nslt,nslv,npbc,
     -  ixyzhex,nneig,ineig,iatnum,molresflag,ixres,atnames,latnam,
     -  resnames,lresnam,iresno,it1,it2,it3,iaskshift,nconf,nosegid,pi,
     -  maxrsd,maxneig,mxat)
      dimension c(3,mxat),ct(3,mxat),atw(mxat),atwtemp(mxat),
     -  edge_gen(3,3),nneig(mxat),iatnum(mxat),ineig(maxneig,mxat),
     -  molresflag(maxrsd),ixres(mxat),rot(3,3),xyzmin(3),it1(mxat),
     -  it2(mxat),it3(mxat),molsltlim(3,maxrsd),ixyzhex(3),crep(3),
     -  c0(3),iresno(mxat)
      character*8 atnames(mxat),resnames(maxrsd)
      character*1 xyz,ans
      common /axislab/ xyz(3)
      dimension edge(3),cell(3,27),cellalt(3,27),x(3),y(3),z(3),cent(3),
     -  c00(3)
      parameter (MAXBRANCH=200)
      dimension ia_branch(MAXBRANCH),ia_predef(3,MAXBRANCH),
     -  r_predef(MAXBRANCH),ba_predef(MAXBRANCH),ta_predef(MAXBRANCH),
     -  loopmem(MAXBRANCH)
c     print *,'TRANSROT maxrsd,mxat=',maxrsd,mxat
      ans=' '
      do while (ans .ne. 'q')
        if (iaskshift .eq. 1) then
          call quiz(ans,itranstyp,' ',' ',0,
     -      '(additional) conformation transformation type',45,0,5,6,
     -      00)
        else
          ans='r'
        end if
        if (ans .eq. 's' .or. ans .eq. 'p' .or. ans .eq. 'a' .or.
     -      ans .eq. 'h') then
c         Shift (first)
          shiftfac=1.0
          if (ans .eq. 's') then
            call getxyz('Value to add to the ',20,' coordinates',12,0.0,
     -        xyzmin,0,0)
          else if (ans .eq. 'p') then
            call getint('Index of atom to be shifted',27,1,1,n,ias,0)
            write (6,1008) 'selected',atnames(ias)(1:latnam),
     -        resnames(ixres(ias))(1:lresnam),iresno(ias)
            if (ias .gt. nslt)
     -        print *,'WARNING: atom ',ias,' is a solvent atom'
            call getxyz('New ',4,' coordinate of the selected atom',
     -        32,0.0,xyzmin,0,0)
            do k=1,3
              xyzmin(k)=xyzmin(k)-c(k,ias)
            end do
          else
            shiftfac=-1.0
            call quiz(ans,ictp,'g',' ',0,'molecular center',16,2,5,6,25)
            if (ans .eq. 'g') then
               call trnsfr(xyzmin,c0,3)
            else
              call cofms(c,xyzmin,n,atw)
            end if
          end if
          call shiftmol(c,n,xyzmin,c,shiftfac)
          write (6,1002) xyzmin
        end if
        if (ans .eq. 'c' .or. ans .eq. 'p') then
c         Now reset into the PBC cell
          call setpbccell(' ',0,edge,edge_gen,cell,ncell,cellalt,
     -      ixyzhex,npbc,ioppbc,iusepbc,vol,nw,rinscr,rcirc,1)
          call setrepats(nmolslt,molsltlim,nslt,nneig,ineig,
     -      molresflag,minresflag,it1,it2,maxneig,maxrsd)
          call setpbcdim(ioppbc,ixyzhex,ixyzexcld,ixyzincld,xyz)
          if (nosegid .eq. 1) print *,'WARNING: PDB file had no ',
     -      'segment ID - solute is considered a single molecule'
          if (ans .eq. 'c') then
c           Centering was requested
c           Set atwtmp to zero for molecular ions
            do ia=1,nslt
              atwtemp(ia)=atw(ia)
              if (molresflag(ixres(ia)) .gt. minresflag) atwtemp(ia)=0.0
            end do
            if (ioppbc .eq. 6 .or. ioppbc .eq. 7)
     -        call askyn('Do you want to try both TO orientations',39,
     -          -1,1,icellalt,137,0)
            call systemcenter(n,nmolslt,nmolsltnoion,molsltlim,c,ct,it1,
     -        atwtemp,cell,ncell,cellalt,icellalt,0,0,nslt,nslv,0,0,
     -        imcent,nconf,mxat,maxrsd)
          else
c           Reset after shift
            do is=1,nmolslt
              natoms=molsltlim(2,is)-molsltlim(1,is)+1
              if (molsltlim(3,is) .eq. 0) then
                call cofms(c(1,molsltlim(1,is)),crep,natoms,atw)
              else
                call trnsfr(crep,c(1,molsltlim(3,is)),3)
              end if
              call pbcreset(c(1,molsltlim(1,is)),natoms,crep,
     -          cell,ncell,0,0,img)
            end do
            nsw=(n-nslt)/nslv
            do iw=1,nsw
              call pbcreset(c(1,nslt+(iw-1)*nslv+1),nslv,
     -          c(1,nslt+(iw-1)*nslv+1),cell,ncell,0,0,img)
            end do
          end if
        else if (ans .eq. 'r') then
c         Rotate by an input angle
          call genrot(rot,pi,iax,angle)
          if (iax .eq. -1) then
            call unitmat(rot)
            ans='q'
          else
            call rotate_c(c,n,rot,c,'TRANSROT',8)
            if (iax .eq. 0) write (6,1000) rot
            if (iax .gt. 0) write (6,1001) xyz(iax),angle
          end if
        else if (ans .eq. 'b') then
          call getint('Atom index to put at the origin',31,1,1,nslt,
     -      iacent,000)
          call getint('Atom index to put on the X axis',31,1,1,nslt,
     -      iaxax,000)
          call getint('Atom index to put in the X-Y plane',34,1,1,nslt,
     -      iaxyplane,000)
          write (6,1008) 'at the origin',atnames(iacent)(1:latnam),
     -        resnames(ixres(iacent))(1:lresnam),iresno(iacent)
          write (6,1008) 'on the X axis',atnames(iaxax)(1:latnam),
     -        resnames(ixres(iaxax))(1:lresnam),iresno(iaxax)
          write (6,1008) 'in the X-Y plane',
     -      atnames(iaxyplane)(1:latnam),
     -      resnames(ixres(iaxyplane))(1:lresnam),iresno(iaxyplane)
          call trnsfr(xyzmin,c(1,iacent),3)
          call shiftmol(c,n,xyzmin,c,-1.0)
          call trnsfr(x,c(1,iaxax),3)
          call norm(x,1.0)
          call trnsfr(y,c(1,iaxyplane),3)
          call vprd(x,y,z)
          call vprd(z,x,y)
          call norm(y,1.0)
          call norm(z,1.0)
          do k=1,3
            rot(1,k)=x(k)
            rot(2,k)=y(k)
            rot(3,k)=z(k)
          end do
          write (6,1000) ((rot(i,j),j=1,3),i=1,3)
          call rotate_c(c,n,rot,c,'TRANSROTb',9)
        else if (ans .eq. 'l') then
          deffac=1.54/0.83
          print *,'Default factor converts ChemAxon SDF file units to ',
     -      'Angstrom'
          call getreal('Factor to multiply the coordinates with',39,
     -      deffac,scfac,1,000)
          do ia=1,n
            do k=1,3
              c(k,ia)=scfac*c(k,ia)
            end do
          end do
        else if (ans .eq. 'e') then
c         Separate molecular complexes
          if (n .gt. nslt) then
            print *,'System includes solvents that will be deleted'
            call askstop(1)
            n=nslt
          end if
          im1=1
          im2=2
          if (nmolslt .lt. 2) then
            print *,'Solute is a single molecule - nothing to separate'
          else
            if (nmolslt .gt. 2) then
              call getint('Index of solute molecule to move',32,1,1,
     -          nmolslt,im1,000)
              call getint('Index of solute molecule to move away from',
     -          42,2,1,nmolslt,im2,000)
            end if
            call getreal('Distance to move',16,5.0,dmove,0,000)
            call setrepats(nmolslt,molsltlim,nslt,nneig,ineig,
     -        molresflag,minresflag,it1,it2,maxneig,maxrsd)
            call zeroit(cent,3)
            do im=1,nmolslt
              if (im .eq. im1 .or. im .eq. im2) then
                if (molsltlim(3,im) .gt. 0) then
                  call trnsfr(ct(1,im),c(1,molsltlim(3,im)),3)
                else if (molsltlim(3,im) .eq. 0) then
                  call extension(c,it1,0,molsltlim(1,im),
     -            molsltlim(2,im),x,y,ct(1,im),0,0,v)
                else
                  call zeroit(c00,3)
                  do ia=molsltlim(1,im),molsltlim(2,im)
                    do k=1,3
                      c00(k)=c00(k)+c(k,ia)
                    end do
                  end do
                  do k=1,3
                    ct(k,im)=c00(k)/(molsltlim(2,im)-molsltlim(1,im)+1)
                  end do
                end if
                do k=1,3
                  cent(k)=cent(k)+ct(k,im)
                end do
              end if
            end do
            do k=1,3
              cent(k)=cent(k)/2.0
            end do
            write (6,1003) cent
            r2=0.0
            do k=1,3
              x(k)=ct(k,im1)-cent(k)
              r2=r2+x(k)**2
            end do
            write (6,1004) im1,x
            r2=sqrt(r2)
            do k=1,3
              y(k)=cent(k)+x(k)*(r2+dmove)/r2
            end do
            write (6,1005) im1,y
            do ia=molsltlim(1,im1),molsltlim(2,im1)
              do k=1,3
                c(k,ia)=c(k,ia)+(y(k)-ct(k,im1))
              end do
            end do
          end if
        else if (ans .eq. 'w') then
c         Swap chirality
          call getint('Index of the chiral center atom',31,999999,1,
     -      nslt,iachir,000)
          if (iatnum(iachir) .ne. 6)
     -      print *,'NOTE atom selected is not a carbon'
          write (6,1008) 'selected',atnames(iachir)(1:latnam),
     -        resnames(ixres(iachir))(1:lresnam),iresno(iachir)
          write (6,1006) (i,ineig(i,iachir),
     -      atnames(ineig(i,iachir))(1:latnam),
     -      resnames(ixres(ineig(i,iachir)))(1:lresnam),
     -      iresno(ineig(i,iachir)),i=1,nneig(iachir))
          if (nneig(iachir) .ne. 4) print *,'Number of neighbors is ',
     -      nneig(iachir),' instead of 4)'
          call getint('Index of the 1st neighbor to swap (1/2/3/4/)',
     -      44,999999,1,nneig(iachir),ixchir1,000)
          ixchir2=ixchir1
          do while (ixchir2 .eq. ixchir1)
            call getint('Index of the 2nd neighbor to swap (1/2/3/4/)',
     -        44,999999,1,nneig(iachir),ixchir2,000)
          end do
          iachir1=ineig(ixchir1,iachir)
          iachir2=ineig(ixchir2,iachir)
          ixchir3=1
          do while (ixchir3 .eq. ixchir1 .or. ixchir3 .eq. ixchir2)
            ixchir3=ixchir3+1
          end do
          iachir3=ineig(ixchir3,iachir)
          rchir1=sqrt(dist2(c(1,iachir),c(1,iachir1)))
          rchir2=sqrt(dist2(c(1,iachir),c(1,iachir2)))
          call trnsfi(it1,nneig,nslt)
          call findbranch(iachir,iachir1,iachir3,c,nbranch1,ia_branch,
     -      ia_predef,r_predef,ba_predef,ta_predef,nneig,it1,ineig,it2,
     -      it3,loopmem,nslt,atnames,latnam,resnames,lresnam,ixres,
     -      iresno,ifail1,maxneig,maxrsd,mxat,maxbranch)
          write (6,1007) iachir1,nbranch1
          call findbranch(iachir,iachir2,iachir3,c,nbranch2,
     -      ia_branch(nbranch1+1),ia_predef(1,nbranch1+1),
     -      r_predef(nbranch1+1),ba_predef(nbranch1+1),
     -      ta_predef(nbranch1+1),nneig,it1,ineig,it2,it3,loopmem,nslt,
     -      atnames,latnam,resnames,lresnam,ixres,iresno,ifail2,
     -      maxneig,maxrsd,mxat,maxbranch)
          write (6,1007) iachir2,nbranch2
          if (ifail1+ifail2 .eq. 0) then
c           Swap the coordinates of iachir1 and iachir2
            write (6,1009) iachir1,
     -        atnames(ineig(ixchir1,iachir))(1:latnam),
     -        resnames(ixres(ineig(ixchir1,iachir)))(1:lresnam),
     -        iresno(ineig(ixchir1,iachir)),iachir2,
     -        atnames(ineig(ixchir2,iachir))(1:latnam),
     -        resnames(ixres(ineig(ixchir2,iachir)))(1:lresnam),
     -        iresno(ineig(ixchir2,iachir))
            call unitvec(c(1,iachir),c(1,iachir1),x)
            call unitvec(c(1,iachir),c(1,iachir2),y)
            do k=1,3
              c(k,iachir1)=c(k,iachir)+y(k)*rchir1
              c(k,iachir2)=c(k,iachir)+x(k)*rchir2
            end do
            if (nbranch1 .gt. 0) write (6,1010) iachir1,(i,ia_branch(i),
     -        atnames(ia_branch(i))(1:latnam),
     -        resnames(ixres(ia_branch(i)))(1:lresnam),
     -        iresno(ia_branch(i)),i=1,nbranch1)
            call buildbranch(c,nbranch1,ia_branch,
     -        ia_predef,r_predef,ba_predef,ta_predef,pi,mxat,maxbranch)
            if (nbranch2 .gt. 0) write (6,1010) iachir1,
     -        (i,ia_branch(nbranch1+i),
     -        atnames(ia_branch(i))(1:latnam),
     -        resnames(ixres(ia_branch(i)))(1:lresnam),
     -        iresno(ia_branch(i)),i=1,nbranch2)
            call buildbranch(c,nbranch2,
     -        ia_branch(nbranch1+1),ia_predef(1,nbranch1+1),
     -        r_predef(nbranch1+1),ba_predef(nbranch1+1),
     -        ta_predef(nbranch1+1),pi,mxat,maxbranch-nbranch1)
          end if
          print *,'NOTE: new conformation may involve clashes!'
        end if
      end do
      return
1000  format(' System was rotated with the rotation matrix:',/,
     -  (10x,3f11.7))
1001  format(' System was rotated arount the ',a1,'-axis by',f8.1,
     -  ' deg')
1002  format(' System was translated by <',2(f9.3,','),f9.3,'> A')
1003  format(' Center of centers;',3f10.5)
1004  format(' im=',i4,' center vector=',3f10.5)
1005  format(' im=',i4,' shifted center vector=',3f10.5)
1006  format(' Neighbor atoms of the chiral center:',/,
     -  (i3,' index=',i6,1x,a,1x,a,i6))
1007  format(' Branch on atom ',i6,' consists of ',i4,
     -  ' additional atoms')
1008  format(' Atom ',a,': ',a,1x,a,i5)
1009  format(' First swapping atoms',i6,1x,a,1x,a,i5,' and',
     -  i6,1x,a,1x,a,i5)
1010  format(' Atoms moved in the branch on atom',i6,':',/,
     -  (i3,1x,i6,1x,a,1x,a,i5))
      end
