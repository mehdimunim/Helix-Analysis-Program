      subroutine molreset(nmolfst,nmolslt,nmolshift,c,ct,molsltlim,it1,
     -  cell,ncell,cellalt,icellalt,imcenter,maxrsd,mxat)
      dimension c(3,mxat),ct(3,mxat),molsltlim(3,maxrsd),it1(mxat),
     -  cell(3,ncell),cellalt(3,ncell)
      dimension c0(3),cmin(3),cmax(3),cmincent(3),cmaxcent(3),
     -  cmin12sav(3),cmax12sav(3),cmin12(3),cmax12(3)
c     print *,'MOLRESET ncell=',ncell
      call extension(c,it1,0,molsltlim(1,imcenter),
     -  molsltlim(2,imcenter),cmincent,cmaxcent,c0,0,0,v)
      nmolshift=0
c     imgw=0
      do im=nmolfst,nmolslt
        ialtopt=0
        if (im .ne. imcenter) then
          ichange=1
          do while (ichange .gt. 0)
            call extension(c,it1,0,molsltlim(1,im),molsltlim(2,im),
     -        cmin,cmax,c0,0,0,v)
            call combinevol(cmin,cmax,cmincent,cmaxcent,cmin12sav,
     -        cmax12sav,volmin)
c             write (6,8733) 1,cmin,cmax,cmincent,cmaxcent,cmin12sav,
c    -          cmax12sav,volmin
            imgopt=1
c           imgw=imgw+1
c           call savepdb(89,'MOLRESET.pdb',12,c,molsltlim(2,nmolslt),
c    -        imgw)
            do img=2,ncell
              call trnsfr(ct,c,3*molsltlim(2,nmolslt))
              call shiftmol(c(1,molsltlim(1,im)),
     -          molsltlim(2,im)-molsltlim(1,im)+1,cell(1,img),
     -          ct(1,molsltlim(1,im)),-1.0)
c             imgw=imgw+1
c             call savepdb(89,'MOLRESET.pdb',12,ct,molsltlim(2,nmolslt),
c    -          imgw)
              call extension(ct,it1,0,molsltlim(1,im),molsltlim(2,im),
     -          cmin,cmax,c0,0,0,v)
              call combinevol(cmin,cmax,cmincent,cmaxcent,cmin12,cmax12,
     -          vol)
c              write (6,8733) img,cmin,cmax,cmincent,cmaxcent,cmin12,
c     -          cmax12,vol
c8733          format(' IMG=',i2,' CMIN=',3f8.2,' CMAX=',3f8.2,/,
c     -               ' CMINCENT=',3f8.2,' CMAXCENT=',3f8.2,/,
c     -               ' CMIN12=',3f8.2,' CMAX12=',3f8.2,' VOL=',f12.2)
              if (vol .lt. volmin) then
                volmin=vol
                call trnsfr(cmin12sav,cmin12,3)
                call trnsfr(cmax12sav,cmax12,3)
                imgopt=img
              end if
c             write (6,9888) img,vol,imgopt,volmin
c9888         format(' IMG=',i2,' V=',f12.3,' IMGOPT=',i2,
c    -        ' VOLMN=',f12.3)
            end do
            do k=1,3
              cmincent(k)=amin1(cmincent(k),cmin12sav(k))
              cmaxcent(k)=amax1(cmaxcent(k),cmax12sav(k))
            end do
            if (icellalt .gt. 0) then
c             Try the alternate (TO) cell orientation
              do img=2,ncell
                call trnsfr(ct,c,3*molsltlim(2,nmolslt))
                call shiftmol(c(1,molsltlim(1,im)),
     -            molsltlim(2,im)-molsltlim(1,im)+1,cell(1,img),
     -            ct(1,molsltlim(1,im)),-1.0)
                call extension(ct,it1,0,molsltlim(1,im),molsltlim(2,im),
     -            cmin,cmax,c0,0,0,v)
                call combinevol(cmin,cmax,cmincent,cmaxcent,cmin12,
     -            cmax12,vol)
                if (vol .lt. volmin) then
                  volmin=vol
                  call trnsfr(cmin12sav,cmin12,3)
                  call trnsfr(cmax12sav,cmax12,3)
                  imgopt=img
                  ialtopt=1
                end if
              end do
            end if
c           print *,'IMGOPT=',imgopt
            if (imgopt .gt. 1) then
              ichange=1
              if (ialtopt .eq. 0) then
                call shiftmol(c(1,molsltlim(1,im)),
     -            molsltlim(2,im)-molsltlim(1,im)+1,cell(1,imgopt),
     -            c(1,molsltlim(1,im)),-1.0)
              else
                call shiftmol(c(1,molsltlim(1,im)),
     -            molsltlim(2,im)-molsltlim(1,im)+1,cellalt(1,imgopt),
     -            c(1,molsltlim(1,im)),-1.0)
              end if
            else
              ichange=0
            end if
          end do
        end if
      end do
c     close (89)
      return
      end
