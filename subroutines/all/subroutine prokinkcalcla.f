      subroutine prokinkcalcla(nrepc,c,nslt,bend,wobble,faceshift,
     -  rmsb,rmsa,
     -  iflatproline,ix5,radtodeg,maxhlx)
      dimension c(3,nslt),ix5(5)
      real*8 calph(3,maxhlx),axisdirb(3),axisinib(3),axisendb(3),
     -  axisdira(3),axisinia(3),axisenda(3),rms
      parameter (MAXHX=50)
      common /prokink/ icab(MAXHX),icaa(MAXHX),icb(MAXHX),ica(MAXHX),
     -  inb(MAXHX),ina(MAXHX),icapr,icpr,inpr,nra,nrb,icbpr,icgpr,icdpr,
     -  iprintpk
      character*60 message
c     ProKink variables
      real*8 ddot,dmag,cosa,pro_alphC,orig,origa,orig3,orig4,
     -  perpvec,perpveca,perpvec3,perpvec4,perpvec34,enda,
     -  bmin3C,bmin4C,zax,xx,prolinering,p0,ringnorm,axfact,rr
      dimension pro_alphC(3),orig(3),origa(3),orig3(3),orig4(3),
     -  perpvec(3),perpveca(3),perpvec3(3),perpvec4(3),perpvec34(3),
     -  bmin3C(3),bmin4C(3),enda(3),zax(3),xx(3),prolinering(3,5),
     -  p0(3),ringnorm(3),axfact(5)
c     print *,'PROKINKCALCLA icpr,inpr,nra,nrb=',icpr,inpr,nra,nrb
      iprintkahn=0
      do k=1,3
        calph(k,1)=c(k,icapr)
      end do
      do ir=1,nra
        do k=1,3
          calph(k,ir+1)=c(k,icaa(ir+1))
        end do
      end do
      message=
     - 'Helix after the kink                                        '
      call kahn(calph,nra+1,.true.,axisdira,axisinia,axisenda,
     -  rms,iprintkahn, message,MAXHX)
      rmsa=rms
      if (iprintpk .gt. 1) then
        write (6,7011) 'after','axis',axisdira
        write (6,7011) 'after','init',axisinia
        write (6,7011) 'after','end ',axisenda
        print *,'rms=',rms
      end if
      do ir=1,nrb
        do k=1,3
c         calph(k,ir+1)=c(k,icab(ir+1))
          calph(k,ir)=c(k,icab(ir))
        end do
      end do
      message=
     -  'Helix before the kink                                       '
      call kahn(calph,nrb,.true.,axisdirb,axisinib,axisendb,rms,
     -  iprintkahn, message,MAXHX)
      rmsb=rms
      if (iprintpk .gt. 1) then
        write (6,7011) 'before','axis',axisdirb
        write (6,7011) 'before','init',axisinib
        write (6,7011) 'before','end ',axisendb
        print *,'RMS=',rms
      end if
      iac=icapr
      im3=icab(3)
      im4=icab(4)
      do k=1,3
        pro_alphC(k)=c(k,icapr)
        bmin3C(k)=c(k,im3)
        bmin4C(k)=c(k,im4)
      end do
      if (iflatproline .eq. 1) then
        do i=1,5
          do k=1,3
            prolinering(k,i)=c(k,ix5(i))
          end do
        end do
        call fitpoints(prolinering,5,3,p0,ringnorm,axfact,iprintpk)
        call dvdif(pro_alphC,p0,xx)
        rr=ddot(xx,ringnorm)
        do k=1,3
          pro_alphC(k)=pro_alphC(k)+rr*ringnorm(k)
        end do
      end if
c     New code for proline kink calculation
c     Bend angle: from the scalar product of the two axis vectors
c     Both axis vectors point away from the proline
      cosa=-ddot(axisdira,axisdirb)
      bend=dacoscheck(cosa,ccc,1,6,'PROK-B')*radtodeg
c     Wobble  angle: from the scalar product of the two normals to the
c     before helix axis (from the C-alpha of Proline and the end of the
c     after helix)
      call calcperp(axisendb,axisdirb,pro_alphC,orig,perpvec,iprintpk)
      call dvsum(orig,axisdira,enda)
      call calcperp(axisinib,axisdirb,enda,origa,perpveca,iprintpk)
      call dvdif(origa,enda,xx)
      cc=dmag(xx)
      if (cc .lt. 0.01) print *,'Wobble is suspect'
      cosa=ddot(perpvec,perpveca)
      wobble=dacoscheck(cosa,ccc,1,6,'PROK-W')*radtodeg
c     Establish sign
      call dcross(axisdirb,perpvec,zax)
      if (ddot(zax,perpveca) .gt. 0.d0) wobble=-wobble
      call calcperp(axisinib,axisdirb,bmin3C,orig3,perpvec3,iprintpk)
      call calcperp(axisinib,axisdirb,bmin4C,orig4,perpvec4,iprintpk)
      call dvnorm(perpvec)
      call dvnorm(perpvec3)
      call dvnorm(perpvec4)
      do k=1,3
        perpvec34(k)=(perpvec3(k)+perpvec4(k))/2.d0
      end do
      call dvnorm(perpvec34)
      cosa=ddot(perpvec,perpvec34)
      faceshift=dacoscheck(cosa,ccc,1,6,'PROK-FS')*radtodeg
c     Establish sign
      cosa3=ddot(perpvec,perpvec3)
      cosa4=ddot(perpvec,perpvec4)
      if (cosa3 .gt. cosa4) faceshift=-faceshift
      if (nrepc .lt. 0) return
      return
7011  format(' Helix ',a,': ',a,'=',3f10.4)
      end
