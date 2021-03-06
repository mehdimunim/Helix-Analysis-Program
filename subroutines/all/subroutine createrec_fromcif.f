      subroutine createrec_fromcif(linew,iocif,iotyp,ctemp,iatnum,
     -  atnam,resnam,segnam1,ires,frocc,charge,iat,blankline)
      dimension ctemp(3)
      character*132 linew,blankline
      character*1 segnam1
      character*4 atnam,segnam,chemnam
      character*6 potnam
      character*8 resnam
      character*2 iatnm2
      common /atnam/ aw(99),nval(99),nvalmax(99),mmofan(99),
     -  mmatno(64),iatnm2(99)
      dimension ibnd(1),ineig(1)
      data chemnam /'    '/,potnam /'      '/,segnam /'    '/
c     print *,'CREATEREC_FROMCIF iat=',iat,' IRES=',ires,' IAN=',iatnum
      chemnam(1:2)=iatnm2(iatnum)
      segnam(1:1)=segnam1
      ihetat=0
      call createrec(linew,iocif,iotyp,ctemp(1),ctemp(2),ctemp(3),' ',
     -  ' ',atnam,resnam,segnam,iat,ires,ires,chemnam,potnam,frocc,
     -  charge,5,1,'  ',1,ineig,1,ibnd,ihetat,blankline)
c     subroutine createrec(line,inpcrdtyp,ioutyp,cx,cy,cz,
c    -  atnam,resnam,segnam,iseqno,iresnum,iresid,chemnam,potnam,
c    -  frocc,q,nqdec,iqspace,mmcgm,nn,in,nconfig,ibnd,ihetat,blankline)
c     character* 132 line,blankline
c     character*1 resnam1
c     character*4 atnam,segnam,chemnam
c     character*6 potnam
c     character*8 resnam
c     character*2 mmcgm
      return
      end
