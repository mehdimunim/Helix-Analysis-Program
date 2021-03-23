      def analyze(nconfig,inpcrdtyp,inpcrdtyporg,ianaltyp,
     -  extnam1,analfile,inpfile,namleni,line,index,n,nslt,nmolslt,
     -  numres,nresslt,c,atw,rprox,cv,charge,iatnum,ifchrg,ncl,nneig,
     -  nneiga,nhbneig,ineig,nhneig,nnneig,ncneig,nsneig,npneig,mmtype,
     -  molsltlim,isegno,altcol,inscol,ninsres,ih,ifgtyp,molresflag,
     -  cres,chn,c2,cdp,cdp2,temp,ibnd,indexn,indexo,indexs,indexa,
     -  indexdel,ifree,icntrl,indexov,indexrmsd,itemp1,itemp2,itemp3,
     -  itemp4,irrix,naslv,islvw,innlist,hblimfac,angmin,iresno,ixres,
     -  ixresno,ixsegno,ifres,ilres,irepav,irepmx,irepeng,irepkm,
     -  irescol1,irescol2,iresncol1,iresncol2,inamcol1,inamcol2,iqcol1,
     -  iqcol2,idcol,segid4,resnamslv,atnames,resnames,irescount1,
     -  irescount2,irescount3,listrefres,listnegres,segnames,ianc_anc,
     -  dssplab,idistdssp,marker,iwhead,nsegslt,isc,reportqmin,ntitlin,
     -  ntitlinw,title,version,iresnrestart,iresidrestart,blankline,
     -  label2d,inptrajtyp,mmctrajtyp,itraj,ireorient,centgra,iconfsel,
     -  numsel,nclstmem,iresshift,lastanal,ifixq,nhbdist,rhbdist,iwt,
     -  rlim,iqspaceask,icharges,iallheavy,keeprem,innlistread,lentest,
     -  pi,radtodeg,asterisk,mflogfile,maxconfsel,maxrepconf,maxng,
     -  maxbox,mxbonds,maxrsd,maxrec,maxrec10,maxrcorr,mx2d)

    """







    """



      dimension nneig(maxrec),ineig(maxng,maxrec),iatnum(maxrec),
     -  c(3,maxrec10),atw(maxrec),cv(maxrec),charge(maxrec),
     -  nhbneig(maxrec),nneiga(maxrec),nhneig(maxrec),nnneig(maxrec),
     -  ncneig(maxrec),nsneig(maxrec),npneig(maxrec),index(maxrec),
     -  ifgtyp(maxrec),molresflag(maxrsd),isegno(maxrec),ifchrg(maxrec),
     -  mmtype(maxrec),iresno(maxrec),ixres(maxrec),ixresno(maxrsd),
     -  ixsegno(maxrsd),ifres(maxrsd),ilres(maxrsd),irepav(mx2d),
     -  irepmx(mx2d),irepeng(mx2d),irepkm(mx2d),irescount1(maxrsd),
     -  irescount2(maxrsd),irescount3(maxrsd),listrefres(maxrsd),
     -  listnegres(maxrsd),ih(maxrec),cres(3,maxrec10),chn(3,maxrec),
     -  c2(3,maxrec),temp(maxrec),ibnd(maxbox,maxrec),ifree(maxrec),
     -  indexov(maxrec),indexrmsd(maxrec),itemp1(maxrec),itemp2(maxrec),
     -  itemp3(maxrec),itemp4(maxrec),irrix(maxrec),indexn(maxrec),
     -  indexo(maxrec),indexs(maxrec),indexa(maxrec),indexdel(maxrec),
     -  iconfsel(maxconfsel),nclstmem(mx2d),idistdssp(9,maxrsd),
     -  icntrl(20),centgra(3),molsltlim(3,maxrsd),ianc_anc(mxbonds),
     -  nhbdist(mxbonds),rhbdist(mxbonds),iwt(mx2d),rlim(maxng),
     -  rprox(maxrec),isc(maxrec)

     
      #real*8 cdp(3,maxrec),cdp2(3,maxrec)
      character*1 asterisk,dssplab(maxrsd),altcol(maxrec),inscol(maxrec)
      character*4 segid4(nsegslt),segnames(maxrsd),extnam1,extnam2,
     -  extnam3
      character*8 resnam,rn,resnamslv,atnames(maxrec),resnames(maxrsd)
      character*6 marker(16)
      character*8 version
      character*24 askcolcode
      character*80 label2d(mx2d)
      character*200 inpfile,analfile,analfile1,analfile2,analfile3,
     -  analfile4,trackfile
      character* 132 line(maxrec),blankline
      character*2 iatnm2
      common /atnam/ aw(99),nval(99),nvalmax(99),mmofan(99),
     -  mmatno(64),iatnm2(99)
      character*11 trajformatname
      common /trajectory/ nmmccheck,iftrajtyp(6),trajformatname(6)
      character*200 trajnam,trajnam2
      common /trajname/ trajnam,trajnam2,ltrajnam,ltrajnam2,ifirsttraj,
     -  ifirsttraj2,ilasttraj,ilasttraj2,incrementtraj,incrementtraj2
      character*80 title,trtitle(32)
      character*4 namfcg
      character*4 tanames
      character*8 tnames
      common /tordat/ ntorn,tanames(4,28),tnames(28)
      common /connatdat/ ramax(99),ramax2(99),hlimfac,ianfg(99),
     -  namfcg(100),nrmw
      #real*8 xtlabc,xtlabc0
      common /boxdat/ xtlabc(6),xtlabc0(6),box(3),box0(3),edge_gen(3,3),
     -  cell0(3,27),cell(3,27),cellalt(3,27),
     -  ncell,ioppbc,noboxinfoar,noboxinfow,noboxrep,
     -  istuner,iboxtypfound,ixcrd(3),ixang(3),ixyzhex(3),
     -  ixyzhextraj(3),isizewarn
      common /colorinfo/ ncolcode,maxcolcode
      common /graphics/ npixx,npixy,maxpixx,maxpixy,idwmain,idwplot,
     -  wx,wy,wz,wxdr
      common /rotmat/ matrot0(4,4),matrot(4,4),nomat0
      character*1 xyz
      common /axislab/ xyz(3)
      common /columnlim/ incol(19),iidcol(19),iialtcol(19),iiinscol(19),
     -  iinamcol(2,19),iirescol(2,19),iiccol(2,19),iiresncol(2,19),
     -  iiseqncol(2,19),iisegcol(2,19),iiresidcol(2,19),iiqcol(2,19),
     -  iipotcol(2,19),iiocccol(2,19),iichemcol(2,19)
      character*5 crdext
      common /iotypes/ iocha,iochaex,iobpdb,iocpdb,ioa3pdb,ioa4pdb,
     -  iommod,iommc,iommc4,iogro,iomol2,iomae,iocif,ioxxx,ioins,
     -  ionxyz,iosxyz,iosxyzrq,iograsp,iofull,lext(19),crdext(19)
      dimension edge(3),evals0(3),evecs0(3,3)
      character*4 pflsv(100)
      character*8 namesv(100)
      character*9 decideb# end
      character*27 qq
      parameter (MAXBRIDGEATOM=8000,MAXBRIDGETYPE=50,MAXBRIDGELEN=4)
c     Maximum number of bridge anchor atoms: MAXBRDIGEATOM
c     Maximum number of bridge destination atoms: MAXBRDIGETYPE
c     Maximum number of H bonds in a bridge:  MAXBRIDGELEN
      common /bridges/ ianchor(MAXBRIDGEATOM),
     -  nbridgetype(MAXBRIDGELEN,MAXBRIDGEATOM),
     -  ibridgetype(MAXBRIDGETYPE,MAXBRIDGELEN,MAXBRIDGEATOM),
     -  lpath(MAXBRIDGETYPE,MAXBRIDGELEN,MAXBRIDGEATOM)
      parameter (MAXDDISTR=200,MAXCDLIST=2000,MAXDDBIN=20)
c     Maximum number of atom pairs to calculate the distance distribution
      dimension listpairdist(2,MAXDDISTR),npairdist(MAXDDBIN,MAXDDISTR),
     -  pairdistsum(MAXDDISTR),pairdistsum2(MAXDDISTR),
     -  pairdistwsum(2,MAXDDISTR),pairdistminmax(2,MAXDDISTR),
     -  iclustermem(MAXCDLIST),ifstclst1(MAXDDISTR),
     -  ifstclst2(MAXDDISTR),ilstclst2(MAXDDISTR)
      character*1 typc
      character*21 ssname
      common /dsspnames/ lssname(9),ssname(9),typc(9)
c     Pseudorotation calculation
      parameter (MAXRING=50)
      #real*8 sinpsrs,cospsrs,qpsrs,qpsr2s,zavs,zsqs
      dimension sinpsrs(MAXRING),cospsrs(MAXRING),qpsr(MAXRING),
     -  qpsrs(MAXRING),qpsr2s(MAXRING)
      dimension psr(MAXRING),ix5(MAXRING),zring(MAXRING),rring(MAXRING)
      dimension iasv(100),qsv(100),ixlist(99),idelse:g(1000),ialist(15),
     -  icatlist(15),psr5(5),itypsse(200),ifsse(200),ilsse(200)
c     All arrays in prokink are of length MAXHX
      parameter (MAXHX=50)
      common /prokink/ icab(MAXHX),icaa(MAXHX),icb(MAXHX),ica(MAXHX),
     -  inb(MAXHX),ina(MAXHX),icapr,icpr,inpr,nra,nrb,icbpr,icgpr,icdpr,
     -  iprintpk
      common /analparm/ nsltref_f,nsltref_l,rcut_cv,icvtyp
c     All arrays are of length maxframe=MAXFRAMES > MAX2D ##
      parameter (MAXFRAMES=50000,MAXCOPY=600)
      parameter (MAXCOPY1=MAXCOPY-1)
      common /analres/ nframe,nframeref,nframetot,maxframe,maxpres,
     -  increment,increment2,res(2,MAXFRAMES,MAXCOPY1),
     -  xyplot(2,MAXFRAMES),x0(MAXCOPY),y0(MAXCOPY),nxselres,
     -  ixselres(MAXCOPY)
      common /logging/ logfile,ipredict
      character*1 separatorchar
      common /filenuminfo/ iaskunderscore,separatorchar
      parameter (MAXNHX=12,MAXNHX2=(MAXNHX*(MAXNHX-1)))
      character*2 ap_pa,in_ex
      common /signs/ itmem,normhx,isg2(MAXNHX2),memdir(MAXNHX),ap_pa(3),
     -  in_ex(2)
c     common /hxtrack/ ihxt1,ihxt2,torsav1(3,4),torsav2(3,4),x1(3),x2(3)
c     Array of lenght maxconfsel (set to MAXFRAMES)
      dimension iconfsel2(MAXFRAMES),rplot(MAXFRAMES)
c     Helix axis analysis arrays
c     Arrays of length involving MAXNHX
      dimension anglechange(MAXHX,MAXNHX),anglechangeref(MAXHX,MAXNHX),
     -  angles(3,MAXNHX),anglesn(3,MAXNHX),helixlen(MAXNHX),
     -  helixlen0(MAXNHX),rcirchx(MAXNHX),turnperres(MAXNHX),
     -  ireshx1(MAXNHX),ireshx2(MAXNHX),iseghx(MAXNHX),nreshx(MAXNHX),
     -  icaahx(MAXHX,MAXNHX),icbahx(MAXNHX),icbreshx(MAXNHX),
     -  indexaxhx(3,MAXNHX),lhxhxlab(MAXNHX2)
      dimension plotdat(2,MAXFRAMES)
      character*24 hxhxlab(MAXNHX2)
      #real*8 calph0(3,MAXHX,MAXNHX),perpvec0(3,MAXHX,MAXNHX),
     -  camod(3,MAXHX,MAXNHX),axfact(MAXHX,MAXNHX),
     -  calph(3,MAXHX,MAXNHX),perpvec(3,MAXHX,MAXNHX),
     -  axisdir(3,MAXNHX),axisini(3,MAXNHX),axis# end(3,MAXNHX),
     -  helixcent(3,MAXNHX),circ(3,MAXNHX),axisdir0(3,MAXNHX),
     -  axisini0(3,MAXNHX),axisen0(3,MAXNHX),helixcent0(3,MAXNHX),
     -  rn0(3,MAXNHX),rnorm(3,MAXNHX),rms,ddd,sinphisum,cosphisum,
     -  sinpsisum,cospsisum,atwsum
      character*1 axdirchar(MAXNHX)
      dimension crmslt0(3),crmslt(3),rot(3,3),trajrot(3,3),nhelixok(3),
     -  dr(3),caref(3),caprev1(3),caprev2(3),
     -  drincr(3),drimg(3),com1(3),com2(3),indexax(3),
     -  dip1(3),sdxyx1(3),lbondname(6)
      parameter (MAXDISTR=1000,MAXDISTRN=10)
      dimension nframes_err(MAXDISTRN),err12(2,MAXDISTR),
     -  rmsfav(MAXDISTR),rmsfavs(2,MAXDISTR),bl(MAXDISTRN)
      #real*8 distrerr(MAXDISTRN,MAXDISTR),rmsfsum(MAXDISTR),
     -  xcum(MAXDISTRN)
      character*1 analtyp,analinp,psrtyp,ansrun,marks(9),unit,ans
      character*1 tmchar,tmchars(2)
      character*4 ext1(45),ext2(45),ext3(45),atnam,bondlab(5)
      character*8 brslv
      character*6 hxoklab(3),crdexti
      character*7 mark0,mark1,mark2
      character*11 xtrajlabs(4),xtrajlab
      character*20 helixcklab
      character*22 bondname(6)
      character*30 shiftlab
      character*80 question,linein,pstitle,system
      character*34 distlab,distlab_cc
      character*36 restitle
      character*22 atomdist
      character*25 prokinklab(5),helixang(6),helixrlab(14),rmsdlab(4),
     -  ramalab(MAXCOPY1),printrlab(10),volumelab(4)
      character*29 resrange
      character*30 talab(MAXCOPY1)
      character*80 plotdescr
      character*100 hostname
      character*200 trajnam1,trajnamr1,trajnamr2
      dimension lprokinklab(5),lhelixang(6),lhelixrlab(14),lrmsdlab(4),
     -  lramalab(MAXCOPY1),lprintrlab(10),lvolumelab(4),
     -  ltalab(MAXCOPY1),ixtor1234(4,MAXCOPY1)
      parameter (MAXPHI=400,MAX2D=5000,MAXBONDS=10000)
      dimension it1(MAXBONDS),it2(MAXBONDS),it3(MAXBONDS),it4(MAXBONDS),
     -  it5(MAXBONDS),ixclst(MAXBONDS),index2d(MAXBONDS),
     -  ifhb2d(MAXBONDS),ilhb2d(MAXBONDS),ifa_s(MAXBONDS),
     -  ila_s(MAXBONDS),rmsdlim(MAXBONDS),value(MAXBONDS),engcl(MAX2D),
     -  ifclst1(MAX2D),ilclst1(MAX2D),index2d1(MAX2D),nsimclst1(MAX2D),
     -  ifclst2(MAX2D),ilclst2(MAX2D),index2d2(MAX2D),nsimclst2(MAX2D),
     -  irepmx1(MAX2D),irepmx2(MAX2D),ixshuffle(MAX2D),
     -  ixshuffleref(MAX2D),xtraj(MAXFRAMES)
      data marks /'*','+','#','|','@','%','x','-','='/
      data bondname /'hydrogen bond','hydrophobic bond','salt bridge',
     -  'heavy atom VdW contact','heavy atom MPX contact',
     -  'hydrogen-bond bridge'/
      data lbondname /13,16,11,22,22,20/
      data bondlab /'D..A','C..C','+..-',2*'A..A'/
      data shiftlab /'Helix shift in the Z direction'/
      data prokinklab /'B# end angle','Wobble','Face shift',
     -  'Wobble-Face shift','Pseudorotation phase'/
      data lprokinklab /10,6,10,22,20/
      data helixang /'Global helix X tilt-angle',
     -  'Global helix Y-tilt angle','Global helix Z-tilt angle',
     -  'Helix rotation','Local helix tilt angle','Turn angle per res.'/
      data lhelixang /25,25,25,14,16,19/
      data helixrlab /'Helix length','Curvature radius',
     -  'Total displacement (cent)','X-displacement (cent)',
     -  'Y-displacement (cent)','Z-displacement (cent)',
     -  'Normal-to-b# end drift (X)','Normal-to-b# end drift (Y)',
     -  'X-displacement (start)','Y-displacement (start)',
     -  'Z-displacement (start)','X-displacement (# end)',
     -  'Y-displacement (# end)','Z-displacement (# end)'/
      data lhelixrlab /12,16,25,3*21,2*24,3*22,3*20/
      data rmsdlab /'RMSD (overlaid)   ','Max dev (overlaid)',
     -              'RMSD              ','Maximum deviation '/
      data lrmsdlab /17,19,4,19/,ltalab /MAXCOPY1*30/,
     -  ramalab/MAXCOPY1*'Residue                  '/,
     -  lramalab /MAXCOPY1*25/
      data volumelab /'V(solvent-excluded shell)',
     -  'V(first solvation shell)','V(solute)','V(interface)'/
      data lvolumelab /25,24,9,12/
      data ext1 /'.nnl','.n14','.fcg','.bnd','.hbn','.hpn','.stb',
     -    '.rsd','.dst','.cnt','.psr','.pkn','.hph','.cvl','.cvr',
     -    '.dss','.hbr','.ram','.dih','.phi','.axd','.rms','.rd2',
     -    '.rdx','.cor','.add','.svl','.pca','.rgh','.sum','.adj',
     -    '.ang','.mld','.ctb',4*'    ','.sdr','.flt','.mpx',4*'    '/
      data ext2 /'    ','    ','.bbn','    ','.ps ','.ps ','.ps ',
     -    '.rsm','.ps ','    ','    ','.ps ','    ','    ','.ps ',
     -    '.ps ','    ','.ps ','.ps ','    ','.ps ','.ps ','.ps ',
     -    '.ps ','.ps ','    ','.ps ','.pdb','.ps ','    ','.ps ',
     -    '.ps ','.mlc','.ps ',4*'    ','.ps ','    ','.ps ',4*'    '/
      data ext3 /7*'    ','.rmp',9*'    ','.rdp',3*'    ','.rmf',
     -  10*'    ','    ',12*'    '/
      data crdexti /'      '/,namlen_root /0/
      data hxoklab /'intact','broken','frayed'/
      data xtrajlabs /' N(frames) ','Picoseconds','Nanoseconds',
     -  'Miliseconds'/,nframes_err /MAXDISTRN*0/
      data lhxhxlab /MAXNHX2*24/
      data tmchars /'i','e'/
      data ifirst /1/, ilast /999999/,ifirst2 /1/, ilast2 /999999/,
     -  nframetotmin /50/,maxdevplot /0/,iunmatchplot /0/,
     -  nhneigmin /0/,numsel2 /0/,noprol /0/
c     ianaltyp=1:  (S) Neighbor, bond, angle and torsion list
c     ianaltyp=2:  (S) 1-4 statistics
c     ianaltyp=3:  (S) defal group and backbone list
c     ianaltyp=4:  (S) Bond length statistics
c     ianaltyp=5:  (S) Hydrogen-bond list
c     ianaltyp=6:  (S) Hydrophobic bond list
c     ianaltyp=7:  (S) Salt bridge list
c     ianaltyp=8:  (S) Calculate residue distances
c     ianaltyp=9:  (S) Calculate a PBC-adjusted distance
c     ianaltyp=10: (S) Check for potentially unphysical contacts
c     ianaltyp=11: (S) Pseudorotation angle calculation
c     ianaltyp=12: (S) Calculate Proline kinks
c     ianaltyp=13: (S) Hydropathy labeling
c     ianaltyp=14: (S) Circular variance labeling of solute and solvent
c     ianaltyp=15: (S) Circular variance residue-residue plot
c     ianaltyp=16: (S) DSSP secondary structure assignment
c     ianaltyp=17: (S) Hydrogen-bond bridge analysis
c     ianaltyp=18: (S) Ramachandran plot
c     ianaltyp=19: (S) Torsion dial plots
c     ianaltyp=20: (S) Delphi map annotation
c     ianaltyp=21: (S) Helix axis directions
c     ianaltyp=22: (T) 1-D RMSD and residue RMS fluctuations
c     ianaltyp=23: (T) 2-D RMSD
c     ianaltyp=24: (T) Cross RMSD
c     ianaltyp=25: (T) Residue correlation matrix calculation
c     ianaltyp=26: (S) Atom-atom distance (distribution) calculation
c     ianaltyp=27: (S) Solvation shell volume calculation
c     ianaltyp=28: (S) Principal axis calculation
c     ianaltyp=29: (S) Radius and dipole calculation
c     ianaltyp=30: (S) Summarize Amber energy partition table
c     ianaltyp=31: (S) Adjacency matrix analysis
c     ianaltyp=32: (S) Angle dial plots
c     ianaltyp=33: (S) Molecule-molecule distance list
c     ianaltyp=34: (T) Heavy atom contact list
c     ianaltyp=35: (S) Calculate eigenvectors from input matrix
c     ianaltyp=36: (S) Compare residue-residue average distance matrices
c     ianaltyp=37: (S) Compare residue-residue bond matrices
c     ianaltyp=38: (S) Compare residue RMSF values
c     ianaltyp=39: (T) Atom-atom SD from trajectory
c     ianaltyp=40: (T) Solvent filtering
c     ianaltyp=41: (T) Mutually proximal contact list
      call testconst(0,1,2,0.0,1.0,2.0,6,nfail,1,'ANST')
      nmc=0
      maxcolcode=8
      iaskcolcode=0
      increment=1
      increment2=1
c     do ir=1,nresslt
c       write (77,*) 'ir=',ir,' iresf,l=',ifres(ir),ilres(ir)
c       do ia=ifres(ir),ilres(ir)
c         write (77,*) line(index(ia))(1:78)
c       
c     
c     stop
      call indexit(ixshuffle,1,MAX2D,0)
      call indexit(ixshuffleref,1,MAX2D,0)
      for i in range(0, MAXFRAMES):
        xtraj(i)=i
      
      nanos=0
      ncolcode=maxcolcode
      write (askcolcode,2152) maxcolcode
      ireseq=-1
      nframesign=0
      nextconfsel=0
c     print *,'numsel,maxconfsel=',numsel,maxconfsel
c     print *,'maxrepconf,maxng,maxbox=',maxrepconf,maxng,maxbox
c     print *,'maxrec,maxrsd=',maxrec,maxrsd
      iseqncol1=iiseqncol(1,inpcrdtyp)
      iseqncol2=iiseqncol(2,inpcrdtyp)
      is1=iisegcol(1,inpcrdtyp)
      is2=iisegcol(2,inpcrdtyp)
      nrescol=irescol2-irescol1+1
      nnamcol=inamcol2-inamcol1+1
      nsegcol=is2-is1+1
      innlistorig=innlist
      call trnsfr(cres,c,3*n)
      ifc=1
      call nextchar(title,ifc,1000)
      call lastchar(title,ltitle,80)
      call blankout(question,1,80)
      question(1:ltitle-ifc+1)=title(ifc:ltitle)
      title=question
      ltitle=ltitle-ifc+1
      ltitle76=min0(ltitle,76)
      natsorig=n
      nsegm=isegno(nslt)
      numanal=0
      iallrama=0
      iconfirmname=0
      iuseinp=0
      limresrange=0
      nhxres=15
      nfravgd=0
      nfravgt=0
9000  idsspcheck=0
      ltrajnam2=0
      noprintconf=0
      iread2d=0
      etot2=0.0
      itypavg=0
      rmsdplotmax=0.0
      nfreqsd=0
      rmsfdplotmaxlotmax=0.0
      isdtyp=1
      ireadtracks=0
      ilastframe=0
      nbfound=0
      nbresfound=0
      nosameseg=0
      iadjust_xtraj=1
      ionefile=1
      numsolv=(n-nslt)/naslv
      ihbondcalc=0
      ibondtype=0
      inputref=0
      nbondavg=1
      nmcmaxbond=0
      maxbondf=0
      ifstrmsd=1
      nocontigrmsd=0
      ilstrmsd=nresslt
      iconndial=1
      mappdf=0
      if (nconfig  ==  1) :
        call getatnumlist(n,iatnum,ifchrg,ialist,icatlist,ixlist,nanos)
c       Get residue name list
        ires=0
        nres=0
        for ia in range(0, n):
          if (iresno(ia)  !=  ires) :
            ires=iresno(ia)
            nres=nres+1
            resnames(nres)(1:irescol2-irescol1+1)=
     -        line(index(ia))(irescol1:irescol2)
          ## end if
        
        namleno=0
        namleno1=0
        namleno2=0
        extnam2='    '
        hbf0=hblimfac
        angm0=angmin
        listbridge=0
      ## end if
9002  if (nconfig  ==  1) :
        numanal=numanal+1
9003    call quiz(analinp,ityp,' ',' ',0,'configuration analysis',22,
     -    0,5,6,0)
        if (analinp  ==  'q') return
        analtyp=analinp
        if (analinp  ==  'g') :
          call quiz(analtyp,ityp,' ',' ',0,
     -      'topology/geometry analysis',26,0,5,6,0)
          ianaltyp=ityp
        elif (analinp  ==  'b') :
          call quiz(analtyp,ityp,' ',' ',0,'bond tracking',13,
     -      0,5,6,0)
          if (analtyp  ==  'h') ianaltyp=5
          if (analtyp  ==  'p') ianaltyp=6
          if (analtyp  ==  's') ianaltyp=7
          if (analtyp  ==  'b') ianaltyp=17
          if (analtyp  ==  'c') ianaltyp=34
          if (analtyp  ==  'm') ianaltyp=37
          if (analtyp  ==  'u') ianaltyp=41
        elif (analinp  ==  'y') :
          call quiz(analtyp,ityp,' ',' ',0,
     -      'atomic property calculation',27,0,5,6,101)
          if (analtyp  ==  'y') ianaltyp=13
          if (analtyp  ==  'c') ianaltyp=14
          if (analtyp  ==  'p') ianaltyp=20
        elif (analinp  ==  'o') :
          call quiz(analtyp,ityp,' ',' ',0,
     -      'molecular property calculation',30,0,5,6,102)
          ianaltyp=26+ityp
          if (analtyp  ==  'b') :
            if (ispdb(inpcrdtyp)  ==  0) :
              print *,'B factors are only used for PDB input'
            else:
              print *,'B-factor segment averages are printed above'
            ## end if
            go to 9003
          ## end if
        elif (analinp  ==  'm') :
          call quiz(analtyp,ityp,' ',' ',0,'RMSD calculation',16,
     -      0,5,6,0)
          ianaltyp=21+ityp
          if (analtyp  ==  'd') ianaltyp=38
        elif (analinp  ==  'u') :
          call quiz(analtyp,ityp,' ',' ',0,'distance analysis',17,
     -      0,5,6,0)
          if (analtyp  ==  'e') ianaltyp=8
          if (analtyp  ==  'l') ianaltyp=33
          if (analtyp  ==  'm') ianaltyp=31
          if (analtyp  ==  'u') ianaltyp=9
          if (analtyp  ==  't') ianaltyp=26
          if (analtyp  ==  's') ianaltyp=39
          if (analtyp  ==  'c') ianaltyp=10
          if (analtyp  ==  'o') ianaltyp=36
        else:
          if (analtyp  ==  'a') ianaltyp=18
          if (analtyp  ==  'i') ianaltyp=32
          if (analtyp  ==  't') ianaltyp=19
          if (analtyp  ==  'k') ianaltyp=12
          if (analtyp  ==  'x') ianaltyp=21
          if (analtyp  ==  'p') ianaltyp=11
          if (analtyp  ==  's') ianaltyp=16
          if (analtyp  ==  'v') ianaltyp=15
          if (analtyp  ==  'f') ianaltyp=40
          ireadcov=0
          if (analtyp  ==  'n') :
            call askyn('Do you have an input covariance matrix',38,1,-1,
     -        ireadcov,000,0)
            if (ireadcov  ==  1) :
              ianaltyp=35
            else:
              ianaltyp=25
            ## end if
          ## end if
          if (analtyp  ==  'd') ianaltyp=30
        ## end if
        call blankout(resrange,1,29)
        if (analtyp  ==  'q') go to 9003
        lastanal=ianaltyp
        if (ianaltyp  <=   7  or  ianaltyp  ==  11  or 
     -    ianaltyp  ==  12  or  ianaltyp  ==  16  or 
     -    ianaltyp  ==  17  or  ianaltyp  ==  19  or 
     -    ianaltyp  ==  34  or  ianaltyp  ==  40  or 
     -    ianaltyp  ==  41) :
          call quiz(ansrun,nntyp,' ',' ',0,'Bond information source',23,
     -      0,5,6,00)
          if (nntyp  >  1) :
            call top_to_bond(nntyp,nneig,nhneig,ineig,iatnum,n,0,
     -        itemp1,itemp2,maxng,maxrec)
            innlist=innlistorig
          else:
            call askyn('Do you want to change bond thresholds',37,1,
     -        -1,newbondlims,0,8)
c           Update bond thresholds
            if (newbondlims  ==  1) :
              for ianl in range(0, nanos):
                ian=ialist(ianl)
                write (linein,2115) iatnm2(ian)
                call get#real(linein,23,ramax(ian),ramaxnew,1,000)
                if (ramaxnew  !=  ramax(ian)) :
                  ramax(ian)=ramaxnew
                  innlist=0
                ## end if
              
            ## end if
          ## end if
        ## end if
        if (ianaltyp  ==   3  or  ianaltyp  ==  15  or 
     -      ianaltyp  ==  30  or 
     -      ianaltyp  ==  35  or  ianaltyp  ==  36  or 
     -      ianaltyp  ==  37  or  ianaltyp  ==  38) :
c         These options do not work on trajectories
          itraj=0
        elif (ianaltyp  ==  22  or  ianaltyp  ==  23  or 
     -           ianaltyp  ==  24  or  ianaltyp  ==  25  or 
     -           ianaltyp  ==  19  or  ianaltyp  ==  32  or 
     -           ianaltyp  ==  39) :
c         These options only work on trajectories
          itraj=1
          if (ianaltyp  ==  22) write (6,2013)
        else:
          itrajdef=-1
          if ((ianaltyp  >=  5  and  ianaltyp  <=  7)  or 
     -         ianaltyp  ==  34  or  ianaltyp  ==  41  or 
     -         ianaltyp  ==  9  or  ianaltyp  ==  16  or 
     -         ianaltyp  ==  18  or  ianaltyp  ==  26) itrajdef=+1
          call askyn('Do you want to analyze a trajectory',35,1,
     -      itrajdef,itraj,0,1)
          if (itraj  ==  0) iaskcolcode=0
        ## end if
        ncl=0
        framefac=1.0
        iframeunit=1
        xtrajlab=xtrajlabs(iframeunit)
        if (itraj  ==  1) :
          if (ianaltyp  ==  05  or  ianaltyp  ==  06  or 
     -        ianaltyp  ==  07  or  ianaltyp  ==  09  or 
     -        ianaltyp  ==  16  or  ianaltyp  ==  18  or 
     -        ianaltyp  ==  21  or  ianaltyp  ==  22  or 
     -        ianaltyp  ==  23  or  ianaltyp  ==  24  or 
     -        ianaltyp  ==  27  or  ianaltyp  ==  28  or 
     -        ianaltyp  ==  29  or  ianaltyp  ==  34  or 
     -        ianaltyp  ==  39  or  ianaltyp  ==  41) :
            call quiz(unit,iframeunit,'f',' ',0,'trajectory unit',15,
     -        0,5,6,0)
            if (iframeunit  >  1) :
              xtrajlab=xtrajlabs(iframeunit)
              linein(1:11)=xtrajlab
              linein(12:45)=' per frames on the trajectory file'
              call get#real(linein,45,1.0,framefac,1,113)
              if (ianaltyp  ==  24) print *,'NOTE: ',
     -          'same time interval is assumed for both trajectories'
              for i in range(0, MAXFRAMES):
                xtraj(i)=i*framefac
              
            ## end if
          ## end if
        ## end if
        for ia in range(0, n):
          if (inpcrdtyp  <=  ioins) :
            atw(ia)=aw(iatnum(ia))
          else:
            atw(ia)=1.0
          ## end if
          rprox(ia)=1.0
        
        if (ianaltyp  ==   5  or  ianaltyp  ==   6   or 
     -      ianaltyp  ==   7  or  ianaltyp  ==  17  or 
     -      ianaltyp  ==  29) :
c         icharges=1: charge in iqcol1:iqcol2; icharges=2: charges read
          if (inpcrdtyp  <  ioa3pdb) :
            call readcharges(nread,nslt,n,charge,iatnum,isv,
     -        icharges,nerr)
            if (icharges  >  0) :
              call checkreschargesum(nslt,iresno,isegno,line,index,
     -          irescol1,irescol2,iresncol1,iresncol2,is1,is2,iqcol1,
     -          iqcol2,charge,molsltlim,'residue',7,resnames,ifres,
     -          ixres,itemp1,itemp2,nreschrg,iallzero,iallnonpos,
     -          iallnonneg,icharges,ifixq,reportqmin,mflogfile,maxrsd,
     -          maxrec)
              if (iallzero+iallnonpos+iallnonneg  >  0) icharges=0
            else:
              write (6,2035)
              call askstop(1)
            ## end if
          ## end if
        ## end if
        if (innlist  <  innlistorig)
     -    call nnlist(n,islvw,naslv,n,iatnum,ifchrg,c,nneig,nneiga,
     -      nhbneig,ineig,nhneig,nnneig,ncneig,nsneig,npneig,line,
     -      irescol1,irescol2,inamcol1,inamcol2,index,nconfig,innlist,
     -      molresflag,hblimfac,angmin,ihbondcalc,ibnd,indexo,isegno,
     -      ixres,maxrepconf,0,0,radtodeg,0,maxbox,maxng,maxrsd,maxrec)
        itrajrot=0
        ifailbond=0
        lpstitle=0
        npspages=0
        ipspage=0
        xm=800.0
        ym=900.0
        nframe=0
        nrep=0
        iw0=0
        iw1=0
        iw2=0
        iw4=0
        iw5=0
        ioutpr=0
        ipdfgrd=2
        iedit=0
        ndprow=0
        ndials=0
        noplotdist=0
        iskip2dplot=0
        ietotsaved=0
        ibondcorr=0
        iresbondcorr=0
        ibondprint=0
        itorcorr=0
        call zeroiti(nhbdist,0,mxbonds)
        call zeroit(rhbdist,mxbonds)
        extnam1=ext1(ianaltyp)
        extnam2=ext2(ianaltyp)
        extnam3=ext3(ianaltyp)
        call indexit(indexs,1,MAX2D,0)
        call indexit(index2d,1,MAX2D,0)
        call indexit(ixclst,1,MAX2D,0)
        if (itraj  ==  0) :
          if (ianaltyp  ==   9  or  ianaltyp  ==  11) extnam1='    '
          if (ianaltyp  ==  16  or  ianaltyp  ==  21  or 
     -        ianaltyp  ==  27  or  ianaltyp  ==  29) extnam2='    '
          if (ianaltyp  ==  20) :
            call askyn('Do you want to write a grid potential file',42,
     -        1,-1,igridfile,0,0)
            if (igridfile  ==  1) extnam2='.grd'
          ## end if
        ## end if
        if (ianaltyp  ==  23  or  ianaltyp  ==  24) :
          if (nslt  >  200  and  n/nslt  >=  2) :
            write (6,2030)
            call askstop(0)
          ## end if
          if (ianaltyp  ==  23) call askyn(
     -      'Do you want to read an existing .rd2 file',41,
     -       1,-1,iread2d,90,0)
          if (ianaltyp  ==  24) call askyn(
     -      'Do you want to read existing .rd2 and .rdx files',48,
     -       1,-1,iread2d,91,0)
          if (iread2d  >  0) :
            iw4=44
            itraj=0
            extnam2='    '
          ## end if
        ## end if
        if (itraj  ==  1) :
          nres=ixres(n)
          if (n  >  nslt  or  inpcrdtyp  ==  iommc) :
c           Ext# end arrays with more solvents than the input has
            nadd=min0(maxrsd-nres-1,(maxrec-index(n))/naslv-1)
c           print *,'nadd=',nadd,' iasv=',iatnum(nslt+1),iatnum(nslt+2),
c    -         iatnum(nslt+3)
            print *,'Generating information for ',nadd,' additional ',
     -         'solvents for possible use'
            if (n  ==  nslt) :
              linein(1:45)=
     -          'Number of atoms in a solvent (     ) molecule'
              linein(31:35)=resnamslv
              if (resnamslv(1:3)  ==  'HOH'  or 
     -            resnamslv(1:3)  ==  'TIP') :
                naslv=3
                iatnum(nslt+1)=8
                iatnum(nslt+2)=1
                iatnum(nslt+3)=1
                charge(nslt+1)=-0.834
                charge(nslt+2)=0.417
                charge(nslt+3)=0.417
              else:
                call getint(linein,45,999999,1,0,naslv,0)
                linein(1:32)='Atomic number of solvent atom 00'
                for ia in range(0, naslv):
                  write (linein(1:32),2091) ia
                  call getint(linein,32,999999,1,99,iatnum(nslt+ia),0)
                  write (linein(1:32),2092) ia
                  if (icharges  >  0)
     -              call get#real(linein,22,0.0,charge(nslt+ia),0,0)
                
              ## end if
            ## end if
            for is in range(0, nadd):
              nres=nres+1
              resnames(nres)=resnamslv
              for ia in range(0, naslv):
                iatnum(n+(is-1)*naslv+ia)=iatnum(nslt+ia)
                charge(n+(is-1)*naslv+ia)=charge(nslt+ia)
                index(n+(is-1)*naslv+ia)=index(n)+(is-1)*naslv+ia
                line(index(n+(is-1)*naslv+ia))=line(index(nslt+ia))
                ixres(n+(is-1)*naslv+ia)=nres
              
            
          ## end if
          if (analtyp  ==  'g'  or  analtyp  ==  'b'  or 
     -        analtyp  ==  'y'  or  analtyp  ==  'o'  or 
     -        analtyp  ==  'm'  or  analtyp  ==  'x') :
            call askyn(
     -        'Do you want to rotate each frame of the trajectory',50,
     -        1,-1,itrajrot,0,0)
            if (itrajrot  >  0) call genrot(trajrot,pi,iax,angle)
          ## end if
        ## end if
        mark0=' '
        lmark0=1
        mark1=' '
        lmark1=1
        mark2=' '
        lmark2=1
        analfile=inpfile
        if (extnam1  !=  '    ') :
c         Open output file(s)
          iw0=40
          if (namleno  !=  0) close (iw0)
          call strip_cext(analfile,namleni,namleno,lenext)
          analfile(namleno+1:namleno+4)=extnam1
          namleno=namleno+4
c         If output is a pdb file, add .pdb extension
          if (extnam1  ==  '.cvl'  or  extnam1  ==  '.hph') :
            analfile(namleno+1:namleno+4)='.pdb'
            namleno=namleno+4
            mark0='REMARK '
            lmark0=7
          ## end if
          if (ianaltyp  ==  14  or  ianaltyp  ==  15) :
c           Select circular variance type
            call quiz(ansrun,icvtyp,'o',' ',0,
     -        'circular variance type',22,0,5,6,0)
              if (icvtyp  ==  2) :
                analfile(namleno:namleno)='w'
                analfile1(namleno1-3:namleno1-3)='w'
              ## end if
          ## end if
          if (iread2d  >  0) :
            analfile4=analfile
            lanalfile4=namleno
            analfile(namleno-3:namleno+2)='_a.rd2'
            namleno=namleno+2
            write (6,2036) analfile(1:namleno)
          ## end if
          call openfile(iw0,0,'analysis',8,'new',analfile,namleno,
     -      notfnd,0,1,1,0,0)
          call datprt(iw0,version,1,mark0,lmark0,hostname,lhostname,
     -      iheadnode,0)
          if (extnam2  !=  '    ') :
            iw1=41
            if (namleno1  !=  0) close (iw1)
            namleno1=namleno
            analfile1=analfile(1:namleno)
            if (extnam2  ==  '.ps ') :
c             Second file is PS - app# end the .ps to the name
              analfile1(namleno+1:namleno+3)='.ps'
              namleno1=namleno+3
              mark1='% '
              lmark1=2
            elif (extnam2  ==  '.pdb') :
c             Second file is PDB - app# end the .pdb to the name
              analfile1(namleno+1:namleno+4)='.pdb'
              namleno1=namleno+4
            else:
c             Replace extension
              analfile1(namleno-3:namleno)=extnam2
              namleno1=namleno
            ## end if
            if (itraj  ==  1  or  (ianaltyp  !=  5  and 
     -          ianaltyp  !=  6  and  ianaltyp  !=  7  and 
     -          ianaltyp  !=  34  and  ianaltyp  !=  41))
     -        call openfile(iw1,0,'analysis',8,'new',analfile1,namleno1,
     -          notfnd,0,1,1,0,0)
            if (extnam2  !=  '.grd'  and  extnam2  !=  '.ps ')
     -        call datprt(iw1,version,1,mark1,lmark1,hostname,lhostname,
     -          iheadnode,0)
          ## end if
          if (extnam3  !=  '    ') :
            if (namleno2  !=  0) close (iw2)
            iw2=42
            namleno2=namleno
            analfile2=analfile(1:namleno)
            analfile2(namleno-3:namleno)=extnam3
            namleno2=namleno
            if (extnam3  ==  '.rdp') :
              analfile2(namleno2+1:namleno2+3)='.ps'
              namleno2=namleno2+3
            ## end if
          ## end if
        ## end if
        call datprt(6,version,1,' ',1,hostname,lhostname,iheadnode,1)
        if (ianaltyp  ==  23  or  ianaltyp  ==  24  or 
     -      ianaltyp  ==  25) :
          if (iheadnode  ==  1  and  iread2d  ==  0) write (6,2037)
        ## end if
        if (ianaltyp  ==  12  or  ianaltyp  ==  18  or  ianaltyp  ==  19
     -     or  ianaltyp  ==  21  or  ianaltyp  ==  32) :
          call askyn(
     -      'Do you want to calculate angle probability distributions',
     -      56,1,1,ipdfc,0,0)
          if (ipdfc  ==  1) :
            ioutpr=iw0
8001        call getint('Angle distribution grid size',28,2,1,30,
     -        ipdfgrd,000)
            if (mod(360,ipdfgrd)  !=  0) :
              print *,'Grid size is not a divisor of 360'
              go to 8001
            ## end if
            call askyn(
     -        'Do you want to map the angle distributions on the dials',
     -        55,1,1,mappdf,0,0)
          ## end if
          call askyn('Do you want to draw arcs in the dial plots',42,1,
     -      1,iconndial,0,0)
        ## end if
        if (ianaltyp  ==   1) :
          write (6,2024) analfile(1:namleno)
          call askyn('Do you want to print bond angles',32,1,-1,iangpr,
     -      0,0)
          call askyn('Do you want to print torsion angles',35,1,-1,
     -      itorpr,0,0)
          call findchiral(nslt,iatnum,nneig,nhneig,ineig,indexa,maxng)
        elif (ianaltyp  ==   2) :
          write (6,2000) analfile(1:namleno)
        elif (ianaltyp  ==  18) :
          if (inpcrdtyp  >  ioins) :
            print *,'ERROR: this input format does not carry atom names'
            return
          ## end if
          if (iresno(nslt)  >  maxframe) :
            write (6,2074) maxframe,'residues',' '
            return
          ## end if
          write (6,2019) 'Psi and phi angles',analfile(1:namleno)
          write (6,2019) 'Ramachandran plot (Postscript)',
     -      analfile1(1:namleno1)
          xm=500.0
          ym=620.0
          if (itraj  ==  0) :
            call openps(iw1,xm,ym,title,ltitle76,'Ramachandran plot',
     -        17,analfile,namleni,analfile,0,1,ipspage)
            write (6,2086) 'residue',' ',' ',nres
          else:
            pstitle='Ramachandran plot'
            lpstitle=17
            write (6,2086) 'frame',' max(frame)'
            nxselres=0
            maxrp=(maxpres-1)/2
            if (nres  <=  maxrp) :
              call askyn('Do you want dial plots for all residues',39,
     -          1,-1,ialldial,0,0)
              if (ialldial  ==  1) :
                call indexit(ixselres,1,nres,0)
                nxselres=nres
              ## end if
            else:
              ialldial=0
              print *,'Maximum number of dial plots can be ',
     -          'tracked=',maxrp
              call getint('Number of residues to track',27,0,1,maxrp,
     -          nxselres,0)
            ## end if
            if (nxselres  >  0) :
              question(1:38)='Phi \ Psi dial pair #   residue number'
              for ix in range(0, nxselres):
                write (question(22:23),2001) ix
8003            isegix=1
                if (ialldial  ==  0) :
                  if (nsegm  >  1) :
                    isegdef=999999
                    question(25:31)='segment'
                    call getint(question,38,isegdef,1,nsegm,isegix,00)
                    isegdef=isegix
                    question(25:31)='residue'
                  ## end if
                  call getint(question,38,999999,1,nres,iresix,00)
                  call findsegres(isegno,iresno,ixres,1,nslt,isegix,
     -              iresix,iaf,ixselres(ix),ifail)
                  if (inpcrdtyp  <=  ioins)
     -              write (6,2111) isegix,iresix,
     -                resnames(ixselres(ix))(1:nrescol)
                ## end if
                irfound=0
                icafound=0
                for ia in range(0, nslt):
                  if (iresno(ia)  ==  ixselres(ix)) :
                    irfound=1
                    atnam=line(index(ia))(inamcol1:inamcol1+3)
                    call leftadjust4(atnam,atnam)
                    if (atnam  ==  'CA  ') icafound=ia
                  ## end if
                
                if (icafound*irfound  ==  0) :
                  if (irfound  ==  0) :
                    print *,'Residue # ',ixselres(ix),' is not found'
                  else:
                    print *,'Residue # ',ixselres(ix),' has no CA atom'
                  ## end if
                  if (ialldial  ==  0) go to 8003
                  ixselres(ix)=0
                else:
                  call ca_to_bb(icafound,iresno,nneig,ineig,index,line,
     -              inamcol1,i1,i2,i4,i5,ires,iprotein,
     -              maxng,maxrec)
                  if (iprotein  ==  0) :
                    write (6,2104) ixselres(ix)
                    if (ialldial  ==  0) go to 8003
                    ixselres(ix)=0
                  ## end if
                ## end if
              
              if (ialldial  ==  1) :
                ndel=0
                for ix in range(0, nxselres):
                  if (ixselres(ix)  ==  0) :
                    ndel=ndel+1
                  else:
                    ixselres(ix-ndel)=ixselres(ix)
                  ## end if
                
              ## end if
              nxselres=nxselres-1
              if (nxselres  >  1)
     -          call getint('Number of dials to draw in a line',33,
     -            min0(4,nxselres),1,12,ndprow,55)
              ix0=8
              if (ndprow  >  4) ix0=4
              if (ndprow  >  8) ix0=2
              for ix in range(0, nxselres):
                write (ramalab(2*ix-1)(ix0:ix0+13),2023)
     -            ixselres(ix),resnames(ixselres(ix)),'Phi'
                write (ramalab(2*ix)(ix0:ix0+13),2023) ixselres(ix),
     -            resnames(ixselres(ix)),'Psi'
                lramalab(2*ix-1)=ix0+13
                lramalab(2*ix)=ix0+13
              
              call openfile(iw2,0,'dial plot',9,'new',analfile2,
     -          namleno2,notfnd,0,1,1,0,0)
              call datprt(iw0,version,1,mark0,lmark0,hostname,lhostname,
     -          iheadnode,0)
              write (6,2019) 'Phi and Psi dial plots',
     -          analfile2(1:namleno2)
              call askyn(
     -         'Do you want the Ramachandran plot with all residues',51,
     -            1,-1,iallrama,0,0)
            else:
              print *,'Only the Ramachandran plot will be prepared ',
     -          'using all residues'
              iallrama=1
            ## end if
            call ramachandran_init(n,ixres)
          ## end if
        elif (ianaltyp  ==  19  or  ianaltyp  ==  32) :
          write (6,2019) 'Selected angles',analfile(1:namleno)
          xm=500.0
          ym=620.0
          if (itraj  ==  1) :
            write (6,2019) 'Dial plots for the selected angles',
     -        analfile1(1:namleno1)
          ## end if
          if (ianaltyp  ==  19) :
c           Torsion angles
            call torslistinp(ixtor1234,talab,ltalab,ntorsel,inpcrdtyp,
     -        ioins,line,index,n,nslt,nresslt,iatnum,nneig,ineig,nhneig,
     -        ixres,ixresno,nsegm,isegno,ifres,ilres,iresno,indexs,
     -        listrefres,irescol1,irescol2,inamcol1,inamcol2,maxng,
     -        maxrsd,maxrec,MAXCOPY1)
            ndials=ntorsel
            if (ndials  ==  0) :
              print *,'ERROR: no torsions were selected'
              go to 9003
            ## end if
          else:
c           Bond angles
            call getint('Number of angles to track',25,0,1,49,
     -        nangsel,0)
            for it in range(0, nangsel):
8004          question(1:35)='Angle   , atomindices ( 3 numbers )'
              write (question(6:8),2043) it
              call getintline(question,35,1,nslt,ixtor1234(1,it),3,0)
              if (inpcrdtyp  <=  ioins)
     -          write (6,2039) 'A',(ixtor1234(k,it),
     -            line(index(ixtor1234(k,it)))(inamcol1:inamcol2),
     -            line(index(ixtor1234(k,it)))(irescol1:irescol2),
     -          k=1,3)
c             print *,' ixtor:', (ixtor1234(k,it),k=1,4)
              nerr=0
              for k in range(2, 3):
                if (isbonded(ixtor1234(k-1,it),ixtor1234(k,it),
     -              nneig,ineig,n,maxng)  ==  0) :
                  write (6,2066) (ixtor1234(kk,it),kk=k-1,k),' NOT '
                  nerr=nerr+1
                ## end if
              
              if (isbonded(ixtor1234(1,it),ixtor1234(3,it),
     -            nneig,ineig,n,maxng)  ==  1) :
                write (6,2066) ixtor1234(1,it),ixtor1234(3,it),' '
                nerr=nerr+1
              ## end if
              if (nerr  >  0) :
                call askyn('Do you want to use this angle',29,1,-1,iok,
     -            0,0)
                if (iok. eq. 0) go to 8004
              ## end if
              call blankout(talab(it),1,30)
              if (inpcrdtyp  >  ioins) :
                write (talab(it),2067) 'A',it,(ixtor1234(k,it),k=1,3)
              else:
                write (talab(it),2108)
     -            (line(index(ixtor1234(kk,it)))(inamcol1:inamcol1+3),
     -            kk=1,2),line(index(ifres(ixres(ixtor1234(2,it)))))
     -             (irescol1:irescol2),
     -            line(index(ixtor1234(3,it)))(inamcol1:inamcol1+3)
              ## end if
              ltalab(it)=25
            
            ndials=nangsel
          ## end if
          if (ndials  >  1) :
            call getint('Number of dials to draw in a line',33,
     -            min0(ndials,3),1,12,ndprow,55)
            call quiz(ans,itorcorr,'j',' ',0,
     -        'Angular (circular) correlation type',35,0,5,6,0)
            itorcorr=itorcorr-1
          else:
            ndprow=1
            itorcorr=0
          ## end if
        elif (ianaltyp  ==   3) :
          write (6,2019) 'defal group and bond lists',
     -      analfile(1:namleno)
          write (6,2019) 'Backbone atoms',analfile1(1:namleno1)
        elif (ianaltyp  ==   4) :
          write (6,2019) 'Bond lengths',analfile(1:namleno)
        elif (ianaltyp  ==  17) :
c         Hydrogen-bond bridge analysis
          ibondtype=6
          write (6,2002) 'ed bridge ',analfile(1:namleno)
          write (iw0,2127) bondname(ibondtype)(1:lbondname(ibondtype))
          call set_hbondlims(hbf0,hblimfac,angm0,angmin,iw0)
          brslv=resnamslv
          nabr=naslv
          if (islvw  ==  1) islvw=2
          if (nslt  ==  n) :
901         call getname(brslv,lbrslv,'Bridge residue name',19,
     -        nrescol,'',0,0,100,0)
            if (brslv(1:nrescol)  !=  resnamslv(1:nrescol)) :
c             Check if residue exists, find # of atoms
              ia=1
              ir=1
              while (resnames(ir)(1:nrescol)  !=  brslv(1:nrescol)
     -            and  ia  <  nslt)
                ia=ia+1
                ir=ixres(ia)
              
              if (resnames(ir)(1:nrescol)  ==  brslv(1:nrescol)) :
                iaf=ia
                while (resnames(ir)(1:nrescol)  ==  brslv(1:nrescol)
     -              and  ia  <  nslt)
                  ia=ia+1
                  ir=ixres(ia)
                
                nabr=ia-iaf
                if (resnames(ir)(1:nrescol)  ==  brslv(1:nrescol))
     -            nabr=nabr+1
              else:
                print *,'ERROR: residue ',brslv(1:nrescol),' not found'
                go to 901
              ## end if
            ## end if
          ## end if
          write (6,2042) brslv(1:nrescol),brslv(1:nrescol),nabr
          write (iw0,2042) brslv(1:nrescol),brslv(1:nrescol),nabr
          call getint('Maximum number of bridge residues in a bridge',
     -      45,MAXBRIDGELEN-1,1,MAXBRIDGELEN-1,maxbridgemem,0)
          call gethbanchordef(line,index,nslt,ixres,iresno,iatnum,
     -      indexa,indexo,nanchor,ianchor,iiqcol(1,inpcrdtyp),
     -      iiqcol(2,inpcrdtyp),inpcrdtyp,iobpdb,iocpdb,icharges,qmin,
     -      iw0,inamcol1,resnames,brslv,nrescol,segid4,molsltlim,
     -      nsegslt,isc,nneig,ineig,'Bridge',6,ianchor2,iselfanc,
     -      nosameseg,iqfsel2,22,ifail,maxng,maxrsd,maxrec,
     -      MAXBRIDGEATOM)
          if (ifail  ==  1) go to 9005
          call zeroiti(nbridgetype,0,nanchor*MAXBRIDGELEN)
          ipb=2*itraj+1
          call askyn('Do you want to print the bridges',32,1,ipb,
     -      listbridge,0,0)
          print *,'Bridge of length 1 is a solute-solute hydrogen bond'
          call getint('Minimum length of bridge to list',32,1,1,
     -      maxbridgemem,minbridgelenprint,000)
          if (itraj  ==  1) :
            call getint('Minimum percent of frames present to list',41,
     -        0,1,100,minbridgepercprint,000)
          else:
            minbridgepercprint=0
          ## end if
          write (iw0,*)
          call zeroiti(lpath,0,MAXBRIDGETYPE*MAXBRIDGELEN*
     -      MAXBRIDGEATOM)
        elif (ianaltyp  >=   5  and  ianaltyp  <=  7  or 
     -           ianaltyp  ==  34  or  ianaltyp  ==  41) :
          npspages=3
          ibondtype=ianaltyp-4
          if (ianaltyp  ==  34) ibondtype=4
          if (ianaltyp  ==  41) ibondtype=5
          if (iallheavy  ==  1  and  ibondtype  ==  2) :
            print *,'There are no hydrogens in the structure'
            call askstop(1)
          ## end if
          if (itraj  ==  1) :
            linein='Do you want to calculate '//
     -        bondname(ibondtype)(1:lbondname(ibondtype))//
     -        ' correlation'
            llinein=37+lbondname(ibondtype)
            call askyn(linein,llinein,1,-1,ibondcorr,0,0)
            if (ibondcorr  >  0) npspages=npspages+1
            linein='Do you want to calculate residue-aggregated '//
     -        'bond correlation'
            llinein=60
            call askyn(linein,llinein,1,-1,iresbondcorr,0,0)
            if (iresbondcorr  >  0) npspages=npspages+1
            linein='Do you want to print every '//
     -        bondname(ibondtype)(1:lbondname(ibondtype))
            llinein=27+lbondname(ibondtype)
            call askyn(linein,llinein,1,-1,ibondprint,0,0)
            if (ibondprint  ==  1  and  ibondtype  ==  1)
     -      write (iw0,2141)
            noprintconf=1-ibondprint
            iaskcolcode=1
            call askyn(
     -        'Do you want to read previously generated bond tracks',52,
     -        1,-1,ireadtracks,0,0)
            if (ireadtracks  >  0) :
              ifail=1
              while (ifail  >  0)
                call getname(trackfile,ltrackfile,'Name of track file',
     -            18,80,'',0,0,0,0)
                iout_track=92
                call openfile(iout_track,0,'previously written track',
     -            24,'old',trackfile,ltrackfile,notfnd,0,1,1,0,0)
                call readtrack(iout_track,iw0,30,nbfound,nbresfound,
     -            nres2d,ixres,trackfile,ltrackfile,ianc_anc,ifail,
     -            maxrec)
              
              xtrajlab=xtrajlabs(iframeunit)
              npspages=2
              go to 3001
            ## end if
            call getint(
     -        'Number of frames to average in the bond number plot',51,
     -        1,1,MAXFRAMES/2,nbondavg,00)
          ## end if
          if (ianaltyp  ==   5) :
c           Set up H bond analysis
            write (6,2002) ' ',analfile(1:namleno)
            write (iw0,2127) bondname(ibondtype)(1:lbondname(ibondtype))
            call set_hbondlims(hbf0,hblimfac,angm0,angmin,iw0)
            nres2d=ixres(n)
            if (n  >  nslt) :
              call askyn('Do you want to include the solvents',35,1,-1,
     -          isolvent,0,0)
              if (isolvent  ==  0) nres2d=ixres(nslt)
            ## end if
            write (6,2140)
            for i in range(0, nanos):
              if (icatlist(i)  ==  0) :
                if (ialist(i)  !=  1  and  ialist(i)  !=  6)
     -            write (6,2142) iatnm2(ialist(i)),
     -              sqrt(ramax2(ialist(i))*hblimfac)
              elif (islvw  >  0) :
                write (6,2143)
     -            iatnm2(ialist(i)),sqrt(ramax2(ialist(i))*hblimfac)
              ## end if
            
            if (itraj  ==  1) :
              call gethbanchordef(line,index,nslt,ixres,iresno,iatnum,
     -          indexa,indexo,nanchor,ianchor,iiqcol(1,inpcrdtyp),
     -          iiqcol(2,inpcrdtyp),inpcrdtyp,iobpdb,iocpdb,icharges,
     -          qmin,iw0,inamcol1,resnames,'      ',nrescol,segid4,
     -          molsltlim,nsegslt,isc,nneig,ineig,'Hydrogen bond',13,
     -          ianchor2,iselfanc,nosameseg,iqfsel2,72,ifail,maxng,
     -          maxrsd,maxrec,MAXBRIDGEATOM)
              if (ifail  ==  1) go to 9005
              call indexit(it3,1,MAXBONDS,0)
              call indexit(it4,1,MAXBONDS,0)
              iaskcolcode=1
            ## end if
            innlist=0
          elif (ianaltyp  ==   6) :
c           Set up hydrophobic contact analysis
            write (6,2019) 'Hydrophobic bond list',analfile(1:namleno)
            write (iw0,2127) bondname(ibondtype)(1:lbondname(ibondtype))
            call get#real('Hydrophopbic bond length limit',30,5.0,
     -        rhphmax,1,76)
            call getint('Minimum number of carbon-bonded hydrogens',41,
     -        1,1,4,nhneigmin,77)
            call gethphanchordef(line,index,nslt,iresno,iatnum,charge,
     -        indexa,indexn,indexo,nneig,ineig,nhneig,nhneigmin,nanchor,
     -        ianchor,ianchor2,iselfanc,nosameseg,iallheavy,iqcol1,
     -        iqcol2,inpcrdtyp,iobpdb,iocpdb,icharges,ibondtype,iw0,
     -        segid4,molsltlim,nsegslt,bondname(ibondtype),
     -        lbondname(ibondtype),ifail,maxng,maxrec,MAXBRIDGEATOM)
            if (ifail  ==  1) go to 9005
            call ext# end_nnlist(nneig,ineig,npneig,nslt,maxng,maxrec)
c            do ia=1,nslt
c              atnam=line(index(ia))(inamcol1:inamcol1+3)
c              write (77,7719) ia,atnam,nneig(ia),npneig(ia),
c     -          (ineig(in,ia),in=1,npneig(ia))
c            
c7719        format(i5,1x,a,' nn=',i3,' nnp=',i6,' in=',(10i5))
            nres2d=ixres(nslt)
            innlist=0
          elif (ianaltyp  ==   7) :
c           Set up salt-bridge contact analysis
            write (iw0,2127) bondname(ibondtype)(1:lbondname(ibondtype))
            write (6,2019) 'Salt-bridge list',analfile(1:namleno)
            call get#real('Salt-bridge length limit',24,5.0,rsltbmax,1,2)
            call getsltbanchordef(line,index,nslt,iresno,iatnum,charge,
     -        indexa,indexn,indexo,temp,nhneig,ineig,nanchor,ianchor,
     -        ianchor2,iselfanc,nosameseg,iqcol1,iqcol2,inpcrdtyp,
     -        iobpdb,iocpdb,icharges,iw0,irescol1,irescol2,inamcol1,
     -        inamcol2,segid4,molsltlim,nsegslt,isegno,ifail,maxng,
     -        maxrec,MAXBRIDGEATOM)
            if (ifail  ==  1) go to 9005
            call ext# end_nnlist(nneig,ineig,npneig,nslt,maxng,maxrec)
            nres2d=ixres(nslt)
            innlist=0
          elif (ianaltyp  ==   34) :
c           Set up heavy-atom VdW contact analysis
            write (iw0,2127) bondname(ibondtype)(1:lbondname(ibondtype))
            write (6,2019) 'Heavy-atom contact list',analfile(1:namleno)
            call get#real('Heavy-atom distance threshold',29,5.0,
     -        rhphmax,1,76)
            call gethphanchordef(line,index,nslt,iresno,iatnum,charge,
     -        indexa,indexn,indexo,nneig,ineig,nhneig,nhneigmin,nanchor,
     -        ianchor,ianchor2,iselfanc,nosameseg,iallheavy,iqcol1,
     -        iqcol2,npcrdtyp,iobpdb,iocpdb,icharges,ibondtype,iw0,
     -        segid4,molsltlim,nsegslt,bondname(ibondtype),
     -        lbondname(ibondtype),ifail,maxng,maxrec,MAXBRIDGEATOM)
            if (ifail  ==  1) go to 9005
            call ext# end_nnlist(nneig,ineig,npneig,nslt,maxng,maxrec)
c            do ia=1,nslt
c              atnam=line(index(ia))(inamcol1:inamcol1+3)
c              write (77,7719) ia,atnam,nneig(ia),npneig(ia),
c     -          (ineig(in,ia),in=1,npneig(ia))
c            
c7719        format(i5,1x,a,' nn=',i3,' nnp=',i6,' in=',(10i5))
            npspages=2
            if (itraj  ==  1) :
              ibondcorr=0
              call askyn('Do you want to print the heavy-atom contacts',
     -          44,1,+1,ihphprint,0,0)
              noprintconf=1-ihphprint
              iaskcolcode=1
            ## end if
            nres2d=ixres(nslt)
            innlist=0
          elif (ianaltyp  ==   41) :
c           Set up mutually proximal contact analysis
            write (iw0,2127) bondname(ibondtype)(1:lbondname(ibondtype))
            call getmpxbdef(nslt,indexa,indexov,indexn,segid4,iresno,
     -        molsltlim,nsegslt,nanchorr,nanchorn,iw0,maxrsd)
            call get#real('Maximum distance for contact',28,99999.0,
     -        rmpxlim,1,000)
            if (rmpxlim  <  99999.0) write (iw0,2128) rmpxlim
c           call ext# end_nnlist(nneig,ineig,npneig,nslt,maxng,maxrec)
            nres2d=ixres(nslt)
          ## end if
        elif (ianaltyp  ==  8) :
c         Set up residue contact list calculation
          if (inpcrdtyp  >  ioins) :
            print *,'Input format does not have residue information'
            go to 9002
          ## end if
          call modrepats
          call get#real('Threshold distance with representative atoms',
     -      44,10.0,resdistlim,1,38)
          call get#real('Threshold distance with closest approach',40,
     -      10.0,resapplim,1,39)
          call askyn(
     -      'Do you want to ignore hydrogens for the closest approach',
     -      56,1,1,ignoreh,106,0)
          if (iresshift  >  0) write (6,2029)
          isegdef1=1
          call getresrange(nsegm,indexs,isegno,ixres,iresno,ifres,
     -      'reference residue to use',24,nresslt,nslt,irefres1,
     -      irefres2,isegdef1,irefseg1,irefseg2,listrefres,nrefres,
     -      nrefrange,0,0,maxrsd,maxrec,111)
          write (6,*)
          if (isegdef1  >  0) isegdef2=nsegm
          if (isegdef1  ==  nsegm) isegdef2=1
          call getresrange(nsegm,indexs,isegno,ixres,iresno,ifres,
     -      'neighbourhood residue to use',28,nresslt,nslt,
     -      inegres1,inegres2,isegdef2,inegseg1,inegseg2,listnegres,
     -      nnegres,nnegrange,0,0,maxrsd,maxrec,111)
          write (6,2050) 'representative',analfile(1:namleno)
          write (6,2050) 'closest',analfile1(1:namleno1)
          if (itraj  >  0) :
            if (irefres2-irefres1  <  MAX2D  and 
     -         inegres2-inegres1  <  MAX2D) :
              print *,'By default, average distances are based on the ',
     -          'representative atoms'
              itypavg=1
              call askyn(
     -          'Do you want to average the contact distances instead',
     -           52,1,-1,iusecont,000,0)
              itypavg=itypavg+iusecont
              npspages=1
              iw4=44
              if (itypavg  ==  1) :
                analfile4=analfile
                namleno4=namleno
              else:
                analfile4=analfile1
                namleno4=namleno1
              ## end if
              analfile4(namleno4+1:namleno4+3)='.ps'
              namleno4=namleno4+3
              call openfile(iw4,0,'average distance matrix',23,'new',
     -          analfile4,namleno4,notfnd,0,1,1,0,0)
              call openps(iw4,500.0,500.0,' ',1,' ',1,inpfile,0,
     -          inpfile,0,1,ipspage)
            else:
              write (6,2112) MAX2D
            ## end if
          ## end if
          increst=0
          incsolvrr=0
          numresrr=numres
          idisjoint=1
          if (inegres1  <=  irefres2  and  inegres2  >=  irefres1) :
            match=0
            for irr in range(0, nrefres):
              for irn in range(0, nnegres):
                if (listrefres(irr)  ==  listnegres(irn)) match=1
              
            
            idistjoint=1-match
          ## end if
          if (idisjoint  ==  1) :
c           The ranges are disjoint
            call openfile(iw2,0,'contact list',12,'new',analfile2,
     -        namleno2,notfnd,0,1,1,0,0)
            call datprt(iw2,version,1,mark0,lmark0,hostname,lhostname,
     -        iheadnode,0)
            write (6,2050) 'contact',analfile2(1:namleno2)
            call askyn(
     -        'Do you want all neighbors of the interface residues',51,
     -        1,1,increst,57,0)
            if (n  >  nslt  and  increst  ==  1) :
              call askyn('Do you want to include solvent neighbors',40,
     -          1,-1,incsolvrr,0,0)
              if (incsolvrr  ==  0) :
                numresrr=nresslt
                incsolvrr=-1
              ## end if
            ## end if
          else:
            print *,'NOTE: residue ranges overlap; no contact list ',
     -        'will be generated'
            iw2=0
          ## end if
          call zeroiti(irescount1,0,nresslt)
          call zeroiti(irescount2,0,nresslt)
          call zeroiti(irescount3,0,nresslt)
          irefresinc=ixresno(ixres(molsltlim(1,irefseg1)))-
     -      ixres(molsltlim(1,irefseg1))
          inegresinc=ixresno(ixres(molsltlim(1,inegseg1)))-
     -      ixres(molsltlim(1,inegseg1))
        elif (ianaltyp  ==  33) :
c         Set up molecule-molecule distance calculation
          if (nmolslt  ==  1) :
            write (6,*) 'Solute is a single molecule'
            go to 9003
          ## end if
          write (6,2120) 'COM-COM distance',analfile(1:namleno)
          write (6,2120) 'closest distance',analfile1(1:namleno1)
          call askyn(
     -      'Do you want to ignore hydrogens for the closest approach',
     -      56,1,1,ignoreh,106,0)
          write (iw0,2121) 'center-of-mass distance'
          write (iw1,2121) 'closest approach'
        elif (ianaltyp  ==  36  or  ianaltyp  ==  37  or 
     -           ianaltyp  ==  38) :
          iwpdb=0
          if (iresshift  >  0) write (6,2156)
          if (ianaltyp  ==  36) :
c           Compare two residue distance matrices
            call compare_rrdist(resnames,nrescol,itemp1,temp,irefres1,
     -        irefres2,analfile4,lanalfile4,maxrsd)
            call askyn(
     -        'Do you want to write a PDB file with average changes',52,
     -         1,0,iwpdb,0,0)
          elif (ianaltyp  ==  37) :
c           Compare two residue bond (HB, SB, or HP) matrices from Simulaid log
            nres=ixres(nslt)
            call compare_bondmat(resnames,nrescol,itemp1,temp,nres,
     -        analfile4,lanalfile4,maxrsd)
            irefres1=1
            irefres2=nres
            call askyn(
     -        'Do you want to write a PDB file with cumulative changes',
     -         55,1,0,iwpdb,0,0)
          else:
c           Compare two residue RMSF lists from Simulaid log
            call askyn(
     -        'Do you want to write a PDB file with the differences',52,
     -        1,0,iwpdb,0,0)
            siglev=0.0
            if (iwpdb  ==  1) call get#real(
     -        'Significance level threshold to set to 0 the difference',
     -        55,0.05,siglev,1,000)
            call compare_rmsf(resnames,nrescol,siglev,temp,analfile4,
     -        lanalfile4,irefres1,irefres2,maxrsd)
          ## end if
          if (iwpdb  >  0) :
            for ir in range(irefres1, irefres2):
              for ia in range(ifres(ir), ilres(ir)):
                cv(ia)=temp(ir)
              
            
            analfile4(lanalfile4:lanalfile4+2)='db'
            lanalfile4=lanalfile4+1
            iw4=44
            call openfile(iw4,0,'Avg difference labeled PDB',26,
     -        'new',analfile4,lanalfile4,notfnd,0,1,1,0,0)
            call writeconf(iw4,inpcrdtyp,iobpdb,inpcrdtyporg,
     -        nslt,nslt,nslt,naslv,islvw,iasv,namesv,qsv,pflsv,1,
     -        iwhead,0,iatnum,ifchrg,nconfig,innlist,c,rprox,cv,ixres,
     -        iresno,atnames,resnames,segnames,charge,isegno,altcol,
     -        inscol,ninsres,marker,ntitlin,ntitlinw,title,ireseq,
     -        iresnrestart,iresidrestart,nneig,nneiga,nhbneig,ineig,
     -        nhneig,nnneig,ncneig,nsneig,npneig,numres,numslv,
     -        resnamslv,line,blankline,mmtype,ibnd,index,indexn,indexo,
     -        1,molresflag,irescount3,itemp1,hblimfac,angmin,0,1,1,1,0,
     -        3,iqspaceask,ianaltyp,0,0.0,0,0,0,keeprem,iwriteatsym,
     -        radtodeg,maxrepconf,maxng,maxrsd,maxrec)
            close (iw4)
          ## end if
        elif (ianaltyp  ==  31) :
c         Set up adjacency matrix-based analysis
          if (inpcrdtyp  >  ioins) :
            print *,'Input format does not have residue information'
            go to 9002
          ## end if
          write (6,2046) analfile(1:namleno)
          call quiz(analtyp,idtyp,' ',' ',0,
     -      'residue distance definition',27,0,5,6,0)
          if (analtyp  ==  'r') :
            resdistdef=7.5
            call modrepats
            ignoreh=0
            irepuse=1
          else:
            call askyn(
     -       'Do you want to ignore hydrogens for the closest approach',
     -        56,1,1,ignoreh,106,0)
            irepuse=0
            resdistdef=5.0
          ## end if
          call get#real('Threshold distance',18,resdistdef,resdistlim,1,
     -      38)
          call quiz(analtyp,iadjtyp,' ',' ',0,'adjacency analysis',18,
     -      0,5,6,108)
          call getint('Highest exponent to raise the adjacency matrix',
     -      46,1,1,100,nexpmax,0)
          call getint('Power interval to plot',22,1,1,100,npint,0)
          call askyn(
     -      'Do you want to scale the column sums to [0,1] range',51,
     -      1,1,iscalesum,000,0)
          isegdef=1
          call getresrange(nsegm,indexs,isegno,ixres,iresno,ifres,
     -      'residue to use',14,nresslt,nslt,ires1,ires2,isegdef,
     -      iseg1,iseg2,listrefres,nresref,nrange,1,1,maxrsd,maxrec,000)
          call askyn('Do you want to mark residues',28,1,0,imarkres,109,
     -      0)
          if (imarkres  ==  1) :
            iw4=44
            namleno4=0
            call openfile(iw4,0,'residue mark',12,'old',analfile4,
     -        namleno4,notfound,3,1,1,0,0)
            read (iw4,2044,# end=771) (itemp4(i-ires1+1),i=ires1,ires2)
 771        close (iw4)
            nmarks=0
            for ii in range(ires1, ires2):
              i=ii-ires1+1
              if (itemp4(i)  >  nmarks) nmarks=itemp4(i)
            
            if (ipredict  ==  0) :
              for imark in range(0, nmarks):
                write (6,2056) imark
                read (5,2095) marks(imark)
              
            else:
              write (6,2045) (i,marks(i),i=1,9)
            ## end if
          ## end if
          npspages=1
          call openps(iw1,xm,ym,title,ltitle76,
     -      'Adjacency matrix analysis',25,inpfile,namleni,inpfile,0,1,
     -      ipspage)
        elif (ianaltyp  ==   9) :
c         Set up distance measuring
          if (itraj  ==  1)
     -       write (6,2019) 'Distances measured',analfile(1:namleno)
          write (6,2047)
          call getint('First  atom number',18,1,1,n,ia1,0)
          call getint('Second atom number',18,1,1,n,ia2,0)
          ifirstref=0
          if (itraj  ==  1) :
            if (ia1  ==  ia2) :
              ifirstref=1
            else:
              call askyn(
     -          'Is the first atom fixed on the input structure',46,
     -          1,-1,ifirstref,0,0)
            ## end if
            print *,'Progression of distance vector will be plotted ',
     -        'projected to a coordinate plane'
            call readax('Axis normal to the plane (1,2,3)',32,3,idax,
     -        indexax)
          ## end if
          if (ia1  ==  ia1) :
            write (atomdist,2052) ia1
            pstitle(1:15)=atomdist
            lpstitle=15
          else:
            write (atomdist,2053) ia1,ia2
            pstitle(1:22)=atomdist
            lpstitle=22
          ## end if
          npspages=4
          if (inpcrdtyp  <=  ioins) :
            ir1=iresno(ia1)
            ir2=iresno(ia2)
            write (6,2057) ia1,ir1,line(index(ia1))(inamcol1:inamcol2),
     -        line(index(ia1))(irescol1:irescol2),
     -        ia2,ir2,line(index(ia2))(inamcol1:inamcol2),
     -        line(index(ia2))(irescol1:irescol2)
            if (itraj  ==  1) write (iw0,2057) ia1,ir1,
     -        line(index(ia1))(inamcol1:inamcol2),
     -        line(index(ia1))(irescol1:irescol2),
     -        ia2,ir2,line(index(ia2))(inamcol1:inamcol2),
     -        line(index(ia2))(irescol1:irescol2)
          else:
            write (6,2058) ia1,ia2
            if (itraj  ==  1) write (iw0,2058) ia1,ia2
          ## end if
          if (ifirstref  ==  1) :
            write (iw0,2051) inpfile(1:namleni)
            write (6,2051) inpfile(1:namleni)
            call trnsfr(caref,cres(1,ia1),3)
          ## end if
          call setpbccell('Do you want to use periodic images',34,
     -      edge,edge_gen,cell,ncell,cellalt,ixyzhex,npbc,
     -      ioppbc,iusepbc,vol,nw,rinscr,rcirc,0)
          if (iusepbc  >=  0) :
            if (ifirstref  ==  1) :
              call pbcdist(caref,c(1,ia2),ia1,ia2,cell,ncell,-iw0,
     -          1,img,drincr,dimg)
              call zeroit(drimg,3)
              call arrdiff(drimg,cell(1,img),drimg,3)
            ## end if
          else:
            ncell=1
          ## end if
          xm=800.0
          ym=900.0
        elif (ianaltyp  ==  10) :
c         Set up check for unphysical features
          write (6,2019) 'List of suspicious contacts',
     -      analfile(1:namleno)
          call get#real('CTFAC',5,1.4,ctfac,1,62)
          call get#real('MINFAC',6,0.4,bondminfac,1,62)
          call getint('MAXDIST',7,min0(nslt,50),1,nslt,maxdist,62)
          isltonly=1
          if (n  >  2000) print *,'NOTE: This check is SLOW'
          if (nslt  <  n)
     -      call askyn('Do you want to include the solvents',
     -        35,0,-1,isltonly,0,0)
          if (nsegslt  >  3  and  ischarmm(inpcrdtyp)  ==  1) :
            call askyn('Do you have LES segments',24,1,-1,iles,0,0)
            if (iles  ==  1) :
              qq='Is segment      a duplicate'
              for is in range(0, nsegslt):
                qq(12:15)=segid4(is)
                call askyn(qq,27,1,-1,idelse:g(is),0,0)
              
            ## end if
           write (6,2114)
          ## end if
          call setmolres(ifres,ilres,isegno,molresflag,
     -      molsltlim,nrescol,irescol1,irescol2,resnames,nresslt,
     -      nmolslt,nsegslt,nmolsltnoion,minresflag,index,indexa,indexs,
     -      line,maxrsd,maxrec)
          call setpbccell('Do you want to use periodic images',34,
     -      edge,edge_gen,cell,ncell,cellalt,ixyzhex,npbc,
     -      ioppbc,iusepbc,vol,nw,rinscr,rcirc,0)
          write (iw0,2061) ctfac,bondminfac,maxdist
          call printbondthres(ialist,nanos,ctfac,bondminfac,
     -      iatnm2,ramax,iw0)
          if (iusepbc  >=  0) :
            write (iw0,*) 'Intermolecular distances include PBC images'
            call prtcell(ioppbc,edge,edge_gen,r,vol,nw,-iw0)
          ## end if
          if (nslt  >  1000) print *,'Wait ...'
        elif (ianaltyp  ==  11) :
c         Psudorot calc
          if (itraj  ==  1) :
            call getring(line,index,ix5,irescol1,irescol2,inamcol1,
     -        inamcol2,nslt,numres,iresring,iresno,ifres,ilres,nmem,
     -        psrtyp,0,incgen,iapex,MAXRING,maxrec)
            write (6,2019) 'Pseudorotation angles',analfile(1:namleno)
          ## end if
        elif (ianaltyp  ==  12) :
          write (6,2071) analfile(1:namleno)
          extnam1='.pkn'
          iprintpk=0
          irespro=0
          isegpkhx=1
          if (iresshift  >  0) write (6,2029)
          if (nsegm  >  1) :
9083        call getint('Segment number of helix',23,1,1,
     -        nsegm,isegpkhx,0)
            call findrange(isegno,1,nslt,isegpkhx,ifseghx,ilseghx,
     -        'segment',7,0,ifail)
            if (ifail  >  0) go to 9083
          else:
            ifseghx=1
            ilseghx=nslt
          ## end if
          while (irespro  ==  0)
            call getint('Residue number of the kink (Proline)',36,
     -        0,1,iresno(nslt),irespro,0)
            call findresnum(iresno,ixres,irespro,ifseghx,ilseghx,ia,
     -        irfound)
          
          ial=ia-1
          irespro=irfound
          resnam(1:nrescol)=line(index(ia))(irescol1:irescol2)
          call leftadjustn(resnam,rn,8)
          noprol=0
          if (rn(1:3)  !=  'PRO') :
            print *,'NOTE: residue ',irespro,
     -        ' is not a proline but a ',resnam(1:nrescol)
            noprol=1
            iflatproline=0
          else:
            call askyn('Do you want to project the proline to a plane',
     -        45,1,-1,iflatproline,96,0)
          ## end if
9082      call getint('Number of helix residues before the kink',40,
     -      7,1,MAXHX-1,nrb,0)
          call getint('Number of helix residues after  the kink',40,
     -      7,1,MAXHX-1,nra,0)
          if (nra  <  4  or  nrb  <  4) :
            print *,'Minimum helix length is four '
            go to 9082
          ## end if
          icaonly=-1
          call findprotbackbone(line,index,iresno,ia,inamcol1,
     -      icapr,icpr,inpr,icbpr,icgpr,icdpr,nslt,icaonly,maxrec)
          write (6,2158) 'Proline',0,
     -      line(index(ia-1))(irescol1:irescol2),iresno(ia-1),
     -      icapr,icpr,inpr
          icaa(1)=icapr
          ica(1)=icpr
          ina(1)=inpr
c         icab(1)=icapr
c         icb(1)=icpr
c         inb(1)=inpr
          for ir in range(0, nra):
            call findprotbackbone(line,index,iresno,ia,inamcol1,
     -        icaa(ir+1),ica(ir+1),ina(ir+1),icbx,icgx,icdx,
     -        nslt,icaonly,maxrec)
            write (6,2158) 'After  ',ir,
     -        line(index(ia-1))(irescol1:irescol2),iresno(ia-1),
     -        icaa(ir),ica(ir),ina(ir)
          
          for ir in range(0, nrb):
            ia=ial
            ires=iresno(ial)
            while (iresno(ia)  ==  ires  and  ia  >  1)
              ia=ia-1
            
            ial=ia
            ia=ia+1
            call findprotbackbone(line,index,iresno,ia,inamcol1,
     -        icab(ir),icb(ir),inb(ir),icbx,icgx,icdx,nslt,
     -        icaonly,maxrec)
            write (6,2158) 'Before ',ir,
     -        line(index(ia-1))(irescol1:irescol2),iresno(ia-1),
     -        icab(ir),icb(ir),inb(ir)
          
          print *,"Kink residue (proline) will be added only to the ",
     -      "helix 'after'"
          if (noprol  ==  0) :
            call getring(line,index,ix5,irescol1,irescol2,inamcol1,
     -        inamcol2,nslt,numres,irespro,iresno,ifres,ilres,nmem,
     -        'p',1,incgen,iapex,MAXRING,maxrec)
            if (nmem  ==  0)
     -        print *,'Pseudorotation angle calculation is canceled'
          ## end if
          write (iw0,2077)
          write (iw0,2040) irespro-nrb,irespro+nra,isegpkhx
          write (resrange(1:22),2087) irespro-nrb,irespro+nra
          write (iw0,2083) irespro
          if (itraj  ==  1) :
c           Open Prokink dial window, initialize trajectory accumulators
            ndials=5
            if (nmem  ==  0) ndials=ndials-1
          ## end if
          nomat0=1
        elif (ianaltyp  ==  21) :
          itmem=0
9194      call getint('Number of helices to analyze',28,1,1,MAXNHX,nhx,
     -      0)
          nresrec=nhxres*nhx+3*nhx*(nhx-1)/2
          npspages=nresrec
          if (nresrec  >  MAXCOPY1) :
            write (6,2161) nresrec,MAXCOPY1
            go to 9194
          ## end if
          for ihx in range(0, nhx):
c           Helix axis directions
            if (nhx  >  1) write (6,2159) ihx
            if (nsegm  >  1) :
9183          call getint('Segment number of helix',23,1,1,
     -          nsegm,iseghx(ihx),0)
              call findrange(isegno,1,nslt,iseghx(ihx),ifseghx,ilseghx,
     -          'segment',7,0,ifail)
              if (ifail  >  0) go to 9183
            else:
              ifseghx=1
              ilseghx=nslt
            ## end if
            if (iresshift  >  0) write (6,2029)
9182        call getrange(ireshx1(ihx),999999,ireshx2(ihx),999999,incr,
     -        0,'helix residue',13,iresno(ilseghx),0)
            if (ireshx1(ihx)  <  iresno(ifseghx)) :
              write (6,2123) 'less',iresno(ifseghx)
              go to 9182
            ## end if
            if (ireshx2(ihx)  >  iresno(ilseghx)) :
              write (6,2123) 'greater',iresno(ifseghx)
              go to 9182
            ## end if
            nreshx(ihx)=ireshx2(ihx)-ireshx1(ihx)+1
            if (nreshx(ihx)  <  4) :
              print *,'Minimum helix length is four '
              go to 9182
            elif (nreshx(ihx)  >  MAXHX-1) :
              print *,'Maximum helix length is ',MAXHX-1
              go to 9182
            ## end if
            print *
            icaonly=0
            for ir in range(ireshx1(ihx), ireshx2(ihx)):
              ir0=ir
              call findresnum(iresno,ixres,ir0,ifseghx,ilseghx,ia,
     -          irfound)
              if (ir0  ==  0) go to 9182
              ir1=ir-ireshx1(ihx)+1
              ia=ifres(irfound)
              call findprotbackbone(line,index,iresno,ia,
     -          inamcol1,icaahx(ir1,ihx),ica(ir1),ina(ir1),icb(ir1),
     -          icgx,icdx,nslt,icaonly,maxrec)
              write (6,2158) ' ',ir,
     -          line(index(icaahx(ir1,ihx)))(irescol1:irescol2),
     -          iresno(icaahx(ir1,ihx)),icaahx(ir1,ihx),ica(ir1),
     -          ina(ir1)
            
            if (nhx  >  1) :
9186          call getint(
     -          'Residue number (rn) that defines the helix rotation',
     -          51,0,1,ireshx2(ihx),iresrot,143)
              if (iresrot  >  0) :
                if (iresrot  <  ireshx1(ihx)) :
                  print *,'Residue number is outside the ',ireshx1(ihx),
     -              ' - ',ireshx2(ihx),'range'
                  go to 9186
                ## end if
                ir1=iresrot-ireshx1(ihx)+1
                icbahx(ihx)=icb(ir1)
                icbreshx(ihx)=ir1
                if (icb(ir1) .eq .0) :
                  print *,'Residue ',iresrot,' has not CB atom'
                else:
                  print *,'Atom index of the CB defining the rotation=',
     -              icbahx(ihx)
                ## end if
              ## end if
            ## end if
9184        call getint(
     -        'Number of residues to ignore for rotation at each # end',
     -        53,0,1,0,incrot,0)
            if (ireshx2(ihx)-ireshx1(ihx)+1-2*incrot  <  4) :
              print *,'ERROR: too few residues remain to form a helix'
              go to 9184
            ## end if
          
          call askyn('Is this a transmembrane protein',31,1,-1,itmem,0,
     -      0)
          if (itmem  ==  1) call getint(
     -      'Coordinate axis normal to the membrane (1/2/3)',46,3,1,3,
     -      normhx,0)
c         call getint('First helix of helix pair to track',34,1,1,nhx,
c    -      ihxt1,0)
c         call getint('Secnd helix of helix pair to track',34,ihxt1+1,1,
c    -      nhx,ihxt2,0)
          ndials=6
          idebughx=0
          call askyn('Do you want debug output',24,1,-1,idebughx,0,0)
          if (idebughx  >  0) :
            call askyn('Do you want to skip b# end analysis',33,1,-1,ibx,
     -        0,0)
            if (ibx  >  0) idebughx=idebughx+1
          ## end if
          if (idebughx  <  2) call get#real(
     -      'Minimum distance (in A) from axis to count as b# end',50,
     -      0.0,axtol,1,0)
          write (iw0,2004) axtol
          if (nhx  >  0) write (iw0,2124)
          write (iw0,2007)
          if (itraj  ==  1) :
            linein(1:45)='Do you want to subtract the solute COM shift '
            linein(46:65)='from the helix shift'
            call askyn(linein,65,1,1,isubcrm,0,0)
            linein(1:41)='Do you want to overlay each frame on the '
            linein(42:60)='reference structure'
            call askyn(linein,60,1,1,ioverlay,0,41)
            if (ioverlay  >  0  and  nsegslt  >  1) :
              linein(1:43)='Do you want to overlay each solute molecule'
              linein(44:54)=' separately'
              call askyn(linein,54,1,1,ioverlaym,0,0)
              if (ioverlaym  >  0) ioverlay=ioverlay+1
            ## end if
            if (idebughx  >  0) :
              linein(1:40)='Do you want to reorient each helix onto '
              linein(41:59)='the reference helix'
              call askyn(linein,59,1,1,ireorienthx,0,0)
            else:
              ireorienthx=1
            ## end if
            idsspcheck=0
            if (icaonly  <  1)
     -        call askyn('Do you want to run DSSP check on each frame',
     -          43,1,1,idsspcheck,0,0)
            call zeroiti(nhelixok,0,3)
            if (idsspcheck  ==  1)
     -        write (iw0,2016) (typc(i),ssname(i)(1:lssname(i)),i=1,9)
          ## end if
          for ihx in range(0, nhx):
            write (iw0,2015) ihx,ireshx1(ihx),ireshx2(ihx),iseghx(ihx),
     -        icbreshx(ihx)
          
        elif (ianaltyp  ==  13) :
c         Set up hydropathy labeling
          write (6,2076) 'Hydropathy scale',analfile(1:namleno)
          ihydtyp=0
        elif (ianaltyp  ==  20) :
          innlist=0
c         Set up Delphi potential labeling
          nl2=0
          call openfile(50,0,'Delphi potential map file',25,'old',
     -      analfile2,nl2,notfnd,0,2,1,0,0)
          call readmap(50,xstart,ystart,zstart,gx,gy,gz,ngx,ngy,ngz,
     -      nconf,c2,maxrec)
          igincr=0
          while (igincr  <  1)
            call getint('Increment to use in reading the potential map',
     -         45,1,1,ngx,igincr,47)
            if (igincr  <  1) print *,'Invalid'
          
          if (igincr  >  1) write (6,2026) igincr
          write (6,2076) 'Delphi potentials',analfile(1:namleno)
          call askyn(
     -      'Do you want to exclude/interpolate grids near atoms',51,
     -      1,-1,iexcl,63,0)
          if (iexcl  >  0) :
            call get#real('Minimum distance between an atom and a grid',
     -        43,1.0,rnear,1,63)
            call askyn(
     -        'Do you want to interpolate the excluded grid values',51,
     -        1,-1,interpol,0,0)
          ## end if
          call askyn('Do you want to query the potential map',38,1,-1,
     -      iquery,64,0)
          iw2=0
          iw3=0
          if (igridfile+iquery  ==  0) :
            print *,'No action requested'
          else:
            if (igridfile  >  0) :
              write (6,2076) 'Delphi grid potentials',
     -          analfile1(1:namleno1)
              write (iw1,2027) 'Delphi grid potentials'
c             Open 2nd grid PDB file
              iw2=iw1+1
              analfile2=analfile1
              analfile2(namleno1:namleno1)='0'
              write (6,2076)
     -          'Delphi grid potentials without excluded grids',
     -          analfile2(1:namleno1)
              call openfile(iw2,0,'grid',4,'new',analfile2,namleno1,
     -          notfnd,0,1,1,0,0)
              write (iw2,2027)
     -          'Delphi grid potentials without excluded grids'
              if (interpol  >  0) :
                iw3=iw2+1
                analfile3=analfile1
                analfile3(namleno1:namleno1)='1'
                write (6,2076)
     -          'Delphi grid potentials, interpolating excluded grids',
     -            analfile3(1:namleno1)
                call openfile(iw3,0,'grid',4,'new',analfile3,namleno1,
     -            notfnd,0,1,1,0,0)
                write (iw3,2027)
     -            'Delphi grid potentials, interpolating excluded grids'
              ## end if
            ## end if
            call delphigrid(iw1,iw2,iw3,c,n,nslt,xstart,ystart,zstart,
     -        gx,gy,gz,ngx,ngy,ngz,igincr,rnear,igridfile,iexcl,
     -        iquery,interpol)
            if (iw1  >  0) close (iw1)
            if (iw2  >  0) close (iw2)
            if (iw3  >  0) close (iw3)
            if (iexcl  >  0) call readmap(50,xstart,ystart,zstart,
     -        gx,gy,gz,ngx,ngy,ngz,nconf,c2,maxrec)
          ## end if
        elif (ianaltyp  ==  14) :
c         Set up circular variance labeling calculation
          write (6,2076) 'Circular variances',analfile(1:namleno)
          nsltref_f=1
          nsltref_f=nslt
          call getrange(nsltref_f,1,nsltref_l,nslt,incr,0,
     -      'solute atom to include in the CV calculation',44,nslt,000)
          call get#real(
     -      'Distance cutoff for the circular variance calculation',
     -      53,10.0,rcut_cv,1,66)
          islvrep=1
          if (n  >  nslt) :
            if (naslv  >  1)
     -        call getint('Index of the representative solvent atom',
     -          40,1,1,naslv,islvrep,25)
            call askyn('Do you want to sort solvents by CV',34,1,+1,
     -        isortslv,0,0)
            if (isortslv  >  0) :
              write (6,2105)
              call get#real('CV threshold to count the # of solvents',
     -          39,0.5,cvlim,1,67)
            ## end if
          ## end if
        elif (ianaltyp  ==  15) :
c         Set up circular variance map calculation
          write (6,2075) 'Circular variance',analfile(1:namleno),
     -      analfile1(1:namleno1)
          xm=500.0
          ym=800.0
          call openps(iw1,xm,ym-20,title,ltitle76,
     -      'Circular variance map',21,inpfile,namleni,inpfile,0,1,
     -      ipspage)
          iaskcolcode=1
          npspages=1
        elif (ianaltyp  ==  16) :
c         Set up DSSP calculation
          write (6,2019)
     -      'Secondary structure assignment according to DSSP',
     -      analfile(1:namleno)
          if (itraj  ==  1) write (6,2019)
     -       'Secondary structure plot',analfile1(1:namleno1)
          write (iw0,2003)
          if (itraj  ==  1) :
            call getrange(ifrdssp,1,ilrdssp,nresslt,incr,0,
     -        'residue sequence number to plot',31,nresslt,0)
            if (iresshift  >  0) :
              write (6,2117) 'first',ixresno(ifrdssp)
              write (6,2117) 'last ',ixresno(ilrdssp)
              call askyn(
     -          'Do you want to use the actual residue # on the Y axis',
     -          53,1,-1,iuseinp,0,0)
            ## end if
            call askyn(
     -        'Do you want to calculate turn information (see GW Rose)',
     -        55,1,-1,irose,0,0)
            iwrose=iw0*irose
          else:
            ifrdssp=1
            ilrdssp=nresslt
            iwrose=0
          ## end if
          if (iuseinp  ==  1) :
            call trnsfi(indexdel,ixresno,nresslt)
            for i in range(nresslt+1, maxrec):
              indexdel(i)=indexdel(nresslt)+(i-nresslt)
            
          ## end if
          if (iuseinp  ==  0) call indexit(indexdel,1,nresslt,0)
          call zeroiti(idistdssp,0,9*maxrsd)
          ncolcode=8
          innlist=0
        elif (ianaltyp  ==  22  or  ianaltyp  ==  23  or 
     -           ianaltyp  ==  24) :
c         Set up RMSD calculations
          iaskcolcode=1
          call zeroiti(indexdel,0,n)
          call indexit(indexov,1,nslt,0)
          nfinalov=nslt
          if (ianaltyp  ==  23  or  ianaltyp  ==  24) :
            absdevmin=100000.0
            absdevmax=0.0
            innlist=0
            xm_2d=600.0
            ym_2d=775.0

            xm=xm_2d
            ym=ym_2d
          ## end if
          noopt2d=0
          if (iread2d  ==  0) :
            if (ianaltyp  ==  23  or  ianaltyp  ==  24) :
              write (iw0,2084)
              call askyn(
     -          'Do you want to superimpose first the frames',43,
     -          0,+1,noopt2d,0,0)
              if (noopt2d  ==  0) write (iw0,2116)
     -          'after obtaining the best fit with the Kabsch method'
              if (noopt2d  ==  1) write (iw0,2116)
     -           'without superimposition'
              write (iw0,*)
              rmsdmin=0.0
              call get#real(
     -          'MAXimum of the RMSD scale (default: actual maximum)',
     -          51,0.0,rmsdmax,1,53)
            ## end if
            if (noopt2d  ==  0) :
              call askyn('Do you want to select atoms for overlay',39,
     -          1,-1,iedit,69,0)
              if (iedit  >  0) :
                nrecdel=0
                call select(line,nrecdel,idcol,asterisk,n,nslt,index,
     -            ixres,is1,is2,iseqncol1,iseqncol2,inamcol1,inamcol2,
     -            irescol1,irescol2,iqcol1,iqcol2,charge,iatnum,
     -            nneig,ineig,indexdel,iw0,maxng,maxrec)
c               indexdel: 0 for atoms to use for the overlay, 1 for the rest
                call masktolist(indexov,indexdel,nslt,nfinalov,0)
              ## end if
            ## end if
c           indexov contains the list of atoms used for overlay
            iqincr=0
            if (nfinalov  ==  nslt) iqincr=1
            call quiz(ans,icalctyp,'a',' ',0,
     -        'atoms for RMSD calculation',26,iqincr,5,6,0)
            if (ans  ==  'a') :
              call indexit(indexrmsd,1,nslt,0)
              nfinalrmsd=nslt
            elif (ans  ==  'o') :
              call trnsfi(indexrmsd,indexov,nfinalov)
              nfinalrmsd=nfinalov
            else:
              nrecdel=0
              call zeroiti(indexdel,0,n)
              call select(line,nrecdel,idcol,asterisk,n,nslt,index,
     -          ixres,is1,is2,iseqncol1,iseqncol2,inamcol1,inamcol2,
     -          irescol1,irescol2,iqcol1,iqcol2,charge,iatnum,nneig,
     -          ineig,indexdel,iw0,maxng,maxrec)
              call masktolist(indexrmsd,indexdel,nslt,nfinalrmsd,0)
              ifstrmsd=indexrmsd(1)
              ilstrmsd=indexrmsd(nfinalrmsd)
              if (ilstrmsd-ifstrmsd+1  !=  nfinalrmsd) :
                nocontigrmsd=1
                write (6,2132)
                write (iw0,2132)
              ## end if
            ## end if
            if (nfinalrmsd  <  nslt) limresrange=1
c           indexrmsd contains the list of atoms used for overlay
            if (nfinalov  <  nslt) :
              write (6,2153) 'overlay',nfinalov,nslt
              write (iw0,2153) 'overlay',nfinalov,nslt,';',
     -          (indexov(i),i=1,nfinalov)
            elif (noopt2d  ==  0) :
              write (6,2154) 'overlay'
              write (iw0,2154) 'overlay'
            ## end if
            if (nfinalrmsd  <  nslt) :
              write (6,2153) 'RMSD',nfinalrmsd,nslt
              write (iw0,2153) 'RMSD',nfinalrmsd,nslt,';',
     -          (indexrmsd(i),i=1,nfinalrmsd)
            else:
              write (6,2154) 'RMSD'
              write (iw0,2154) 'RMSD'
            ## end if
            atwsum=0.0
            for iaa in range(0, nfinalov):
              ia=indexov(iaa)
              atw(ia)=aw(iatnum(ia))
              atwsum=atwsum+atw(ia)
              
          ## end if
          if (ianaltyp  ==  22) :
            write (6,2098)
            call askyn('Do you want to use the input structure instead',
     -        46,1,-1,inputref,0,0)
            call askyn('Do you want to plot also the maximum deviation',
     -        46,1,-1,maxdevplot,0,0)
            call askyn('Do you want to plot without overlay',35,1,-1,
     -        iunmatchplot,0,0)
            write (iw0,2082)
            npspages=1+iunmatchplot
            if (inptrajtyp  ==  3) npspages=npspages+1
            call get#real(
     -        'MAXimum of the RMSD scale (default: actual maximum)',
     -        51,0.0,rmsdplotmax,1,0)
            call get#real(
     -        'MAXimum of the max dev scale (default: actual maximum)',
     -        54,0.0,rmaxdevplotmax,1,0)
            if (nresslt  <=  MAXDISTR) :
              call zeroitd(rmsfsum,nresslt)
              call get#real(
     -          'MAXimum of the RMSF scale (default: actual maximum)',
     -          51,0.0,rmsfplotmax,1,0)
              call getint('Frequency of plotting error bars',32,
     -          max0(1,nresslt/10),1,max0(1,nresslt/2),nfreqsd,0)
            else:
              write (6,2062) nresslt,MAXDISTR
              npdpages=npspages-1
            ## end if
            call zeroitd(cdp,3*nfinalrmsd)
            call zeroitd(cdp2,3*nfinalrmsd)
          ## end if
          if (ianaltyp  ==  23) :
c           2D RMSD
            call askyn(
     -        'Do you want to also plot the RMSD distributions',47,
     -        0,+1,noplotdist,93,0)
            ireplot=0
            pstitle='2-D RMSD map'
            lpstitle=12
            if (iread2d  ==  0) :
              npspages=3
              call askyn('Do you want RMSD-based clustering',33,
     -          1,+1,irmsdclust,0,0)
              call askyn('Do you want to skip plotting the matrix',39,
     -          1,-1,iskip2dplot,121,0)
              if (irmsdclust  >  0)
     -          call askyn(
     -            'Do you want an RMSD map rearranged by clusters',46,
     -            1,1,ireplot,0,0)
            else:
              irmsdclust=0
c             Read 2-d RMSD map, run clustering
              lq=32+lanalfile4
              question(1:lq)='Input rd2 file name: '
     -          //analfile4(1:lanalfile4)//' - is it OK'
              call askyn(question,lq,1,1,ians,00,0)
              if (ians  ==  0) lanalfile4=0
              call openfile(iw4,0,'previously written .rd2',23,
     -          'old',analfile4,lanalfile4,notfnd,3,1,1,0,0)
              extnam2='    '
              write (6,2019) 'Clustering results',analfile(1:namleno)
              write (iw0,2064) analfile4(1:namlen4)
              call read_2drmsd(iw4,system,lsystem,trajnam,ltrajnam,
     -          trajnam2,ltrajnam2,ifirsttraj,ilasttraj,incrementtraj,
     -          ifirsttraj2,ilasttraj2,incrementtraj2,1,rmsdmn,rmsdmx,
     -          nframe,nframex,nframey,ietotsaved,0,limresrange,ierr)
              if (ierr  >  0) go to 9003
              nframeref=nframe
              close(iw4)
c             Run clustering
              call indexit(iconfsel,1,MAX2D,0)
              call indexit(ixclst,1,MAX2D,0)
              call askyn(
     -          'Do you want an RMSD map rearranged by clusters',46,
     -          1,1,ireplot,0,0)
              call clusterdistr(nframe,iw0,rmsdlim,rmsdmn,rmsdmx,
     -          nhbdist,it1,it2,it3,itemp4,indexn,indexo,ncl,indexa,
     -          iconfsel,ixclst,it4,value,ifa_s,ila_s,ih,cv,0.0,rdclust,
     -          res(1,1,11),ietotsaved,'RMSD',4,1,1,irepav,irepmx,
     -          irepeng,irepkm,engcl,c,chn,c2,1,27,iclstyp,iwt,0,
     -          label2d,80,0,1,1,1,mx2d,maxframe)
              for i in range(0, ncl):
                nclstmem(i)=indexo(i)-indexn(i)+1
              
              call trnsfi(ixshuffle,iconfsel,MAX2D)
c             iconfsel contains the sorted cluster members
c             indexn, indexo contain the cluster limits
              call countsim(indexn,indexo,iconfsel,ncl,rdclust,rmsdsim,
     -          nsimclst1,iw0,mx2d)
              call askyn('Do you want to replot the input RMSD map',40,
     -          1,-1,ireplotinp,0,0)
              if (ireplot  ==  0) call askyn(
     -          'Do you want an RMSD map rearranged by clusters',46,
     -          1,1,ireplot,0,0)
              if (ireplot+ireplotinp  >  0) :
                call getint(askcolcode,24,maxcolcode,1,maxcolcode,
     -            ncolcode,95)
c               Generate plot file name
                npspages=2
                analfile2=analfile
                analfile2(namleni-2:namleni+6)='rd2_sh.ps'
                namleno2=namleni+6
                write (6,2019)
     -            pstitle(1:lpstitle)//' and membership plot',
     -            analfile2(1:namleno2)
                iw1=iw0+1
                call openfile(iw1,0,'(shuffled) RMSD map',19,
     -            'new',analfile2,namleno2,notfnd,0,1,1,0,0)
                call openps(iw1,xm_2d,ym_2d,title(1:76),76,pstitle,
     -            lpstitle,' ',0,analfile4(1:namlen4),
     -            namelen4,npspages,ipspage)
                if (ireplot  ==  1) :
                  call indexit(it3,1,MAX2d,0)
                  call adjust_xtraj(xtraj,ifirsttraj,ilasttraj,
     -              incrementtraj,framefac,iadjust_xtraj)
                  call plot2drmsd(nrep,iw0,iw1,xtraj,maxrec,title,
     -              'Frames sorted by clusters',25,0,xtrajlab,11,
     -              ncolcode,maxcolcode,iedit,noopt2d,limresrange,0,
     -              rmsdmin,rmsdmax,absdevmin,absdevmax,rmsdmn,rmsdmx,
     -              indexa,it3,ixshuffle,ixshuffle,indexo,ncl,indexo,
     -              ncl,ym_2d,it1,it2,temp,1,1,noplotdist,0,1,1,
     -              ipspage,0)
                ## end if
                if (ncl  >  0)
     -            call clusterplot(iw1,xtraj,value,indexn,indexo,ncl,
     -              ixclst,nframe,xtrajlab,11,ipspage,ireplotinp,mx2d)
                call indexit(ixshuffle,1,MAX2D,0)
                if (ireplotinp  ==  1) :
                  call indexit(it3,1,MAX2d,0)
                  call plothead(iw1,xm_2d,ym_2d,title(1:76),76,
     -              'RMSD file read:'//analfile4(1:namlen4),
     -              namlen4+15,trajnam,ltrajnam,'',0)
                  call plot2drmsd(nrep,iw0,iw1,xtraj,maxrec,title,
     -              'Input RMSD matrix',17,0,xtrajlab,11,ncolcode,
     -              maxcolcode,iedit,
     -              noopt2d,limresrange,0,rmsdmin,rmsdmax,absdevmin,
     -              absdevmax,rmsdmn,rmsdmx,indexa,it3,ixshuffle,
     -              ixshuffle,indexo,0,indexo,0,ym_2d,it1,it2,temp,1,1,
     -              noplotdist,0,1,1,ipspage,0)
                ## end if
              ## end if
              go to 9005
            ## end if
            call indexit(iconfsel,1,MAX2D,0)
          elif (ianaltyp  ==  24) :
c           Cross RMSD
            call get#real(
     -        'MINimum of the RMSD scale (default: actual minimum)',
     -        51,0.0,rmsdmin,1,53)
            if (iread2d  ==  0) :
              call askyn('Do you want to find matching structures',39,
     -          1,-1,matchconf,0,0)
              call get#real('MAXimum RMSD for similarity statistics',38,
     -          2.5,rmsdsim,1,86)
              if (rmsdsim  >  0) matchconf=1
              pstitle='Cross RMSD map'
              lpstitle=14
              npspages=2
              call indexit(iconfsel,1,MAX2D,0)
            else:
c             Read RMSD maps and try to match them
              notfnd=1
              while (notfnd  >  0)
                call getname(analfile4,namlen4,
     -            'Name of .rd2 file from the first trajectory',43,
     -             200,'',0,0,0,0)
                call openfile(iw4,0,'previously written .rd2',23,
     -            'old',analfile4,namlen4,notfnd,0,1,1,0,0)
              
              write (iw0,2064) analfile4(1:namlen4)
              call read_2drmsd(iw4,system,lsystem,trajnam1,ltrajnam1,
     -          trajnam1,ltrajnam2x,if,il,inc,if2,il2,inc2,1,rmsdmn,
     -          rmsdmx,nframe1,nframex,nframey,ietotsaved,0,limresrange,
     -          ierr)
              close (iw4)
c             Run clustering
              write (6,2019) 'Clustering results',analfile(1:namleno)
              call indexit(index2d1,1,MAX2D,0)
              call indexit(ixclst,1,MAX2D,0)
              call clusterdistr(nframe1,iw0,rmsdlim,rmsdmn,rmsdmx,
     -          nhbdist,it1,it2,it3,itemp4,ifclst1,ilclst1,ncl1,indexa,
     -          index2d1,ixclst,it4,value,ifa_s,ila_s,ih,cv,0.0,
     -          rdclust,res(1,1,11),ietotsaved,'RMSD',4,1,1,irepav,
     -          irepmx,irepeng,irepkm,engcl,c,chn,c2,1,27,iclstyp,iwt,0,
     -          label2d,80,0,1,1,1,mx2d,maxframe)
              call trnsfi(irepmx1,irepmx,MAX2D)
              call trnsfi(ixshuffle,index2d1,MAX2D)
              call countsim(ifclst1,ilclst1,index2d1,ncl1,rdclust,
     -          rmsdsim1,nsimclst1,iw0,mx2d)
              namlen4=0
              notfnd=1
              while (notfnd  >  0)
                call getname(analfile4,namlen4,
     -            'Name of .rd2 file from the second trajectory',44,
     -             200,'',0,0,0,0)
                call openfile(iw4,0,'previously written .rd2',23,
     -            'old',analfile4,namlen4,notfnd,0,1,1,0,0)
              
              write (iw0,2064) analfile4(1:namlen4)
              call read_2drmsd(iw4,system,lsystem,trajnam2,ltrajnam2,
     -          trajnam1,ltrajnam2x,if,il,inc,if2,il2,inc2,1,rmsdmn,
     -          rmsdmx,nframe2,nframex,nframey,ietotsaved2,0,
     -          limresrange,ierr)
              close (iw4)
c             Run clustering
              call indexit(index2d2,1,MAX2D,0)
              call indexit(ixclst,1,MAX2D,0)
              call clusterdistr(nframe2,iw0,rmsdlim,rmsdmn,rmsdmx,
     -          nhbdist,it1,it2,it3,itemp4,ifclst2,ilclst2,ncl2,indexa,
     -          index2d2,ixclst,it4,value,ifa_s,ila_s,ih,cv,0.0,
     -          rdclust,res(1,1,11),ietotsaved,'RMSD',4,1,1,irepav,
     -          irepmx,irepeng,irepkm,engcl,c,chn,c2,1,27,iclstyp,iwt,0,
     -          label2d,80,0,1,1,1,mx2d,maxframe)
              call trnsfi(irepmx2,irepmx,MAX2D)
              call trnsfi(ixshuffleref,index2d2,MAX2D)
              call countsim(ifclst2,ilclst2,index2d2,ncl2,rmsdsim1,
     -          rmsdsim2,nsimclst2,iw0,mx2d)
              namlen4=0
              notfnd=1
              while (notfnd  >  0)
                call getname(analfile4,namlen4,
     -            'Name of the previously written .rdx file',40,200,
     -             '',0,0,0,0)
                analfile=inpfile
                analfile(namleni+1:namleni+4)='.rdx'
                call openfile(iw4,0,'previously written .rdx',23,
     -            'old',analfile4,namlen4,notfnd,0,1,1,0,0)
                extnam2='    '
              
              call read_2drmsd(iw4,system,lsystem,trajnamr1,ltrajnamr1,
     -          trajnamr2,ltrajnamr2,if,il,inc,if2,il2,inc2,0,rmsdmn,
     -          rmsdmx,nframe2x,nframex,nframey,ietotsaved,0,
     -          limresrange,ierr)
              nframe=nframex
              nframeref=nframey
              write (iw0,2064) analfile4(1:namlen4)
              nerr=0
              if (trajnamr1(1:ltrajnamr1)  !=  trajnam1(1:ltrajnam1))
     -            :
                write (6,2063) 'first',trajnamr1(1:ltrajnamr1),
     -            trajnam1(1:ltrajnam1)
                nerr=nerr+1
              ## end if
              if (trajnamr2(1:ltrajnamr2)  !=  trajnam2(1:ltrajnam2))
     -            :
                write (6,2063) 'second',trajnamr2(1:ltrajnamr2),
     -            trajnam2(1:ltrajnam2)
                nerr=nerr+1
              ## end if
              if (nframex  !=  nframe1  or  nframey  !=  nframe2) :
                write (6,2106) nframex,nframey,nframe1,nframe2
                nerr=nerr+1
              ## end if
              if (nerr  >  0) call askstop(1)
              call askyn(
     -          'Do you want a cross-RMSD map rearranged by clusters',
     -            51,1,1,ireplot,0,0)
              if (ireplot  ==  1) :
c               Generate plot file name
                call getint(askcolcode,24,maxcolcode,1,maxcolcode,
     -            ncolcode,95)
                pstitle='Cluster-ordered cross RMSD map'
                lpstitle=30
                npspages=2
                analfile2=analfile
                analfile2(namleni-2:namleni+6)='rdx_sh.ps'
                namleno2=namleni+6
                write (6,2019) pstitle(1:lpstitle),analfile2(1:namleno2)
                iw1=iw0+1
                call openfile(iw1,0,'shuffled cross-RMSD map',23,
     -            'new',analfile2,namleno2,notfnd,0,1,1,0,0)
                call openps(iw1,xm_2d,ym_2d,title(1:76),76,pstitle,
     -            lpstitle,' ',2,' ',2,npspages,ipspage)
                call adjust_xtraj(xtraj,ifirst,ilast,increment,
     -            framefac,iadjust_xtraj)
                call plot2drmsd(nrep,iw0,iw1,xtraj,maxrec,title,
     -            'Frames sorted by clusters',25,ltrajnamr2,xtrajlab,11,
     -            ncolcode,maxcolcode,iedit,noopt2d,limresrange,0,
     -            rmsdmin,rmsdmax,absdevmin,absdevmax,rmsdmn,rmsdmx,
     -            indexa,iconfsel,ixshuffle,ixshuffleref,ilclst1,ncl1,
     -            ilclst2,ncl2,ym_2d,it1,it2,temp,1,0,noplotdist,0,1,
     -            1,ipspage,0)
                call indexit(ixshuffle,1,MAX2D,0)
              ## end if
              call countsimx(ifclst1,ilclst1,index2d1,ncl1,
     -          ifclst2,ilclst2,index2d2,ncl2,rmsdsim1,iw0,mx2d)
              write (iw0,*)
              call get#real('Maximum RMSD for mapping',24,rmsdsim1,
     -          rmsdmapmax,1,87)
              call mapclustx(2,ifclst1,ilclst1,irepmx1,ncl1,
     -          ifclst2,ilclst2,index2d2,ncl2,it1,it2,
     -          rmsdmapmax,nframe2,trajnam1,ltrajnam1,
     -          trajnam2,ltrajnam2,iw0,mx2d)
              write (iw0,*)
              call mapclustx(1,ifclst2,ilclst2,irepmx2,ncl2,
     -          ifclst1,ilclst1,index2d1,ncl1,it1,it2,
     -          rmsdmapmax,nframe1,trajnam1,ltrajnam1,
     -          trajnam2,ltrajnam2,iw0,mx2d)
              go to 9005
            ## end if
          ## end if
          if (iread2d  ==  0  and  ianaltyp  >  22) write (6,2075)
     -        'RMSD',analfile(1:namleno),analfile1(1:namleno1)
        elif (ianaltyp  ==  25) :
c         Residue correlation matrix calculation
9185      call getrange(ifrcorr,1,ilrcorr,ixres(nslt),incrres,1,
     -      'residue to calculate correlation and covariance map',51,
     -       numres,0)
          if ((ilrcorr-ifrcorr+1)/incrres+1  >  maxrcorr) :
            write (6,2096) maxrcorr,maxrcorr
            go to 9185
          ## end if
          iucorrmat=0
          call askyn(
     -      'Do you also want to calculate eigenvalues, eigenvectors',
     -      55,1,-1,ieig,125,0)
          if (ieig  ==  1) iucorrmat=46
          call askyn('Do you also want to plot the covariance matrix',
     -      46,1,-1,icovmatplot,0,0)
          call modrepats
          ncorr=0
          for ir in range(ifrcorr, ilrcorr, incrres):
            ncorr=ncorr+1
            indexs(ncorr)=ir
            call findat(indexa(ncorr),ifres(ir),ilres(ir),line,index,
     -        irescol1,inamcol1,MAXREC)
          
          ncolcode=0
          pstitle=' '
          lpstitle=0
          npspages=2
          ym=800.0
          write (iw0,*) 'Residue covariance and correlation matrix ',
     -     'calculation over a trajectory'
        elif (ianaltyp  ==  35) :
c         Input covariance matrix to get eigenvectors, eigenvalues
          call askyn('Is the covariance matrix file binary',36,1,-1,
     -      ibin,000,0)
          if (ibin  ==  0) :
            call askyn('Is the matrix broken into 5 numbers/line',40,
     -        1,-1,i5num,000,0)
            iformcov=2+i5num
          else:
            iformcov=1
          ## end if
          iw0=40
          iw1=iw0+1
          notfound=1
          while (notfound  >  0)
            call openfile(iw1,0,'input',5,'old',analfile1,lanalfile1,
     -        notfound,0,ibin+1,1,0,0)
            if (notfound  ==  1) lanalfile1=0
          
          call openfile(iw0,0,'output',6,'new',analfile,lanalfile,
     -      notfoundout,0,1,1,0,0)
          if (notfoundout  >  0) return
          call askyn('Do you want annotated output',28,1,1,iannout,128,
     -      0)
          call normalmodes(ncorr,iw1,0,iformcov,iw0,
     -      iannout,ierr,index2d,value,ifa_s,ila_s,it1,rmsdlim,MAXBONDS)
          if (ierr  >  0) return
          go to 9005
        elif (ianaltyp  ==  26) :
c         Atom-atom distance distribution calculation
          write (iw0,2028)
          call askyn(
     -      'Do you want to calculate distances between atom clusters',
     -      56,1,-1,iclusterdist,105,0)
          if (iclusterdist  ==  0) :
            write (6,2028) ' Enter atom pair indices'
            call getlist(listpairdist,npairs,1,nslt,2,MAXDDISTR)
          else:
            call getclusterpairs(npairs,iclustermem,ifstclst1,ifstclst2,
     -        ilstclst2,nslt,0,MAXDDISTR,MAXCDLIST)
          ## end if
          call zeroit(pairdistsum,npairs)
          call zeroit(pairdistwsum,2*npairs)
          call zeroit(pairdistsum2,npairs)
          call zeroiti(npairdist,0,10*npairs)
          call get#real('Maximum distance to calculate distribution',42,
     -      10.0,rmaxpair,1,65)
          pairgrid=rmaxpair/float(MAXDDBIN)
          for ip in range(0, npairs):
            pairdistminmax(1,ip)=10000.0
            pairdistminmax(2,ip)=-pairdistminmax(1,ip)
          
        elif (ianaltyp  ==  39) :
c         Atom-atom distance SD calculation
          write (iw0,2119)
          write (6,2019) 'Atom-atom distances and SDs',
     -      analfile(1:namleno)
          call select(line,nrecdel,idcol,asterisk,n,nslt,index,
     -      ixres,is1,is2,iseqncol1,iseqncol2,inamcol1,inamcol2,
     -      irescol1,irescol2,iqcol1,iqcol2,charge,iatnum,
     -      nneig,ineig,indexdel,6,maxng,maxrec)
c           indexdel: 0 for atoms to use, 1 for the rest
          call masktolist(ianchor,indexdel,nslt,nanchor,0)
          call quiz(ans,isdtyp,'n',' ',0,'SD normalization',16,0,5,6,0)
          xm_2d=600.0
          ym_2d=775.0
          pstitle='Atom-Atom distance SD plot'
          lpstitle=26
          call openps(iw1,xm_2d,ym_2d,title(1:76),76,pstitle,
     -      lpstitle,trajnam,ltrajnam,' ',0,npspages,ipspage)
        elif (ianaltyp  ==  27) :
c         Solvation shell volume calculation
          call get#real('Solvent radius',14,1.4,rsolv,1,0)
          nranddef=20000*(rsolv/1.4)**2
          call getint('Number of random points for volume calculation',
     -      46,nranddef,1,0,nrand,33)
          write (iw0,2025) rsolv,nrand
          levout=0
          call getint('Print level',11,0,0,0,levout,0)
          nframetotmin=10
          write (6,2019) 'Calculated volumes',analfile(1:namleno)
          if (itraj  ==  1) write (6,2019) analfile1(1:namleno1)
        elif (ianaltyp  ==  28) :
c         Principal axis calculation
          print *,'Select atoms for principal axis calculation'
          call select(line,nrecdel,idcol,asterisk,n,nslt,index,ixres,
     -      is1,is2,iseqncol1,iseqncol2,inamcol1,inamcol2,irescol1,
     -      irescol2,iqcol1,iqcol2,charge,iatnum,nneig,ineig,indexdel,
     -      0,maxng,maxrec)
          call masktolist(indexa,indexdel,n,natspax,0)
          write (6,2099)
          ic=44
          call condenselist(indexa,natspax,ic,6)
          write (iw0,2099)
          ic=44
          call condenselist(indexa,natspax,ic,iw0)
          write (6,2097) analfile(1:namleno)
          if (itraj  ==  1) :
            write (6,2098)
            call askyn('Do you want to use the input structure instead',
     -        46,1,-1,inputref,0,0)
            xm=800.0
            ym=500.0
            write (6,2019) 'Evolution and distribution plots',
     -        analfile1(1:namleno1)
c           pstitle='Principal axis evolution'
c           lpstitle=24
          ## end if
        elif (ianaltyp  ==  29) :
c         Radius of gyration, hydrodynamic radius, com, dipole moment
          print *,'Select the atoms for the calculation'
          call select(line,nrecdel,idcol,asterisk,n,nslt,index,ixres,
     -      is1,is2,iseqncol1,iseqncol2,inamcol1,inamcol2,
     -      irescol1,irescol2,iqcol1,iqcol2,charge,iatnum,nneig,ineig,
     -      indexdel,0,MAXNEIG,MAXREC)
          call masktolist(indexa,indexdel,n,ngyrats,0)
          write (6,2100)
          ic=12
          call condenselist(indexa,ngyrats,ic,6)
          write (6,2107) analfile(1:namleno)
          if (itraj  ==  1) :
            xm=800.0
            ym=500.0
            write (6,2019) 'Evolution and distribution plots',
     -        analfile1(1:namleno1)
c           pstitle='Gyration and hydrodynamic radii evolution'
c           lpstitle=41
          ## end if
          write (iw0,2100)
          ic=12
          call condenselist(indexa,ngyrats,ic,iw0)
        elif (ianaltyp  ==  30) :
          print *,'Amber analysis is now under the Miscellaneous menu'
          go to 9005
        elif (ianaltyp  ==  40) :
c         Filter solvents
          write (6,2019) 'Filtering results',analfile(1:namleno)
          ierr=0
          call getint('Representative atom (center) of the solvent',
     -      43,1,1,naslv,iarepslv,0)
          call get#real('Solvent space grid spacing',26,8.0,spacing,1,
     -      116)
          rsltmax=0.0
          rcvmax=0.0
          cvmin=0.0
          call quiz(ans,ifilttyp,'c',' ',0,'solvent filtering option',
     -      24,0,5,6,000)
          if (ifilttyp  ==  1  or  ifilttyp  ==  3) call get#real(
     -      'Maximum distance from the nearest solute heavy atom',51,
     -      5.0,rsltmax,1,0)
          if (ifilttyp  >  1  and  ifilttyp  <  5) call get#real(
     -      'CV calculation cutoff distance',30,15.0,rcvmax,1,0)
          if (ifilttyp  ==  2  or  ifilttyp  ==  3) call get#real(
     -      'Minimum CV to keep a solvent',28,0.65,cvmin,1,0)
          if (ifilttyp  ==  4) :
            if (nsegslt  <  2) :
              print *,'There is only one solute molecule'
              return
            elif (nsegslt  ==  2) :
              mol11=1
              mol12=1
              mol21=2
              mol22=2
            else:
              print *,'There are ',nsegslt,' solute molecules'
              call quiz(ans,intftyp,' ',' ',0,'interface partners',18,
     -          0,5,6,0)
              if (ans  ==  'p') :
                call getint('First molecule of the interface',31,1,1,
     -            nsegslt,mol1,000)
                call getint('Second molecule of the interface',32,
     -            nsegslt,1,nsegslt,mol2,000)
                mol11=mol1
                mol12=mol1
                mol21=mol2
                mol22=mol2
              elif (ans  ==  's') :
                call getint('Selected interface molecule',27,1,1,
     -            nsegslt,mol1,000)
                mol11=mol1
                mol12=mol1
                mol21=1
                mol22=nsegslt
              else:
                mol11=1
                mol12=nsegslt
                mol21=1
                mol22=nsegslt
c               All pairs

              ## end if
            ## end if
            ang12min1=110.0
            ang12min2=145.0
            r12max=15.0
            cvrlim=3.5
            cvmin=0.60
            rneigh=5.0
            write (6,2125) cvmin,r12max,ang12min1,ang12min2,cvrlim,
     -        rneigh
            call askyn('Do you want to change these limits',34,1,-1,
     -        inp,0,0)
            if (inp  ==  1) :
              call get#real('Minimum CV of an interface solvent',34,
     -        cvmin,cvmin,1,0)
              call get#real(
     -          'Maximum R1-R2 distance of an interface solvent',48,
     -          r12max,r12max,1,0)
              call get#real(
     -          'Minimum R1-SLV-R2 angle when d(R1-R2) <= 6 A',44,
     -          ang12min1,ang12min1,1,0)
              call get#real(
     -          'Minimum R1-SLV-R2 angle when d(R1-R2) > 6 A',43,
     -          ang12min2,ang12min2,1,0)
              call get#real(
     -          'Maximum CV (wrt the two solute molecules) ratio',47,
     -          cvrlim,cvrlim,1,0)
              call get#real('Radius of sphere that can not be empty',
     -          38,rneigh,rneigh,1,0)
              write (6,2125)cvmin,r12max,ang12min1,ang12min2,cvrlim,
     -          rneigh
            ## end if
          elif (ifilttyp  ==  5) :
c           Hydrogen bond bridging solvent
            call set_hbondlims(hbf0,hblimfac,angm0,angmin,iw0)
            innlist=0
            ihbondcalc=1
          ## end if
          analfile1=inpfile
          call strip_cext(analfile1,namleni,namleno,lenext)
          crdexti(1:lenext)=analfile1(namleni-lenext+1:namleni)
c         print *,'NAMLENI,NAMLENO,LENEXT=',namleni,namleno,lenext
          analfile1(namleno+1:namleno+5)='_filt'
          namleno=namleno+5
          namlen_root=namleno
          if (itraj  ==  1) :
            call askyn(
     -        'Do you want the filtered structures in a single file',52,
     -        1,1,ionefile,0,0)
            write (6,2110) analfile1(1:namlen_root)
            write (iw0,2167)
            if (ifilttyp  ==  1) write (iw0,2164)
     -        'solute solvent distance'
            if (ifilttyp  ==  2) write (iw0,2164)
     -        'CV w.r.t. the solute'
            if (ifilttyp  ==  3) write (iw0,2164)
     -        'solute solvent distance and CV w.r.t. the solute'
            if (ifilttyp  ==  4) write (iw0,2164) 'interface'
            if (ifilttyp  ==  1  or  ifilttyp  ==  3)
     -        write (iw0,2165) rsltmax
            if (ifilttyp  ==  2  or  ifilttyp  ==  3)
     -        write (iw0,2166) cvmin,rcvmax
          ## end if
          if (ionefile  ==  0  and  (itraj  ==  1  or  nconfig  ==  2))
     -        write (6,2102) analfile1(1:namlen_root),
     -        crdexti(1:lenext)//'.pdb'
        else:
          write (6,2041) analtyp
          go to 9002
        ## end if
c         pstitle='Principal axis calculation'
c         lpstitle=26
c##       Early openps and landscape conflict#
        if (iw0  >  0  and  ianaltyp  !=  35) :
          if (title(1:4)  !=  '@#$%') write (iw0,*) 'System:',
     -      title(1:ltitle)
          if (itraj  ==  0) write (iw0,2006)
     -      mark0(1:lmark0),'Structure',inpfile(1:namleni)
          if (itraj  ==  1) write (iw0,2118) inpfile(1:namleni)
        ## end if
        if (ianaltyp  !=  29  and  ianaltyp  !=  35) :
          if (iw1  >  0  and   extnam2  !=  '.ps') :
            if (title(1:4)  !=  '@#$%') write (iw1,2014)mark1(1:lmark1),
     -        mark1(1:lmark1),mark1(1:lmark1),title(1:ltitle)
            if (itraj  ==  0) write (iw1,2006)
     -        mark1(1:lmark1),'Structure',inpfile(1:namleni)
            if (itraj  ==  1) write (iw1,2118) inpfile(1:namleni)
          ## end if
          if (iw2  >  0  and  ianaltyp  !=  22) :
            if( analfile2(namleno2-2:namleno2)  !=  '.ps ') :
              if (title(1:4)  !=  '@#$%')
     -          write (iw2,2014)mark2(1:lmark2), mark2(1:lmark2),
     -            mark2(1:lmark2),title(1:ltitle)
              if (itraj  ==  0) write (iw2,2006)
     -          mark2(1:lmark2),'Structure',inpfile(1:namleni)
              if (itraj  ==  1) write (iw2,2118) inpfile(1:namleni)
            ## end if
          ## end if
        ## end if
      else:
c       nconfig > 1
        if (iw0  >  0) write (iw0,2173) nconfig,nconfig
        if (iw1  >  0) write (iw1,2173) nconfi,nconfigg
        itraj=0
      ## end if
      call cofms(cres,crmslt0,nslt,atw)
      if (ianaltyp  ==  21  or  ianaltyp  ==  22  or  ianaltyp  ==  23
     -   or  ianaltyp  ==  27  or  ianaltyp  ==  28  or  ianaltyp  == 29
     -   or  ianaltyp  ==  33) :
        call askyn('Do you want mass-weighting',26,1,1,masswt,0,0)
        if (masswt  ==  0) :
          write (iw0,2022)
          for ia in range(0, n):
            atw(ia)=1.0
          
        ## end if
      ## end if
      if (ibondtype  >  0  and  ibondtype  !=  1  and  itraj  ==  0)
     -  call printanchorlist(bondname,lbondname,ibondtype,ianchor,
     -    nanchor,indexa,indexov,nanchorr,nanchorn,line,index,iresno,
     -    isegno,segid4,inamcol1,inamcol2,irescol1,irescol2,iw0,maxrsd,
     -    maxrec)
      if (iaskcolcode  ==  1)
     -  call getint(askcolcode,24,maxcolcode,1,maxcolcode,ncolcode,95)
3001  if (itraj  ==  0) :
c       Analysis of single conformations
        if (ianaltyp  ==   1) :
c---------(S) Neighbor, bond, angle and torsion list
          call printbondlist(iangpr,itorpr,1,nslt,nslt,c,nneig,ineig,
     -      line,inamcol1,inamcol2,irescol1,irescol2,iresncol1,
     -      iresncol2,index,cv,indexa,maxng,radtodeg,iw0,maxrec)
        elif (ianaltyp  ==   2) :
c---------(S) 1-4 statistics
          call stat14(c,n,nneig,ineig,index,line,nconfig,pi,iatnm2,
     -      irescol1,irescol2,inamcol1,inamcol2,inpcrdtyp,ioins,iw0,
     -      maxrepconf,maxng,maxrec)
        elif (ianaltyp  ==   3) :
C---------(S) defal group analysis, backbone plot
          call findfg(1,n,iatnum,mmtype,indexo,indexn,nhbneig,ih,ifgtyp,
     -      ibnd,nfg,6,nneig,nhneig,nneiga,ineig,iw0,inpcrdtyp,ioins,
     -      inamcol1,inamcol2,irescol1,irescol2,iresncol1,iresncol2,
     -      line,index,maxng,maxrec)
          call findbackbone(1,nslt,iw1,nneig,ineig,nneiga,
     -      ibnd,inamcol1,inamcol2,irescol1,irescol2,line,blankline,
     -      index,indexn,indexo,indexs,ifree,maxng,maxbox,maxrec)
        elif (ianaltyp  ==   4) :
c---------(S) Bond length statistics
          call bondlenstat(c,n,iatnum,nneig,ineig,iatnm2,nconfig,
     -      nanos,ixlist,ialist,iw0,inpcrdtyp,ioins,inamcol1,inamcol2,
     -      irescol1,irescol2,line,index,radtodeg,maxrepconf,maxng,
     -      maxrec)
        elif (ianaltyp  ==   5) :
c---------(S) Hydrogen-bond list
          call hblist(c,n,nslt,islvw,naslv,iatnum,ifchrg,isegno,nneig,
     -      nneiga,nhbneig,ineig,nhneig,nnneig,ncneig,nsneig,npneig,
     -      ixres,nconfig,hbf0,angm0,molresflag,hblimfac,angmin,
     -      iw0,inamcol1,inamcol2,irescol1,irescol2,iresncol1,
     -      iresncol2,blankline,line,index,indexn,ibnd,indexo,isolvent,
     -      1,0,nosameseg,0,radtodeg,maxrepconf,maxng,maxbox,
     -      maxrsd,maxrec)
          if (ireorient  >  0) call zeroit(centgra,3)
        elif (ianaltyp  ==   6  or  ianaltyp  ==  34) :
c---------(S) Hydrophobic bond or heavy atom VdW contact list
          call nnlisthph_sltb(nslt,ianchor2,iselfanc,indexa,
     -      nneig,npneig,nhbneig,ineig,c,rhphmax,0,isegno,nosameseg,
     -      bondlab(ibondtype),4,ibnd,indexo,ifail,maxng,maxng,
     -      maxrec)
          call hph_sltblist(ibondtype,c,nslt,nhbneig,ineig,iw0,is1,
     -      is2,inamcol1,inamcol2,irescol1,irescol2,iresncol1,
     -      iresncol2,line,index,bondname,lbondname,bondlab(ibondtype),
     -      4,maxng)
        elif (ianaltyp  ==  41) :
c---------(S) Mutually proximal heavy atom contact list
          call nnlistmpx(n,nanchorr,nanchorn,indexa,indexov,iatnum,
     -      nneig,nhbneig,c,isegno,ibnd,indexo,rmpxlim,itemp1,itemp2,
     -      itemp3,itemp4,indexn,temp,cv,ifail,maxng,maxrec)
          call mpxblist(n,itemp1,itemp2,cv,itemp3,itemp4,is1,is2,
     -      inamcol1,inamcol2,irescol1,irescol2,iresncol1,iresncol2,
     -      line,index,ixres,nhbdist,rhbdist,nbfound,nbresfound,nbonds,
     -      c,rmpxlim,ifailbond,iw0,MAXBONDS,maxrec)
        elif (ianaltyp  ==   7) :
c---------(S) Salt bridge list
          call nnlisthph_sltb(nslt,ianchor2,iselfanc,indexa,
     -      nneig,npneig,nhbneig,ineig,c,rsltbmax,1,isegno,nosameseg,
     -      bondlab(ibondtype),4,ibnd,indexo,ifail,maxng,maxng,
     -      maxrec)
          call hph_sltblist(ibondtype,c,nslt,nhbneig,ineig,iw0,
     -      is1,is2,inamcol1,inamcol2,irescol1,irescol2,
     -      iresncol1,iresncol2,line,index,bondname,lbondname,
     -      bondlab(3),4,maxng)
        elif (ianaltyp  ==   8) :
c---------(S) Calculate residue-residue distances
          call rrdist(c,n,nres,ixresno,iresno,ixres,resnames,ifres,
     -      ilres,increst,incsolvrr,irefres1,irefres2,irefseg1,irefseg2,
     -      irefresinc,resdistlim,inegres1,inegres2,inegseg1,inegseg2,
     -      inegresinc,resapplim,listrefres,nrefres,nrefrange,l
     -      istnegres,nnegres,nnegrange,iatnum,ignoreh,0,0,iw0,iw1,iw2,
     -      irescol1,irescol2,inamcol1,inamcol2,is1,is2,irescount1,
     -      irescount2,irescount3,line,index,indexn,indexo,indexs,
     -      indexa,itemp1,itemp2,temp,cv,maxrec,maxrsd)
        elif (ianaltyp  ==   9) :
c---------(S) Calculate a distance
          call pbcdist(c(1,ia1),c(1,ia2),ia1,ia2,cell,ncell,0,0,
     -      img,dr,d)
        elif (ianaltyp  ==  10) :
c---------(S) Check for potentially unphysical contacts
          call checkunphys(c,nslt,n,naslv,islvw,iatnum,ifchrg,isegno,
     -      idelse:g,indexo,nneig,nneiga,nhbneig,ineig,nhneig,nnneig,
     -      ncneig,nsneig,npneig,ixres,line,irescol1,irescol2,inamcol1,
     -      inamcol2,index,nconfig,innlist,molresflag,ioppbc,cell,ncell,
     -      edge,ixyzhex,molsltlim,nmolslt,hblimfac,angmin,ctfac,
     -      bondminfac,maxdist,iles,iw0,isltonly,ibnd,indexn,rlim,
     -      innlistread,0,radtodeg,nerrunphys,maxrepconf,maxng,maxbox,
     -      maxrsd,maxrec)
          if (nerrunphys  ==  0) write (6,*) 'Nothing untoward detected'
        elif (ianaltyp  ==  11) :
c---------(S) Pseudorotation angle calculation
          write (6,2078)
          nmem=5
          while ( True )
            call getring(line,index,ix5,irescol1,irescol2,inamcol1,
     -        inamcol2,nslt,numres,iresring,iresno,ifres,ilres,nmem,
     -        psrtyp,0,incgen,iapex,MAXRING,maxrec)
            if (nmem  ==  0) go to 9001
            call pseudorot(c,n,ix5,nmem,nneig,ineig,psr5,psr,npsr,
     -        psrtyp,incgen,zring,rring,sinpsrs,cospsrs,qpsr,qpsrs,
     -        qpsr2s,zav,zsq,zavs,zsqs,1,6,0,radtodeg,maxng,maxrec)
            print *
          
        elif (ianaltyp  ==  12) :
c---------(S) Calculate Proline kinks
          if (noprol  ==  0)
     -      call pseudorot(c,n,ix5,nmem,nneig,ineig,psr5,psr,npsr,'p',
     -        incgen,zring,rring,sinpsrs,cospsrs,qpsr,qpsrs,qpsr2s,
     -        zav,zsq,zavs,zsqs,1,6,0,radtodeg,maxng,maxrec)
          if (nconfig  >  2) iprintpk=0
          call prokinkcalcla(1,c,nslt,b# end,wobble,faceshift,
     -      rmsb,rmsa,
     -      iflatproline,ix5,radtodeg,MAXHX)
          if (noprol  ==  0)
     -      call pseudorot(c,n,ix5,nmem,nneig,ineig,psr5,psr,npsr,'p',
     -        incgen,zring,rring,sinpsrs,cospsrs,qpsr,qpsrs,qpsr2s,
     -        zav,zsq,zavs,zsqs,0,6,0,radtodeg,maxng,maxrec)
          write (iw0,2072) nconfig,b# end,wobble,faceshift,
     -      wobble-faceshift,psr(2),rmsb,rmsa
        elif (ianaltyp  ==  13) :
c---------(S) Hydropathy labeling
          call hydropathylist(n,nslt,ixres,resnames,cv,ihydtyp,
     -      maxrsd,maxrec)
          call writeconf(iw0,inpcrdtyp,inpcrdtyp,inpcrdtyporg,n,n,nslt,
     -      naslv,islvw,iasv,namesv,qsv,pflsv,1,iwhead,0,iatnum,ifchrg,
     -      nconfig,innlist,c,rprox,cv,ixres,iresno,atnames,resnames,
     -      segnames,charge,isegno,altcol,inscol,ninsres,marker,ntitlin,
     -      ntitlinw,title,ireseq,iresnrestart,iresidrestart,nneig,
     -      nneiga,nhbneig,ineig,nhneig,nnneig,ncneig,nsneig,npneig,
     -      numres,numslv,resnamslv,line,blankline,mmtype,ibnd,index,
     -      indexn,indexo,1,molresflag,irescount3,itemp1,hblimfac,
     -      angmin,0,1,1,1,0,2,iqspaceask,13,0,0.0,0,0,0,keeprem,
     -      iwriteatsym,radtodeg,maxrepconf,maxng,maxrsd,maxrec)
        elif (ianaltyp  ==  14) :
c---------(S) Circular variance labeling of solute atoms and solvent molecules
          call cvlist(c,n,nslt,nsltref_f,nsltref_l,naslv,islvrep,icvtyp,
     -      rcut_cv,rprox,cv,itemp1,c2,temp,indexn,indexo,indexs,chn,
     -      isortslv,cvlim)
          call writeconf(iw0,inpcrdtyp,inpcrdtyp,inpcrdtyporg,n,n,nslt,
     -      naslv,islvw,iasv,namesv,qsv,pflsv,1,iwhead,0,iatnum,ifchrg,
     -      nconfig,innlist,c,rprox,cv,ixres,iresno,atnames,resnames,
     -      segnames,charge,isegno,altcol,inscol,ninsres,marker,ntitlin,
     -      ntitlinw,title,ireseq,iresnrestart,iresidrestart,nneig,
     -      nneiga,nhbneig,ineig,nhneig,nnneig,ncneig,nsneig,npneig,
     -      numres,numslv,resnamslv,line,blankline,mmtype,ibnd,index,
     -      indexn,indexo,1,molresflag,irescount3,itemp1,hblimfac,
     -      angmin,0,1,1,1,0,3,iqspaceask,14,0,0.0,0,0,0,keeprem,
     -      iwriteatsym,radtodeg,maxrepconf,maxng,maxrsd,maxrec)
        elif (ianaltyp  ==  15) :
c---------(S) Circular variance residue-residue plot
          call cvplot(c,n,nslt,icvtyp,line,index,indexs,inamcol1,
     -      inamcol2,iresncol1,iresncol2,'CA      ',title,ltitle,
     -      ncolcode,maxcolcode,iw0,iw1,1,maxrec,ipspage)
        elif (ianaltyp  ==  16) :
c---------(S) DSSP secondary structure assignment
          call dssp(c,1,n,nslt,line,index,inamcol1,inamcol2,iresncol1,
     -      iresncol2,nneig,ineig,nneiga,ibnd,indexn,indexo,npneig,
     -      nsneig,nnneig,dssplab,idistdssp,chn,c2,cv,indexs,indexa,
     -      nsse,itypsse,ifsse,ilsse,0,nconfig,iw0,0,1,ifail,radtodeg,
     -      maxrepconf,maxng*10,maxrec/10,maxng*20,maxrec/20,200,
     -      maxrsd,maxrec)
        elif (ianaltyp  ==  17) :
c---------(S) Hydrogen-bond bridge analysis
          call hblist(c,n,nslt,islvw,naslv,iatnum,ifchrg,isegno,nneig,
     -      nneiga,nhbneig,ineig,nhneig,nnneig,ncneig,nsneig,npneig,
     -      ixres,nframe,hbf0,angm0,molresflag,hblimfac,angmin,
     -      iw0,inamcol1,inamcol2,irescol1,irescol2,iresncol1,
     -      iresncol2,blankline,line,index,indexn,ibnd,indexo,1,0,1,
     -      nosameseg,0,radtodeg,maxrepconf,maxng,maxbox,maxrsd,maxrec)
          call hbbridge(nanchor,ianchor,indexa,ianchor2,iselfanc,lpath,
     -      nbridgetype,ibridgetype,maxbridgemem,n,nhbneig,ineig,indexn,
     -      indexo,indexs,ixres,resnames,brslv,nabr,nrescol,0,ifail,
     -      listbridge,iw0,maxng,MAXBRIDGELEN,MAXBRIDGETYPE,
     -      minbridgelenprint,maxrsd,maxrec)
          call hbbridgeprint(nanchor,ianchor,lpath,nbridgetype,
     -      ibridgetype,maxbridgemem,line,index,iresno,inamcol1,
     -      inamcol2,irescol1,irescol2,iw0,maxbondcount,
     -      maxhbtype,minbridgelenprint,0,1,MAXBRIDGELEN,MAXBRIDGETYPE,
     -      maxrec)
        elif (ianaltyp  ==  18) :
c---------(S) Ramachandran plot
          call ramachandran(c,nslt,index,line,nconfig,pi,nresfound,
     -      iresno,ixres,irescol1,irescol2,inamcol1,iw0,
     -      maxng,maxrec)
          call ramachandranplot(nresfound,iw1,xm,iallrama)
        elif (ianaltyp  ==  32) :
c---------(S) Angle dial plots
          call angledials(c,nslt,nangsel,ixtor1234,index,line,pi,
     -      irescol1,irescol2,inamcol1,inamcol2,iw0,maxrec)
        elif (ianaltyp  ==  19) :
c---------(S) Torsion dial plots
          call torsiondials(c,nslt,ntorsel,ixtor1234,index,line,pi,
     -      irescol1,irescol2,inamcol1,inamcol2,iw0,maxrec)
        elif (ianaltyp  ==  20) :
c---------(S) Delphi map annotation
          call delphilabel(c,n,nslt,xstart,ystart,zstart,gx,gy,gz,cv)
          call writeconf(iw0,inpcrdtyp,inpcrdtyp,inpcrdtyporg,n,n,nslt,
     -      naslv,islvw,iasv,namesv,qsv,pflsv,1,iwhead,0,iatnum,ifchrg,
     -      nconfig,innlist,c,rprox,cv,ixres,iresno,atnames,resnames,
     -      segnames,charge,isegno,altcol,inscol,ninsres,marker,ntitlin,
     -      ntitlinw,title,ireseq,iresnrestart,iresidrestart,nneig,
     -      nneiga,nhbneig,ineig,nhneig,nnneig,ncneig,nsneig,npneig,
     -      numres,numslv,resnamslv,line,blankline,mmtype,ibnd,index,
     -      indexn,indexo,1,molresflag,irescount3,itemp1,hblimfac,
     -      angmin,0,1,1,1,0,0,iqspaceask,20,0,0.0,0,0,0,keeprem,
     -      iwriteatsym,radtodeg,maxrepconf,maxng,maxrsd,maxrec)
        elif (ianaltyp  ==  21) :
c---------(S) Helix axis directions
          for ihx in range(0, nhx):
            call helixaxis(c,nslt,iw0,calph(1,1,ihx),axisdir(1,ihx),
     -        axisini(1,ihx),axis# end(1,ihx),helixcent(1,ihx),
     -        perpvec(1,1,ihx),camod(1,1,ihx),anglechangeref(1,ihx),
     -        circ(1,ihx),rnorm(1,ihx),axfact(1,ihx),axtol,rot,rms,
     -        helixlen(ihx),angles(1,ihx),decideb# end,nup,ndown,nrun,
     -        nnear,rcirchx(ihx),turnperres(ihx),anglesn(1,ihx),0,
     -        incrot,nrep,0,nreshx(ihx),icaahx(1,ihx),ihx,nhxres,
     -        idebughx,radtodeg,pi,MAXHX)
            iadssp1=ifres(ixres(icaahx(1,ihx)))
            iadssp2=ilres(ixres(icaahx(nreshx(ihx),ihx)))
            dssp(c,iadssp1,n,iadssp2,line,index,inamcol1,inamcol2,iresncol1,iresncol2,nneig,ineig,nneiga,ibnd,indexn,indexo,npneig,nsneig,nnneig,dssplab,idistdssp,chn,c2,cv,indexs,indexa,nsse,itypsse,ifsse,ilsse,ireshx1(ihx)-1,nframe,0,0,0,ifail,radtodeg,maxrepconf,maxng*10,maxrec/10,maxng*20,maxrec/20,200,maxrsd,maxrec)
            lhelixcklab=9+1
            write (helixcklab,2160) 'H',ihx
            checkforhelix(ireshx2(ihx)-ireshx1(ihx)+1,dssplab,indexn,iw0,hxoklab,ihxok,helixcklab,lhelixcklab,maxrsd,maxrec)
            nhelixok(ihxok)=nhelixok(ihxok)+1
          
          if (nhx  >  1): 
            multihelix(iw0,nhx,nhxres,radtodeg,c,icaahx,icbahx,icbreshx,maxrec,MAXHX,MAXNHX)
        elif (ianaltyp  ==  26) :
c---------(S) Distance calculation
          if (iclusterdist  ==  0) :
            call pairdistcalc(c,nslt,npairs,listpairdist,pairdistsum,
     -        pairdistsum2,pairdistwsum,npairdist,pairdistminmax,
     -        pairgrid,iw0,MAXDDBIN,MAXDDISTR)
          else:
            call clusterdistcalc(c,nslt,npairs,iclustermem,
     -        ifstclst1,ifstclst2,ilstclst2,pairdistsum,pairdistsum2,
     -        pairdistwsum,npairdist,pairdistminmax,pairgrid,
     -        iw0,MAXDDBIN,MAXDDISTR,MAXCDLIST)
          ## end if
          call pairdistprint(nframe,npairs,listpairdist,iclusterdist,
     -      iclustermem,ifstclst1,ifstclst2,ilstclst2,pairdistsum,
     -      pairdistsum2,pairdistwsum,npairdist,pairdistminmax,
     -      pairgrid,rmaxpair,line,index,inamcol1,inamcol2,irescol1,
     -      irescol2,iresncol1,iresncol2,inpcrdtyp,ioins,iw0,nslt,
     -      MAXDDBIN,MAXDDISTR,MAXCDLIST,maxrec)
        elif (ianaltyp  ==  27) :
c---------(S) Solvation shell volume calculation
          call volcalc(nrand,c,nslt,isegno,iatnum,nsegslt,molsltlim,
     -      itemp1,cres,chn,c2,atw,vmac,vshell,vint,vmacsd,vshellsd,
     -      vintsd,vfrstsh,vfrstshsd,rsolv,ih,indexn,indexo,indexs,
     -      indexa,indexdel,1,iw0,levout,maxrsd,maxrec)
        elif (ianaltyp  ==  28) :
c---------(S) Principal axis calculation
          call princax(c,c2,atw,temp,n,indexdel,evecs0,evals0,iw0,0,1,
     -      radtodeg,0)
          call askyn('Do you want to align the structure to an axis',
     -      45,1,0,ialign,0,0)
          if (ialign  ==  1) :
            call getint('Axis index (1/2/3) to align',27,3,1,3,iax,000)
            iemax=1
            for k in range(2, 3):
              if (evals0(k)  >  evals0(iemax)) iemax=k
            
            call cofms(c,crmslt,nslt,aw)
            call shiftmol(c,nslt,crmslt,c2,-1.0)
            for k in range(0, 3):
              rot(iax,k)=evecs0(iemax,k)
              rot(mod(iax,3)+1,k)=evecs0(mod(iemax,3)+1,k)
              rot(mod(iax+1,3)+1,k)=evecs0(mod(iemax+1,3)+1,k)
            
            call rotate_c(c2,n,rot,c2,'PRINCAX',7)
            call writeconf(iw1,inpcrdtyp,inpcrdtyp,inpcrdtyporg,n,n,
     -        nslt,naslv,islvw,iasv,namesv,qsv,pflsv,1,iwhead,0,iatnum,
     -        ifchrg,nconfig,innlist,c2,rprox,cv,ixres,iresno,atnames,
     -        resnames,segnames,charge,isegno,altcol,inscol,ninsres,
     -        marker,ntitlin,ntitlinw,title,ireseq,iresnrestart,
     -        iresidrestart,nneig,nneiga,nhbneig,ineig,nhneig,nnneig,
     -        ncneig,nsneig,npneig,numres,numslv,resnamslv,line,
     -        blankline,mmtype,ibnd,index,indexn,indexo,1,molresflag,
     -        irescount3,itemp1,hblimfac,angmin,0,1,1,1,0,2,iqspaceask,
     -        13,0,0.0,0,0,0,keeprem,iwriteatsym,radtodeg,maxrepconf,
     -        maxng,maxrsd,maxrec)
          else:
            close (iw1,status='delete')
          ## end if
        elif (ianaltyp  ==  29) :
c---------(S) Radius and dipole calculation
            call molrad(c,indexa,ngyrats,iw0,MAXREC)
            call celldipole(c,n,nslt,indexa,ngyrats,charge,icharges,atw,
     -        iw0,1)
        elif (ianaltyp  ==  31) :
c---------(S) Calculate adjacency-matrix based analysis
          call rrconn(c,n,ires1,ires2,ifres,ilres,iatnum,ignoreh,
     -      irepuse,iadjtyp,iscalesum,resdistlim,nexpmax,npint,ipspage,
     -      npspages,iw0,iw1,line,index,irescol1,inamcol1,
     -      imarkres,marks,itemp4,itemp1,itemp2,itemp3,
     -      temp,inpfile,namleni,analfile4,namleno4,maxrec)
        elif (ianaltyp  ==   33) :
c---------(S) Calculate molecule-molecule distance matrix
            call mmdist(c,n,atw,iatnum,nmolslt,molsltlim,c2,temp,
     -        ignoreh,iw0,iw1)
        elif (ianaltyp  ==   40) :
c---------(S) Filter solvents
          if (ifilttyp  ==  5)
     -      call nnlist(n,islvw,naslv,n,iatnum,ifchrg,c,nneig,nneiga,
     -        nhbneig,ineig,nhneig,nnneig,ncneig,nsneig,npneig,line,
     -        irescol1,irescol2,inamcol1,inamcol2,index,nconfig,innlist,
     -        molresflag,hblimfac,angmin,ihbondcalc,ibnd,indexo,isegno,
     -        ixres,maxrepconf,0,0,radtodeg,0,maxbox,maxng,maxrsd,
     -        maxrec)
          call filterslv(c,iatnum,nslt,n,naslv,numsolv,molsltlim,
     -      nsegslt,index,nconfig,ifilttyp,intftyp,rsltmax,rcvmax,cvmin,
     -      ang12min1,ang12min2,r12max,cvrlim,rneigh,mol11,mol12,mol21,
     -      mol22,spacing,itemp1,iarepslv,numsolvleft,nhbneig,ineig,
     -      ixres,nrecdel,ierr,6,maxrec,maxng)
          write (iw0,2126) numsolvleft
          if (nconfig  ==  2)
     -      call askyn(
     -      'Do you want all filtered structures in a single file',52,
     -        1,1,ionefile,0,0)
          if (numsolvleft  >  0) :
            if (nconfig  >  1  and  ionefile  ==  0) :
              namleno=namlen_root+1
              analfile1(namleno:namleno)='.'
              call writeint(analfile1,namleno,nconfig,lenc)
              namleno=namleno+lenc
            ## end if
            if (nconfig  ==  1  or  ionefile  ==  0) :
              analfile1(namleno+1:namleno+lenext)=crdexti(1:lenext)
              namleno=namleno+lenext
              call openfile(45,0,'FILT',4,'new',analfile1,namleno,
     -          notfnd,0,1,1,0,0)
            ## end if
            call writeout(45,inpcrdtyp,inpcrdtyp,line,index,isegno,
     -        n-nrecdel,marker,1,0,ntitlin,ntitlinw,title,blankline,
     -        1,1,0,0.0,0,0,0,keeprem,iwriteatsym,iatnum,maxrec)
            close (45)
          else:
            write (6,2070) nconfig
          ## end if
        ## end if
9001    if (nconfig  ==  1) go to 9005
      else:
c       Trajectory scan
c       Open, initialize trajectory
        if (ireadtracks  ==  0) call asktrajform(inptrajtyp,ioutrajtyp,
     -    mmctrajtyp,resnamslv,-1,1)
        inpt=70
        iverbose=1
        notitprint=0
        iopen_rep=0
      ## end if
      numsolvmin=999999
      numsolvmax=0
      numsolvsum=0
      irepscan=1
      nrepinc=0
      nframeread=0
      if (itraj  ==  1  and  ianaltyp  ==  23) :
c       2D RMSD calculation in blocks
        irepscan=0
        ifail=0
        maxconf=999999
        while (maxconf  >  mx2d)
          call opentraj(c,1,inpt,inptrajtyp,n,ntitltr,trtitle,
     -      inpcrdtyp,ifirst,ilast,increment,maxconf,ninconf,noutconf,
     -      natom,nfreeat,ifree,icntrl,0,mmctrajtyp,trajnam,ltrajnam,
     -      'trajectory',10,iconfsel,numsel,iverbose,notitprint,0,0,
     -      icellfound,notfnd,iconfirmname,iopen_rep,lentest,iw0,
     -      ianaltyp,mx2d,maxconfsel,maxrec)
        
        call save_traj_lim(ifirst,ilast,increment,1)
        if (iw0  >  0) write (iw0,2006) mark0(1:lmark0),
     -    'Trajectory',trajnam(1:ltrajnam)
        write (iw0,2060) 'Trajectory',ifirst,ilast,increment
        nframetot=(ilast-ifirst)/increment+1
        nframetot=nframetot*(nframetot-1)/2
        nframesign=nframetot/10
        ncmem=(maxrec10-n)/nslt+1
        call adjust_xtraj(xtraj,ifirst,ilast,increment,framefac,
     -    iadjust_xtraj)
        write (6,2157) ncmem
        write (iw0,2157) ncmem
        nframe1=0
        nframe2d=0
        ic1=1
        ncop=0
        iehead=1
        while (ninconf  <  ilast  and  ifail  ==  0)
c         Read a conformation
          call readtraj(inpt,inptrajtyp,mmctrajtyp,n,naslv,nwatr,
     -      nslt,trajformatname(inptrajtyp),ntitltr,trtitle,trajnam,
     -      ltrajnam,natom,nfreeat,ifree,icntrl,c(1,ic1),ninconf,
     -      noutconf,increment,inpcrdtyp,ietotsaved,etot,ifail,ifirst,
     -      ilast,iconfsel,numsel,maxrepconf*irepscan,nmc,lentest,1.0,
     -      maxconf,maxconfsel,maxrec)
          if (ietotsaved  ==  1  and  iw0  >  0) :
            if (iehead  ==  1) write (6,*)
     -        'Energies read will be listed on the log file'
            iehead=0
            write (iw0,2175) ninconf,etot
          ## end if
          icsel=0
          if (ifail  ==  0) :
            nframeread=nframeread+1
            call selectconf(numsel,ninconf,ifirst,increment,
     -        iconfsel,nextconfsel,nframeref,icsel,ifail,maxconfsel)
          else:
            ninconf=ilast
            nframetot=nframe
            if (nframeread  ==  0) go to 9002
          ## end if
          if (icsel  ==  1) :
            ic1=ic1+nslt
            ncop=ncop+1
            it1(ncop)=ninconf
            nframe1=nframe1+1
            res(2,nframe1,11)=etot2
            res(1,nframe1,11)=etot
          ## end if
          if (ninconf  ==  ilast  or  ncop  ==  ncmem) :
c           Calculate RMSD between the ncop structures in c
            for nc1 in range(0, ncop):
              nframe=nframe1-ncop+nc1
              ic1=(nc1-1)*nslt+1
              for nc2 in range(nc1+1, ncop):
                ic2=(nc2-1)*nslt+1
                nframeref=nframe1-ncop+nc2
                write (iw0,2173) nframe,it1(nc1),' ',nframeref,it1(nc2)
                call rmsd(c(1,ic2),c(1,ic1),nslt,nfinalov,nfinalrmsd,
     -            atw,atwsum,temp,chn,c2,indexov,indexov,noopt2d,
     -            1,indexrmsd,indexrmsd,rot,com1,com2,etot,etot2,iw0,
     -            devmax,devmaxnoopt,maxrec)
                if (absdevmin  >  devmax) absdevmin=devmax
                if (absdevmax  <  devmax) absdevmax=devmax
                nframe2d=nframe2d+1
                call progress_rep(0,nframe2d,nframesign)
              
            
c           print *,'Self block done nframe1=',nframe1
c           Open same trajectory and skip to nframeread+increment
            ifirst=nframeread+increment
            nframe2=nframe1
            rewind inpt
            call opentraj(c2,0,inpt,inptrajtyp,n,ntitltr,trtitle,
     -        inpcrdtyp,ifirst,ilast,increment,maxconf,ninconf2,
     -        noutconf,natom,nfreeat,ifree,icntrl,0,mmctrajtyp,trajnam,
     -        ltrajnam,'trajectory',10,iconfsel,numsel2,0,1,0,0,
     -        icellfound,notfnd,iconfirmname,iopen_rep,lentest,0,
     -        ianaltyp,mx2d,maxconfsel,maxrec)
            notitprint=1
            iopen_rep=1
c           Read blocks into cref and calculate RMSD with the c block
            ic2=1
            ncop2=0
            nframeread2=0
            ifail2=0
            while (ninconf2  <  ilast  and  ifail2  ==  0)
c             Read a conformation
              call readtraj(inpt,inptrajtyp,mmctrajtyp,n,naslv,nwatr,
     -          nslt,trajformatname(inptrajtyp),ntitltr,trtitle,trajnam,
     -          ltrajnam,natom,nfreeat,ifree,icntrl,cres(1,ic2),
     -          ninconf2,noutconf,increment,inpcrdtyp,ietotsaved,etot,
     -          ifail2,ifirst,ilast,iconfsel,numsel2,
     -          maxrepconf*irepscan,nmc,lentest,1.0,maxconf,maxconfsel,
     -          maxrec)
              icsel=0
              if (ifail2  ==  0) :
                nframeread2=nframeread2+1
                call selectconf(numsel,ninconf2,ifirst,increment,
     -            iconfsel,nextconfsel,nframeref,icsel,ifail2,
     -            maxconfsel)
              ## end if
              if (icsel  ==  1) :
                ic2=ic2+nslt
                ncop2=ncop2+1
                nframe2=nframe2+1
                it2(ncop2)=ninconf2
              ## end if
              if (ninconf2  ==  ilast  or  ncop2  ==  ncmem) :
c               RMSD calc between c and cres
                for nc1 in range(0, ncop):
                  ic1=(nc1-1)*nslt+1
                  nframe=nframe1-ncop+nc1
                  for nc2 in range(0, ncop2):
                    ic2=(nc2-1)*nslt+1
                    nframeref=nframe2-ncop2+nc2
                    write (iw0,2173) nframe,it1(nc1),' ',nframeref,
     -                it2(nc2)
                    call rmsd(c(1,ic1),cres(1,ic2),nslt,nfinalov,
     -                nfinalrmsd,atw,atwsum,temp,chn,c2,indexov,indexov,
     -                noopt2d,1,indexrmsd,indexrmsd,rot,com1,com2,etot,
     -                etot2,iw0,devmax,devmaxnoopt,maxrec)
                    if (absdevmin  >  devmax) absdevmin=devmax
                    if (absdevmax  <  devmax) absdevmax=devmax
                    nframe2d=nframe2d+1
                    call progress_rep(0,nframe2d,nframesign)
                  
                
c               print *,'Block pair done nframe1,2=',nframe1,nframe2
                ic2=1
                ncop2=0
              ## end if
            
            numsel2=1
            ic1=1
            ncop=0
            if (ninconf  <  ilast) :
              nskip=ninconf
              rewind inpt
              call opentraj(c,1,inpt,inptrajtyp,n,ntitltr,trtitle,
     -          inpcrdtyp,ifirst,ilast,increment,maxconf,ninconf,
     -          noutconf,natom,nfreeat,ifree,icntrl,0,mmctrajtyp,
     -          trajnam,ltrajnam,'trajectory',10,iconfsel,numsel,
     -          0,notitprint,1,0,icellfound,notfnd,iconfirmname,
     -          iopen_rep,lentest,iw0,ianaltyp,mx2d,maxconfsel,maxrec)
              while (ninconf  <  nskip)
                call readtraj(inpt,inptrajtyp,mmctrajtyp,n,naslv,nwatr,
     -            nslt,trajformatname(inptrajtyp),ntitltr,trtitle,
     -            trajnam,ltrajnam,natom,nfreeat,ifree,icntrl,c(1,ic1),
     -            ninconf,noutconf,increment,inpcrdtyp,ietotsaved,etot,
     -            ifail,ifirst,ilast,iconfsel,numsel,
     -            maxrepconf*irepscan,nmc,lentest,1.0,maxconf,
     -            maxconfsel,maxrec)
              
            ## end if
          ## end if
        
        nframe=nframe1
        nframeref=nframe2+1
        close (inpt)
        write (iw0,2085)
        write (iw0,2020) trajnam(:ltrajnam),nframeread,nframe,' ',
     -    ifirst-1
        if (lpstitle  >  0) :
          call openps(iw1,xm,ym,title(1:76),76,pstitle,lpstitle,trajnam,
     -      ltrajnam,trajnam2,0,npspages,ipspage)
        ## end if
      elif (itraj  ==  1  and  ianaltyp  ==  24) :
c       Cross-RMSD in blocks
        irepscan=0
        ifail1=0
        maxconf=999999
        ifirst1=1
        ifirst2=1
        while (maxconf  >  mx2d)
          call opentraj(c,1,inpt,inptrajtyp,n,ntitltr,trtitle,
     -      inpcrdtyp,ifirst1,ilast1,increment1,maxconf,ninconf1,
     -      noutconf,natom,nfreeat,ifree,icntrl,0,mmctrajtyp,trajnam,
     -      ltrajnam,'trajectory',10,iconfsel,numsel,iverbose,
     -      notitprint,0,0,icellfound,notfnd,iconfirmname,iopen_rep,
     -      lentest,iw0,ianaltyp,mx2d,maxconfsel,maxrec)
        
        call save_traj_lim(ifirst1,ilast1,increment1,1)
        if (increment1  >   2  and  nslt  >  200) :
          write (6,2090)
          call askstop(0)
        ## end if
        if (iw0  >  0) write (iw0,2006) mark0(1:lmark0),
     -    'Trajectory',trajnam(1:ltrajnam)
        write (iw0,2060) 'Trajectory',ifirst1,ilast1,increment1
        nframetot1=(ilast1-ifirst1)/increment1+1
        ncmem=(maxrec10-n)/nslt+1
        call adjust_xtraj(xtraj,ifirst1,ilast1,increment1,framefac,
     -    iadjust_xtraj)
        write (6,2157) ncmem
        write (iw0,2157) ncmem
        nframe1=0
        ic11=1
        ncop1=0
        nframe2d=0
        nframeread1=0
        nframesign=1
        inpt2=80
        while (ninconf1  <  ilast1  and  ifail1  ==  0)
c         Read a conformation
          call readtraj(inpt,inptrajtyp,mmctrajtyp,n,naslv,nwatr,
     -      nslt,trajformatname(inptrajtyp),ntitltr,trtitle,trajnam,
     -      ltrajnam,natom,nfreeat,ifree,icntrl,c(1,ic11),ninconf1,
     -      noutconf,increment,inpcrdtyp,ietotsaved,etot,ifail1,ifirst1,
     -      ilast1,iconfsel,numsel,maxrepconf*irepscan,nmc,lentest,1.0,
     -      maxconf,maxconfsel,maxrec)
          icsel=0
          if (ifail1  ==  0) :
            nframeread1=nframeread1+1
            call selectconf(numsel,ninconf1,ifirst1,increment1,
     -        iconfsel,nextconfsel,nframeref,icsel,ifail1,maxconfsel)
          else:
            ninconf1=ilast
            nframetot1=nframe
            if (nframeread1  ==  0) go to 9002
          ## end if
          if (icsel  ==  1) :
            ic11=ic11+nslt
            ncop1=ncop1+1
            it1(ncop1)=ninconf1
            nframe1=nframe1+1
            res(2,nframe1,11)=etot2
            res(1,nframe1,11)=etot
          ## end if
          if (ninconf1  ==  ilast1  or  ncop1  ==  ncmem) :
c           Read blocks of traj2
            call opentraj(cres,0,inpt2,inptrajtyp,n,ntitltr,trtitle,
     -        inpcrdtyp,ifirst2,ilast2,increment2,maxconf,ninconf2,
     -        noutconf,natom,nfreeat,ifree,icntrl,1,mmctrajtyp,
     -        trajnam2,ltrajnam2,'second trajectory',17,iconfsel2,
     -        numsel2,iverbose,1,0,0,icellfound,notfnd,0,iopen_rep2,
     -        lentest,0,ianaltyp,mx2d,maxconfsel,maxrec)
            call save_traj_lim(ifirst2,ilast2,increment2,2)
            if (increment2  >   2  and  nslt  >  200) :
              write (6,2090)
              call askstop(0)
            ## end if
            if (nframe2d  ==  0) :
              write (iw0,2006) mark0(1:lmark0),'Second trajectory',
     -          trajnam2(1:ltrajnam2)
              write (iw0,2060) 'Second trajectory',ifirst2,ilast2,
     -          increment2
              nframetot2=(ilast2-ifirst2)/increment2+1
              nframesign=nframetot1*nframetot2/10
            ## end if
            nframe2=0
            nframeread2=0
            ninconf2=0
            ic21=1
            ncop2=0
            ifail2=0
            while (ninconf2  <  ilast2  and  ifail2  ==  0)
c             Read a conformation
              call readtraj(inpt2,inptrajtyp,mmctrajtyp,n,naslv,nwatr,
     -          nslt,trajformatname(inptrajtyp),ntitltr,trtitle,trajnam,
     -          ltrajnam,natom,nfreeat,ifree,icntrl,cres(1,ic21),
     -          ninconf2,noutconf,increment,inpcrdtyp,ietotsaved,etot,
     -          ifail2,ifirst2,ilast2,iconfsel,numsel2,
     -          maxrepconf*irepscan,nmc,lentest,1.0,maxconf,maxconfsel,
     -          maxrec)
              icsel=0
              if (ifail2  ==  0) :
                nframeread2=nframeread2+1
                call selectconf(0,ninconf2,ifirst2,increment2,
     -            iconfsel,nextconfsel,nframeref,icsel,ifail2,
     -            maxconfsel)
              else:
                print *,'TRAJ 2 # endED NINCONF2=',ninconf2
                ninconf2=ilast
                nframetot=nframe2
                if (nframeread2  ==  0) go to 9002
              ## end if
              if (icsel  ==  1) :
                ic21=ic21+nslt
                ncop2=ncop2+1
                it2(ncop2)=ninconf2
                nframe2=nframe2+1
              ## end if
              if (ninconf2  ==  ilast2  or  ncop2  ==  ncmem) :
c               RMSD calc between c and cres
c               print *,'NFRAME1,NCOP1,NFRAME2,NCOP2=',
c    -            NFRAME1,NCOP1,NFRAME2,NCOP2
                for nc1 in range(0, ncop1):
                  ic1=(nc1-1)*nslt+1
                  nframe=nframe1-ncop1+nc1
                  for nc2 in range(0, ncop2):
                    ic2=(nc2-1)*nslt+1
                    nframeref=nframe2-ncop2+nc2
                    write (iw0,2173) nframe,it1(nc1),' ',nframeref,
     -                it2(nc2)
                    call rmsd(c(1,ic1),cres(1,ic2),nslt,nfinalov,
     -                nfinalrmsd,atw,atwsum,temp,chn,c2,indexov,indexov,
     -                noopt2d,0,indexrmsd,indexrmsd,rot,com1,com2,etot,
     -                etot2,iw0,devmax,devmaxnoopt,maxrec)
                    if (absdevmin  >  devmax) absdevmin=devmax
                    if (absdevmax  <  devmax) absdevmax=devmax
                    nframe2d=nframe2d+1
                    call progress_rep(0,nframe2d,nframesign)
                  
                
                ic21=1
                ncop2=0
c               print *,'NINCONF1,2=',ninconf1,ninconf2,
c    -            ' ILAST1,2=',ilast1,ilast2
              ## end if
            
            close (inpt2)
            numsel2=1
            ic11=1
            ncop1=0
            iverbose=0
          ## end if
        
        write (iw0,2020) trajnam(:ltrajnam),nframeread,nframe
        if (lpstitle  >  0) :
          call openps(iw1,xm,ym,title(1:76),76,pstitle,lpstitle,trajnam,
     -      ltrajnam,trajnam2,ltrajnam2,npspages,ipspage)
        ## end if
      elif (itraj  ==  1) :
        if (ireadtracks  >  0) go to 8006
        iopen_rep2=0
        nframeref=1
        nextconfsel2=2
        nframe2d=0
        nrep=nrep+1
        nframe=0
        nframeread=0
        ifail=0
        if (nrep  >  1  or  nframeref  >  1) :
          iverbose=0
          nrepinc=1
        ## end if
        ninp=n
        separatorchar='_'
        call opentraj(c,nrep+nrepinc,inpt,inptrajtyp,n,ntitltr,trtitle,
     -    inpcrdtyp,ifirst,ilast,increment,maxconf,ninconf,noutconf,
     -    natom,nfreeat,ifree,icntrl,0,mmctrajtyp,trajnam,ltrajnam,
     -    'trajectory',10,iconfsel,numsel,iverbose,notitprint,0,0,
     -    icellfound,notfnd,iconfirmname,iopen_rep,lentest,iw0,
     -    ianaltyp,mx2d,maxconfsel,maxrec)
        call save_traj_lim(ifirst,ilast,increment,1)
        call adjust_xtraj(xtraj,ifirst,ilast,increment,framefac,
     -    iadjust_xtraj)
        iopen_rep=1
        nerr_int=max0(1,maxconf/MAXDISTRN)
        if (numsel  >  0) nextconfsel=1
        if (ianaltyp  ==  16) :
c         Round of the max # of frames for plots
          maxtrajplot=ilast/increment
          call roundlimint(maxtrajplot,idiv,ndivdssp)
          maxtrajplot=idiv*ndivdssp
          print *,'maxtrajplot,idiv,ndivdssp=',
     -      maxtrajplot,idiv,ndivdssp
        ## end if
        notitprint=1
        if (ntitltr  >  0  and  title(1:4)  ==  '@#$%') :
          title=trtitle(1)
          if (trtitle(1)(1:21)  ==  'Created by DCD plugin')
     -      title=trtitle(2)
          call lastchar(title,ltitle,80)
        ## end if
        if (iw0  >  0  and  nframeref  ==  1) write (iw0,2006)
     -    mark0(1:lmark0),'Trajectory',trajnam(1:ltrajnam)
        if (n  !=  ninp) :
          write (6,2034) ninp,n
          write (iw0,2034) ninp,n
          if (nrep  ==  1  and  nframe  ==  0  and  nframeref  ==  1)
     -      call askstop(1)
        ## end if
        if (ibondtype  >  0)
     -    call printanchorlist(bondname,lbondname,ibondtype,ianchor,
     -      nanchor,indexa,indexov,nanchorr,nanchorn,line,index,iresno,
     -      isegno,segid4,inamcol1,inamcol2,irescol1,irescol2,iw0,
     -      maxrsd,maxrec)
        if (nrep  ==  1  and  nframe  ==  0  and  nframeref  ==  1):
c         PS files with trajectory name
          if (lpstitle  >  0) call openps(iw1,xm,ym,title(1:76),76,
     -      pstitle,lpstitle,trajnam,ltrajnam,trajnam2,ltrajnam2,
     -      npspages,ipspage)
        ## end if
        if (nframeref  ==  1) :
          if (ianaltyp  ==  16  and  ilast  ==  999999) :
            call getint('Last configuration in the trajectory',36,
     -        999999,1,0,ilast,0)
          ## end if
          nframetot=999999
          if (ilast  !=  999999) :
            nframetot=(ilast-ifirst)/increment+1
            nframesign=nframetot/10
          ## end if
          if (nframetot  >  maxframe  and 
     -      (ianaltyp  ==  12  or  ianaltyp  ==  21)) :
            if (ilast  !=  999999) :
              write (6,2074) maxframe,'structures',' '
              close (inpt)
              return
            else:
              write (6,2074) maxframe,'structures',' you may have to '
            ## end if
          ## end if
          if (ianaltyp  ==  23  and  increment  >  10  and 
     -        ilast  >  1000) :
            write (6,2021)
            call askstop(0)
          ## end if
          if (nrep  <=  1) :
            if (ianaltyp  ==  21  or  ianaltyp  ==  22)
     -         write (iw0,2012) inpfile(1:namleni)
            if (ianaltyp  ==  21) :
              if (isubcrm+ioverlay  >  0) :
                if (ioverlay  ==  0) write (iw0,2017)
                if (ioverlay  ==  1) write (iw0,2018) ' '
                if (ioverlay  ==  2)
     -            write (iw0,2018) ' by molecule'
                atwsum=0.0
                for ia in range(0, nslt):
                  atw(ia)=aw(iatnum(ia))
                  atwsum=atwsum+atw(ia)
                
              ## end if
              write (iw0,2085)
              iextra=0
              if (itmem  >  0)
     -          call askyn('Does the normal point to extracellular',
     -               38,1,1,iextra,0,0)
              for ihx in range(0, nhx):
                write (iw0,2010)
                call helixaxis(cres,nslt,iw0,calph0(1,1,ihx),
     -            axisdir0(1,ihx),axisini0(1,ihx),axisen0(1,ihx),
     -            helixcent0(1,ihx),perpvec0(1,1,ihx),camod(1,1,ihx),
     -            anglechangeref(1,ihx),circ(1,ihx),rn0(1,ihx),
     -            axfact(1,ihx),axtol,rot,rms,helixlen0(ihx),
     -            angles(1,ihx),decideb# end,nup,ndown,nrun,nnear,
     -            rcirchx(ihx),turnperres(ihx),anglesn(1,ihx),0,incrot,
     -            nrep,0,nreshx(ihx),icaahx(1,ihx),ihx,nhxres,
     -            idebughx,radtodeg,pi,MAXHX)
c               Establish helix direction(s)
                ihax=0
                projmax=0.0
                for k in range(0, 3):
                  if (abs(axisdir0(k,ihx))  >  projmax) :
                    projmax=abs(axisdir0(k,ihx))
                    ihax=k
                  ## end if
                
                axdirchar(ihx)=xyz(ihax)
                write (6,2129) ihx,axdirchar(ihx)
                indexaxhx(1,ihx)=ihax
                indexaxhx(2,ihx)=mod(ihax,3)+1
                indexaxhx(3,ihx)=mod(ihax+1,3)+1
                if (itmem  >  0) :
c                 Establish intra/extra cellular helix # end labeling
                  if (normhx  ==  ihax) :
                    memdir(ihx)=1
                    if (axisdir0(ihax,ihx)  <  0) memdir(ihx)=2
                    if (iextra  ==  0) memdir(ihx)=3-memdir(ihx)
                    tmchar=tmchars(memdir(ihx))
                  else:
                    memdir(ihx)=0
                    tmchar='n'
                  ## end if
                  anghx=acos(abs(axisdir0(normhx,ihx)))*radtodeg
                  write (6,2130) xyz(normhx),anghx
                  call quiz(ans,ihxtyp,tmchar,' ',0,'Helix type',10,0,
     -              5,6,0)
                  memdir(ihx)=ihxtyp-1
                ## end if
                if (memdir(ihx)  ==  0) :
                   write (iw0,2129) ihx,axdirchar(ihx)
                else:
                  write (iw0,2129) ihx,axdirchar(ihx),' (TM)'
                ## end if
c               Calculate reference turn angles
                for ir in range(2, (ireshx2(ihx)-ireshx1(ihx)+1)):
                  call angcomp(perpvec0(1,ir-1,ihx),axisdir0(1,ihx),
     -              perpvec0(1,ir,ihx),anglechangeref(ir,ihx))
                  if (anglechangeref(ir,ihx)  <  0.0)
     -              anglechangeref(ir,ihx)=anglechangeref(ir,ihx)+2.0*pi
                
c               Shift the helix so that the start is at the origin
                for ir in range(0, (ireshx2(ihx)-ireshx1(ihx)+1)):
                  call dvdif(calph0(1,ir,ihx),axisini0(1,ihx),
     -            calph0(1,ir,ihx))
                
              
              if (nhx  >  1) call multihelix(iw0,nhx,nhxres,radtodeg,c,
     -          icaahx,icbahx,icbreshx,maxrec,MAXHX,MAXNHX)
            ## end if
          ## end if
          write (iw0,2005)
          write (6,2081) maxrepconf
        ## end if
c       Open, read and analyze trajectory
        nntest=0
cx      print *,'Start scan nmc=',nmc
        ilastframe=0
        while (ninconf  <  ilast  and  ifail  ==  0)
c         Read a conformation
          call readtraj(inpt,inptrajtyp,mmctrajtyp,n,naslv,nwatr,
     -      nslt,trajformatname(inptrajtyp),ntitltr,trtitle,trajnam,
     -      ltrajnam,natom,nfreeat,ifree,icntrl,c,ninconf,
     -      noutconf,increment,inpcrdtyp,ietotsaved,etot,ifail,ifirst,
     -      ilast,iconfsel,numsel,maxrepconf*irepscan,nmc,lentest,1.0,
     -      maxconf,maxconfsel,maxrec)
          if (mod(n-nslt,naslv)  !=  0) :
            write (6,2168) n,nslt,naslv
            call askstop(-1)
          ## end if
          numsolv=(n-nslt)/naslv
c         if (ianaltyp  ==  23) irepscan=0
          icsel=0
          if (ifail  ==  0) :
            nframeread=nframeread+1
            call selectconf(numsel,ninconf,ifirst,increment,
     -        iconfsel,nextconfsel,nframeref,icsel,ifail,maxconfsel)
          else:
            ninconf=ilast
            nframetot=nframe
            if (nframeread  ==  0) go to 9002
          ## end if
          if (nframe  >  MAXFRAMES) write (6,2059) MAXFRAMES
          if (icsel  >  0) :
            nframe=nframe+1
c           print *,'Selected nmc=',nmc
            if (ianaltyp  ==   9  and  ioppbc  >=  0) :
              if (nframe  ==  1) :
                if (ifirstref  ==  0  and  ioppbc  >=  0) :
                  call pbcdist(c(1,ia1),c(1,ia2),ia1,ia2,cell,ncell,
     -              -iw0,1,img,drincr,dimg)
                  call zeroit(drimg,3)
                  call arrdiff(drimg,cell(1,img),drimg,3)
                ## end if
              elif (noboxinfoar  ==  0) :
c               Update box info (assuming isotropic box fluctuation)
                call updatecell(inptrajtyp,edge)
              ## end if
            ## end if
            if (nntest  ==  0  and  innlist  >  0) :
              call checknnlist(1,n,ineig,nneig,nerr,maxng)
              call comparetop(c,n,nneig,ineig,iatnum,innlist,nslt,naslv,
     -          cell,ncell,ioppbc,maxng,maxrec)
              nntest=1
            ## end if
            if (noprintconf  ==  0  and  ianaltyp  !=  12  and 
     -        (ianaltyp  !=  17  or  listbridge  ==  1)) :
              if (iw0  >  0) write (iw0,2173) nframe,ninconf
              if (iw1  >  0  and  extnam2  !=  '.ps ') :
                write (iw1,2173) nframe,ninconf
              elif (iw2  >  0  and  ianaltyp  !=  22) :
                if (analfile2(namleno2-2:namleno2)  !=  '.ps ')
     -            write (iw2,2173) nframe,ninconf
              ## end if
            ## end if
            if (itrajrot  >  0) call rotate_c(c,n,trajrot,c,
     -        'TRAJROT',7)
c           Configuration is ready, choose the analysis
            if (ianaltyp  ==   1) :
c-------------(T) Neighbor, bond, angle and torsion list
              call printbondlist(iangpr,itorpr,1,nslt,nslt,c,nneig,
     -          ineig,line,inamcol1,inamcol2,irescol1,irescol2,
     -          iresncol1,iresncol2,index,cv,indexa,maxng,radtodeg,
     -          iw0,maxrec)
            elif (ianaltyp  ==   2) :
c-------------(T) 1-4 statistics
              call stat14(c,n,nneig,ineig,index,line,nframe,pi,iatnm2,
     -          irescol1,irescol2,inamcol1,inamcol2,inpcrdtyp,ioins,iw0,
     -          maxrepconf,maxng,maxrec)
            elif (ianaltyp  ==   3) :
C-------------(T) defal group analysis
              call findfg(1,n,iatnum,mmtype,indexo,indexn,nhbneig,ih,
     -          ifgtyp,ibnd,nfg,6,nneig,nhneig,nneiga,ineig,iw0,
     -          inpcrdtyp,ioins,inamcol1,inamcol2,irescol1,irescol2,
     -          iresncol1,iresncol2,line,index,maxng,maxrec)
            elif (ianaltyp  ==   4) :
c-------------(T) Bond length statistics
              call bondlenstat(c,n,iatnum,nneig,ineig,iatnm2,nframe,
     -          nanos,ixlist,ialist,iw0,inpcrdtyp,ioins,inamcol1,
     -          inamcol2,irescol1,irescol2,line,index,radtodeg,
     -          maxrepconf,maxng,maxrec)
            elif (ianaltyp  ==   5) :
c-------------(T) Hydrogen-bond list
    write (77,*) 'ninconf,nframe=',ninconf,nframe,
     -            ' if,il=',ifirst,ilast
              call hblist(c,n,nslt,islvw,naslv,iatnum,ifchrg,isegno,
     -          nneig,nneiga,nhbneig,ineig,nhneig,nnneig,ncneig,nsneig,
     -          npneig,ixres,nframe,hbf0,angm0,molresflag,hblimfac,
     -          angmin,iw0,inamcol1,inamcol2,irescol1,irescol2,
     -          iresncol1,iresncol2,blankline,line,index,indexn,ibnd,
     -          indexo,isolvent,ibondprint,1,nosameseg,nframe,
     -          radtodeg,maxrepconf,maxng,maxbox,maxrsd,maxrec)
              if (ifailbond  ==  0) :
                call selectbond(ixres,nanchor,ianchor,indexa,ianchor2,
     -            iselfanc,iqfsel2,nhbneig,ineig,nbfound,nbresfound,
     -            nmc,nhbdist,rhbdist,'hydrogen',8,ianc_anc,nosameseg,
     -            isegno,c,it1,it2,nbonds,ifailbond,iw0,maxng,maxrec,
     -            mxbonds)
                if (maxbondf  <  nbonds) :
                  maxbondf=nbonds
                  nmcmaxbond=nframe
                ## end if
              else:
                go to 9005
              ## end if
            elif (ianaltyp  ==   6  or  ianaltyp  ==  34) :
c-------------(T) Hydrophobic bond or heavy atom contact list
              call nnlisthph_sltb(nslt,ianchor2,iselfanc,indexa,
     -          nneig, npneig,nhbneig,ineig,c,rhphmax,0,isegno,
     -          nosameseg,bondlab(ibondtype),4,ibnd,indexo,ifail,
     -          maxng,maxng,maxrec)
              if (ihphprint  ==  1)
     -          call hph_sltblist(ibondtype,c,nslt,nhbneig,ineig,iw0,
     -            is1,is2,inamcol1,inamcol2,irescol1,irescol2,
     -            iresncol1,iresncol2,line,index,bondname,lbondname,
     -            bondlab(ibondtype),4,maxng)
              if (ifailbond  ==  0) :
                call selectbond(ixres,nanchor,ianchor,indexa,ianchor2,
     -            iselfanc,0,nhbneig,ineig,nbfound,nbresfound,nmc,
     -            nhbdist,rhbdist,'hydrophobic',11,ianc_anc,nosameseg,
     -            isegno,c,it1,it2,nbonds,ifailbond,iw0,maxng,maxrec,
     -            mxbonds)
                if (maxbondf  <  nbonds) :
                  maxbondf=nbonds
                  nmcmaxbond=nframe
                ## end if
              else:
                go to 9005
              ## end if
          elif (ianaltyp  ==  41) :
c-----------(T) Mutually proximal heavy atom contact list
            call nnlistmpx(n,nanchorr,nanchorn,indexa,indexov,iatnum,
     -        nneig,nhbneig,c,isegno,ibnd,indexo,rmpxlim,itemp1,itemp2,
     -        itemp3,itemp4,indexn,temp,cv,ifail,maxng,maxrec)
            call mpxblist(n,itemp1,itemp2,cv,itemp3,itemp4,is1,is2,
     -        inamcol1,inamcol2,irescol1,irescol2,iresncol1,iresncol2,
     -        line,index,ixres,nhbdist,rhbdist,nbfound,nbresfound,
     -        nbonds,c,rmpxlim,ifailbond,iw0,MAXBONDS,maxrec)
              if (maxbondf  <  nbonds) :
                maxbondf=nbonds
                nmcmaxbond=nframe
              ## end if
              if (ifailbond  >  0) go to 9005
            elif (ianaltyp  ==   7) :
c-------------(T) Salt bridge list
              call nnlisthph_sltb(nslt,ianchor2,iselfanc,indexa,
     -          nneig,npneig,nhbneig,ineig,c,rsltbmax,1,isegno,
     -          nosameseg,bondlab(ibondtype),4,ibnd,indexo,ifail,
     -          maxng,maxng,maxrec)
              if (ibondprint  ==  1)
     -          call hph_sltblist(ibondtype,c,nslt,nhbneig,ineig,iw0,
     -            is1,is2,inamcol1,inamcol2,irescol1,irescol2,
     -            iresncol1,iresncol2,line,index,bondname,lbondname,
     -            bondlab(3),4,maxng)
              if (ifailbond  ==  0) :
                call selectbond(ixres,nanchor,ianchor,indexa,ianchor2,
     -            iselfanc,0,nhbneig,ineig,nbfound,nbresfound,nmc,
     -            nhbdist,rhbdist,'salt-bridge',11,ianc_anc,nosameseg,
     -            isegno,c,it1,it2,nbonds,ifailbond,iw0,maxng,maxrec,
     -            mxbonds)
                if (maxbondf  <  nbonds) :
                  maxbondf=nbonds
                  nmcmaxbond=nframe
                ## end if
              else:
                go to 9005
              ## end if
            elif (ianaltyp  ==   8) :
c-------------(T) Calculate residue-residue distances
              call rrdist(c,n,nres,ixresno,iresno,ixres,resnames,ifres,
     -          ilres,increst,incsolvrr,irefres1,irefres2,irefseg1,
     -          irefseg2,irefresinc,resdistlim,inegres1,inegres2,
     -          inegseg1,inegseg2,inegresinc,resapplim,listrefres,
     -          nrefres,nrefrange,listnegres,nnegres,nnegrange,iatnum,
     -          ignoreh,nframe,itypavg,iw0,iw1,iw2,irescol1,irescol2,
     -          inamcol1,inamcol2,is1,is2,irescount1,irescount2,
     -          irescount3,line,index,indexn,indexo,indexs,indexa,
     -          itemp1,itemp2,temp,cv,maxrec,maxrsd)
            elif (ianaltyp  ==   9) :
c-------------(T) Calculate a distance
              if (ifirstref  ==  0) :
                if (nframe  >  1) :
c                 See if c(1,ia1) has not been reset to the cell
                  call pbcdist(caprev1,c(1,ia1),ia1,ia2,cell,ncell,-iw0,
     -              nframe,img,drincr,dimg)
                  if (img  >  1) :
c                   PBC change found
                    call arrsum(drimg,cell(1,img),drimg,3)
                  ## end if
c                 See if c(1,ia2) has not been reset to the cell
                  call pbcdist(caprev2,c(1,ia2),ia1,ia2,cell,ncell,-iw0,
     -              nframe,img,drincr,dimg)
                  if (img  >  1) :
c                   PBC change found
                    call arrdiff(drimg,cell(1,img),drimg,3)
                  ## end if
                ## end if
                call trnsfr(caprev1,c(1,ia1),3)
                call trnsfr(caprev2,c(1,ia2),3)
                call arrdiff(c(1,ia2),c(1,ia1),dr,3)
                call arrsum(dr,drimg,dr,3)
                d=sqrt(scprod(dr,dr))
                write (iw0,2049) nframe,
     -            ia1,(c(k,ia1),k=1,3),ia2,(c(k,ia2),k=1,3)
              else:
                if (nframe  >  1) :
c                 See if c(1,ia2) has not been reset to the cell
                  call pbcdist(caprev2,c(1,ia2),ia1,ia2,cell,ncell,-iw0,
     -              nframe,img,drincr,dimg)
                  if (img  >  1) :
c                   PBC change found
                    call arrdiff(drimg,cell(1,img),drimg,3)
                  ## end if
                ## end if
                call trnsfr(caprev2,c(1,ia2),3)
                call arrdiff(c(1,ia2),caref,dr,3)
                call arrsum(dr,drimg,dr,3)
                dd=scprod(dr,dr)
                d=sqrt(dd)
                write (iw0,2048) nframe,ia2,(c(k,ia2),k=1,3)
              ## end if
              write (iw0,2054) nframe,d,dd,dr
              call trajlimtest(nframe,MAXFRAMES)
              res(1,nframe,1)=d
              res(2,nframe,1)=dr(idax)
              res(1,nframe,2)=dr(indexax(2))
              res(2,nframe,2)=dr(indexax(3))
              res(1,nframe,3)=dd
            elif (ianaltyp  ==  10) :
c-------------(T) Check for potentially unphysical contacts
              call checkunphys(c,nslt,n,naslv,islvw,iatnum,ifchrg,
     -          isegno,idelse:g,indexo,nneig,nneiga,nhbneig,ineig,nhneig,
     -          nnneig,ncneig,nsneig,npneig,ixres,line,irescol1,
     -          irescol2,inamcol1,inamcol2,index,nframe,innlist,
     -          molresflag,ioppbc,cell,ncell,edge,ixyzhex,molsltlim,
     -          nmolslt,hblimfac,angmin,ctfac,bondminfac,maxdist,iles,
     -          iw0,isltonly,ibnd,indexn,rlim,innlistread,nframe,
     -          radtodeg,nerrunphys,maxrepconf,maxng,maxbox,maxrsd,
     -          maxrec)
            elif (ianaltyp  ==  11) :
c-------------(T) Pseudorotation angle calculation
              call pseudorot(c,n,ix5,nmem,nneig,ineig,psr5,psr,npsr,
     -          psrtyp,incgen,zring,rring,sinpsrs,cospsrs,qpsr,qpsrs,
     -          qpsr2s,zav,zsq,zavs,zsqs,1,iw0,nframe,radtodeg,maxng,
     -          maxrec)
            elif (ianaltyp  ==  12) :
c-------------(T) Calculate Proline kinks
              if (nframe  >  maxframe) :
                write (6,2074) maxframe,'structures',' '
                go to 8002
              ## end if
              call prokinkcalcla(nrep,c,nslt,b# end,wobble,
     -          faceshift,rmsb,rmsa,
     -          iflatproline,ix5,radtodeg,MAXHX)
              call pseudorot(c,n,ix5,nmem,nneig,ineig,psr5,psr,npsr,
     -          'p',incgen,zring,rring,sinpsrs,cospsrs,qpsr,qpsrs,
     -          qpsr2s,zav,zsq,zavs,zsqs,0,6,nframe,radtodeg,maxng,
     -          maxrec)
              if (nrep  ==  1) call savekinkdat(nmem,b# end,wobble,
     -          faceshift,psr(2),radtodeg)
              if (nrep  <=  1) write (iw0,2072) ninconf,b# end,wobble,
     -          faceshift,wobble-faceshift,psr(2),rmsb,rmsa
            elif (ianaltyp  ==  13) :
c-------------(T) Hydropathy labeling
              call hydropathylist(n,nslt,ixres,resnames,cv,ihydtyp,
     -          maxrsd,maxrec)
              call writeconf(iw0,inpcrdtyp,inpcrdtyp,inpcrdtyporg,n,n,
     -          nslt,naslv,islvw,iasv,namesv,qsv,pflsv,1,iwhead,0,
     -          iatnum,ifchrg,nframe,innlist,c,rprox,cv,ixres,iresno,
     -          atnames,resnames,segnames,charge,isegno,altcol,inscol,
     -          ninsres,marker,ntitlin,ntitlinw,title,ireseq,
     -          iresnrestart,iresidrestart,nneig,nneiga,nhbneig,ineig,
     -          nhneig,nnneig,ncneig,nsneig,npneig,numres,numslv,
     -          resnamslv,line,blankline,mmtype,ibnd,index,indexn,
     -          indexo,2,molresflag,irescount3,itemp1,hblimfac,angmin,0,
     -          1,1,1,0,1,iqspaceask,13,1,0.0,0,0,0,keeprem,iwriteatsym,
     -          radtodeg,maxrepconf,maxng,maxrsd,maxrec)
            elif (ianaltyp  ==  14) :
c-------------(T) Circular variance labeling of solute ats and solvent molecs
              call cvlist(c,n,nslt,nsltref_f,nsltref_l,naslv,islvrep,
     -          icvtyp,rcut_cv,rprox,cv,itemp1,c2,temp,indexn,indexo,
     -          indexs,chn,isortslv,cvlim)
              call writeconf(iw0,inpcrdtyp,inpcrdtyp,inpcrdtyporg,n,n,
     -          nslt,naslv,islvw,iasv,namesv,qsv,pflsv,1,iwhead,0,
     -          iatnum,ifchrg,nframe,innlist,c,rprox,cv,ixres,iresno,
     -          atnames,resnames,segnames,charge,isegno,altcol,inscol,
     -          ninsres,marker,ntitlin,ntitlinw,title,ireseq,
     -          iresnrestart,iresidrestart,nneig,nneiga,nhbneig,ineig,
     -          nhneig,nnneig,ncneig,nsneig,npneig,numres,numslv,
     -          resnamslv,line,blankline,mmtype,ibnd,index,indexn,
     -          indexo,1,molresflag,irescount3,itemp1,hblimfac,angmin,0,
     -          1,1,1,0,3,iqspaceask,14,1,0.0,0,0,0,keeprem,iwriteatsym,
     -          radtodeg,maxrepconf,maxng,maxrsd,maxrec)
            elif (ianaltyp  ==  16) :
c-------------(T) DSSP secondary structure assignment
              call dssp(c,1,n,nslt,line,index,inamcol1,inamcol2,
     -          iresncol1,iresncol2,nneig,ineig,nneiga,ibnd,indexn,
     -          indexo,npneig,nsneig,nnneig,dssplab,idistdssp,chn,c2,cv,
     -          indexs,indexa,nsse,itypsse,ifsse,ilsse,0,nframe,iw0,
     -          iwrose,nframe,ifail,dtodeg,maxrepconf,maxng*10,
     -          maxrec/10,maxng*20,maxrec/20,200,maxrsd,maxrec)
              if (ifail  >  0) go to 8002
              call plotdssp(iw1,ifsse,ilsse,itypsse,nsse,ninconf,
     -          maxtrajplot,framefac,ndivdssp,10,title(1:76),76,
     -          xtrajlab,11,'N(res)',6,nframe,ifrdssp,ilrdssp,indexdel)
            elif (ianaltyp  ==  17) :
c-------------(T) Hydrogen-bond bridge analysis
              call hblist(c,n,nslt,islvw,naslv,iatnum,ifchrg,isegno,
     -          nneig,nneiga,nhbneig,ineig,nhneig,nnneig,ncneig,nsneig,
     -          npneig,ixres,nframe,hbf0,angm0,molresflag,hblimfac,
     -          angmin,iw0,inamcol1,inamcol2,irescol1,irescol2,
     -          iresncol1,iresncol2,blankline,line,index,indexn,ibnd,
     -          indexo,1,0,1,nosameseg,nframe,radtodeg,
     -          maxrepconf,maxng,maxbox,maxrsd,maxrec)
              call hbbridge(nanchor,ianchor,indexa,ianchor2,iselfanc,
     -          lpath,nbridgetype,ibridgetype,maxbridgemem,n,nhbneig,
     -          ineig,indexn,indexo,indexs,ixres,resnames,brslv,nabr,
     -          nrescol,nmc,ifailhbr,listbridge,iw0,maxng,MAXBRIDGELEN,
     -          MAXBRIDGETYPE,minbridgelenprint,maxrsd,maxrec)
              if (ifailhbr  ==  1) go to 8002
            elif (ianaltyp  ==  18) :
c-------------(T) Ramachandran plot
              call ramachandran(c,nslt,index,line,nframe,pi,nresfound,
     -          iresno,ixres,irescol1,irescol2,inamcol1,iw0,
     -          maxng,maxrec)
              call ramachandranplot(nresfound,iw1,xm,iallrama)
            elif (ianaltyp  ==  19) :
c-------------(T) Torsion dial plots
            call torsiondials(c,nslt,ntorsel,ixtor1234,index,line,pi,
     -        irescol1,irescol2,inamcol1,inamcol2,iw0,maxrec)
            elif (ianaltyp  ==  20) :
c-------------(T) Delphi map labeling
              call delphilabel(c,n,nslt,xstart,ystart,zstart,gx,gy,gz,
     -          cv)
              call writeconf(iw0,inpcrdtyp,inpcrdtyp,inpcrdtyporg,n,n,
     -          nslt,naslv,islvw,iasv,namesv,qsv,pflsv,1,iwhead,0,
     -          iatnum,ifchrg,nframe,innlist,c,rprox,cv,ixres,iresno,
     -          atnames,resnames,segnames,charge,isegno,altcol,inscol,
     -          ninsres,marker,ntitlin,ntitlinw,title,ireseq,
     -          iresnrestart,iresidrestart,nneig,nneiga,nhbneig,ineig,
     -          nhneig,nnneig,ncneig,nsneig,npneig,numres,numslv,
     -          resnamslv,line,blankline,mmtype,ibnd,index,indexn,
     -          indexo,1,molresflag,irescount3,itemp1,hblimfac,angmin,0,
     -          1,1,1,0,0,iqspaceask,20,1,0.0,0,0,0,keeprem,iwriteatsym,
     -          radtodeg,maxrepconf,maxng,maxrsd,maxrec)
            elif (ianaltyp  ==  21) :
c-------------(T) Helix directions
cd77              write (77,*) ' ------------ Nframe=',nframe,' -----------'
              if (nframe  ==  1) :
                for ihx in range(0, nhx):
                  iadssp1=ifres(ixres(icaahx(1,ihx)))
                  iadssp2=ilres(ixres(icaahx(nreshx(ihx),ihx)))
                  call dssp(cres,iadssp1,n,iadssp2,line,index,
     -              inamcol1,inamcol2,iresncol1,iresncol2,nneig,ineig,
     -              nneiga,ibnd,indexn,indexo,npneig,nsneig,nnneig,
     -              dssplab,idistdssp,chn,c2,cv,indexs,indexa,nsse,
     -              itypsse,ifsse,ilsse,ireshx1(ihx)-1,nframe,0,0,0,
     -              ifail,radtodeg,maxrepconf,maxng*10,maxrec/10,
     -              maxng*20,maxrec/20,200,maxrsd,maxrec)
                  lhelixcklab=9+11
                  write (helixcklab,2160) 'Reference h',ihx
                  call checkforhelix(ireshx2(ihx)-ireshx1(ihx)+1,
     -              dssplab,indexn,6,hxoklab,ihxok,helixcklab,
     -              lhelixcklab,maxrsd,maxrec)
                
              ## end if
              if (nrep  ==  1)
     -          call soluteoverlay(isubcrm,ioverlay,nslt,nsegslt,c,cres,
     -            chn,c2,crmslt0,crmslt,atw,temp,cv,rprox,molsltlim,
     -            indexs,idebughx,iw0,maxrsd,maxrec)
              for ihx in range(0, nhx):
                call helixcomp(nslt,nreshx(ihx),calph0(1,1,ihx),
     -            axisdir0(1,ihx),perpvec0(1,1,ihx),axisini0(1,ihx),
     -            axisen0(1,ihx),helixcent0(1,ihx),rn0(1,ihx),
     -            camod(1,1,ihx),axfact,calph(1,1,ihx),axisdir(1,ihx),
     -            perpvec(1,1,ihx),axisini(1,ihx),axis# end(1,ihx),
     -            helixcent(1,ihx),anglechange(1,ihx),
     -            anglechangeref(1,ihx),axtol,iw0,torsion,
     -            rotation,nrep,indexaxhx(1,ihx),incrot,ireorienthx,
     -            idebughx,radtodeg,pi,c2,nreshx(ihx),icaahx(1,ihx),ihx,
     -            nhxres,MAXHX,maxrec)
                if (idsspcheck  ==  1) :
                  iadssp1=ifres(ixres(icaahx(1,ihx)))
                  iadssp2=ilres(ixres(icaahx(nreshx(ihx),ihx)))
                  call dssp(c,iadssp1,n,iadssp2,line,index,inamcol1,
     -              inamcol2,iresncol1,iresncol2,nneig,ineig,nneiga,
     -              ibnd,indexn,indexo,npneig,nsneig,nnneig,dssplab,
     -              idistdssp,chn,c2,cv,indexs,indexa,nsse,itypsse,
     -              ifsse,ilsse,ireshx1(ihx)-1,nframe,0,0,0,ifail,
     -              radtodeg,maxrepconf,maxng*10,maxrec/10,maxng*20,
     -              maxrec/20,200,maxrsd,maxrec)
                  lhelixcklab=9+1
                  write (helixcklab,2160) 'H',ihx
                  call checkforhelix(ireshx2(ihx)-ireshx1(ihx)+1,
     -              dssplab,indexn,iw0,hxoklab,ihxok,helixcklab,
     -              lhelixcklab,maxrsd,maxrec)
                  nhelixok(ihxok)=nhelixok(ihxok)+1
                ## end if
              
              if (nhx  >  1) call multihelix(iw0,nhx,nhxres,radtodeg,c,
     -          icaahx,icbahx,icbreshx,maxrec,MAXHX,MAXNHX)
            elif (ianaltyp  ==  22) :
c-------------(T) 1-D RMSD and residue RMS fluctuations
              if (nframe  ==  1  and  inputref  ==  0)
     -          call trnsfr(cres,c,3*nslt)
              nframeref=0
              call rmsd(cres,c,nslt,nfinalov,nfinalrmsd,atw,atwsum,temp,
     -          chn,c2,indexov,indexov,0,1,indexrmsd,indexrmsd,rot,com1,
     -          com2,etot,etot2,iw0,devmax,devmaxnoopt,maxrec)
              for ia in range(0, nfinalrmsd):
                for k in range(0, 3):
                  cdp(k,ia)=cdp(k,ia)+c2(k,indexrmsd(ia))
                  cdp2(k,ia)=cdp2(k,ia)+c2(k,indexrmsd(ia))**2
                
c               write (iw0,*) ia,' indexrmsd(ia)=',indexrmsd(ia)
              
              if (nresslt  <=  MAXDISTR  and  nocontigrmsd  ==  0) :
                call rmsf(chn,c2,atw,nfinalrmsd,ixres,indexrmsd,
     -            rmsfsum,nresslt,nslt)
c               write (77,8791) nframe,(rmsfsum(i),i=1,nresslt)
c8791           format(' NFRAME=',i4,' rmsfsum:',/,(10f8.2))
                if (mod(nframe,nerr_int)  ==  0  or 
     -              nframe  ==  maxconf) :
                  ns=nframe/nerr_int
                  nframes_err(ns)=nframe
                  for ir in range(0, nresslt):
                    distrerr(ns,ir)=rmsfsum(ir)
                  
                ## end if
              ## end if
c-------------2D RMSD and cross RMSD calculations done separately (in blocks)
            elif (ianaltyp  ==  25) :
c-------------(T) Residue correlation calculation
              call residcorr(c,c,nslt,indexa,ncorr,nframe)
            elif (ianaltyp  ==  26) :
c-------------(T) Atom-atom distribution calculation
              if (iclusterdist  ==  0) :
                call pairdistcalc(c,nslt,npairs,listpairdist,
     -            pairdistsum,pairdistsum2,pairdistwsum,npairdist,
     -            pairdistminmax,pairgrid,iw0,MAXDDBIN,MAXDDISTR)
              else:
                call clusterdistcalc(c,nslt,npairs,iclustermem,
     -            ifstclst1,ifstclst2,ilstclst2,pairdistsum,
     -            pairdistsum2,pairdistwsum,npairdist,pairdistminmax,
     -            pairgrid,iw0,MAXDDBIN,MAXDDISTR,MAXCDLIST)
              ## end if
            elif (ianaltyp  ==  27) :
c-------------(T) Solvation shell calculation
              call volcalc(nrand,c,nslt,isegno,iatnum,nsegslt,molsltlim,
     -          itemp1,cres,chn,c2,atw,vmac,vshell,vint,vmacsd,vshellsd,
     -          vintsd,vfrstsh,vfrstshsd,rsolv,ih,indexn,indexo,indexs,
     -          indexa,indexdel,nframe,iw0,levout,maxrsd,maxrec)
              call trajlimtest(nframe,MAXFRAMES)
              res(1,nframe,1)=vshell
              res(2,nframe,1)=vfrstsh
              res(1,nframe,2)=vmac
              res(2,nframe,2)=vint
              res(1,nframe,3)=vshellsd
              res(2,nframe,3)=vfrstshsd
              res(1,nframe,4)=vmacsd
              res(2,nframe,4)=vintsd
            elif (ianaltyp  ==  28) :
c-------------(T) Principal axis calculation
              if (nframe  ==  1  and  inputref  ==  1)
     -          call princax(cres,c2,atw,temp,n,indexdel,evecs0,evals0,
     -            iw0,inputref,1,radtodeg,0)
              call princax(c,c2,atw,temp,n,indexdel,evecs0,evals0,iw0,
     -          inputref,0,radtodeg,0)
            elif (ianaltyp  ==  29) :
c-------------(T) Radius and dipole calculation
              call molrad(c,indexa,ngyrats,iw0,MAXREC)
              call celldipole(c,n,nslt,indexa,ngyrats,charge,icharges,
     -          atw,iw0,1)
            elif (ianaltyp  ==  32) :
c-------------(T) Angle dial plots
              call angledials(c,nslt,nangsel,ixtor1234,index,line,pi,
     -          irescol1,irescol2,inamcol1,inamcol2,iw0,maxrec)
            elif (ianaltyp  ==   33) :
c-------------(T) Calculate molecule-molecule distance matrix
                call mmdist(c,n,atw,iatnum,nmolslt,molsltlim,c2,temp,
     -            ignoreh,iw0,iw1)
            elif (ianaltyp  ==  39) :
c-------------(T) Atom-atom distance SD calculation
              call atomdist_sd(c,nslt,ianchor,nanchor,nframe)
            elif (ianaltyp  ==  40) :
c-------------(T) Solvent filtering
              if (ifilttyp  ==  5)
     -          call nnlist(n,islvw,naslv,n,iatnum,ifchrg,c,nneig,
     -          nneiga,nhbneig,ineig,nhneig,nnneig,ncneig,nsneig,
     -          npneig,line,irescol1,irescol2,inamcol1,inamcol2,index,
     -          nconfig,innlist,molresflag,hblimfac,angmin,ihbondcalc,
     -          ibnd,indexo,isegno,ixres,maxrepconf,0,nframe,radtodeg,
     -          0,maxbox,maxng,maxrsd,maxrec)
              call filterslv(c,iatnum,nslt,n,naslv,numsolv,molsltlim,
     -          nsegslt,index,nframe,ifilttyp,intftyp,rsltmax,rcvmax,
     -          cvmin,ang12min1,ang12min2,r12max,cvrlim,rneigh,mol11,
     -          mol12,mol21,mol22,spacing,itemp1,iarepslv,numsolvleft,
     -          nhbneig,ineig,ixres,nrecdel,ierr,iw0,maxrec,maxng)
              numsolvsum=numsolvsum+numsolvleft
              if (numsolvmin  >  numsolvleft) numsolvmin=numsolvleft
              if (numsolvmax  <  numsolvleft) numsolvmax=numsolvleft
              if (numsolvleft  >  0) :
                namleno=namlen_root
                if (nframe  >  1  and  ionefile  ==  0) :
                  namleno=namlen_root+1
                  analfile(namleno:namleno)='.'
                  namleno=namleno+1
                  call writeint(analfile,namleno,nframe,lenc)
                  namleno=namleno-1
                ## end if
                if (nframe  ==  1  or  ionefile  ==  0) :
                  analfile1(namleno+1:namleno+lenext)=crdexti(1:lenext)
                  namleno=namleno+lenext
                  call openfile(45,0,'FILT',4,'new',analfile1,namleno,
     -              notfnd,0,1,1,0,0)
                ## end if
                call writeout(45,inpcrdtyp,inpcrdtyp,line,index,isegno,
     -            n-nrecdel,marker,1,0,ntitlin,ntitlinw,title,
     -            blankline,1,1,0,0.0,0,0,0,keeprem,iwriteatsym,iatnum,
     -            maxrec)
                if (ionefile  ==  0) close (45)
              else:
                write (6,2070) nframe
              ## end if
            ## end if
            if (nframetot  >  nframetotmin) :
c             Print progress report
              if (nframesign  >  0) :
                call progress_rep(nframe,nframe2d,nframesign)
              else:
                nframetot=(ilast-ifirst+1)/increment
c               if (ianaltyp  ==  23)
c    -            nframetot=nframetot*(nframetot-1)/2
                nframesign=nframetot/10
              ## end if
            ## end if
          ## end if
          if ((ninconf-ifirst)/increment  >  2) iprintpk=0
          ifirstscan=0
        
c       Trajectory scan done
8002    close (inpt)
        nmw=0
        if (iw0  >  0) :
          write (iw0,2085)
          write (iw0,2020) trajnam(:ltrajnam),nframeread,nframe,' ',
     -      ifirst-1
          if (ianaltyp  ==  24) write (iw0,2065) nframeref
        ## end if
      ## end if
8006  if (itraj  ==  1) :
        call testconst(0,1,2,0.0,1.0,2.0,6,nfail,1,'ANPO')
c       Post scan activities
        if (ntitlin  ==  0  and  ltitle  ==  0) :
c         Use trajectory title
          title=trtitle(1)
        ## end if
        if (iw2  >  0  and  ianaltyp  !=  22) :
          if (analfile2(namleno2-2:namleno2)  !=  '.ps ')
     -      write (iw2,2085)
        ## end if
        if (ianaltyp  ==  12) :
          write (title,2069) irespro,nra,nrb
          call dialps(iw1,prokinklab,lprokinklab,title,ltitle,resrange,
     -      22,5,ndials,ndprow,nframe,pi,ipspage,0,nfravgd,0,0,0,
     -      iconndial,ioutpr,mappdf,ipdfgrd)
          write (6,2073) analfile1(1:namleno1)
          call trajstat(iw0,ndials,5,prokinklab,lprokinklab,0,1,0,
     -      rmsdlab,lrmsdlab,corr12,0,0,0,' ',1,radtodeg)
        elif (ianaltyp  ==   9  or  ianaltyp  ==  19  or 
     -           ianaltyp  ==  21  or  ianaltyp  ==  22  or 
     -           ianaltyp  ==  23  or  ianaltyp  ==  24  or 
     -           ianaltyp  ==  26  or  ianaltyp  ==  27  or 
     -           ianaltyp  ==  28  or  ianaltyp  ==  29  or 
     -           ianaltyp  ==  32  or  ianaltyp  ==  39) :
          if (ianaltyp  ==  9) :
            irhx0=3+indexax(1)
            irhx1=3+indexax(2)
            irhx2=3+indexax(3)
            call arminmax2(res(1,1,1),1,nframe,2,armin1,armax1,
     -        armin2,armax2,0,2)
            call roundlim(armax1,y1div,ny1div)
            call plot2fun(iw1,2,xtraj,res(1,1,1),res(1,1,1),nframe,
     -        0.0,0.0,0,0.0,y1div,ny1div, 0.0,y2div,ny2div,title,80,
     -        atomdist,22,xtrajlab,11 ,'D',1,helixrlab(irhx0),
     -        lhelixrlab(irhx0),1,0,6,2,1,0,1,0,0,
     -        ipspage,1,1,0)
            call arminmax2(res(1,1,3),1,nframe,1,armin1,armax1,
     -        armin2,armax2,0,2)
            call roundlim(armax1,y1div,ny1div)
            call plot2fun(iw1,1,xtraj,res(1,1,3),res(1,1,3),nframe,
     -        0.0,0.0,0,0.0,y1div,ny1div, 0.0,y2div,ny2div,title,80,
     -        atomdist,22,xtrajlab,11 ,'D^2',3,
     -        'Total displacement square',25,1,0,6,
     -        2,1,0,1,0,0,ipspage,1,1,0)
            call arminmax2(res(1,1,2),1,nframe,2,armin1,armax1,
     -        armin2,armax2,0,2)
            ix=4.0*armin1
            if (armin1  <  0.0) ix=ix-1
            xmn=float(ix)/4.0
            call roundlim(armax1-xmn,xdv,nxdv)
            iy=4.0*armin2
            if (armin2  <  0.0) iy=iy-1
            ymn=float(iy)/4.0
            call roundlim(armax2-ymn,ydv,nydv)
            call plot2fun(iw1,2,xtraj,res(1,1,2),res(1,1,2),nframe,
     -        0.0,0.0,00,xmn,xdv,nxdv,ymn,ydv,nydv,title,80,atomdist,22,
     -        xtrajlab,11 ,helixrlab(irhx1),lhelixrlab(irhx1),
     -        helixrlab(irhx2),lhelixrlab(irhx2),1,
     -        0,6,2,1,0,1,0,0,ipspage,1,1,0)
            plotdescr='Progression of distance vector projected to the '
            plotdescr=plotdescr(1:48)//xyz(indexax(2))//'-'
            plotdescr=plotdescr(1:50)//xyz(indexax(3))//' plane'
            call plot2d(iw1,res(1,1,2),nframe,nfravgt,xmn,xdv,nxdv,
     -        ymn,ydv,nydv,title,80,plotdescr,56,atomdist,22,
     -        xyz(1),1,xyz(2),1,xtrajlab,11,0,0,6,1,ipspage,0,MAXFRAMES)
          elif (ianaltyp  ==  21) :
            noopen=0
            noclose=1
            now6=0
            lresrange=22
            if (nhx  >  1) :
              now6=1
              lresrange=29
            ## end if
            for ihx in range(0, nhx):
              write (resrange(1:22),2087) ireshx1(ihx),ireshx2(ihx)
              shiftlab(20:20)=axdirchar(ihx)
              if (nhx  >  1) write (resrange(23:29),2163) ihx
              incrhx=(ihx-1)*nhxres
              call dialps(iw1,helixang,lhelixang,title,ltitle,resrange,
     -          lresrange,5,ndials,ndprow,nframe,pi,ipspage,9,nfravgd,
     -          incrhx,noopen,1,iconndial,ioutpr,mappdf,ipdfgrd)
              noopen=1
              if (idsspcheck  >  0) :
                write (6,*)
                write (iw0,*)
                for k in range(0, 3):
                  if (nhx  ==  1) write (6,2089) hxoklab(k),nhelixok(k)
                  write (iw0,2089) hxoklab(k),nhelixok(k)
                
              ## end if
              for k in range(0, 3):
                printrlab(k)=helixrlab(k)
                lprintrlab(k)=lhelixrlab(k)
              
              printrlab(4)=helixrlab(3+indexaxhx(1,ihx))
              lprintrlab(4)=lhelixrlab(3+indexaxhx(1,ihx))
              for k in range(0, 2):
                printrlab(4+k)=helixrlab(3+indexaxhx(k+1,ihx))
                lprintrlab(4+k)=lhelixrlab(3+indexaxhx(k+1,ihx))
                printrlab(6+k)=helixrlab(8+indexaxhx(k+1,ihx))
                lprintrlab(6+k)=lhelixrlab(8+indexaxhx(k+1,ihx))
                printrlab(8+k)=helixrlab(11+indexaxhx(k+1,ihx))
                lprintrlab(8+k)=lhelixrlab(11+indexaxhx(k+1,ihx))
              
              write (iw0,2172) ihx
              call trajstat(iw0,ndials,6,helixang,lhelixang,10,10,7,
     -          printrlab,lprintrlab,corr12,incrhx,itorcorr,now6,' ',1,
     -          radtodeg)
              call arminmax2(res(1,1,incrhx+7),1,nframe,2,armin1,armax1,
     -          armin2,armax2,0,2)
              iy1=armin1
              y1mn=iy1
              id1=armax1-iy1+1
              if (mod(id1-iy1,2)  ==  1  or  id1. eq. 0) id1=id1+1
              y1div=float(id1)/10.0
              iy2=armin2
              if (armin2  <  0.0) iy2=iy2-1
              y2mn=iy2
              call roundlim(armax2-y2mn,y2div,ny2div)
              call plot2fun(iw1,2,xtraj,res(1,1,incrhx+7),
     -          res(1,1,incrhx+7),nframe,0.0,0.0,00,y1mn,y1div,10,y2mn,
     -          y2div,ny2div,title,80,resrange,lresrange,xtrajlab,11,
     -          'Helix length',12,'Radius of fitted circle',23,1,
     -          0,6,2,1,0,1,0,0,ipspage,1,1,0)
              call arminmax2(res(1,1,incrhx+8),1,nframe,2,armin1,armax1,
     -          armin2,armax2,0,2)
              iy1=armin1
              y1mn=iy1
              call roundlim(armax1-y1mn,y1div,ny1div)
              iy2=4.0*armin2
              if (armin2  <  0.0) iy2=iy2-1
              y2mn=float(iy2)/4.0
              call roundlim(armax2-y2mn,y2div,ny2div)
              call plot2fun(iw1,2,xtraj,res(1,1,incrhx+8),
     -          res(1,1,incrhx+8),nframe,0.0,0.0,00,y1mn,y1div,ny1div,
     -          y2mn,y2div,ny2div,title,80,resrange,lresrange,xtrajlab,
     -          11,'Total helix displacement',24,shiftlab,30,1,
     -          0,6,2,1,0,1,0,0,ipspage,1,1,0)
c             Plot the progress of the helix center in the plane of the other ax
              irhx1=3+indexaxhx(2,ihx)
              irhx2=3+indexaxhx(3,ihx)
              call arminmax2(res(1,1,incrhx+9),1,nframe,2,armin1,armax1,
     -          armin2,armax2,0,2)
              ix=4.0*armin1
              if (armin1  <  0.0) ix=ix-1
              xmn=float(ix)/4.0
              call roundlim(armax1-xmn,xdv,nxdv)
              iy=4.0*armin2
              if (armin2  <  0.0) iy=iy-1
              ymn=float(iy)/4.0
              call roundlim(armax2-ymn,ydv,nydv)
              call plot2fun(iw1,2,xtraj,res(1,1,incrhx+9),
     -          res(1,1,incrhx+9),nframe,0.0,0.0,00,xmn,xdv,nxdv,ymn,
     -          ydv,nydv,title,80,resrange,lresrange,xtrajlab,11,
     -          helixrlab(irhx1),lhelixrlab(irhx1),helixrlab(irhx2),
     -          lhelixrlab(irhx2),1,0,6,2,1,0,1,0,0,
     -          ipspage,1,1,0)
              plotdescr='Helix center track in the '//
     -          xyz(indexaxhx(2,ihx))//'-'
              plotdescr=plotdescr(1:28)//xyz(indexaxhx(3,ihx))//' plane'
              call plot2d(iw1,res(1,1,incrhx+9),nframe,nfravgt,xmn,
     -          xdv,nxdv,ymn,ydv,nydv,title,ltitle,plotdescr,35,
     -          resrange,lresrange,xyz(indexaxhx(2,ihx)),1,
     -          xyz(indexaxhx(3,ihx)),1,xtrajlab,11,0,0,6,0,ipspage,1,
     -          MAXFRAMES)
c             Plot the move of the helix start in the plane of the other axes
              irhx1=8+indexaxhx(2,ihx)
              irhx2=8+indexaxhx(3,ihx)
              call arminmax2(res(1,1,incrhx+10),1,nframe,2,armin1,
     -          armax1,armin2,armax2,0,2)
              ix=4.0*armin1
              if (armin1  <  0.0) ix=ix-1
              xmn=float(ix)/4.0
              call roundlim(armax1-xmn,xdv,nxdv)
              iy=4.0*armin2
              if (armin2  <  0.0) iy=iy-1
              ymn=float(iy)/4.0
              call roundlim(armax2-ymn,ydv,nydv)
              call plot2fun(iw1,2,xtraj,res(1,1,incrhx+10),
     -          res(1,1,incrhx+10),nframe,0.0,0.0,00,xmn,xdv,nxdv,ymn,
     -          ydv,nydv,title,80,resrange,lresrange,xtrajlab,11,
     -          helixrlab(irhx1),lhelixrlab(irhx1),helixrlab(irhx2),
     -          lhelixrlab(irhx2),1,0,6,2,1,0,1,0,0,
     -          ipspage,1,1,0)
              plotdescr(7:12)='start '
              call plot2d(iw1,res(1,1,incrhx+10),nframe,nfravgt,
     -          xmn,xdv,nxdv,ymn,ydv,nydv,title,ltitle,plotdescr,35,
     -          resrange,lresrange,xyz(indexaxhx(2,ihx)),1,
     -          xyz(indexaxhx(3,ihx)),1,xtrajlab,11,0,0,6,0,ipspage,1,
     -          MAXFRAMES)
c             Plot the progress of the helix # end in the plane of the other axes
              irhx1=11+indexaxhx(2,ihx)
              irhx2=11+indexaxhx(3,ihx)
              call arminmax2(res(1,1,incrhx+11),1,nframe,2,armin1,
     -          armax1,armin2,armax2,0,2)
              ix=4.0*armin1
              if (armin1  <  0.0) ix=ix-1
              xmn=float(ix)/4.0
              call roundlim(armax1-xmn,xdv,nxdv)
              iy=4.0*armin2
              if (armin2  <  0.0) iy=iy-1
              ymn=float(iy)/4.0
              call roundlim(armax2-ymn,ydv,nydv)
              call plot2fun(iw1,2,xtraj,res(1,1,incrhx+11),
     -          res(1,1,incrhx+11),nframe,0.0,0.0,00,xmn,xdv,nxdv,ymn,
     -          ydv,nydv,title,80,resrange,lresrange,xtrajlab,11 ,
     -          helixrlab(irhx1),lhelixrlab(irhx1),helixrlab(irhx2),
     -          lhelixrlab(irhx2),1,0,6,2,1,0,1,0,0,
     -          ipspage,1,1,0)
              plotdescr(7:12)=' # end  '
              call plot2d(iw1,res(1,1,incrhx+11),nframe,nfravgt,
     -          xmn,xdv,nxdv,ymn,ydv,nydv,title,ltitle,plotdescr,35,
     -          resrange,lresrange,xyz(indexaxhx(2,ihx)),1,
     -          xyz(indexaxhx(3,ihx)),1,xtrajlab,11,0,0,6,0,ipspage,1,
     -          MAXFRAMES)
c             Plot the progress of the helix plane normal's projection
              call arminmax2(res(1,1,incrhx+12),1,nframe,2,armin1,
     -          armax1,armin2,armax2,0,2)
              xmin=amin1(armin1,armin2,-armax1,-armax2)
              call roundlim(-10.0*xmin,xdv,nxdv)
              xdv=xdv/5.0
              nxdv=nxdv/2
              xmn=-xdv*nxdv
              nxdv=2*nxdv
              linein(1:lresrange)=resrange
              linein(lresrange+1:lresrange+33)=
     -          ' Helix axis is in the X direction'
              plotdescr='Track of the projection of the normal to the '
     -          //xyz(indexaxhx(2,ihx))//'-'//xyz(indexaxhx(3,ihx))//
     -          ' plane fitted to the helix b# end'
              call plot2d(iw1,res(1,1,incrhx+12),nframe,nfravgt,
     -          xmn,xdv,nxdv,xmn,xdv,nxdv,title,ltitle,plotdescr,79,
     -          linein,lresrange+33,xyz(1),1,xyz(2),1,xtrajlab,11,0,0,
     -          6,0,ipspage,1,MAXFRAMES)
              for i in range(0, nframe):
                scp_ax=res(1,1,incrhx+13)*res(1,i,incrhx+13)+
     -            res(2,1,incrhx+13)*res(2,i,incrhx+13)+
     -            res(1,1,incrhx+14)*res(1,i,incrhx+14)
                cv(i)=(180./3.141592)*dacoscheck(ddd,scp_ax,0,6,
     -            'HELIXAN')
              
              call arminmax2(cv,1,nframe,1,armin1,armax1,armin2,armax2,
     -          0,1)
              iy1=armin1
              y1mn=iy1
              call roundlim(armax1-y1mn,y1div,ny1div)
              iprt=0
              call plot2fun(iw1,1,xtraj,cv,cv,nframe,0.0,0.0,00,y1mn,
     -          y1div,ny1div,y2mn,y2div,ny2div,title,80,resrange,
     -          lresrange,xtrajlab,11,
     -          'Angle between start and current axes',36,shiftlab,30,
     -          1,iprt,iw0,1,1,0,1,0,0,ipspage,noclose,1,0)
            
            if (nhx  >  1) :
              for i in range(0, 4):
                call blankout(printrlab(i),1,25)
                lprintrlab(i)=25
              
              incrhx=nhx*nhxres+1
              incrhx0=1
              for ihx in range(0, nhx):
                for jhx in range(ihx+1, nhx):
c                 Plot the two helix # end-# end distances
c                 res(2,nframes,ixres)=10000*idistss+idisthh
                  for i in range(0, nframe):
                    ii=res(2,i,incrhx)
                    plotdat(1,i)=float(ii/10000)/50.0
                    plotdat(2,i)=float(mod(ii,10000))/50.0
                    call trnsfr(res(1,i,1),plotdat(1,i),2)
                  
                  call arminmax2(plotdat,1,nframe,2,armin1,
     -              armax1,armin2,armax2,0,2)
                  iy1=armin1
                  y1mn=iy1
                  id1=armax1-iy1+1
                  if (mod(id1-iy1,2)  ==  1  or  id1. eq. 0) id1=id1+1
                  y1div=float(id1)/10.0
                  iy2=armin2
                  if (armin2  <  0.0) iy2=iy2-1
                  y2mn=iy2
                  call roundlim(armax2-y2mn,y2div,ny2div)
                  write (distlab,2162) ihx,jhx,'1st # end-# end',' ',
     -              ap_pa(isg2(incrhx0)+2)
                  write (distlab_cc,2162) ihx,jhx,'2nd # end-# end',' ',
     -              ap_pa(isg2(incrhx0)+2)
                  if (itmem  >  0  and  memdir(ihx)+memdir(jhx)  >  0)
     -                :
                    distlab(11:12)=in_ex(memdir(ihx))
                    distlab_cc(11:12)=in_ex(3-memdir(ihx))
                    distlab(13:13)=' '
                    distlab_cc(13:13)=' '
                  ## end if
                  printrlab(1)(1:22)=distlab(11:32)
                  printrlab(2)(1:22)=distlab_cc(11:32)
                  call plot2fun(iw1,2,xtraj,plotdat,plotdat,nframe,0.0,
     -              0.0,00,y1mn,y1div,10,y2mn,y2div,ny2div,title,80,' ',
     -              1,xtrajlab,11,distlab_cc,32,distlab,32,1,0,6,2,1,0,
     -              1,0,0,ipspage,noclose,1,0)
c                 Plot the two kinds of helix-helix distances
c                 res(1,nframes,ixres)=10000*idist+idist_cc
                  for i in range(0, nframe):
                    ii=res(1,i,incrhx)
                    plotdat(1,i)=float(ii/10000)/50.0
                    plotdat(2,i)=float(mod(ii,10000))/50.0
                    call trnsfr(res(1,i,2),plotdat(1,i),2)
                  
                  call arminmax2(plotdat,1,nframe,2,armin1,armax1,
     -              armin2,armax2,0,2)
                  iy1=armin1
                  y1mn=iy1
                  id1=armax1-iy1+1
                  if (mod(id1-iy1,2)  ==  1  or  id1. eq. 0) id1=id1+1
                  y1div=float(id1)/10.0
                  iy2=armin2
                  if (armin2  <  0.0) iy2=iy2-1
                  y2mn=iy2
                  call roundlim(armax2-y2mn,y2div,ny2div)
                  write (distlab,2162) ihx,jhx,'closest'
                  write (distlab_cc,2162) ihx,jhx,'center-center'
                  printrlab(3)(1:21)=distlab(11:31)
                  printrlab(4)(1:24)=distlab_cc(11:34)
                  call plot2fun(iw1,2,xtraj,plotdat,plotdat,nframe,0.0,
     -              0.0,00,y1mn,y1div,10,y2mn,y2div,ny2div,title,80,' ',
     -              1,xtrajlab,11,distlab,31,distlab_cc,34,1,0,6,2,1,0,
     -              1,0,0,ipspage,noclose,1,0)
                  write (distlab,2174) ihx,jhx,ap_pa(isg2(incrhx0)+2)
                  call trajstat(iw0,0,1,helixang,lhelixang,4,4,
     -              1,printrlab,lprintrlab,corr12,itorcorr,1,1,distlab,
     -              15,radtodeg)
                  incrhx=incrhx+1
                  incrhx0=incrhx0+1
                
              
              nhx2tot=(nhx*(nhx-1))/2
c             Plot helix axis angles
              incrhx=nhx*nhxres+1
              incrlab=1
              for ihx in range(0, nhx):
                for jhx in range(ihx+1, nhx):
c                 Convert angles back to sin/cos
c                 res(1,nframes,nhx2tot+ixres)=10000*iang+idhang
                  for i in range(0, nframe):
                    iangs=res(1,i,nhx2tot+incrhx)
                    ang=(float(iangs/10000)/10.0)/radtodeg
                    res(1,i,incrhx)=cos(ang)
                    res(2,i,incrhx)=sin(ang)
                  
                  write (hxhxlab(incrlab),2169) ihx,jhx
                  incrlab=incrlab+1
                  incrhx=incrhx+1
                
              
              call dialps(iw1,hxhxlab,lhxhxlab,title,ltitle,resrange,
     -          lresrange,5,nhx2tot,ndprow,nframe,pi,ipspage,9,
     -          nfravgd,nhx*nhxres,1,1,iconndial,ioutpr,mappdf,ipdfgrd)
              call trajstat(iw0,nhx2tot,nhx2tot,hxhxlab,lhxhxlab,0,1,
     -          1,printrlab,lprintrlab,corr12,nhx*nhxres,itorcorr,1,
     -          ' ',1,radtodeg)
c             Plot helix axis dihedral angles
              incrhx=nhx*nhxres+1
              incrlab=1
              for ihx in range(0, nhx):
                for jhx in range(ihx+1, nhx):
c                 Convert angles back to sin/cos
c                 res(1,nframes,nhx2tot+ixres)=10000*iang+idhang
                  for i in range(0, nframe):
                    iangs=res(1,i,nhx2tot+incrhx)
                    dhang=(float(mod(iangs,10000))/10.0)/radtodeg
                    res(1,i,incrhx)=cos(dhang)
                    res(2,i,incrhx)=sin(dhang)
                  
                  write (hxhxlab(incrlab),2170) ihx,jhx
                  incrhx=incrhx+1
                  incrlab=incrlab+1
                
              
              call dialps(iw1,hxhxlab,lhxhxlab,title,ltitle,resrange,
     -          lresrange,5,nhx2tot,ndprow,nframe,pi,ipspage,9,
     -          nfravgd,nhx*nhxres,1,1,iconndial,ioutpr,mappdf,ipdfgrd)
              call trajstat(iw0,nhx2tot,nhx2tot,hxhxlab,lhxhxlab,0,1,
     -          1,printrlab,lprintrlab,corr12,nhx*nhxres,1,itorcorr,
     -          ' ',1,radtodeg)
              nhx2tot2=2*nhx2tot
c             res(1,nframes,2*nhx2tot+ixres)=10000*idist_ar1+idist_ar2
c             Plot helix rotations
              incrhx=nhx*nhxres+1
              incrhx1=nhx*nhxres+1
              incrlab=1
              for ihx in range(0, nhx):
                for jhx in range(ihx+1, nhx):
                  for i in range(0, nframe):
                    iangs=res(1,i,nhx2tot2+incrhx1)
                    angrot1=(float(iangs/10000)/10.0)/radtodeg
                    angrot2=(float(mod(iangs,10000))/10.0)/radtodeg
                    res(1,i,incrhx)=cos(angrot1)
                    res(2,i,incrhx)=sin(angrot1)
                    res(1,i,incrhx+1)=cos(angrot2)
                    res(2,i,incrhx+1)=sin(angrot2)
                  
                  write (hxhxlab(incrlab),2171) ihx,jhx
                  write (hxhxlab(incrlab+1),2171) jhx,ihx
                  incrhx1=incrhx1+1
                  incrhx=incrhx+2
                  incrlab=incrlab+2
                
              
              call dialps(iw1,hxhxlab,lhxhxlab,title,ltitle,resrange,
     -          lresrange,5,nhx2tot2,ndprow,nframe,pi,ipspage,9,
     -          nfravgd,nhx*nhxres,1,0,iconndial,ioutpr,mappdf,ipdfgrd)
              call trajstat(iw0,nhx2tot2,nhx2tot2,hxhxlab,lhxhxlab,0,1,
     -          1,printrlab,lprintrlab,corr12,nhx*nhxres,1,itorcorr,
     -          ' ',1,radtodeg)
            ## end if
          elif (ianaltyp  ==  19) :
c           Torsion angle dial plots
            call dialps(iw1,talab,ltalab,title,ltitle,
     -        'Torsion angle dial plots',24,5,ntorsel,ndprow,nframe,
     -        pi,ipspage,0,nfravgd,0,0,1,iconndial,ioutpr,mappdf,
     -        ipdfgrd)
            call trajstat(iw0,ntorsel,MAXCOPY1,talab,ltalab,0,1,0,
     -        rmsdlab,lrmsdlab,corr12,0,itorcorr,0,' ',1,radtodeg)
            if (itorcorr  >  0) :
              call askyn(
     -          'Do you want to plot the circular correlation matrix',
     -          51,1,-1,iplotccc,0,0)
              if (iplotccc  ==  1) call plot_ccc(ntorsel,ioutpr,iw1,
     -          ipspage,talab,ltalab,inpfile,namleni)
              call askyn(
     -          'Do you want to calculate the eigenvectors/values',48,
     -          1,-1,ieigccc,0,0)
              if (ieigccc  ==  1) :
                call indexit(index2d,1,ndials,0)
                iout_ccc=47
                call read_write_ccc(ndials,iout_ccc,+1)
                call normalmodes(ndials,iout_ccc,nframe,1,iw0,1,ierr,
     -            index2d,temp,ifa_s,ila_s,it1,cv,MAX2D)
                call read_write_ccc(ndials,iout_ccc,-1)
                close (iout_ccc,status='delete')
              ## end if
              call askyn(
     -          'Do you want to cluster the angles by the correlation',
     -          52,1,-1,iclstccc,0,0)
              if (iclstccc  ==  1) call cluster_ccc(ndials,
     -          talab,ltalab,ixclst,ifclst1,ilclst1,index2d,it1,it2,
     -          it3,it4,' ',1,iw0,inpfile,namleni)
            ## end if
          elif (ianaltyp  ==  22) :
            call arminmax2(res(1,1,7),1,nframe,2,armin1,armax1,
     -        armin2,armax2,0,2)
            call roundlim(armax1,y1div,ny1div)
            rmsdmin=armin1
            rmsdmax=armax1
            if (rmsdplotmax  >  0.0) :
              ny1div=10
              y1div=rmsdplotmax/ny1div
            ## end if
            call arminmax2(res(2,1,7),1,nframe,2,armin1,armax1,
     -        armin2,armax2,0,2)
            call roundlim(armax1,y2div,ny2div)
            if (rmaxdevplotmax  >  0.0) :
              ny2div=10
              y2div=rmsdplotmax/ny2div
            ## end if
            call plot2fun(iw1,maxdevplot+1,xtraj,res(1,1,7),res(1,1,7),
     -        nframe,0.0,0.0,00,0.0,y1div,ny1div, 0.0,y2div,ny2div,
     -        title,80,'Best overlap',12,xtrajlab,11 ,'RMSD',4,
     -        'Maximum deviation',17,1,0,6,2,1,0,1,0,3,
     -        ipspage,1,1,0)
            if (iunmatchplot  ==  1) :
              call arminmax2(res(1,1,8),1,nframe,2,armin1,armax1,
     -          armin2,armax2,0,2)
              call roundlim(armax1,y1div,ny1div)
              call roundlim(armax2,y2div,ny2div)
              call plot2fun(iw1,maxdevplot+1,xtraj,res(1,1,8),
     -          res(1,1,8),nframe,0.0,0.0,00,0.0,y1div,ny1div,
     -          0.0,y2div,ny2div,title,80,
     -          'Comparison is done without overlay',34,xtrajlab,11 ,
     -          'RMSD',4,'Maximum deviation',17,1,
     -          0,6,2,1,0,1,0,0,ipspage,1,1,0)
            ## end if
            call trajstat(iw0,0,6,helixang,lhelixang,4,4,7,rmsdlab,
     -        lrmsdlab,corr12,itorcorr,0,0,' ',1,radtodeg)
            if (inptrajtyp  ==  3) :
              call correl(iw0,res(1,1,7),1,'RMSD',4,av1,sd1,res(1,1,9),
     -          1,'Energy',6,av2,sd2,corr,nframe,1,now6)
              call arminmax2(res(1,1,9),1,nframe,2,armin1,armax1,
     -          armin2,armax2,0,2)
              write (linein(1:48),2033) corr
              call scatterps(iw1,rmsdmin,rmsdmax,575.0,armin1,armax1,
     -          720.0,res(1,1,7),res(1,1,9),1,1,nframe,linein,48,1)
            ## end if
            if (nresslt  <=  MAXDISTR  and  nocontigrmsd  ==  0) :
              nosd=0
              for i in range(0, MAXDISTRN):
                xcum(i)=nframes_err(i)
c               print *,i,' nframes_err=',nframes_err(i)
              
              for ir in range(0, nresslt):
                rmsfavs(2,ir)=sqrt(rmsfsum(ir)/float(nframe))
              
              call rmsf_av(cdp,cdp2,atw,nfinalrmsd,indexrmsd,ixres,
     -           rmsfsum,rmsfav,nframe,nslt,nresslt)
              nresplot=0
              ifstres=1
              while (nresplot  ==  0)
                call getrange(ifstplot,1,ilstplot,nresslt,incr,0,
     -            'residue index in the RMSF plot',30,nresslt,000)
                if (ifstplot  <  ixres(ifstrmsd)  or 
     -              ilstplot  >  ixres(ilstrmsd)) :
                  write (6,2133) ixres(ifstrmsd),ixres(ilstrmsd)
                else:
                  nresplot=ilstplot-ifstplot+1
                  ifstres=ifstplot
                ## end if
              
              for ir in range(ifstplot, ilstplot):
                rmsfavs(1,ir)=rmsfav(ir)
                call blockfromcum(bl,distrerr(1,ir),xcum,MAXDISTRN)
                call batchmean(MAXDISTRN,0,bl,'RMSF error',10,iw0,1,
     -            avg,sd,ci)
                err12(1,ir)=sd
                err12(2,ir)=sqrt(ci)
              
              for ir in range(ifstplot, ilstplot):
                ci=err12(2,ir)
                sd=err12(1,ir)
                if (ci  !=  999.0) :
                  write (iw0,2155) ir,(rmsfavs(k,ir),k=1,2),sd,' ',ci
                else:
                  write (iw0,2155) ir,(rmsfavs(k,ir),k=1,2),sd
                ## end if
                err12(1,ir)=0.0
              
              call arminmax2(rmsfavs,ifstplot,ilstplot,2,armin1,
     -          armax1,armin2,armax2,0,2)
              if (rmsfplotmax  ==  0.0) :
                call arminmax2(rmsfavs,ifstplot,ilstplot,2,armin1,
     -            armax1,armin2,armax2,0,2)
                call roundlim(amax1(armax1,armax2),y1div,ny1div)
              else:
                ny1div=10
                y1div=rmsfplotmax/ny1div
              ## end if
              if (ixresno(ifstplot)  !=  ifstplot) :
                print *,'The residue number of the ',ifstplot,
     -            '-th residue is ',ixresno(ifstplot)
                call askyn('Do you want to use the residue number',37,
     -            1,+1,iusern,0,0)
                if (iusern  ==  1) ifstres=ixresno(ifstplot)
              ## end if
              ixmin=ifstres-mod(ifstres,10)
              ilstres=ifstres+nresplot-1
              ixmax=ilstres-mod(ilstres,10)
              if (mod(ilstres,10)  >  0) ixmax=ixmax+10.0
              xdiv=float(ixmax-ixmin)/10.0
              xmin=ixmin
              for ir in range(0, nresplot):
                cv(ir)=ifstres-1+ir
              
              plotdescr=
     -        'Residue RMSF wrt average (Cav) and input (Cref) position'
              call askyn('Do you want error bars plotted on RMSF(Cref)',
     -           44,0,+1,nosd,0,0)
              armax=amax1(armax1,armax2)
              if (nosd  ==  0  and  rmsfplotmax  ==  0.0) :
                armax=armax+1.0
                call roundlim(armax,y1div,ny1div)
              ## end if
              if (rmsfplotmax  >  0.0  and  rmsfplotmax  <  armax)
     -          print *,'NOTE part of the RMSF plots will be outside ',
     -            'the plot bounding box'
              ny2div=ny1div
              y2div=y1div
              print *
              call plot2fun(iw1,2,cv,rmsfavs(1,ifstplot),
     -          err12(1,ifstplot),nresplot,xmin,xdiv,10,0.0,y1div,
     -          ny1div,0.0,y2div,ny2div,title,80,plotdescr,56,'N(res)',
     -          6,'RMSF(Cav)',9,'RMSF(Cref)',10,1,0,6,2,nosd,nfreqsd,1,
     -          0,0,ipspage,0,1,0)
              call askyn('Do you want to write a PDB file with RMSF',41,
     -          1,0,iwpdb,0,0)
              if (iwpdb  >  0) :
                for ir in range(nresslt, 1, -1):
                  for ia in range(ifres(ir), ilres(ir)):
                    cv(ia)=rmsfav(ir)
                  
                
                analfile2(namleno2:namleno2)='f'
                analfile2(namleno2+1:namleno2+4)='.pdb'
                call openfile(iw2,0,'RMSF',4,'new',analfile2,namleno2+4,
     -            notfnd,0,1,1,0,0)
                write (6,2031) analfile2(1:namleno2+4)
                call writeconf(iw2,inpcrdtyp,inpcrdtyp,inpcrdtyporg,n,n,
     -            nslt,naslv,islvw,iasv,namesv,qsv,pflsv,1,iwhead,0,
     -            iatnum,ifchrg,nconfig,innlist,c,rprox,cv,ixres,iresno,
     -            atnames,resnames,segnames,charge,isegno,altcol,inscol,
     -            ninsres,marker,ntitlin,ntitlinw,title,ireseq,
     -            iresnrestart,iresidrestart,nneig,nneiga,nhbneig,ineig,
     -            nhneig,nnneig,ncneig,nsneig,npneig,numres,numslv,
     -            resnamslv,line,blankline,mmtype,ibnd,index,indexn,
     -            indexo,1,molresflag,irescount3,itemp1,hblimfac,angmin,
     -            0,1,1,1,0,3,iqspaceask,ianaltyp,0,0.0,0,0,0,keeprem,
     -            iwriteatsym,radtodeg,maxrepconf,maxng,maxrsd,maxrec)
              ## end if
            ## end if
          elif (ianaltyp  ==  23) :
c           Prepare 2-D RMSD plot
            noclose=0
            if (irmsdclust  >  0) noclose=1
            nframeref=nframeref-1
            call plot2drmsd(nrep,iw0,iw1,xtraj,maxrec,title,'',0,
     -        0,xtrajlab,11,
     -        ncolcode,maxcolcode,iedit,noopt2d,limresrange,
     -        0,rmsdmin,rmsdmax,absdevmin,absdevmax,rmsdmn,rmsdmx,
     -        indexa,indexs,ixshuffle,ixshuffle,indexa,0,indexa,0,
     -        ym_2d,it1,it2,temp,0,1,noplotdist,iskip2dplot,0,
     -        noclose,ipspage,1)
            if (irmsdclust  >  0) :
              inspect=0
              write (6,2068)
8010          call trnsfi(iconfsel,indexs,nframe)
              call clusterdistr(nframe,iw0,rmsdlim,rmsdmn,rmsdmx,
     -          nhbdist,it1,it2,it3,itemp4,indexn,indexo,ncl,
     -          indexa,iconfsel,ixclst,it4,value,ifa_s,ila_s,ih,cv,0.0,
     -          rdclust,res(1,1,11),ietotsaved,'RMSD',4,1,1,irepav,
     -          irepmx,irepeng,irepkm,engcl,c,chn,c2,1,26,iclstyp,iwt,0,
     -          label2d,80,0,1,1,1,mx2d,maxframe)
c             iconfsel contains the sorted cluster members
c             indexn, indexo contain the cluster limits
              for i in range(0, ncl):
                nclstmem(i)=indexo(i)-indexn(i)+1
              
              call countsim(indexn,indexo,iconfsel,ncl,rdclust,rmsdsim,
     -          nsimclst1,iw0,mx2d)
              if (ireplot  ==  1) :
                call trnsfi(ixshuffle,iconfsel,MAX2D)
                call indexit(it3,1,MAX2D,0)
                if (inspect  ==  1) :
                  call openfile(iw1,0,'analysis',8,'old',analfile1,
     -              namleno1,notfnd,0,1,1,0,0)
                  while ( True )
                    read (iw1,2005,# end=8011) ans
                  
8011              continue
                ## end if
                call plot2drmsd(nrep,iw0,iw1,xtraj,maxrec,title,
     -            'Frames sorted by clusters',25,
     -            0,xtrajlab,11,ncolcode,
     -            maxcolcode,iedit,noopt2d,limresrange,0,rmsdmin,
     -            rmsdmax,absdevmin,absdevmax,rmsdmn,rmsdmx,indexa,
     -            it3,ixshuffle,ixshuffle,indexo,ncl,indexo,ncl,
     -            ym_2d,it1,it2,temp,1,1,noplotdist,iskip2dplot,1,1,
     -            ipspage,0)
                call askyn('Do you want to inspect the clustered map'//
     -            ' and redo the clustering',64,1,-1,inspect,000,0)
                if (inspect  >  0) :
                  call clusterplot(iw1,xtraj,value,indexn,indexo,ncl,
     -              ixclst,nframe,xtrajlab,11,ipspage,0,mx2d)
                  close (iw1)
                  go to 8010
                ## end if
              ## end if
              call clusterplot(iw1,xtraj,value,indexn,indexo,ncl,ixclst,
     -          nframe,xtrajlab,11,ipspage,0,mx2d)
            ## end if
          elif (ianaltyp  ==  24) :
c           Prepare cross RMSD plot
            call plot2drmsd(nrep,iw0,iw1,xtraj,maxrec,title,'',0,
     -        ltrajnam2,xtrajlab,11,
     -        ncolcode,maxcolcode,iedit,noopt2d,limresrange,
     -        1,rmsdmin,rmsdmax,absdevmin,absdevmax,rmsdmn,rmsdmx,
     -        indexa,indexs,ixshuffle,ixshuffle,indexa,0,indexa,0,
     -        ym_2d,it1,it2,temp,0,0,noplotdist,0,0,0,ipspage,0)
            if (matchconf  ==  1)
     -        call matchtraj(rmsdsim,iw0)
            write (6,2038)
            write (iw0,2038)
          elif (ianaltyp  ==  26) :
c           Prepare distance distribution output
            call pairdistprint(nframe,npairs,listpairdist,iclusterdist,
     -        iclustermem,ifstclst1,ifstclst2,ilstclst2,pairdistsum,
     -        pairdistsum2,pairdistwsum,npairdist,pairdistminmax,
     -        pairgrid,rmaxpair,line,index,inamcol1,inamcol2,irescol1,
     -        irescol2,iresncol1,iresncol2,inpcrdtyp,ioins,iw0,nslt,
     -        MAXDDBIN,MAXDDISTR,MAXCDLIST,maxrec)
            call askyn(
     -        'Do you want to a PDB file with bonds for the distances',
     -         54,1,0,iwpdb,0,0)
            if (iwpdb  ==  1) :
              analfile4=analfile
              analfile4(namleno-2:namleno+4)='dsm.pdb'
              iw4=44
              call openfile(iw4,0,'PDB with distances marked',25,
     -          'new',analfile4,namleno+4,notfnd,0,1,1,0,0)
              call writeconf(iw4,inpcrdtyp,iobpdb,inpcrdtyporg,
     -          nslt,nslt,nslt,naslv,islvw,iasv,namesv,qsv,pflsv,1,
     -          iwhead,0,iatnum,ifchrg,nconfig,innlist,c,rprox,cv,ixres,
     -          iresno,atnames,resnames,segnames,charge,isegno,altcol,
     -          inscol,ninsres,marker,ntitlin,ntitlinw,title,ireseq,
     -          iresnrestart,iresidrestart,nneig,nneiga,nhbneig,ineig,
     -          nhneig,nnneig,ncneig,nsneig,npneig,numres,numslv,
     -          resnamslv,line,blankline,mmtype,ibnd,index,indexn,
     -          indexo,1,molresflag,irescount3,itemp1,hblimfac,angmin,
     -          0,1,1,1,0,3,iqspaceask,90+icontyp,0,0.0,0,0,1,keeprem,
     -          iwriteatsym,radtodeg,maxrepconf,maxng,maxrsd,maxrec)
              for i in range(0, npairs):
                write (iw4,2131) (listpairdist(j,i),j=1,2)
              
              write (iw4,2005) '# end'
              close (iw4)
              write (6,2019) 'PDB file with distances marked',
     -          analfile4(1:namleno+4)
            ## end if
          elif (ianaltyp  ==  39) :
c           Prepare distance SD matrix plot
            call plot_atomdist_sd(nslt,line,index,inamcol1,inamcol2,
     -        irescol1,irescol2,iresncol1,iresncol2,ianchor,nanchor,
     -        nframe,indexa,ixshuffle,xtraj,title,trajnam,ltrajnam,
     -        isdtyp,iw0,iw1,ipspage,maxrec)
          elif (ianaltyp  ==  27) :
            call trajstat(iw0,0,1,helixang,lhelixang,4,4,1,volumelab,
     -        lvolumelab,corr12,0,0,0,' ',1,radtodeg)
            call arminmax2(res(1,1,1),1,nframe,2,armin1,armax1,
     -        armin2,armax2,0,2)
            call setdivxy(armin1,armax1,ny1div,y1div,y1min)
            call setdivxy(armin2,armax2,ny2div,y2div,y2min)
            write (linein,2148) nrand,corr12
            call plot2fun(iw1,2,xtraj,res(1,1,1),res(1,1,3),nframe,
     -        0.0,0.0,0,y1min,y1div,ny1div,y2min,y2div,ny2div,title,80,
     -        linein,70,xtrajlab,11 ,'V(solvent-excluded shell)',25,
     -        'V(first solvation shell)',24,1,0,6,2,00,0,
     -        1,0,4,ipspage,1,1,0)
            call arminmax2(res(1,1,2),1,nframe,2,armin21,armax21,
     -        armin22,armax22,0,2)
            linein(1:1)=' '
            call setdivxy(armin21,armax21,ny1div,y1div,y1min)
            call setdivxy(armin22,armax22,ny2div,y2div,y2min)
            call plot2fun(iw1,2,xtraj,res(1,1,2),res(1,1,4),nframe,
     -        0.0,0.0,0,y1min,y1div,ny1div,y2min,y2div,ny2div,title,80,
     -        linein,01,xtrajlab,11 ,'V(macromolecule)',16,
     -        'V(interface)',12,1,0,6,2,01,1,1,0,0,
     -        ipspage,1,1,0)
            nbin=max0(10,nframe/100)
            call zeroiti(indexo,0,nbin)
            range1=armax1-armin1
            call roundlim(range1,xdiv,nxdiv)
            nx1=armin1/xdiv
            xmin=nx1*xdiv
            facnorm=nbin/amax1(1.0,range1)
            for i in range(0, nframe):
              ix=min0(nbin,int(facnorm*(res(1,i,1)-armin1))+1)
              indexo(ix)=indexo(ix)+1
            
            histmax=0.0
            for i in range(0, nbin):
              cv(i)=xmin+(i-0.5)*range1/nbin
              charge(i)=indexo(i)
              if (charge(i)  >  histmax) histmax=charge(i)
            
            call roundlim(histmax,y1div,ny1div)
            call plot2fun(iw1,1,cv,charge,res(1,1,3),nbin,xmin,xdiv,
     -        nxdiv,0.0,y1div,ny1div, 0.0,y2div,ny2div,title,80,linein,
     -        36,'V(solvent-excluded shell)',25,xtrajlab,11 ,' ',1,
     -        1,0,6,1,01,0,1,0,0,ipspage,1,1,0)
            call zeroiti(indexo,0,nbin)
            call setdivxy(armin2,armax2,nxdiv,xdiv,xmin)
            facnorm=nbin/amax1(1.0,armax2-armin2)
            for i in range(0, nframe):
              ix=min0(nbin,int(facnorm*(res(2,i,1)-armin2))+1)
              indexo(ix)=indexo(ix)+1
            
            histmax=0.0
            for i in range(0, nbin):
              cv(i)=xmin+(i-0.5)*(armax2-armin2)/nbin
              charge(i)=indexo(i)
              if (charge(i)  >  histmax) histmax=charge(i)
            
            call roundlim(histmax,y2div,ny2div)
            call plot2fun(iw1,1,cv,charge,res(1,1,3),nbin,xmin,xdiv,
     -        nxdiv,0.0,y1div,ny1div, 0.0,y2div,ny2div,title,80,linein,
     -        36,'V(first solvation shell)',24,xtrajlab,11 ,' ',1,
     -        1,0,6,1,01,0,1,0,0,ipspage,0,1,0)
          elif (ianaltyp  ==  28) :
c           Principal axis calculation
            call arminmax2(res(1,1,6),1,nframe,2,armin1,armax1,
     -        armin2,armax2,0,2)
            call setdivxy(armin1,armax1,ny1div,y1div,y1min)
            call setdivxy(armin2,armax2,ny2div,y2div,y2min)
            call setdivxy(0.0,xtraj(nframe),nxdiv,xdiv,xmin)
            write (linein,2101) 'first and second'
            call plot2fun(iw1,2,xtraj,res(1,1,6),res(1,1,6),nframe,xmin,
     -        xdiv,nxdiv,y1min,y1div,ny1div,y2min,y2div,ny2div,title,80,
     -        linein,66,xtrajlab,11 ,'First principal axis',20,
     -        'Second principal axis',21,1,0,6,2,1,0,1,0,
     -        2,ipspage,1,1,0)
            call arminmax2(res(1,1,7),1,nframe,1,armin1,armax1,
     -        armin2,armax2,0,2)
            call setdivxy(armin1,armax1,ny1div,y1div,y1min)
            write (linein,2101) 'third'
            call plot2fun(iw1,1,xtraj,res(1,1,7),res(1,1,7),nframe,xmin,
     -        xdiv,nxdiv,y1min,y1div,ny1div,y2min,y2div,ny2div,title,80,
     -        linein,55,xtrajlab,11 ,'Third principal axis',20,' ',1,
     -        1,0,6,2,1,0,1,0,0,ipspage,0,1,0)
          elif (ianaltyp  ==  29) :
c           Molecular radii, moments of inertia and dipole moment calculation
            call averageres(nframe,res,1,1,MAXFRAMES,MAXCOPY1,rgav,rgsd)
            call averageres(nframe,res,1,2,MAXFRAMES,MAXCOPY1,rhav,rhsd)
            call averageres(nframe,res,2,2,MAXFRAMES,MAXCOPY1,rmxmn,
     -        rmxmnsd)
            write (iw0,2079) rgav,rgsd,1.0/rhav,rhsd,rmxmn
            call arminmax2(res(1,1,1),1,nframe,2,armin1,armax1,
     -        armin2,armax2,0,2)
            call setdivxy(armin1,armax1,ny1div,y1div,y1min)
            call setdivxy(armin2,armax2,ny2div,y2div,y2min)
            call setdivxy(0.0,xtraj(nframe),nxdiv,xdiv,xmin)
            noclose=0
            if (icharges  >  0) noclose=1
            call plot2fun(iw1,2,xtraj,res(1,1,1),res(1,1,1),nframe,xmin,
     -        xdiv,nxdiv,y1min,y1div,ny1div,y2min,y2div,ny2div,title,80,
     -        linein,00,xtrajlab,11 ,'Radius of gyration',18,
     -        'Hydrodynamic radius',19,1,0,6,2,1,0,1,0,1,
     -        ipspage,noclose,1,0)
            if (icharges  >  0) :
              call averageres(nframe,res,1,10,MAXFRAMES,MAXCOPY1,avcell,
     -          sdcell)
              call averageres(nframe,res,2,10,MAXFRAMES,MAXCOPY1,avslt,
     -          sdslt)
              write (iw0,2080) avcell,sdcell,avslt,sdslt
              call arminmax2(res(1,1,10),1,nframe,2,armin1,armax1,
     -          armin2,armax2,0,2)
              call setdivxy(armin1,armax1,ny1div,y1div,y1min)
              call setdivxy(armin2,armax2,ny2div,y2div,y2min)
              call averageres(nframe,res,1,11,MAXFRAMES,MAXCOPY1,
     -          dip1(1),sdxyx1(1))
              call averageres(nframe,res,2,11,MAXFRAMES,MAXCOPY1,
     -          dip1(2),sdxyx1(2))
              call averageres(nframe,res,1,12,MAXFRAMES,MAXCOPY1,
     -          dip1(3),sdxyx1(3))
              write (iw0,2088)
     -           'Total cell',(xyz(k),dip1(k),sdxyx1(k),k=1,3)
              nplot=1
              if (n  >  nslt) :
                nplot=2
                call averageres(nframe,res,1,13,MAXFRAMES,MAXCOPY1,
     -          dip1(1),
     -            sdxyx1(1))
                call averageres(nframe,res,2,13,MAXFRAMES,MAXCOPY1,
     -          dip1(2),sdxyx1(2))
                call averageres(nframe,res,1,14,MAXFRAMES,MAXCOPY1,
     -          dip1(3),sdxyx1(3))
                write (iw0,2088)
     -            'Solute',(xyz(k),dip1(k),sdxyx1(k),k=1,3)
              ## end if
              call plot2fun(iw1,nplot,xtraj,res(1,1,10),res(1,1,10),
     -          nframe,xmin,xdiv,nxdiv,y1min,y1div,ny1div,y2min,y2div,
     -          ny2div,title,80,linein,00,xtrajlab,11,'Total dipole',
     -          12,'Solute dipole',13,1,0,6,2,1,0,1,0,0,
     -          ipspage,0,1,0)
              call correl(iw0,res(1,1,1),1,'Radius of gyration',18,
     -          av1,sd1,res(1,1,10),1,'Solute dipole',13,av2,sd2,corr,
     -          nframe,1,0)
            ## end if
          elif (ianaltyp  ==  32) :
c           Angle dial plots
            call dialps(iw1,talab,ltalab,title,ltitle,
     -        'Angle dial plots',16,5,nangsel,ndprow,nframe,pi,ipspage,
     -        0,nfravgd,0,0,1,iconndial,ioutpr,mappdf,ipdfgrd)
            call trajstat(iw0,nangsel,MAXCOPY1,talab,ltalab,0,1,0,
     -        rmsdlab,lrmsdlab,corr12,0,itorcorr,0,' ',1,radtodeg)
          ## end if
          if (iw1  >  0) write (6,2073) analfile1(1:namleno1)
        elif (ianaltyp  ==  8) :
c         Print contact count
          write (iw0,2113)
          write (iw1,2113)
          if (iw2  >  0) write (iw2,2113)
          call blankout(linein,1,80)
c         print *,'NSEGCOL=',nsegcol
          for ir in range(0, nresslt):
            if (irescount1(ir)+irescount2(ir)+irescount3(ir)  >  0)
     -          :
              write (linein,2114) ir,iresno(ifres(ir)),
     -          resnames(ir)(1:nrescol),
     -          segid4(isegno(ifres(ir)))(1:nsegcol),
     -          irescount1(ir),irescount2(ir),irescount3(ir)
              call lastchar(linein,lc,80)
              write (iw0,2005) linein(1:lc)
              write (iw1,2005) linein(1:lc)
              if (iw2  >  0) write (iw2,2005) linein(1:lc)
            ## end if
          
          call print_rrdist(itypavg,nframe,irefres1,irefres2,inegres1,
     -      inegres2,listrefres,nrefres,listnegres,nnegres,nrescol,iw0,
     -      iw1,iw4,ipspage,resnames,inpfile,namleni,itemp1,maxrsd)
          call askyn(
     -      'Do you want to write a PDB file with contact counts',51,1,
     -      0,iwpdb,0,0)
          if (iwpdb  >  0) :
            call quiz(ans,icontyp,' ',' ',0,'Contact type',12,0,5,6,
     -        0)
            for ir in range(0, nresslt):
              for ia in range(ifres(ir), ilres(ir)):
                if (icontyp  ==  1) :
                  cv(ia)=irescount1(ir)
                elif (icontyp  ==  2) :
                  cv(ia)=irescount2(ir)
                elif (icontyp  ==  3) :
                  cv(ia)=irescount3(ir)
                elif (icontyp  ==  4) :
                  cv(ia)=irescount1(ir)+irescount2(ir)+irescount3(ir)
                ## end if
                rprox(ia)=1.0
              
            
            analfile4=analfile1
            analfile4(namleno1-2:namleno1+4)='cnt.pdb'
            iw4=44
            call openfile(iw4,0,'Contact-labeled PDB',19,
     -        'new',analfile4,namleno1+4,notfnd,0,1,1,0,0)
            call writeconf(iw4,inpcrdtyp,iobpdb,inpcrdtyporg,
     -        nslt,nslt,nslt,naslv,islvw,iasv,namesv,qsv,pflsv,1,
     -        iwhead,0,iatnum,ifchrg,nconfig,innlist,c,rprox,cv,ixres,
     -        iresno,atnames,resnames,segnames,charge,isegno,altcol,
     -        inscol,ninsres,marker,ntitlin,ntitlinw,title,ireseq,
     -        iresnrestart,iresidrestart,nneig,nneiga,nhbneig,ineig,
     -        nhneig,nnneig,ncneig,nsneig,npneig,numres,numslv,
     -        resnamslv,line,blankline,mmtype,ibnd,index,indexn,indexo,
     -        1,molresflag,irescount3,itemp1,hblimfac,angmin,0,1,1,1,0,
     -        3,iqspaceask,90+icontyp,0,0.0,0,0,0,keeprem,iwriteatsym,
     -        radtodeg,maxrepconf,maxng,maxrsd,maxrec)
            close (iw4)
          ## end if
        elif (ianaltyp  ==  12) :
          write (title,2069) irespro,nra,nrb
          call dialps(iw1,prokinklab,lprokinklab,title,ltitle,resrange,
     -      22,5,ndials,ndprow,nframe,pi,ipspage,0,nfravgd,0,0,0,
     -      iconndial,ioutpr,mappdf,ipdfgrd)
          write (6,2073) analfile1(1:namleno1)
          call trajstat(iw0,ndials,5,prokinklab,lprokinklab,0,1,0,
     -      rmsdlab,lrmsdlab,corr12,0,0,0,' ',1,radtodeg)
        elif (ianaltyp  ==  25) :
          call openps(iw1,xm,ym-30,title,ltitle,
     -      ' ',1,trajnam,0,trajnam2,ltrajnam2,npspages,ipspage)
          call plotresidcorr(ncorr,nframe,indexa,indexs,
     -      ncolcode,maxcolcode,nrep,iw0,iw1,xm,ym,title,ltitle,
     -      ipspage,iucorrmat,icovmatplot,cv,maxrec)
          if (iucorrmat  >  0) :
c           Calculate eigenvalues/eigenvectors
            close (iw0)
            analfile(namleno+1:namleno+4)='.eig'
            namleno=namleno+4
            write (6,2019) 'Eigenvalues and eigenvectors',
     -        analfile(1:namleno)
            call openfile(iw0,0,'eigenvalues \ eigenvectors',26,'new',
     -        analfile,namleno,notfnd,0,1,1,0,0)
            call normalmodes(ncorr,iucorrmat,nframe,
     -        0,iw0,1,ierr,index2d,value,ifa_s,ila_s,it1,rmsdlim,
     -        MAXBONDS)
            close (iucorrmat,status='delete')
          ## end if
        elif (ianaltyp  ==  16  or  ianaltyp  ==  18) :
          if (ianaltyp  ==  18) :
            call rainbowscale(iw1,50,450,25,nframe,xtraj(nframe),0.0,
     -        0.0,xtrajlab,11)
            call dialps(iw2,ramalab,lramalab,title,ltitle,' ',0,5,
     -        2*nxselres,ndprow,nframe,pi,ipspage,0,nfravgd,0,0,0,
     -        iconndial,ioutpr,mappdf,ipdfgrd)
            close (iw2)
            call ramachandran_hist(nresfound,resnames,nrescol,iw0,1,0)
            call ramachandran_hist(nresfound,resnames,nrescol,iw0,0,1)
            iw3=iw2+1
            analfile3=analfile
            analfile3(namleno+1:namleno+7)='.trc.ps'
            namleno3=namleno+7
            write (6,2019) '2D traces on the psi-phi map',
     -        analfile3(1:namleno3)
            call openfile(iw3,0,'2D_Rama_trace',13,'new',analfile3,
     -         namleno3,notfnd,0,1,1,0,0)
            analfile2=analfile
            analfile2(namleno+1:namleno+7)='.auc.ps'
            namleno2=namleno+7
            write (6,2019) 'Psi-Phi autocorrelation defs',
     -        analfile2(1:namleno2)
            call openfile(iw2,0,'Psi-Phi autoc',13,'new',analfile2,
     -         namleno2,notfnd,0,1,1,0,0)
            restitle='Residue 00000 (     ) nconf=        '
            incrementac=1
            if (nframe  >  500) call getint(
     -        'Frame increment for autocorrelation calculation',47,1,1,
     -        nframe/2,incrementac,0)
            for i in range(0, nframe/(2*incrementac)):
              rplot(i)=i*incrementac
            
            npspages=nxselres
            for ir in range(0, nxselres):
              if (ir  ==  nxselres) noclose=0
              iresdr=ixselres(ir)
              write (restitle(9:13),2103) iresdr
              restitle(16:15+nrescol)=resnames(iresdr)(1:nrescol)
              sinphisum=0.0
              cosphisum=0.0
              sinpsisum=0.0
              cospsisum=0.0
              for i in range(0, nframe):
                xyplot(1,i)=dacoscheck(ddd,res(1,i,2*ir-1),0,6,'RAMA')*
     -            radtodeg
                if (res(2,i,2*ir-1)  <  0) xyplot(1,i)=-xyplot(1,i)
                xyplot(2,i)=dacoscheck(ddd,res(1,i,2*ir),0,6,'RAMA')*
     -            radtodeg
                if (res(2,i,2*ir)  <  0) xyplot(2,i)=-xyplot(2,i)
                cosphisum=cosphisum+res(1,i,2*ir-1)
                sinphisum=sinphisum+res(2,i,2*ir-1)
                cospsisum=cospsisum+res(1,i,2*ir)
                sinpsisum=sinpsisum+res(2,i,2*ir)
              
              cosphisum=cosphisum/nframe
              sinphisum=sinphisum/nframe
              sqsum=sqrt(cosphisum**2+sinphisum**2)
              if (sqsum  ==  0.0) sqsum=1.0
              cosphisum=cosphisum/sqsum
              cvphi=1.0-sqsum
              phiav=dacoscheck(cosphisum,cps,1,6,'PHIAV')*radtodeg
              if (sinphisum  <  0.0) phiav=-phiav
              cospsisum=cospsisum/nframe
              sinpsisum=sinpsisum/nframe
              sqsum=sqrt(cospsisum**2+sinpsisum**2)
              if (sqsum  ==  0.0) sqsum=1.0
              cospsisum=cospsisum/sqsum
              psiav=dacoscheck(cospsisum,cps,1,6,'PSIAV')*radtodeg
              if (sinpsisum  <  0.0) psiav=-psiav
              cvpsi=1.0-sqsum
              write (6,2093) iresdr,resnames(iresdr)(1:nrescol),
     -          phiav,cvphi,psiav,cvpsi
              write (iw0,2093) iresdr,resnames(iresdr)(1:nrescol),
     -          phiav,cvphi,psiav,cvpsi
              call plot2d(iw3,xyplot,nframe,1,-180.0,30.0,12,
     -          -180.0,30.0,12,title,80,'2D trace on the psi-phi map',
     -          27,restitle,21,'Phi',3,'Psi',3,xtrajlab,11,1,0,6,
     -          npspages,ipspage3,1,MAXFRAMES)
              call resautocorr(ir,incrementac,ncf,xyplot,MAXFRAMES)
              write (restitle(29:36),2094) ncf
              call plot2fun(iw2,1,rplot,xyplot,xyplot,ncf,0.0,0.0,0,
     -          0.0,0.1,10,0.0,0.1,10,
     -          'Psi-Phi autocorrelation def',32,restitle,36,
     -          xtrajlab,11,'Autocorrelation',15,' ',1,
     -          1,0,77,2,1,0,1,0,npspages,ipspage2,1,1,0)
              npspages=0
            
          ## end if
          write (iw1,*) 'showpage'
          close (iw1)
          if (ianaltyp  ==  16) :
c           Print statistics
            write (iw0,2149) (typc(it),ssname(it)(1:lssname(it)),it=1,9)
            write (iw0,2150) (typc(it),it=1,9)
            write (iw0,2151)
     -         (ir,(idistdssp(it,ir),it=1,9),ir=ifrdssp,ilrdssp)
          ## end if
          close (iw2)
          close (iw3)
        elif (ianaltyp  ==  17) :
          call hbbridgeprint(nanchor,ianchor,lpath,nbridgetype,
     -      ibridgetype,maxbridgemem,line,index,iresno,inamcol1,
     -      inamcol2,irescol1,irescol2,iw0,maxbondcount,
     -      maxhbtype,minbridgelenprint,minbridgepercprint,nframe,
     -      MAXBRIDGELEN,MAXBRIDGETYPE,maxrec)
        elif (ianaltyp  ==   5) :
          call finalizebonds(n,nbfound,nbfoundorig,nbresfound,iw0,iw1,
     -      numres,nresslt,npspages,ipspage,nres2d,ibondcorr,
     -      iresbondcorr,nhneigmin,hblimfac,angmin,0.0,inamcol1,
     -      inamcol2,irescol1,irescol2,ifhb2d,ilhb2d,nhbdist,rhbdist,
     -      iframeunit,framefac,title,ltitle,xtrajlab,11,xtraj,value,
     -      ifa_s,ila_s,ih,cv,temp,index,iresno,ifres,isegno,atnames,
     -      resnames,ixres,ixresno,ixsegno,indexa,indexs,index2d,
     -      ianc_anc,isc,ixclst,irepav,irepmx,irepeng,irepkm,rmsdlim,
     -      engcl,it1,it2,it3,it4,it5,irrix,itemp1,itemp2,itemp3,itemp4,
     -      line,'hydrogen',8,ibondtype,label2d,iselfanc,ianchor2,
     -      iresshift,ifailbond,nbondavg,inpfile,namleni,maxbondf,
     -      nmcmaxbond,ncolcode,maxcolcode,maxbondcount,MAXBONDS,maxrsd,
     -      MAXFRAMES,maxrec,mx2d)
        elif (ianaltyp  ==   6) :
          call finalizebonds(n,nbfound,nbfoundorig,nbresfound,iw0,iw1,
     -      numres,nresslt,npspages,ipspage,nres2d,ibondcorr,
     -      iresbondcorr,nhneigmin,hblimfac,angmin,rhphmax,inamcol1,
     -      inamcol2,irescol1,irescol2,ifhb2d,ilhb2d,nhbdist,rhbdist,
     -      iframeunit,framefac,title,ltitle,xtrajlab,11,xtraj,value,
     -      ifa_s,ila_s,ih,cv,temp,index,iresno,ifres,isegno,atnames,
     -      resnames,ixres,ixresno,ixsegno,indexa,indexs,index2d,
     -      ianc_anc,isc,ixclst,irepav,irepmx,irepeng,irepkm,rmsdlim,
     -      engcl,it1,it2,it3,it4,it5,irrix,itemp1,itemp2,itemp3,itemp4,
     -      line,'hydrophobic',11,ibondtype,label2d,iselfanc,ianchor2,
     -      iresshift,ifailbond,nbondavg,inpfile,namleni,maxbondf,
     -      nmcmaxbond,ncolcode,maxcolcode,maxbondcount,MAXBONDS,maxrsd,
     -      MAXFRAMES,maxrec,mx2d)
        elif (ianaltyp  ==   34) :
          call finalizebonds(n,nbfound,nbfoundorig,nbresfound,iw0,iw1,
     -      numres,nresslt,npspages,ipspage,nres2d,ibondcorr,
     -      iresbondcorr,nhneigmin,hblimfac,angmin,rhphmax,inamcol1,
     -      inamcol2,irescol1,irescol2,ifhb2d,ilhb2d,nhbdist,rhbdist,
     -      iframeunit,framefac,title,ltitle,xtrajlab,11,xtraj,value,
     -      ifa_s,ila_s,ih,cv,temp,index,iresno,ifres,isegno,atnames,
     -      resnames,ixres,ixresno,ixsegno,indexa,indexs,index2d,
     -      ianc_anc,isc,ixclst,irepav,irepmx,irepeng,irepkm,rmsdlim,
     -      engcl,it1,it2,it3,it4,it5,irrix,itemp1,itemp2,itemp3,itemp4,
     -      line,'heavy-atom contact',18,ibondtype,label2d,iselfanc,
     -      ianchor2,iresshift,ifailbond,nbondavg,inpfile,namleni,
     -      maxbondf,nmcmaxbond,ncolcode,maxcolcode,maxbondcount,
     -      MAXBONDS,maxrsd,MAXFRAMES,maxrec,mx2d)
        elif (ianaltyp  ==   41) :
          call finalizebonds(n,nbfound,nbfoundorig,nbresfound,iw0,iw1,
     -      numres,nresslt,npspages,ipspage,nres2d,ibondcorr,
     -      iresbondcorr,nhneigmin,hblimfac,angmin,rhphmax,inamcol1,
     -      inamcol2,irescol1,irescol2,ifhb2d,ilhb2d,nhbdist,rhbdist,
     -      iframeunit,framefac,title,ltitle,xtrajlab,11,xtraj,value,
     -      ifa_s,ila_s,ih,cv,temp,index,iresno,ifres,isegno,atnames,
     -      resnames,ixres,ixresno,ixsegno,indexa,indexs,index2d,
     -      ianc_anc,isc,ixclst,irepav,irepmx,irepeng,irepkm,rmsdlim,
     -      engcl,it1,it2,it3,it4,it5,irrix,itemp1,itemp2,itemp3,itemp4,
     -      line,'mutually proximal atom pairs',28,ibondtype,label2d,
     -      iselfanc,ianchor2,iresshift,ifailbond,nbondavg,inpfile,
     -      namleni,maxbondf,nmcmaxbond,ncolcode,maxcolcode,
     -      maxbondcount,MAXBONDS,maxrsd,MAXFRAMES,maxrec,mx2d)
        elif (ianaltyp  ==   7) :
          call finalizebonds(n,nbfound,nbfoundorig,nbresfound,iw0,iw1,
     -      numres,nresslt,npspages,ipspage,nres2d,ibondcorr,
     -      iresbondcorr,nhneigmin,hblimfac,angmin,rsltbmax,inamcol1,
     -      inamcol2,irescol1,irescol2,ifhb2d,ilhb2d,nhbdist,rhbdist,
     -      iframeunit,framefac,title,ltitle,xtrajlab,11,xtraj,value,
     -      ifa_s,ila_s,ih,cv,temp,index,iresno,ifres,isegno,atnames,
     -      resnames,ixres,ixresno,ixsegno,indexa,indexs,index2d,
     -      ianc_anc,isc,ixclst,irepav,irepmx,irepeng,irepkm,rmsdlim,
     -      engcl,it1,it2,it3,it4,it5,irrix,itemp1,itemp2,itemp3,itemp4,
     -      line,'salt-bridge',11,ibondtype,label2d,iselfanc,ianchor2,
     -      iresshift,ifailbond,nbondavg,inpfile,namleni,maxbondf,
     -      nmcmaxbond,ncolcode,maxcolcode,maxbondcount,MAXBONDS,maxrsd,
     -      MAXFRAMES,maxrec,mx2d)
        elif (ianaltyp  ==  11) :
          zav=zavs/dfloat(nframe)
          zsd=math.sqrt(abs(zsqs/dfloat(nframe)-zav**2))
          write (6,2008) zav,zsd
          write (iw0,2008) zav,zsd
          pi=atan(1.0)*4.0
          write (6,2009)
          write (iw0,2009)
          for m in range(2, npsr):
            sn=sinpsrs(m)/dfloat(nframe)
            cs=cospsrs(m)/dfloat(nframe)
            psr(m)=atansc(sn,cs,1,radtodeg)
            cvm=1.0-math.sqrt(sinpsrs(m)**2+cospsrs(m)**2)/dfloat(nframe)
            qpsrav=qpsrs(m)/dfloat(nframe)
            qsd=math.sqrt(abs(qpsr2s(m)/dfloat(nframe)-qpsrav**2))
            write (6,2011) m,psr(m),cvm,qpsrav,qsd
            write (iw0,2011) m,psr(m),cvm,qpsrav,qsd
          
        elif (ianaltyp  ==  40) :
          avnumsolv=float(numsolvsum)/float(nframe)
          write (6,2109) avnumsolv,numsolvmin,numsolvmax
          write (iw0,2109) avnumsolv,numsolvmin,numsolvmax
        ## end if
        iout_track=0
        if (ianaltyp  ==  5  or  ianaltyp  ==  6  or 
     -      ianaltyp  ==  7  or  ianaltyp  ==  34) :
          if (ireadtracks  ==  0) :
            call askyn(
     -        'Do you want to write the bond tracks',36,1,-1,itrackw,0,
     -        2)
            if (itrackw  >  0) :
              call getname(trackfile,ltrackfile,
     -          'Name of the bond track file',27,80,'',0,0,000,0)
              iout_track=92
              call openfile(92,0,' ',1,'new',trackfile,ltrackfile,
     -          notfnd,0,1,1,0,0)
              call writetrack(iout_track,iw0,30,nbfoundorig,nres2d,
     -          trackfile,ltrackfile,ianc_anc)
            ## end if
          ## end if
          iwconn=0
          if (ilastframe  ==  0) write (6,2122) 'input structure'
          if (ilastframe  ==  1) write (6,2122)
     -      'last structure read from the trajectory file'
          question='Do you want a PDB file of the current structure '//
     -      'with % bonds'
          call askyn(question(1:60),60,1,-1,iwpdb,0,0)
          if (iwpdb  ==  0) :
            call askyn('Do you want a PDB file with lines between the'//
     -        ' residue pairs',59,1,-1,iwconn,140,0)
          else:
            call askyn(
     -        'Do you want to add lines between the residue pairs',50,
     -        1,-1,iwconn,140,0)
            iwpdb=1
          ## end if
          if (iwpdb  >  0) :
            call askyn(
     -        'Do you want percentages added to all atoms in a residue',
     -        55,1,-1,icperc,138,0)
            if (icperc  ==  1) :
              for ir in range(0, numres):
                sum=0.0
                for ia in range(ifres(ir), ilres(ir)):
                  sum=sum+cv(ia)
                
                if (sum  >  0.0) :
                  for ia in range(ifres(ir), ilres(ir)):
                    cv(ia)=sum
                  
                ## end if
              
            ## end if
            analfile4=analfile
            lanalfile4=namleno
            analfile4(lanalfile4-2:lanalfile4+4)='pcb.pdb'
            lanalfile4=lanalfile4+4
            iw4=44
            call openfile(iw4,0,'Percent bond labeled PDB',24,
     -        'new',analfile4,lanalfile4,notfnd,0,1,1,0,0)
            call writeconf(iw4,inpcrdtyp,iobpdb,inpcrdtyporg,
     -        nslt,nslt,nslt,naslv,islvw,iasv,namesv,qsv,pflsv,1,
     -        iwhead,0,iatnum,ifchrg,nconfig,innlist,c,rprox,cv,ixres,
     -        iresno,atnames,resnames,segnames,charge,isegno,altcol,
     -        inscol,ninsres,marker,ntitlin,ntitlinw,title,ireseq,
     -        iresnrestart,iresidrestart,nneig,nneiga,nhbneig,ineig,
     -        nhneig,nnneig,ncneig,nsneig,npneig,numres,numslv,
     -        resnamslv,line,blankline,mmtype,ibnd,index,indexn,indexo,
     -        1,molresflag,irescount3,itemp1,hblimfac,angmin,0,1,1,1,0,
     -        3,iqspaceask,ianaltyp,0,0.0,0,0,1,keeprem,iwriteatsym,
     -        radtodeg,maxrepconf,maxng,maxrsd,maxrec)
            if (iwconn  >  0) :
              call writeconn(iw4,ifres,ilres,line,index,inamcol1,
     -          inamcol2,c,n,iresbondcorr,itemp4,nframe,iw0,mxbonds,
     -          maxrsd,maxrec)
            ## end if
            write (iw4,2005) '# end'
            close (iw4)
          ## end if
        ## end if
      ## end if
9005  write (iw0,*)
      call datprt(iw0,version,0,mark0,lmark0,hostname,lhostname,
     -  iheadnode,0)
      close (iw0)
      n=natsorig
      iconfirmname=1
      call testconst(0,1,2,0.0,1.0,2.0,6,nfail,1,'ANTD')
      go to 9000
      return
2000  format[,' List of 1-4 neighbors, their distances and the ',
     -  'torsion angles will be',/,' written to ',a)
2001  format(i2)
2002  format(' Hydrogen-bond',a,'list will be written to file',/,
     -  5x,a)
2003  format[,20x,' DSSP (Kabsch-Sander) ANALYSIS',]
2004  format[,20x,' HELIX ORIENTATION ANALYSIS',//,
     -  ' S and E: coordinates of the helix axis initial and ',
     -  'final points (in A)',/,
     -  ' RMS: diagnostic of the irregularities in the helix',/,
     -  ' Len: length of the helix (in A)',/,
     -  ' D: unit vector in the helix direction',/,
     -  ' D-X,D-Y,D-Z angles: the tilt angles of the helix w.r.t. the',
     -  ' laboratory frame',/,
     -  ' Shape: Bent, Random or Oscillating, ',
     -  'based on Nup/dn (nup, ndown), and Ncross',/,
     -  ' Tolerance for being on the axis=',f6.2,
     -  ' (Nax: number of CAs within tolerance)'/,
     -  ' Rc: radius of circle fitted to the alpha carbons',/
     -  ' TPR: turn angle/residue',/,
     -  ' C: coordinates of the center of mass of the helix',/,
     -  ' N-X,N-Y,N-Z angles: the tilt angles of the normal to the ',
     -  ' plane fitting',/,
     -  '    the alpha carbons w.r.t. the laboratory frame',/,
     -  ' Rotation: the angle of revolution of helix around ',
     -  'its axis from the start',/,
     -  ' Local tilt: angle between the first and current ',
     -  'helix axes',/,
     -  ' SD: fluctuation of the rotation angles calculated from ',
     -  'each alpha carbon',/,
     -  ' N/Nr angle: angle between the normals N in the ',
     -  'current and reference state',]
2005  format(a)
2006  format(a,a,' file analyzed:',a)
2007  format(' Helix axis calculation Copyright 1996 Jon A.Christopher',
     -  ' and Thomas O. Baldwin.',/,
     -  25x,'Computers in Chemistry Vol. 20, pp 338-349 (1989)',/,
     -  ' Axes generated with the algorithm of Kahn ',/,
     -  25x,'Computers in Chemistry Vol 13, pp 185-189 (1989)')
2008  format(' Average mean distance from the ring plane=',f6.2,
     -  ' A SD=',f6.2)
2009  format(' Average general puckering coordinates ',
     -  '(Cramer \ Pople)=',/,7x,'m   angle    CV      qm     SD')
2010  format(' Reference structure:')
2011  format(' SUM ',i3,f8.2,f7.4,2f8.4)
2012  format[,' Reference structure file:',a)
2013  format(' Note: RMSD calculation is working with trajectory ',
     -  'scan only',/,' but a set of PDB or Charmm CRD files can be ',
     -  'read as MMC trajectories')
2014  format(a,/,a,'System:',/,a,a)
2015  format(' Residue range of the helix ',i3,':',i5,' - ',i5,
     -  ' (segment ',i3,') Center residue:',i6)
2016  format[,' Secondary structure element type:',/,
     -  7(9x,a1,': ',a,],9x,'?',': Unrecognized',]
2017  format(' The COM of each trajectory frame will be shifted to the',
     -  ' reference frame COM')
2018  format(' Each trajectory frame will be overlaid on the reference',
     -  ' structure',a)
2019  format(1x,a,' will be written to file ',/,5x,a)
2020  format(' Scan of trajectory ',a,' finished',/,
     -  ' Number of frames read=',i6,/,' Number of frames analyzed=',i6,
     -  a,/,' Number of frames skipped=',i6)
2021  format[,' NOTE: map preparation would be much faster if you ',
     - 'used a trajectory',/,' that has only the frames for which ',
     - 'calculations are required.',/,' Use the trajectory conversion ',
     - 'option of Simulaid to do this')
2022  format(' NOTE: Mass-weighting is turned off')
2023  format(i5,1x,a4,1x,a3)
2024  format[,' List of neighbors, bondlengths, angles and torsions ',
     -  '(if requested) will be',/,' written to file ',a)
2025  format(' Solvent radius=',f7.2,' A',/,
     -  ' Number of random point generated=',i10)
2026  format(' Only every ',i2,'-th grid point will be used')
2027  format('REMARK ',a)
2028  format(' Calculation of distance distribution over selected atom',
     -  ' pairs',/,a)
2029  format(' NOTE: residue numbers refer to the actual residue ',
     -  'number in the input structure')
2030  format(' WARNING: The presence of a large number of solvents ',
     -  'slows down the calculation',/,10x,'- you may want to ',
     -  'eliminate the solvents using the Edit option')
2031  format(' PDB file with RMSF as the B factor:',a)
2033  format('RMSD (X) - Energy (Y) scatterplot. Corr=',f8.5)
2034  format(' Number of atoms changed from ',i6,' to ',i6,
     -  ' (trajectory)')
2035  format(' NOTE: lack of charge information may hamper the bond ',
     -  'definitions',/,7x,'- you may want to consider other structure',
     -  ' file formats')
2036  format(' Simulaid will try to open file ',a,/,
     -  10x,'to write the results on. If this is the file ',
     -  'that you want to read,',/,
     -  10x,'make sure not to overwrite it')
2037  format(' ##### WARNING: This calculation may be quite time ',
     -  'consuming',/,
     -  7x,'Please, log on to a compute node and restart')
2038  format(' To match clusters of the two trajectory, generate the ',
     -  '2D RMSD map for both',/,' and rerun the cross-RMSD ',
     -  'calculation using the already calculated',/,
     -  ' 2D RMSD and cross-RMSD maps')
2039  format(1x,a1,':',i5,' (',a,1x,a,')',3(' -',i5,' (',a,1x,a,')'))
2040  format(' Residue range of the helix:',i5,' - ',i5,' (segment ',i3,
     -  ')')
2041  format(1x,a,': Invalid selection, try again')
2042  format(' Bridge residue name=',a,' Number of atoms in residue ',a,
     -  '=',i3)
2043  format(i3)
2044  format(100i1)
2045  format(' Default marks will be used:',/,9(i2,': ',a1))
2046  format(' Adjacency-matrix analysis result will be written to ',
     -  'file',/,1x,a)
2047  format(' Enter the indices of the two atoms whose distance is ',
     -  'to be calculated/tracked',/,
     -  ' To track the diffusion of an atom, type the same index',
     -  ' for both atoms')
2048  format(' N=',i6,' c(',i5,')=',3f11.2,' A')
2049  format(' N=',i6,' c(',i5,')=',3f11.2,' c(',i5,')=',3f11.2,' A')
2050  format[,' Residue neighbour list based on ',a,' atoms will be ',
     -  'written to',/,1x,a)
2051  format(' First atom is on the reference structure ',a)
2052  format(i5,' diffusion')
2053  format(i5,' - ',i5,' distance')
2054  format(' N=',i6,' D=',f11.2,' D^2=',f11.2,' r=',3f11.2,' A')

2056  format(' Symbol for residue type ',i1,'=',$)
2057  format(' Calculating distance between atom',i6,' residue',i5,
     -  ' (',a5,a5,') and',/,' atom',i6,' residue',i5,' (',a5,a5,')')
2058  format(' Calculating distance between atoms',i6,' and',i6)
2059  format(' WARNING: number of frames exceeded the parameter ',
     -  'MAXFRAMES (',i8,')',/,10x,'Several analysis options fail when',
     -  ' this limit is passed')
2060  format(1x,a,' first and last frame and increment:',3i7)
2061  format(' Check for unphysical distances',/,
     -  ' Clash threshold = Bondlength threshold * sqrt(',f8.4,')',/,
     -  ' Short bond threshold = Bondlength threshold * sqrt(',f8.4,')',
     -  /,' Index distance threshold for non-S atoms = ',i4)
2062  format(' NOTE: The number of solute residues (',i6,') exceeds ',
     -  'the limit (',i5,')',/,' Increase the parameter MAXDISTR and',
     -  ' recompile to calculate RMSF')
2063  format(' WARNING: ',a,' trajectory name (',a,')',/,8x,
     -  'is different from the trajectory name of the corresponding ',
     -  '2D-RMSD map',/,8x,'(',a,')')
2064  format[,' RMSD values were read from file ',a)
2065  format(' Number of frames used from the second trajectory=',i6)
2066  format(' WARNING: atoms',i6,' and',i6,' are',a,'bonded',/,
     -  ' You can create bonds with the Edit option')
2067  format(a1,i2,' i:',i4,3('-',i4))
2068  format(' NOTE: after quitting clustering you will have the ',
     - 'option of inspecting',/,6x,' the map and repeat the clustering')
2069  format('Proline kink dials for residue',i5,
     -  ' Pre and post proline helix length=',2i3)
2070  format(' Configuration ',i6,': all solvents were filtered out')
2071  format[,' Calculation of proline kink angles defined by ',
     -  'Visiers et al.',/,
     -  ' Note: kink residue does not have to be a proline',/,
     -  ' Axes generated with the algorithm of Kahn ',/,
     -  ' using the defs of J.A. Christopher \ T.O. Baldwin',/,
     -  ' Angles are calculated with vector-algebra based code',/,
     -  ' B# end, angle, wobble angle and phase shift will',
     -  ' be written to file',/,1x,a)
2072  format(i6,' B=',f6.1,' W=',f6.1,' FS=',f6.1,' W-FS=',f6.1,
     -  ' PR=',f6.1,' axis RMS b,a=',f4.2,f5.2)
2073  format(' The Postscript plots will be written to the file',/,
     -  1x,a)
2074  format(' Plotting is limited to ',i6,1x,a,/,
     -  ' -',a,'increase the the arrays (maxframe)',/
     -  ,' in the common block /analres/')
2075  format(1x,a,' map values will be written to file',/,
     -  1x,a,/,' Postscript plot will be written to file',/,1x,a)
2076  format(1x,a,' will be written into the data column ',
     -  'in file',/,1x,a)
2077  format[,20x,' PROLINE KINK CALCULATION',/,
     -  ' b# end (B), wobble (W), phase shift (FS)',/,
     -  ' pseudorotation angle of the proline (PR)',/,' Goodness of ',
     -  ' the helices before and after the proline (axis RMS b,f)')
2078  format[,' Calculation of pseudorotation angles')
2079  format[,' Average radius of gyration=',f10.3,' S.D.=',f10.3,' A',
     -  /,' Average hydrodynamic radius=',f10.3,' S.D.(1/rHD)=',f10.3,/,
     -  ' Average ratio of the largest and smallest moments of ',
     -  'inertia=',f10.3)
2080  format[,' Average cell dipole moment=',f10.3,' S.D.=',f10.3,/
     -  ' Average solute dipole moment=',f10.3,' S.D.=',f10.3,' au*A')
2081  format[,' NOTE: warnings, summaries (if any) will be turned off',
     -  ' after the',i3,'-th frame')
2082  format(' Calculation of the RMSD from the input structure',/,
     -  ' both after obtaining the best fit with the Kabsch ',
     -  'method and without fitting',]
2083  format(' Kink residue:',i5)
2084  format(' Calculation of the RMSD matrix between frames of a ',
     -  'trajectory')
2085  format[,1x,78('-'))
2086  format(' Plot entries color coded by the ',a,' number:',/,
     -  ' Red-Yellow-Green-Cyan-Blue <--> 1 -',a,a,i5,')')
2087  format('Residues',i6,' -',i6)
2088  format(1x,a,' dipole components:',/,
     -  3(' dip',a,'=',f8.2,' +/-',f6.2))
2089  format(' Number of ',a,' helices: ',i7)
2090  format(' WARNING: The calculation would be significantly faster',
     -  ' if a trajectory were first',/,10x,'extracted with only the ',
     -  'frames needed.',/,10x,'- it can be done with the trajectory ',
     -  'conversion option')
2091  format('Atomic number of solvent atom ',i2)
2092  format('Partal charge of solvent atom ',i2)
2093  format(' Residue',i5,' (',a,') <Phi>=',f8.2,' CV=',f6.4,
     -  ' <Psi>=',f8.2,' CV=',f6.4)
2094  format(i8)
2095  format(a1)
2096  format(' Current dimensioning allows only for',i6,' * ',i6,
     -  ' matrix')
2097  format(' Principal axes calculated will be written to file',/,
     -  1x,a)
2098  format(' The reference conformation may be the first ',
     -  'conformation of the trajectory'/,' (this is the default) or ',
     -  'the structure inputted at the start')
2099  format(' Calculation of principal axes. Atoms used:',$)
2100  format(' Calculation of radii of gyration, hydrodynamic radius ',
     -  ' moments of inertia and',/,
     -  ' (when charges are available) and dipole moment',/
     -  ' Atoms used:',$)
2101  format(' Angle between inital and current ',a,' principal axis')
2102  format(' Subsequent structures will be saved into files',/,
     -  1x,a,'*N',a)
2103  format(i5)
2104  format(' Residue # ',i6,' has no full backbone - a bond may be ',
     -  'missing',/,' You can add bonds usind the Edit option')
2105  format(' Solvents will be sorted by increasing CV (out -> in)')
2106  format(' ERROR: dimensions of the cross-MSD matrix (',i5,' x ',i5,
     -  ')',/,8x,'differ from the lengths of the trajectories read:',
     -  2i6)
2107  format(' Radius of gyration, hydrodynamic radius moments of ',
     -  'inertia',/,' and (if charges are available) dipole moment ',
     -  'will be written to file'/,1x,a)
2108  format('A ',a4,'-',a4,'(',a4,') -',a4)
2109  format(' Average number of filtered solvents=',f8.1,/,
     -  ' Range: [',i6,',',i6,']')
2110  format(' Filtering summaries are written to file ',a,'.flt')
2111  format(' Segment',i4,' residue #',i5,' is ',a)
2112  format(' At least one of the residue ranges exceeds MAX2D (',i5,
     -  ')',/,' - no averages will be calculated')
2113  format[,'List of the counts of different tyes of contacts',]
2114  format(i5,' resid=',i5,' (',a,1x,a,') Nrepr=',i6,' Nclose=',i6,
     -  ' Ncontact=',i6)
2115  format(' New value of ramax(',a2,')')
2116  format(' RMSDs were calculated ',a)
2117  format(' Actual residue number of the ',a,' residue=',i6)
2118  format(' Reference structure file=',a)
2119  format(' Sigma-r plot (SD of interatomic distances'/,
     -  ' Hao Zhou, to be published')
2120  format(' Molecule-molecule distance matrix based on ',a,/,
     -  ' will be written to ',a)
2121  format(' Molecule-molecule distance matrix based on ',a,]
2122  format(' Current structure is the ',a)
2123  format(' ERROR: residue number can not be ',a,' than ',i5)
2124  format(' dist: closest approach distance of two helices',/,
     -  ' cc_dist: distance between two helix centers',/,
     -  ' CAs: the atom numbers of the alpha carbons closest ',
     -  ' to the',/,4x,'point on the helix axis nearest to the other ',
     -  'helix',/,' ang: The angle between two helix axes',/,' dhang: ',
     -  'The torsion angle formed by two helix axes around the line ',
     -  'perp# endicular to them',/,' TM: transmembrane helix',/,
     -  'IN: intracellular, EX: extracellular')
2125  format(' A solvent is considered to be the interface if',/,
     -  ' (a) its CV w.r.t the complex is > ',f6.4,/,
     -  ' (b) Dist(R1,R2) <',f6.1,' where R1, R2 are the two nearest ',
     -  'solute heavy atoms',/,
     -  ' (c) When Dist(R1,R2) <= 6 A, the angle R1-SLV-R2 > ',f10.2,/,
     -  ' (d) When Dist(R1,R2) >  6 A, the angle R1-SLV-R2 > ',f10.2,/,
     -  ' (e) The ratio of the CVs w.r.t the two solute molecules < ',
     -  f6.3,/,' (e) There is at least one solvent or solute atom ',
     -  'within ',f5.1,' A')
2126  format(' Number of filtered solvents found=',i7)
2127  format(' Analysis of ',a,'s')
2128  format(' Mutually proximal atom pairs are dropped if their ',
     -  'distance is > ',f4.1,' A')
2129  format(' Helix',i3,' direction established from the reference ',
     -  'structure as the ',a,' axis',a)
2130  format(' Angle between the ',a,'-axis and the helix=',f5.1)
2131  format('CONECT',2i5)
2132  format(' Residue list selected is not contiguous - RMSF ',
     -  'calculations will be skipped')
2133  format(' ERROR: plot range is not within the residue range [',
     -  i6,',',i6,']')

2140  format(' Potentially a hydrogen bond is formed between')
2141  format(10x,'Calculation of all hydrogen bonds',//,
     -  ' Each line describes a hydrogen bond, formed between ',
     -  ' acceptor A',/,' and donor H (bonded to the heavy atom D) ',
     -  'as follows:',/,5x,'Atom number, atom name, residue number, ',
     -  'residue name for donor atom H',/,
     -  5x,'Atom number, atom name, residue number, ',
     -  'residue name for acceptor atom A',/,
     -  5x,'Hydrogen bond length rHA; hydrogen bond angle D-H-A; and ',
     -  'distance rDA;',/,5x,'Hydrogen bond type X-Y, where X and Y ',
     -  'can be B, S, V, ?',/,5x,'indicating backbone, sidechain, ',
     -  'solvent or undetermined origin',/,
     -  ' For each hydrogen, only the two shortest bond will be ',
     -  'retained',]
2142  format('   all ',a2,' atoms within',f6.2,' A of a hydrogen')
2143  format('   all ',a2,' atoms within',f6.2,' A of a water oxygen')
c2147  format(' WARNING: x and z cell change factors differ:',2f10.7,/,
c     -  10x,'Cell change may be anisotropic')
2148  format('Number of random points/cell=',i7,
     -  '    Correlation coefficient=',f6.3)
2149  format[,' Distribution of the different SS elements:',/,
     -  (1x,a1,' : ',a))
2150  format(1x,'Residue #',9(4x,a1))
2151  format(i10,9i5)
2152  format('Number of colors ( <',i2,')')
2153  format(' Number of atoms used for ',a,' calculation:',i6,
     -  ' (out of',i6,')',a,' atom list:'/,(15i5))
2154  format(' All atoms are used for ',a,' calculation')
2155  format(i5,' RMSF(Cav)=',f6.2,' RMSF(Cref)=',f6.2,' SD=',f6.2,a,
     -  'CI=',f6.2)
2156  format(' NOTE: residue indices on the plot are residue SEQUENCE',
     -  ' numbers',/,7x,'- see the map to the reside and segment ',
     -  'numbers on the',/,7x,'output of the bond history calculation')
2157  format(' Number of frames read in a block=',i5)
2158  format(' Backbone ',a,' rn=',i5,' res=',a,i5,
     -  ' Atnos for CA, C, N=',3i6)
2159  format(' Specification of helix # ',i2)
2160  format(1x,a,'elix #',i2)
2161  format(' ERROR: number of data items (',i5,') exceeds limit (',
     -  i5,')',/,' - reduce the number of helices or increase the ',
     -  'parameter MAXCOPY1')
2162  format('HX',i2,'-HX',i2,':',a,' dist.',a,'(',a2,')')
2163  format(' HX#',i3)
2164  format(' Filtering criterion: ',a)
2165  format(' Solute solvent distance threshold=',f6.1)
2166  format(' CV threshold=',f6.1,' CV calculation radius=',f6.1)
2167  format[,' Filtering solvents')
2168  format(' WARNINGN: number of solvent atoms (',i7,'-',i6,') is ',
     -  'not divisible',/,' by the number of atoms/solvent molecule (',
     -  i3,')')
2169  format('HX',i2,'-HX',i2,' angle         ')
2170  format('HX',i2,'-HX',i2,' torsion angle ')
2171  format('HX',i2,' rot vs HX',i2)
2172  format[' HELIX',i3,' statistics')
2173  format[,' Config #',i6,' (inp #: ',i10,')',a,
     -  'Ref #',i6,' (inp #:',i10,')')
2174  format('HX',i2,'-HX',i2,' (',a2,'):')
2175  format(' Trajectory frame ',i6,' Energy=',e12.5)
      # end
