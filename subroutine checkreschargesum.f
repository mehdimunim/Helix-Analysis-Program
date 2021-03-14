      subroutine checkreschargesum(nslt,iresno,isegno,line,index,
     -  irescol1,irescol2,iresncol1,iresncol2,isegcol1,isegcol2,iqcol1,
     -  iqcol2,charge,molsltlim,label,llabel,resnames,ifres,
     -  ixres,ixchrg,ixnochrg,nreschrg,iallzero,iallnonpos,iallnonneg,
     -  ichecksavedq,ifixq,reportqmin,iout,maxrsd,maxrec)
      dimension index(maxrec),iresno(maxrec),ifres(maxrec),
     -  isegno(maxrec),ixchrg(maxrec),ixnochrg(maxrec),charge(maxrec),
     -  molsltlim(3,maxrsd),ixres(maxrec)
      character*(*) label
      character*8 resnames(maxrsd)
      character* 132 line(maxrec)
      character*8 rnprev
      character*6 rnuprev
      real*8 qsum
      character*19 segline
      data segline /' Segment     (    )'/
c     Check if residue charges are integers, print charge sum
c     print *,'CHECKRESCHARGESUM isegcol1,2=',isegcol1,isegcol2
      nrescol=irescol2-irescol1+1
      nresncol=iresncol2-iresncol1+1
      iresprev=iresno(1)
      ixresprev=ixres(1)
      isegprev=isegno(1)
      ndigits=3
      nsegm=isegno(nslt)
      isegnoprev=0
      if (reportqmin .lt. 0.0) call getreal(
     -  'Minimum deviation from integral residue charge to print',55,
     -  0.0,reportqmin,1,000)
      do isg=1,nsegm
        iallzero=1
        iallnonpos=1
        iallnonneg=1
        nreschrg=0
        nresnochrg=0
        nresck=0
        nonintegral=0
        iallzerores=1
        qrsum=0.0
        qsum=0.d0
        qdevmax=0.d0
c       print *,'isg=',isg,' mlim=',molsltlim(1,isg),molsltlim(2,isg)
        do ia=molsltlim(1,isg),molsltlim(2,isg)
          ir=iresno(ia)
          is=isegno(ia)
          if (nsegm .gt. 1 .and. is .gt. isegnoprev) then
            write (segline(9:12),1006) isg
            ncol=12
            if (isegcol2 .gt. isegcol1) then
              segline(15:15+isegcol2-isegcol1)=
     -          line(index(ia))(isegcol1:isegcol2)
              ncol=19
            end if
            write (6,1008) segline(1:ncol)
          end if
          isegnoprev=is
          q=charge(ia)
          if (q .lt. 0.0) iallnonneg=0
          if (q .gt. 0.0) iallnonpos=0
          if (q .ne. 0.0) iallzero=0
          if (q .ne. 0.0) iallzerores=0
          if (ir .ne. iresprev .or. is .ne. isegprev .or.
     -        ia .eq. molsltlim(2,isg)) then
            rnprev(1:nrescol)=line(index(ia-1))(irescol1:irescol2)
            rnuprev(1:nresncol)=line(index(ia-1))(iresncol1:iresncol2)
            qscheck=qrsum
c           resnames(iresprev)=rnprev
            if (ia .eq. nslt .or. ia .eq. molsltlim(2,isg))
     -        qscheck=qscheck+q
            if (isinteger(qscheck,ndigits,qoffset) .eq. 0) then
              if (abs(qoffset) .ge. reportqmin)
     -          write (6,1000) label(1:llabel),iresprev,
     -            rnuprev(1:nresncol),rnprev(1:nrescol),qscheck
              nonintegral=nonintegral+1
              if (abs(qoffset) .gt. qdevmax) qdevmax=abs(qoffset)
            else if (abs(qscheck) .gt. 0.1) then
              nreschrg=nreschrg+1
              ixchrg(nreschrg)=ixresprev
            end if
            if (iallzerores .eq. 1) then
              nresnochrg=nresnochrg+1
              ixnochrg(nresnochrg)=ixresprev
            end if
            qrsum=0.0
            iresprev=ir
            ixresprev=ixres(ia)
            isegprev=is
            nresck=nresck+1
            iallzerores=1
          end if
          qrsum=qrsum+q
          qsum=qsum+q
        end do
        if (iallzero .eq. 1) then
          write (6,1005) ' ',isg
        else
          if (iallnonneg .eq. 1) write (6,1005) ' positive or ',isg
          if (iallnonpos .eq. 1) write (6,1005) ' negative or ',isg
        end if
        if (iallzero+iallnonpos+iallnonneg .eq. 0 .or. nslt .le. 10)then
          nerr=0
          if (ichecksavedq .eq. 1) then
            do ia=molsltlim(1,isg),molsltlim(2,isg)
              read(line(index(ia))(iqcol1:iqcol2),*,err=555,end=555) q
              if (abs(q-charge(ia)) .gt. 0.0001 .and. nerr .lt. 10) then
                write (6,1004) ia,charge(ia),q
                nerr=nerr+1
              end if
              go to 556
555           nerr=nerr+1
              if (nerr .le. 10) print *,'ERROR: Invalid charge on atom',
     -          ia,':',line(index(ia))(iqcol1:iqcol2)
              if (nerr .eq. 10)
     -            print *,'Further error messages are suppressed'
556           continue
            end do
          end if
          write (iout,1001) nresck,label(1:llabel),qsum
          if (nonintegral .gt. 0) write (iout,1009) nonintegral
          if (nerr .gt. 0)
     -      print *, 'WARNING: ',nerr,' atom records had errors'
          qsum4=qsum
          if (isinteger(qsum4,ndigits,qoffset) .eq. 0)
     -      print *,'WARNING: Total charge is not integer'
          if (nresnochrg .gt. 0 .and. iallzero .eq. 0) write (iout,1007)
     -      (ixnochrg(ir),resnames(ixnochrg(ir))(1:nrescol),
     -      ir=1,nresnochrg)
          if (nreschrg .gt. 0) write (iout,1002)
     -      (iresno(ifres(ixchrg(ir))),resnames(ixchrg(ir))(1:nrescol),
     -       ir=1,nreschrg)
          if (nonintegral .gt. 0) then
            if (nresck .gt. 1) write (iout,1003) qdevmax
            if (ifixq .eq. -1) call askyn(
     -      'Do you want to redistribute charges to integral resid sums'
     -        ,58,1,-1,ifixq,126,0)
            if (ifixq .gt. 0) then
              qmin=10.0**(-ndigits)
              qround=10.0**(-ndigits-2)
              qrsum=0.0
              qcorrsum=0.0
              iresf=1
              iresprev=iresno(1)
              isegprev=isegno(1)
              do ia=molsltlim(1,isg),molsltlim(2,isg)
                ir=iresno(ia)
                is=isegno(ia)
                read(line(index(ia))(iqcol1:iqcol2),*) q
                if (ir .ne. iresprev .or. is .ne. isegprev .or.
     -              ia .eq. molsltlim(2,isg)) then
                  qscheck=qrsum
                  if (ia .eq. molsltlim(2,isg)) qscheck=qscheck+q
                  ii=isinteger(qscheck,ndigits,qoffset)
                  qcorrsum=qcorrsum-qoffset
                  if (ii .eq. 0) then
c                   Divide extra charge (qdistr) into small chunks
                    ntotchunks=abs(qoffset)/qmin
                    qmins=qmin
                    if (qoffset .gt. 0.0) qmins=-qmin
                    qdistr=ntotchunks*qmins
                    qdistr1=(ntotchunks+1)*qmins
                    if (abs(qoffset+qdistr) .gt.
     -                  abs(qoffset+qdistr1))then
                      qdistr=qdistr1
                      ntotchunks=ntotchunks+1
                    end if
                    natsres=ia-iresf
                    nchunks=ntotchunks/natsres
                    nchunkext=ntotchunks-natsres*nchunks
                    do ja=iresf,ia-1
                      qadd=nchunks*qmins
                      if (nchunkext .gt. 0 .and. ja-iresf .lt.nchunkext)
     -                  qadd=qadd+qmins
                      qfix=charge(ja)+qadd
                      call putreal(line(index(ja))(iqcol1:iqcol2),
     -                  iqcol2-iqcol1+1,qfix,ndigits)
                      charge(ja)=qfix
                    end do
                  end if
                  qrsum=0.0
                  iresprev=ir
                  isegprev=is
                  iresf=ia
                end if
                qrsum=qrsum+q
              end do
              if (abs(qcorrsum) .gt. abs(qmin)) print *,
     -          'WARNING: total charge distributed is not zero:',
     -          qcorrsum
            end if
          end if
        end if
      end do
      return
1000  format(' Charge sum on ',a,i6,' (',a,1x,a,') is not integer:',
     -  f10.4)
1001  format(' Checked ',i6,1x,a,'s for charge sum. Total charge=',
     -  f10.5,' e')
1002  format(' Charged residues (residues with nonintegral charges ',
     -  'omitted):',/,5(i5,'(',a,')'))
1003  format(' Largest absolute deviation from integral charge=',f8.5)
1004  format(' Program ERROR: solute atom',i6,' charge saved=',f8.5,
     -  ' charge read=',f8.5)
1005  format(' WARNING: all charges are',a,'zero in segment',i4)
1006  format(i4)
1007  format(' WARNING: the following residues have no charges:',/,
     -  5(i5,'(',a,')'))
1008  format(a)
1009  format(' Number of residues with nonintegral charge sum=',i4)
      end
