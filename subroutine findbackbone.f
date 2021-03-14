      subroutine findbackbone(n0,n,iwout,nneig,ineig,
     -  nneiga,ineiga,in1,in2,ir1,ir2,line,blankline,index,iback,
     -  isideback,iparent,list,maxng,maxbox,maxrec)
c#    Based on MMC findcent
c*****Finds the backbone of a molecule
      dimension nneig(n),nneiga(n),ineig(maxng,n),ineiga(maxbox,n),
     -  index(n),iback(n),isideback(n),iparent(n),list(n)
      character* 132 line(maxrec),blankline,linep
c     First grow neighbour list from n0
c     print *,'FINDBACKBONE n0,n,maxbox=',n0,n,maxbox
      ifirst=n0
      ilast=0
      ibb=0
      do while (ilast .lt. n)
        write (iwout,1005)
        call growchain(n0,n,0,ifirst,ial,nsteps,nneig,ineig,
     -    list,iparent,iback,ilast,maxng)
c       Now grow from ial
        call growchain(n0,n,0,ial,iaf,lbackb,nneig,ineig,
     -    list,iparent,iback,iii,maxng)
c       Save nn list for manipulations
        do ia=n0,n
          nneiga(ia)=nneig(ia)
          call trnsfi(ineiga(1,ia),ineig(1,ia),nneig(ia))
        end do
c       Drop backbone atoms from nn list to avoid loopbacks
        do ia=1,lbackb
          iba=iback(ia)
          if (ia .gt. 1) then
            do in=1,nneiga(iba)
              if (ineiga(in,iba) .eq. iback(ia-1)) then
                ineiga(in,iba)=ineiga(nneiga(iba),iba)
                nneiga(iba)=nneiga(iba)-1
                go to 201
              end if
            end do
          end if
201       if (ia .lt. lbackb) then
            do in=1,nneiga(iba)
              if (ineiga(in,iba) .eq. iback(ia+1)) then
                ineiga(in,iba)=ineiga(nneiga(iba),iba)
                nneiga(iba)=nneiga(iba)-1
                go to 202
              end if
            end do
          end if
202       continue
        end do
        do ia=1,lbackb
          iba=iback(ia)
          do in=1,nneiga(iba)
            ina=ineiga(in,iba)
            do ja=1,nneiga(ina)
              if (ineiga(ja,ina) .eq. iba) then
                ineiga(ja,ina)=ineiga(nneiga(ina),ina)
                nneiga(ina)=nneiga(ina)-1
                go to 200
              end if
            end do
200         continue
          end do
        end do
c        do ia=n0,n
c          write (77,6733) ia,(ineiga(i,ia),i=1,nneiga(ia))
c6733      format(i5,' reduced nn=',30i4)
c        end do
        do ia=1,lbackb
c         If there is a sidechain, get its backbone too
          iba=iback(ia)
          ibn1=0
          lchainmax=0
          do in=1,nneiga(iba)
            ibna=ineiga(in,iba)
            if (nneiga(ibna) .gt. 0) then
c             write (77,*) 'ia,in,ib,ibna=',ia,in,ib,ibna
              call growchain(n0,n,iba,ibna,ial,nsteps,nneiga,ineiga,
     -          list,iparent,isideback,iii,maxbox)
              if (lchainmax .lt. nsteps) then
                lchainmax=nsteps
                ichainmax=ibna
              end if
            else if (ibn1 .eq. 0) then
c             (First) monovalent neighbor
              ibn1=ibna
            end if
          end do
          if (lchainmax .gt. 0) then
c           Get the chain from ichainmax
c           write (77,*) 'Longest sidechain:'
            call growchain(n0,n,iba,ichainmax,ial,nsteps,nneiga,ineiga,
     -        list,iparent,isideback,iii,maxbox)
          end if
c         Print line(s)
          linep=blankline
          ibb=ibb+1
          write (linep(20:26),1003) ibb
          write (linep(27:44),1000) iba,
     -      line(index(iba))(ir1:ir2),line(index(iba))(in1:in2)
          if (ibn1 .gt. 0) then
            write (linep(1:18),1000) ibn1,
     -        line(index(ibn1))(ir1:ir2),line(index(ibn1))(in1:in2)
            linep(18:18)='-'
          end if
          if (lchainmax .gt. 0) then
            ic0=44
            do ic=1,lchainmax
              icis=isideback(lchainmax-ic+1)
c             write (linep(ic0+1:ic0+18),1001) icis,
              write (linep(ic0+1:ic0+18),1000) icis,
     -          line(index(icis))(ir1:ir2),line(index(icis))(in1:in2)
              linep(ic0:ic0)='-'
              ic0=ic0+18
              if (mod(ic,2) .eq. 0 .and. ic .ne. lchainmax) then
                write (iwout,1002) linep(1:79)
                linep=blankline
                ic0=44
              end if
            end do
          end if
          write (iwout,1002) linep(1:79)
        end do
        if (ilast .ne. n) write (iwout,1004)
        ifirst=ilast+1
      end do
      return
1000  format(i5,' (',a4,a5,') ')
c1001  format(' -',i4,' (',a4,a5,') ')
1002  format(a)
1003  format('BB',i4,':')
1004  format(/,'New backbone start',/)
1005  format(/,' Backbone search ',/,
     -  ' BB  seqno: atomindex (residue name atom name)')
      end
