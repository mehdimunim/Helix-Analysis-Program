      subroutine buildbranch(c,nbranch,ia_branch,
     -  ia_predef,r_predef,ba_predef,ta_predef,pi,maxat,maxbranch)
      dimension c(3,maxat),ia_branch(maxbranch),ia_predef(3,maxbranch),
     -  r_predef(maxbranch),ba_predef(maxbranch),ta_predef(maxbranch)
      do ia=1,nbranch
c       print *,'IA,IA1,2,3=',ia,ia_predef(1,ia),ia_predef(2,ia),
c    -    ia_predef(3,ia)
        call addatom(1,c(1,ia_predef(1,ia)),c(1,ia_predef(3,ia)),
     -    c(1,ia_predef(2,ia)),c(1,ia_predef(1,ia)),c(1,ia_branch(ia)),
     -    r_predef(ia),ba_predef(ia),ta_predef(ia),0.0,ia_branch(ia),
     -    pi,0,ifail)
      end do
      return
      end
