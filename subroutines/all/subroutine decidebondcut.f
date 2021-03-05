      subroutine decidebondcut(ian1,ian2,rlim)
      common /connatdat/ ramax(99),ramax2(99),hlimfac,ianfg(99),
     -  namfcg(100),nrmw
      rlim=amax1(ramax2(ian1),ramax2(ian2))
      if (ian1 .eq. 1 .or. ian2 .eq. 1) rlim=rlim*hlimfac
      return
      end
