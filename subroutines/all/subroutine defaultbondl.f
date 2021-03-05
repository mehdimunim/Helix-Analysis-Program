      subroutine defaultbondl(nni1,nni2,ian1,ian2,r12)
c     Establish default bondlength based on atomic numbers and valences
c     To be extended
      if (ian1 .lt. ian2) then
        nngi1=nni2
        nngi2=nni1
        iatn1=ian2
        iatn2=ian1
      else
        nngi1=nni1
        nngi2=nni2
        iatn1=ian1
        iatn2=ian2
      end if
      r12=999999.0
      if (iatn2 .eq. 1) r12=1.08
      if (iatn1 .eq. 6) then
c       Carbon
c       if (nngi1 .eq. 4) then
          if (iatn2 .eq. 1) r12=1.08
          if (iatn2 .eq. 6) r12=1.53
c       else if (nngi1 .eq. 3) then
c         if (iatn2 .eq. 1) r12=1.03
c         if (iatn2 .eq. 6) r12=1.03
        else if (nngi1 .eq. 2) then
c         if (iatn2 .eq. 1) r12=1.03
c         if (iatn2 .eq. 6) r12=1.03
c       end if
      else if (iatn1 .eq. 7) then
c       Nitrogen
c       if (nngi1 .ge. 3) then
          if (iatn2 .eq. 1) r12=1.01
          if (iatn2 .eq. 6) r12=1.47
          if (iatn2 .eq. 7) r12=1.25
c       end if
      else if (iatn1 .eq. 8) then
c       Oxigen
        if (iatn2 .eq. 1) r12=0.96
        if (nngi1 .eq. 2) then
          if (iatn2 .eq. 6) r12=1.43
c         if (iatn2 .eq. 7) r12=1.14
c         if (iatn2 .eq. 8) r12=1.03
        else if (nngi1 .eq. 1) then
          if (iatn2 .eq. 6) r12=1.23
c         if (iatn2 .eq. 7) r12=1.03
c         if (iatn2 .eq. 8) r12=1.03
        end if
      else if (iatn1 .eq. 15) then
c       Phosphorus
c       if (iatn2 .eq. 1) r12=1.03
c       if (iatn2 .eq. 6) r12=1.03
c       if (iatn2 .eq. 7) r12=1.03
c       if (iatn2 .eq. 8) r12=1.03
      else if (iatn1 .eq. 16) then
c       Sulphur
        if (iatn2 .eq. 1) r12=0.96
        if (iatn2 .eq. 6) r12=1.81
c       if (iatn2 .eq. 7) r12=1.03
c       if (iatn2 .eq. 8) r12=1.03
      end if
      return
      end
