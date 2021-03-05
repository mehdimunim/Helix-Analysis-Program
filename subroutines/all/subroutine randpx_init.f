      subroutine randpx_init(iseed)
c*****Set the seed of the random number generator (for reproducibility)
      integer*4 ixo
      common /rangen/ ixo
      ixo=iseed
      return
      end
