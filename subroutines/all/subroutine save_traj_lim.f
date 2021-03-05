      subroutine save_traj_lim(ifirst,ilast,incr,i12)
      character*200 trajnam12
      common /trajname/ trajnam12(2),ltrajnam12(2),ifirsttraj12(2),
     -  ilasttraj12(2),incrementtraj12(2)
      ifirsttraj12(i12)=ifirst
      ilasttraj12(i12)=ilast
      incrementtraj12(i12)=incr
      return
      end
