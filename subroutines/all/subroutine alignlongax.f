      subroutine alignlongax(c,n,edge,rot,longax)
      dimension c(3,n),edge(3),rot(3,3)
      dimension permute(3,3)
      character*1 xyz(3)
      data xyz /'X','Y','Z'/
      edgelong=amax1(edge(1),edge(2),edge(3))
      if (edge(longax) .eq. edgelong) return
      call zeroit(permute,9)
      k=1
      longedge=0
      do while (longedge .eq. 0 .and. k .le. 3)
        if (edge(k) .eq. edgelong) longedge=k
        k=k+1
      end do
      permute(longedge,longax)=1.0
      print *,'longedge,longax=',longedge,longax
      kx=mod(longedge,3)+1
      ky=mod(longax,3)+1
      permute(kx,ky)=1.0
      print *,'kx,ky=',kx,ky
      kx=mod(kx,3)+1
      ky=mod(ky,3)+1
      permute(kx,ky)=1.0
      print *,'kxx,kyy=',kx,ky
      write (6,2001) ((permute(i,j),j=1,3),i=1,3)
      write (6,2000) xyz(longax),xyz(longedge)
      call rotate_c(c,n,permute,c,'ALIGNLONG',9)
      call rotate_c(edge,1,permute,edge,'ALIGNLONG1',10)
      call matprod(permute,rot,rot)
      return
2000  format(' System is rotated to make the ',a,' axis the longest ',
     -  '(current longest: ',a,' axis)')
2001  format(' Permute:',/,(3f4.1))
      end
