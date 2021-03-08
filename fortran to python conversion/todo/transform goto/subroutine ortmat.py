import numpy as np

def ortmat(a, b, natoms, iverb):
	"""
	3-atom part originally written by P.K. Mehrotra.
	Compute the rotation matrix orm  by solving the following equation
						b=orm *a
	
	It is assumed that the coordinates of the molecule are given
	columnwise, i.e., the coordinates of the atom 1 constitute
	the first column, coordinates of the atom 2 constitute the
	the second column, etc.
	the elements of a and b are preserved.

	***
	Parameters:
	a: coordinates of the molecule before rotation (local system)
	b: coordinates of the molecule after rotation (global system).

	Returns:
	orm: rotation matrix (only rotation not translation)
	
	"""
	if (natoms .eq. 2) go to 20
      do i=1,3
        do k=1,3
          d(k,i)=a(k,i)-a(k,1)
          e(k,i)=b(k,i)-b(k,1)
        end do
      end do
      call vprod(d,2,3,1)
      call vprod(d,1,2,3)
      call vprod(e,2,3,1)
      call vprod(e,1,2,3)
c     Check for colinearity of the atoms
      dsum=0.0
      do k=1,3
        dsum=dsum+abs(d(k,1))
      end do
      if (dsum .gt. 1.e-6) go to 30
c     switch to two-atom algorithm
      if (iverb .gt. 0) write (iout,1000)
c     Diatomic
20    dsum=0.0
      do k=1,3
        d(k,1)=a(k,2)-a(k,1)
        dsum=dsum+abs(d(k,1))
        e(k,1)=b(k,2)-b(k,1)
        d(k,2)=0.0
        e(k,2)=0.0
      end do
C     Create the second vector as perpendicular to the bond
c     If the two atoms coincide,generate unit matrix
      if (dsum .gt. 1.e-6) go to 25
      if (iverb .gt. 0) write (iout,1002)
      go to 40
25    do k=1,3
        if (abs(d(k,1)) .le. 1.e-7) then
c         Zero component found
          d(k,2)=1.0
          go to 23
        end if
      end do
c     No zero component
      d(2,2)=1.0
      d(3,2)=-d(2,1)/d(3,1)
23    call vprod(d,1,2,3)
      do k=1,3
        if (abs(e(k,1)) .le. 1.e-7) then
c         Zero component found
          e(k,2)=1.0
          go to 26
        end if
      end do
c     No zero component
      e(2,2)=1.0
      e(3,2)=-e(2,1)/e(3,1)
26    call vprod(e,1,2,3)
30    call mnorm(d)
      call mnorm(e)
c     Now, d=orm*e, thus orm=d*inv(e) and
c     The inverse of an orthonormal matrix is its transpose
      do i=1,3
        do j=1,3
          sum=0.0
          do k=1,3
            sum=sum+e(i,k)*d(j,k)
          end do
          orm(i,j)=sum
        end do
      end do
      return
40    call unitmat(orm)
      return
1000  format(' NOTE: the first three solute atoms are on',
     -  ' the same line - two-atom algorithm will be used')
1002  format(1x,20('-'),' NOTE: the two atoms of a diatomic solute ',
     -  'coincide - unit matrix will be used as orientation matrix')
      end
