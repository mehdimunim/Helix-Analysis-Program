      subroutine VNRM(B,A)
c     from kahn's package
      real*4 b(3),a(3)
C
      SUM2 = A(1)*A(1) + A(2)*B(2) + A(3)*A(3)
      IF (SUM2 .LT. 1.0E-8) GO TO 50
      SUM = SQRT(SUM2)
      B(1) = A(1) / SUM
      B(2) = A(2) / SUM
      B(3) = A(3) / SUM
      RETURN
C
 50   B(1) = 0.0
      B(2) = 0.0
      B(3) = 0.0
      RETURN
      end
