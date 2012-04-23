logical function ask_equiv(r,XJ,iCase,tiny_equiv)
!........................................................................
!    It checks whether vector r() is equivalent to lattice translations
!  reciprocal to vectors XJ (up to 2*pi):
!
!  if XJ - reciprocal lattice vectors (without 2*pi), then r() must be 
!          in the direct space;
!  if XJ - direct lattice vectors, then r() must be in the reciprocal
!          space;
!........................................................................
!  iCase = 3 - if 3D space;
!  iCase = 2 - if 2D space
!........................................................................
!  .true. - if r() is equivalent to latice translations
!  .false.- if not.
!........................................................................
  implicit none
  character answer
  real*8 r(3),XJ(3,3),x,tiny_equiv
  integer iCase,i
  ask_equiv=.false.
!
!____ multiply on all direct lattice vectors: the result must be integer
! if r() is an arbitrary sum of reciprocal lattice vectors (note, r() is
! defined up to 2*pi because of the choice of XJ used to generate it)
!
  do i=1,iCase
     x=r(1)*XJ(i,1)+r(2)*XJ(i,2)+r(3)*XJ(i,3)
     if( abs(x-nint(x)) .gt. tiny_equiv ) return
  end do
  ask_equiv=.true.
end function ask_equiv

double precision function pos_num(x)
!..........................................................................
! returns x = 0 if x<0 or if x is very close to zero.
! This is used to print positive numbers to be plotted by gnuplot as the
! latter gets confused with very small positive numbers and shows a large 
! noise
!..........................................................................
real*8, parameter :: err=1.0e-30
real*8 x
if(x.gt.err) then
   pos_num=x
else
   pos_num=0.0d0
end if
end function pos_num

