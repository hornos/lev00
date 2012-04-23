module explore
!.............................................................
!    This parameters are needed for the exploration program
! which explores the density inside your cell.
!
!        Mpeak - maximal possible number of peaks;
! (2*nresol+1) - number of points along 1-dimensional plots of the
!                density (along three lattice vectors)
!
!.............................................................
integer, parameter::  Mpeak=20,nresol=50
real*8 :: x(3,Mpeak),Rad(Mpeak),High(Mpeak),rCharg(Mpeak)
end module explore
