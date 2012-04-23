module kpoints
!.... max. number of points in the DOS plot
!       parameter (Npoints=10000)
!       parameter (N0_Gauss=30)
!
!.........................................................................
! VKPT(3,NKPTS) = the coordinates of the special k points used for
!                     Brillouin zone averages
! WTKPT(NKPTS) = the weight to be attached to each special k point
!.........................................................................
  real*8, dimension(:,:),allocatable :: VKPT
  real*8, dimension(:),allocatable   :: WTKPT
!  yeskp - if k-points are available
!  yesen - if KS eigenvalues are available
  logical yeskp,yesen
end module kpoints


