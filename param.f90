module param
!
! grid               = NGX x NGY x NGZ
!       NPLWV        = NGX*NGY*NGZ
! number of species  = NSPEC
! number of ions     = NIONS
! number of k points = NKPTS
! number of bands    = NBANDS
! if spin-polarised  = ispin
!
  integer :: NGX,NGY,NGZ,NIONS,ispin,NSPEC,NKPTS,NBANDS,NPLWV
end module param
