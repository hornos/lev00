module dos
!...........................................................................
!  Ntask         - the number of jobs (spheres or layers)
!  Ntetra        - number of tetrahedra generated
!  Broad_Band    - a number of dispersions to be taken on both sides
!  jspin         - the actual spin direction (1 = up, 2 = down) to be used
!                  in the DOS calculation
!  k_ref(tetr,k) - k-point numbers (k=1,..,4) for tetrahedra tetr
!
!.......... projected DOS 
!
!  Point(xyz,it) - atom position for the projected DOS job number=it
!  Radius(it)    - radius of the sphere for job number=it
!  method(it)    - projected DOS method for job number=it
!  v_atm(it)     - atom number
!  v_spec(it)    - species
!  which_task(it)  - .true. if chosen to be summed up
!...........................................................................
integer Ntask,Ntetra,jspin

real*8, dimension(:), allocatable          :: Radius
real*8, dimension(:,:), allocatable        :: Point
integer, dimension (:,:), allocatable      :: k_ref
character(Len=1), dimension(:),allocatable :: band_yes
real*8,dimension(:),allocatable            :: RWIGS 
integer, dimension (:), allocatable        :: v_atm,v_spec
logical, dimension(:),allocatable          :: which_task
real*8  :: Dispers=0.3
real*8  :: Broad_Band=1.5

!.....................................................................
! PSI2(NKP,NB,iTask) - integrated values of |psi|**2 for every k-point
!           NKP and every band NB (it is used for the projected DOS)
!.....................................................................
real*8,dimension(:,:,:), allocatable       ::  PSI2
character(len=3),dimension(:),allocatable  ::  Flag
character(len=10),dimension(:),allocatable ::  method
character CaseDos*6
end module dos
