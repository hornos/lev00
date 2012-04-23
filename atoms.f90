module atoms

!     DIRC(#,xyz)          - current basis vectors in real space for the super-cell
!     RECC(#,xyz)          - current basis vectors in reciprocal space (with 2*pi)
!     BCELL(#,xyz)         - current basis vectors in reciprocal space (without 2*pi)
!     VOLC                 - current volume of the super-cell
!     TI(component,number) - atoms in the cell (Cartesian)
!     NspN(nsp)            - number of atoms in every species
!     Species(nsp)         - chemical symbol for every species
!     Mendel(nsp)          - position in the Mendeleev's table for every species
!     QradI(ion)           - atomic radius for the ion used to eliminate
!                            a singularity in the Madelung potential within the
!                            sphere of this radius near the ion 
!     atom_species(ion)    - atomic species (1,...,NSPEC) for every atom 
!     yesrlx               - if relaxation flags are available
!     old_to_new           - reference between atomic number before and after the species ordering (only SIESTA)

  real*8, parameter                      :: tiny_equiv=0.00001
  real*8, dimension(:,:), allocatable    :: TI
  character, dimension(:,:), allocatable :: relax_flag*1
  real*8, dimension(:), allocatable      :: QradI
  integer, dimension(:), allocatable     :: Z_atom
  integer, dimension(:), allocatable     :: NspN, Mendel,old_to_new
  integer, dimension(:), allocatable     :: atom_species
  character(Len=2), dimension(:), allocatable :: Species
  real*8, dimension(3,3)                 :: DIRC,RECC,BCELL
  real*8                                 :: VOLC
  logical                                :: yesrlx
end module atoms
