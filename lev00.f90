program lev00
!.........................................................................
!..................  This is a post-processing routine    ................
!..................      for DFT codes, such as           ................
!..................           VASP, SIESTA               ................
!.........................................................................
!..................  What it can do:                      ................
!..................  ^^^^^^^^^^^^^^^                      ................
!..................  - total and projected DOS;           ................
!..................  - total or partial charge density    ................
!..................    along lines or in planes           ................
!..................  - integration of charge density      ................
!..................  - manipulation with charge densities ................
!..................  - and a lot, lot more                ................
!.........................................................................
!.................    Started around:         19.08.94    ................
!.................    F90 version started:    30.07.07    ................
!.........................................................................
!.................    Lev.Kantorovitch@kcl.ac.uk          ................
!.........................................................................
!
! the code (SIESTA, VASP, ...) is chosen by user and the dimensions
! for arrays are determined from the corresponding DFT output files 
use param
use atoms
use code
use kpoints
use dos
implicit none
logical Err 
!
call hat()
!
!............... ask about the code
!
Err=.true. ; yeskp=.false.  ; yesen=.false. ; yesrlx=.false.
do while(Err) 
   call which_PW_code()
   write(*,*) 'The code is '//Which_Code
!
!................read parameters from output files of the DFT codes
!                and do eigenvalues
!
   if(Which_Code=='  VASP') then
      call do_param(Err)
   else if(Which_Code=='SIESTA') then
      call get_param_siesta(Err)
   end if
end do
NPLWV=NGX*NGY*NGZ

write(*,*)'=========> Parameters found <=========='
write(*,*) 'Grid               = ',NGX,'x',NGY,'x',NGZ
write(*,*) 'Number of atoms    = ',NIONS
write(*,*) 'Number of species  = ',NSPEC
write(*,*) 'Number of K-points = ',NKPTS
write(*,*) 'Number of bands    = ',NBANDS
if(ispin == 1) then
   write(*,*)'This is NOT spin-polarised calculation'
else if(ispin == 2) then
   write(*,*)'This is spin-polarised calculation'
else
   write(*,*)'WARNING: spin NOT recognised; set to no-spin'
   ispin=0
end if
NPLWV=NGX*NGY*NGZ
!
!............. allocate necessary intital memory
!
allocate(TI(3,NIONS))
allocate(QradI(NIONS))
allocate(atom_species(NIONS))
allocate(NspN(NSPEC))
allocate(relax_flag(3,NIONS))
if(Which_Code=='  VASP') allocate(RWIGS(NSPEC))
allocate(VKPT(3,NKPTS))
allocate(WTKPT(NKPTS))
!
!.............. read geometry
!
Err=.true.
do while(Err) 
   if(Which_Code=='  VASP') then
      call read_vasp_geom(Err)
      call get_rwigs(RWIGS)
   else if(Which_Code=='SIESTA') then
      call read_siesta_geom(Err)
   end if
end do
!
!.............. main menu 
!
call mainmenu()
stop "I owe you nothing!"
end program lev00


