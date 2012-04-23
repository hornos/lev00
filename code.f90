module code
!........................................................................
!........... Common statements ..........................................
!........................................................................
  implicit none
! seed     - job name (e.g. SIESTA)
! seedfdf  - name of SIESTA fdf file (e.g. INPUT_FILE.fdf)
! nameout  - name of the DFT output file
! nameg    - name of the DFT geometry file
! rho_form = .true. if siesta RHO file is formatted; otherwise = .false.
  character    ::  seed*100,seedfdf*100,Which_Code*6
  character    ::  nameout*100,nameg*100
  logical rho_form

CONTAINS

  subroutine which_PW_code()
!
!................ ask about the code
!    
    implicit none
    character item*10
10  write(*,*)'..................................................'
    write(*,*)'Choose the input format:'
    write(*,*) 
    write(*,*)'=========== ALL features =========================' 
    write(*,*) 
    write(*,*)'       V.   VASP'
    write(*,*)'       S. SIESTA'
    write(*,*) 
    write(*,*)'====== ONLY density MANIPULATION =================' 
    write(*,*) 
    write(*,*)'      VD.   VASP'
    write(*,*)'      SD. SIESTA'
    write(*,*) 
    write(*,*)'       Q. Quit'
    write(*,*)
    write(*,*)'------------>'
    read(*,*) item
  
    if(trim(item)=='V') then
       Which_Code='  VASP'
       return

    else if(trim(item)=='S') then
       Which_Code='SIESTA'
       return

    else if(trim(item)=='SD') then
       Which_Code='SIESTA'
       call manipulate()

    else if(trim(item)=='VD') then
       Which_Code='  VASP'
       call manipulate()

    else if(trim(item)=='Q') then
       stop
    else
       write(*,*) "Incorrect! Please, try again!"
    end if
    go to 10
  end subroutine which_PW_code

end module code
