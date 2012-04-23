subroutine mainmenu()
use atoms
use code
use param
use kpoints
use dos
use siesta_EIG
implicit none
character problem*3,item*2
logical Err

!................ ask about the problem to solve
!
10 write(*,*)'.....................................................'
   write(*,*)'Choose the problem to solve:'
   write(*,*)
   write(*,*)'============== DOS options ======================='
   if(yeskp .and. yesen) then
      write(*,*)'=====> 1. Method of Tetrahedra (Blochl) <=========='
      write(*,*)' { You must have [band.out, brill.dat] from TETR }'   
      write(*,*)' Cd. Choose dispersion for the DOS smearing'
      write(*,*)' TD. Total DOS'   
      if(Which_Code.eq.'  VASP') write(*,*)' PD. Projected DOS calculation'
   end if
   if(Which_Code.eq.'SIESTA') then
      write(*,*)'============> 2. Sum of Gaussians <================'
      write(*,*)' DO. Total DOS: you need [seed].EIG'   
      write(*,*)' PD. Projected DOS: you need [seed].PDOS file'   
   end if
   write(*,*)
   write(*,*)'============== DENSITY options ==================='
   if(Which_Code.eq.'  VASP') then
      write(*,*) ' D. total/partial electron density, electrostatic potential'
   else
      write(*,*) ' D. total or partial electron density; electrostatic potential'
   end if
   if((Which_Code.eq.'  VASP' .and. ispin.eq.2) .or. Which_Code.eq.'SIESTA') &
      write(*,*) 'SD. total or partial spin density; electrostatic potential'
   
   write(*,*)
   write(*,*)'============== ELECTR. potential options ========='
   write(*,*)'CP. Coulomb potential (point,line,plane)'
   if(Which_Code.eq.'SIESTA') then
      write(*,*)'============== OTHER ============================='
      if(yesrlx) write(*,*)'Ck  Check if SIESTA relaxation finished properly'
   end if
   write(*,*)
   write(*,*)' Q. Quit'
   write(*,*)
   write(*,*)'------------>'
   read(*,'(a)',err=11) item

![DSm]_____ total DOS for various smearing parameters
   if(trim(item).eq.'Cd' .and. yeskp .and. yesen) then
      problem='DSm'
      call ask_spin(jspin)
      call prep_disp()

![PRO]_____ projected DOS 
   else if(trim(item).eq.'PD') then
      if(yeskp .and. yesen .and. Which_Code.eq.'  VASP') then
         problem='PRO'
         call ask_spin(jspin)
         call get_psi2_from_vasp(Err)
         call prep_dos(1)
      else if(Which_Code.eq.'SIESTA') then
         call ask_spin(jspin)
         call read_siesta_pdos(jspin)
      end if
      
![tDS]_____ total DOS: run an interactive 
!  program to sum contributions from every tetrahedra (those were built 
!  earlier by the 1st part of the code.

   else if(trim(item).eq.'TD' .and. yeskp .and. yesen) then
      problem='tDS'
      call ask_spin(jspin)
      call prep_dos(0)

   else if(Which_Code.eq.'SIESTA' .and. trim(item).eq.'DO') then
      problem='tDS'
      call ask_spin(jspin)
      call dos_from_EIG(jspin)
      
![MAP]_____ Prepare output files for the densities in the format 
!  understandable by plotters (i.e. when the coordinates are in 
!  Angstroms (instead of being through the grid numbers) and the 
!  density is along a line or in a plane: run an interactive program

![D]_________ plot the charge density in the case of DNS and quit
   else if(trim(item).eq.'D') then
      problem='DNS'
      call density(.false.)

![SD]_________ plot the spin density in the case of DNS and quit
   else if(trim(item).eq.'SD' .and. ispin.eq.2 ) then
      problem='DNS'
      call density(.true.)
      
![CP]__________ calculate Coulomb potential for point charges 
!
   else if(trim(item).eq.'CP') then
      problem='Mad'
      call lev_coulmb()

![Ck]__________ check siesta forces and if the relaxation finished
!
   else if(trim(item).eq.'Ck' .and. Which_Code.eq.'SIESTA' .and. yesrlx) then
      problem='Chk'
      call check_siesta_forces()
      
   else if(trim(item).eq.'Q') then
      return
   end if
   go to 10
11 write(*,*) "Incorrect item number! Try again!"
   go to 10
end subroutine mainmenu
 
subroutine ask_spin(jspin)
!.................................................................
!  in the case of spin-polarised calculation (ispin=2)
!                ask for which spin: 1 or 2
!.................................................................
! jspin = 1  - non spin-polarised calculation
!         1  - spin polarised, spin UP
!         2  - spin polarised, spin DOWN
!.................................................................
use param
implicit none
integer jspin
1 if(ispin.eq.1) then
     jspin=1
  else if(ispin.eq.2) then
     write(*,*)'Specify which spin: up (1) or down (2)?'
     read(*,*,err=1) jspin
     if(jspin.eq.1) then
        write(*,*) 'DOS will be calculated for spin UP'
     else if(jspin.eq.2) then
        write(*,*) 'DOS will be calculated for spin DOWN'
     else
        go to 1
     end if
  end if
end subroutine ask_spin

subroutine check_siesta_forces()
use param
use atoms
use code
implicit none
character line*200
integer LinEnd(100),LinPos(100),NumLin,iErr,nat,i,at,i1
real*8 f(3),ff,ftol

   write(*,'(a)')'... Reading forces from the main SIESTA output '//trim(nameout)//'...'
   open(19,file=trim(nameout),status='old',form='formatted',err=10)

!___________ check if relaxation was successful   

   call find_string('outcoor: Relaxed atomic coordinates',35,line,19,.true.,iErr)
   if(iErr.eq.0) then
      write(*,*)'GOOD: this calculation HAS FINISHED!'
   else
      write(*,*)'WARNING: this calculation has NOT finished!'
   end if

!..................... check individual forces

20 write(*,*)'Specify force tolerance for checking:'
   read(*,*,err=20) ftol ; if(ftol.lt.0.0) ftol=abs(ftol)
   
!___________ find the last forces
   rewind(19) ; i=0
1  call find_string('siesta: Atomic forces (eV/Ang):',31,line,19,.true.,iErr)
   if(iErr.eq.0) then
      i=i+1 ; go to 1
   end if
   write(*,*)'... Found the forces ',i,' times ...'
   rewind(19) ; i1=0
2  call find_string('siesta: Atomic forces (eV/Ang):',31,line,19,.true.,iErr)
   i1=i1+1 ; if(i1.lt.i) go to 2

!___________ read and check the last forces

   do at=1,NIONS
      call CutStr(Line,NumLin,LinPos,LinEnd,19,0,iErr)
      if(NumLin.eq.5) then
         read(line(LinPos(2):LinEnd(5)),*,err=40) nat,(f(i),i=1,3)
      else if(NumLin.eq.4) then
         read(line(LinPos(1):LinEnd(4)),*,err=40) nat,(f(i),i=1,3)
      end if
      write(*,*) nat,(f(i),i=1,3)
      if(at.ne.nat) go to 40
      ff=sqrt(f(1)*f(1)+f(2)*f(2)+f(3)*f(3))
!      write(*,'(i5,5x,3a)') at,(relax_flag(j,nat),j=1,3)
      if( (relax_flag(1,nat).eq.'T' .and. abs(f(1)).gt.ftol) .or. &
            (relax_flag(2,nat).eq.'T' .and. abs(f(2)).gt.ftol) .or. &
              (relax_flag(3,nat).eq.'T' .and. abs(f(3)).gt.ftol) ) &
          write(*,'(a,3(x,f10.5),a,i5,a,f10.5)') &
              'WARNING: large force ( ',f,' ) on atom =',nat,' abs= ',ff
   end do
   write(*,*)'Finished checking individual atomic forces ...'
   write(*,*)'Found ',nat,' atoms in the system'
   call CutStr(Line,NumLin,LinPos,LinEnd,19,0,iErr)
   if(NumLin.eq.5) then
      read(line(LinPos(3):LinEnd(5)),*,err=40) (f(i),i=1,3)
   else if(NumLin.eq.4) then
      read(line(LinPos(2):LinEnd(4)),*,err=40) (f(i),i=1,3)
   end if
   ff=sqrt(f(1)*f(1)+f(2)*f(2)+f(3)*f(3))
   write(*,'(a,3(f10.5,x),a,f10.5)') '... Checking the net force => ',f,' abs= ',ff
   close(19)  ; return
10 write(*,*)'FATAL! Cannot open SIESTA on-screen output file'
   close(19) ; return
40 write(*,*)'ERROR reading SIESTA forces from on-screen output file'
   close(19)
end subroutine check_siesta_forces
