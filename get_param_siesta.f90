subroutine get_param_siesta(Err)
use param
use atoms
use code
use kpoints
implicit none
character Line*200
integer LinEnd(100),LinPos(100)
integer iErr,nk,nb,isp,j,j1,i,nspin,NumLin,ii,ispin1,NKPTS1,nk1
real*8 DIR(3,3),a,b,c
real*8, parameter                   :: au_A=1.889725989,tiny=0.001
real*8, dimension(:,:,:), allocatable :: E
logical Err,alloc,know_grid,bandsf
!
  Err=.false.
  alloc=.false.
!.................. which siesta job name
  write(*,'(/a/)') '... Reading from SIESTA input files (fdf-files)...'
  write(*,'(/a)')'***** fdf-files in the current directory: *****'
  call system('ls *fdf* ')
  write(*,'(a)') '***********************************************'
  write(*,'(a)') '********   Give the full file name!   *********'
  write(*,'(a/)')'***********************************************'
  read(*,'(a)') seedfdf
  write(*,'(/a/)') '... Reading from '//trim(seedfdf)//' ...'
  open(14,file=trim(seedfdf),form='formatted',status='old',err=29)
 
!_______ system label

  call find_string('SystemLabel',11,line,14,.false.,iErr)
  if(iErr.eq.1) then
     write(*,*)'WARNING! System label NOT found => specify manually:'
     write(*,*)'Specify SYSTEM LABEL (char*20) for your SIESTA run:'
     write(*,'(/a)')'*** e.g. XV-files in the current directory: **'
     call system('ls *.XV ')
     write(*,'(a)') '**********************************************'
     write(*,'(a)') '*   Extention .XV should NOT be specified!   *'
     write(*,'(a/)')'**********************************************'
     read(*,'(a)') seed
  else
     do i=1,100
        seed(i:i)=' '
     end do
     call CutStr(Line,NumLin,LinPos,LinEnd,0,0,iErr)
     read(line(LinPos(2):LinEnd(2)),'(a)') seed
  end if
  write(*,'(a)') '... SIESTA job seed found is '//trim(seed)
!
!========================================================================
!  (a) read number of atoms NIONS
!  (b) number of species NSPEC
!  (c) grid NGX, NGY, NGZ
!  (d) k points number NKPTS
!========================================================================
!
![1-2]............... using fdf file: NIONS, NSPEC, k points and spin
!
   ii=0

!________ number of species
   rewind(14)
   write(*,'(/a)')'[1] checking number of species ...'
   call find_string('NumberOfSpecies',15,line,14,.false.,iErr)
   if(iErr.eq.1) then
      write(*,'(a)')'    ERROR: string "NumberOfSpecies" not found! Exiting ...'
      go to 29
   end if
   call CutStr(Line,NumLin,LinPos,LinEnd,0,0,iErr)
   read(Line(LinPos(2):LinEnd(2)),*) NSPEC
   write(*,'(a,i1)') ' ===> SIESTA: found NSPEC =', NSPEC

!_________ number of atoms 
   rewind(14)
   write(*,'(/a)')'[2] checking number of ions ...'
   call find_string('NumberOfAtoms',13,line,14,.false.,iErr)
   if(iErr.eq.1) then
      write(*,'(a)')'    ERROR: string "NumberOfAtoms" not found! Exiting ...'
      go to 29
   end if
   call CutStr(Line,NumLin,LinPos,LinEnd,0,0,iErr)
   read(Line(LinPos(2):LinEnd(2)),*) NIONS
   write(*,'(a,i5)')' ===>  SIESTA: found NIONS = ',NIONS

!_________ spin 
   rewind(14)
   ISPIN=1
   write(*,'(/a)')'[3] checking spin polarisation ...'
   call find_string('SpinPolarized',13,line,14,.false.,iErr)
   if(iErr.eq.1) then
      write(*,'(a)')'    WARNING: string "SpinPolarized" not found!'
      write(*,'(a)')'             Assuming there is none ...'
      go to 21
   end if
   if(index(line,'true').ne.0 .or. index(line,'TRUE').ne.0 .or. &
                                         index(line,'True').ne.0 ) ISPIN=2
21 write(*,'(a,i5)') ' ===>  SIESTA: found ISPIN = ',ISPIN

!_________ k points
   rewind (14)
   write(*,'(/a)')'[4] checking k-points ...'
   NKPTS=1
   call find_string('BandPoints',10,line,14,.false.,iErr)
   if(iErr.eq.1) then
      write(*,'(a)')'    WARNING: string "BandPoints" not found!'
      write(*,'(a)')' ===>  METHOD of TETRAHEDRA for DOS is NOT available!  <==='
      go to 23
   end if
22 read(14,'(a)',end=29) line
   if(index(line,'endblock').eq.0) then
      ii=ii+1
      go to 22
   end if
   close (14)
   NKPTS=ii
   if(NKPTS.eq.0) NKPTS=1
23 write(*,'(a,i5)')' ===>  SIESTA: setting NKPTS = ',NKPTS
   write(*,*) 'O.K.! This '//trim(seedfdf)//' is fine! Let us proceed!' 
!
![3]............... using bands file (for NBANDS and energies)
!
   NBANDS=1
   write(*,'(a)')'[5] checking electronic bands ...'
   write(*,'(/a/)') '... Reading from '//trim(seed)//'.bands  ...'
   open(14,file=trim(seed)//'.bands',form='formatted',status='old',err=37)
   read(14,*,err=37) a,b,c, NBANDS,ispin1,NKPTS1
   if(ispin1.ne.ISPIN) then
      write(*,'(a)')'FATAL! '//trim(seed)//'.bands and '//&
                               trim(seed)//'.fdf are inconsistent in spin'
      go to 37
   end if
   if(NKPTS1.ne.NKPTS) then
      write(*,'(a)')'FATAL! '//trim(seed)//'.bands and '//&
                               trim(seed)//'.fdf are inconsistent in NKPTS'
      go to 37
   end if
!_________________ reading energies
   write(*,'(a,i5)')' ===>  SIESTA: setting NBANDS = ',NBANDS
   allocate(E(NKPTS,NBANDS,ISPIN)) ; alloc=.true.
   do nk=1,NKPTS
      read(14,*,err=37) a,b,c,((E(nk,nb,isp),nb=1,NBANDS),isp=1,ISPIN)
   end do
   close (14)
   write(*,'(4x,a)')'KS eigenvalues read in successfully!'

!____________________write energies into band.out files

   do isp=1,ISPIN
      if(isp.eq.1) then
         if(ISPIN.eq.1) then
            open(3,file='band.out',form='formatted')
            write(*,'(a)')' ===> Producing file <band.out> ...'
         else
            open(3,file='band.out.1',form='formatted')
            write(*,'(a)')'===> Producing file <band.out.1> ...'
         end if
      else
         open(3,file='band.out.2',form='formatted')
         write(*,'(a)')' ===> Producing file <band.out.2> ...'
      end if
      do nk=1,NKPTS
         do j=1,NBANDS,5
            j1=j+4
            if(j1.gt.NBANDS) j1=NBANDS
            write(3,'(i3,3x,5(1x,f10.5))') nk,(E(nk,nb,isp),nb=j,j1)
         end do
      end do
      close (3)
   end do
   write(*,'(a)')' ===> Done band.out[1,2] ! <=== '
   yesen=.true.
   go to 38
!
37 write(*,'(4x,a)')'WARNING! Bad or absent SIESTA '//trim(seed)//'.bands file!'
   write(*,'(a)')' ===>  DOS option will not be available!  <==='
38 if(alloc) deallocate(E)
!
!___________ getting the grid for the electron/spin density
!
   know_grid=.false.
!__________________ RHO: try formatted first

   write(*,'(/a)')'[6] checking electronic grid ...'
   write(*,'(4x,a)')'... Trying FORMATTED density ...'
   open(19,file=trim(seed)//'.RHO',status='old',form='formatted',err=10)
   write(*,'(4x,a)')'... Opened as FORMATTED'
   read(19,*,err=10,end=10) ((DIR(i,j),j=1,3),i=1,3)
   write(*,'(4x,a)')'... reading ngx ...'
   read(19,*,err=10,end=10) NGX,NGY,NGZ,nspin
   rho_form=.true.
   write(*,'(4x,a)') 'YES! Found FORMATTED density file '//trim(seed)//'.RHO'
   go to 20

!__________________ RHO: try unformatted then

10 close (19)
   write(*,'(4x,a)')'... Trying UNFORMATTED density ..'
   open(19,file=trim(seed)//'.RHO',status='old',form='unformatted',err=50)
   write(*,'(4x,a)')'... Opened as UNFORMATTED'
   write(*,'(4x,a)')'... reading lattice vectors ...'
   read(19,err=50,end=50) ((DIR(i,j),j=1,3),i=1,3)
   write(*,'(4x,a)')'... reading ngx, ngy, ngz ...'
   read(19,err=50,end=50) NGX,NGY,NGZ,nspin
   write(*,'(4x,a)')'YES! Found UNFORMATTED density file '//trim(seed)//'.RHO'
   rho_form=.false.

!__________________________________ proceed from here

20 write(*,'(3(a,i5))') ' ===> SIESTA: found NGX=',NGX,',NGY=',NGY,',NGZ=',NGZ
   if(nspin.eq.1) then
      write(*,'(4x,a)')' ===> SIESTA: file '//trim(seed)// &
           '.RHO does not contain SPIN density'
   else if(nspin.eq.2) then
      write(*,'(4x,a)')' ===> SIESTA: file '//trim(seed)// &
           '.RHO does contain SPIN density'
   else
      go to 50
   end if
   close (19)
   know_grid=.true.
   go to 51

!__________________ OUT: if attempt with RHO was unsuccessful, try the main output

50 write(*,'(4x,a)')'WARNING! File '//trim(seed)//'.RHO not found/bad!'
   write(*,'(/4x,a)')'... Trying main SIESTA output then ...'
51 write(*,'(4x,a)')'Name of the on-screen SIESTA output file (char*100):'
   write(*,'(a)') '**********************************************'
   write(*,'(a)') '********   Give the full file name!   *********'
   write(*,'(a/)')'**********************************************'
   read(*,'(a)') nameout
   open(19,file=trim(nameout),status='old',form='formatted',err=100)

   if(.not.know_grid) then
!_____________________________ read out spin

      call find_string('Number of spin components',25,line,19,.false.,iErr)
      if(iErr.eq.1) go to 100
      call CutStr(Line,NumLin,LinPos,LinEnd,0,0,iErr)
      if(iErr.ne.0) go to 100
      read(Line(LinPos(7):LinEnd(7)),*) nspin
      write(*,'(a,i1)') ' ===> SIESTA: found nspin =', nspin

!_____________________________ read out grid

      call find_string('InitMesh: MESH',14,line,19,.false.,iErr)
      if(iErr.eq.1) go to 100
      call CutStr(Line,NumLin,LinPos,LinEnd,0,0,iErr)
      if(iErr.ne.0) go to 100
      read(Line(LinPos(4):LinEnd(4)),*) NGX
      read(Line(LinPos(6):LinEnd(6)),*) NGY
      read(Line(LinPos(8):LinEnd(8)),*) NGZ
      write(*,'(3(a,i5))') ' ===> SIESTA: found NGX=',NGX,',NGY=',NGY,',NGZ=',NGZ
   end if
   close (19)
   return
!
!......... errors
19 write(*,*)'FATAL! Bad or absent SIESTA '//trim(seed)//'.XV file!'
   go to 200
29 write(*,*)'FATAL! Bad or absent SIESTA '//trim(seedfdf)//' file!'
   go to 200
100 write(*,*)'FATAL! File '//trim(nameout)//' not found or bad!'
200 Err=.true.
end subroutine get_param_siesta
