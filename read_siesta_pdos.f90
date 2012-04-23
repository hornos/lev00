subroutine read_siesta_pdos(jspin)
!.................................................................
!  reading SIESTA projected DOS information from seed.PDOS file
!.................................................................
!  coordinates are given in a.u. in this file; not actually used
!.................................................................
! Up to 20 PDOS curves can be previewed (togteher with the total 
! DOS) simultaneously.
!.................................................................
use param
use code
use atoms
use siesta_EIG
implicit none
real*8,dimension(:),allocatable    :: energy,totdos
real*8,dimension(:),allocatable    :: pdos
integer,dimension(:,:),allocatable :: index
logical,dimension(:,:),allocatable :: listat
character(len=2),dimension(:),allocatable :: Labels
integer spin,iErr,norbitals,i,ii,ien,iorb,nat,ngroup,lenght(20)
character line*200,m*2,item*2,list(20)*200
integer at,n,l,iz,len(20),jspin
real*8 :: x(3),en1,en2,Up_lim=0.0,Up_lim_E=0.0
logical TagYes(20),DoneDos,TagY,ErrF
character Title*50,title_pl(20)*8,filen(20)*12,cha1,cha2*2,ShowTotDos
data Title/'                                                  '/

!........... species by atom number
   allocate(Labels(NIONS))
   nat=0
   do i=1,NSPEC
      do ii=1,NspN(i)
         nat=nat+1 ; Labels(nat)=Species(i)
      end do
   end do

!............ reading HOMO/LUMO info and Fermi energy (just for printing)

   call read_siesta_fermi(ErrF,jspin)

!........... reading and anlysing siesta PDOS file: making up
!    - max numbers of n,l for each atom and overall; 
!    - number of energies
!    - number of spins

   write(*,*)'Reading from '//trim(seed)//'.PDOS ...'
   open(14,file=trim(seed)//'.PDOS',form='formatted',status='old',err=19) 

   read(14,'(/7x,i1)') spin
   if(spin.ne.ispin) go to 19
   read(14,'(11x,i4/)') norbitals
   write(*,*)'... found ',norbitals,'  orbitals ...'
      
!________ reading ENERGIES: ien - number of energies and hence
!                                 of the pdos values for each orbital
   ien=0
1  read(14,'(a)') line
   if(line(1:16).ne.'</energy_values>') then
      ien=ien+1 ; go to 1
   else
      go to 2
   end if
2  write(*,*)'... found ',ien,'  energies ...'
   rewind(14)
   call find_string('<energy_values',14,line,14,.false.,iErr)
   allocate(energy(ien)) ; allocate(pdos(ien))
   read(14,*) (energy(i),i=1,ien)
   en1=energy(1) ; en2=energy(ien)

!___________ start analysing ORBITALS STRUCTURE
!            iorb - orbital number
   allocate(index(NIONS,2))
   nat=0
5  call read_orb_signature(iorb,at,x,n,l,m,iz,14,iErr)
   if(iErr.eq.1) go to 6
   if(nat.ne.at) then
      index(at,1)=iorb
      if(nat.ne.0) index(at-1,2)=iorb-1
      nat=at
   end if
   go to 5
6  index(nat,2)=iorb
   write(*,*)'... scanned ',iorb,'  orbitals ...'
   write(*,'(/a)') '.... found DISTRIBUTION of ORBITALS in atoms: ....'
   do i=1,NIONS
      write(*,'(3(a,i10))') 'atom= ',i,' => between ',index(i,1),' and ',index(i,2) 
   end do
   rewind(14)

!____________ do the TOTAL DOS _____________________________________

   allocate(listat(NIONS,20))
   listat=.true.
   allocate(totdos(ien))
   call do_pdos(listat,index,totdos,ien,14,jspin,iErr)
   if(iErr.eq.1) then
      write(*,*)'ERROR in calculating the total DOS.'
      write(*,*)'Returning ...'
      return
   end if

!____________ start the MAIN MENU __________________________________   

   TagYes=.false.
   listat=.false.
   DoneDos=.false.
   ShowTotDos='Y'
   ngroup=0
21 write(*,*)'.....................................................'
   if(.not.ErrF) write(*,'(a,f10.5,a)')'   Fermi energy = ',Fermi,' eV'
   write(*,*)'Choose the necessary options:'
   write(*,*)
   write(*,'(a,i5)')' N. Number of groups of atoms: ',ngroup
   Tagy=.false.
   if(ngroup.ne.0) then
      do i=1,ngroup
         if(i.lt.10) then
            write(cha1,'(i1)') i
            cha2=' '//cha1
         else
            write(cha2,'(i2)') i
         end if
         if(TagYes(i)) then
            write(*,'(2a)') cha2,'. Choose atoms for the PDOS in this group'
            write(*,'(a)')'     >> Current set of chosen atoms: '//list(i)(1:len(i))
         else
            write(*,'(2a)') cha2,'. Choose atoms for the PDOS in this group <== underfined!'
         end if
         if(i.eq.1) then
            TagY=TagYes(1)
         else
            TagY=TagY.and.TagYes(i)
         end if
      End do
      if(TagY) then
         if(DoneDos) then
            write(*,'(a)')'CD. Calculate the PDOS for the chosen atoms <= Done!'
            if(ngroup.eq.1) then
               write(*,'(a)')'    Data file: pdos_1.dat'
            else if(ngroup.eq.2) then
               write(*,'(a)')'    Data files: pdos_1.dat, pdos_2.dat'
            else
               if(ngroup.lt.10) then
                  write(cha1,'(i1)') ngroup
                  write(*,'(a)')'    Data files: pdos_1.dat,...,pdos_'//cha1//'.dat'
               else
                  write(cha2,'(i2)') ngroup
                  write(*,'(a)')'    Data files: pdos_1.dat,...,pdos_'//cha2//'.dat'
               end if
            end if
            write(*,'(a)')' P. Plot the PDOS'
         else
            write(*,'(a)')'CD. Calculate the PDOS for the chosen atoms'
         end if
      end if
   end if

   write(*,'(/a)')'=========== Settings for the plot ================'
   write(*,'(a,2(f10.5,a))')' E. Calculation: change energy interval [',en1,',',en2,']'
   if(Up_lim.eq.0.0) then
      write(*,'(a)')' Y. Upper limit to the dos-range of the plot: DISABLED'
   else
      write(*,'(a46,f10.5)') ' Y. Upper limit to the dos-range of the plot: ',Up_lim
   end if
   if(Up_lim_E.eq.0.0) then
      write(*,'(a)')' X. Upper limit to the energy-range of the plot: DISABLED'
   else
      write(*,'(a46,f10.5)') &
           ' X. Upper limit to the energy-range of the plot: ',Up_lim_E
   end if
   if(ShowTotDos.eq.'Y') then
      write(*,'(a)')' TD. Total DOS to be shown as well <== YES'
   else
      write(*,'(a)')' TD. Total DOS to be shown as well <== NO'
   end if
   write(*,'(/a)')'================= General help ==================='
   write(*,'(a)')' Co. Show current atomic positions in fractional/Cartesian'
   write(*,'(a)')' Q. Quit'
   write(*,*)
   write(*,*)'------------>'
   read(*,'(a)') item

![T]............ choose atoms for the PDOS for each group

   if(ngroup.ne.0) then
      read(item,'(i2)',err=11) i
      if(i.ge.1 .and. i.le.20) then
         call TTag(TagYes(i),NIONS,TI,Labels,listat(1,i),list(i),len(i))
         DoneDos=.false.
      end if
      go to 21
   end if

![N]............ number of groups of atoms

11 IF(trim(item).eq.'N') then

13    write(*,*)'Specify number of PDOS plots to display simultaneously'
      write(*,*)'=========> this should be no more thatn 20 <=========='
      read(*,*,err=13) ngroup
      if(ngroup.lt.0 .or. ngroup.gt.20) go to 13

![CD].............. calculate PDOS

   ELSE IF(trim(item).eq.'CD' .and. TagY) then

      do i=1,ngroup
         call do_pdos(listat(1,i),index,pdos,ien,14,jspin,iErr)
         if(iErr.eq.1) then
            write(*,*)'ERROR in calculating the PDOS'
            DoneDos=.false.
            go to 21
         end if
         if(i.lt.10) then
            lenght(i)=10
            write(cha1,'(i1)') i
            filen(i)='pdos_'//cha1//'.dat'
            title_pl(i)='group '//trim(cha1)
         else
            lenght(i)=11
            write(cha2,'(i2)') i
            filen(i)='pdos_'//cha2//'.dat'
            title_pl(i)='group '//trim(cha2)
         end if
         open(24,file=trim(filen(i)))
         write(24,'(a)') '#energy   total-DOS  PDOS'
         do ii=1,ien
            if(energy(ii).ge.en1.and.energy(ii).le.en2) then
               write(24,'(3(f10.5,x))') energy(ii),totdos(ii),pdos(ii)
            end if
         end do
         close (24)
      end do
      DoneDos=.true.

![P].............. plot PDOS

   ELSE IF(trim(item).eq.'P' .and. DoneDos) then

      call Plot2(ngroup,filen,lenght,Title,title_pl, &
                    'Energy (eV)         ', &
                    'DOS/PDOS (arb.units)', 'Screen', 33, &
                     Up_lim,Up_lim_E,ShowTotDos)

![Y]............ assign some upper limit on the y-axis for the plot

     ELSE IF(trim(item).eq.'Y') THEN
78      write(*,*)'Specify the upper limit for y-axis on your plot:'
        read(*,*,err=78) Up_lim

![X]............ assign some upper limit on the x-axis for the plot

     ELSE IF(trim(item).eq.'X') THEN
79      write(*,*)'Specify the upper limit for x-axis on your plot:'
        read(*,*,err=79) Up_lim_E

![E]............. choose energy interval

   ELSE IF(trim(item).eq.'E') then

30    write(*,'(2(a,f10.5),a)') 'Specify energy interval for the PDOS between ',&
           energy(1),' and ',energy(ien),' eV'
      read(*,*,err=30) en1,en2
      if(en1.ge.en2) go to 30
      if(en1.lt.energy(1)) en1=energy(1)
      if(en2.gt.energy(ien)) en2=energy(ien)
      if(en2.gt.Up_lim_E) Up_lim_E=0.0
      DoneDos=.false.

![TD].............. show total DOS or not

   ELSE IF(trim(item).eq.'TD') then
      
      if(ShowTotDos.eq.'Y') then
         ShowTotDos='N'
      else
         ShowTotDos='Y'
      end if

![Co].............. show current atoms and the choice

   ELSE IF(trim(item).eq.'Co') then

      call show_atoms_tagged(listat,ngroup)

![Q]................ quit

   ELSE IF(trim(item).eq.'Q') then
      go to 50
   ELSE
      go to 21
   END IF
   go to 21

!......... errors
19 write(*,*)'FATAL! Bad or absent SIESTA '//trim(seed)//'.PDOS file!'
   close (14)
   return
59 write(*,*)'Error in reading energies in SIESTA '//trim(seed)//'.PDOS file!'
50 close (14)
   deallocate(Labels)
   deallocate(listat)
   deallocate(energy)
   deallocate(pdos)
   deallocate(totdos)
   deallocate(index)
end subroutine read_siesta_pdos

subroutine read_orb_signature(ind,at,x,n,l,m,iz,NFIL,iErr)
!..................................................................
! it reads the next signature from the PDOS siesta file containing:
!    ind - orbital number
!    at  - atom number
!    x   - atom coordinates
!    n   - main quantum number
!    l   - orbital quantum number
!    m   - magnetic quantum number
!    iz  - polarisation number 
! Then it jumps to the data line so that PDOS could be read in 
! immediately.
!..................................................................
implicit none
integer ind,at,n,l,iz,NFIL,iErr,pos(2)
character line*200,m*2
real*8 x(3)

  call find_string('<orbital',8,line,NFIL,.false.,iErr)
  if(iErr.eq.1) return
  call position_quote(line,pos,NFIL,iErr)
  read(line(pos(1):pos(2)),*) ind
  call position_quote(line,pos,NFIL,iErr)
  read(line(pos(1):pos(2)),*) at
  read(NFIL,*)
  call position_quote(line,pos,NFIL,iErr)
  read(line(pos(1):pos(2)),*) x
  call position_quote(line,pos,NFIL,iErr)
  read(line(pos(1):pos(2)),*) n
  call position_quote(line,pos,NFIL,iErr)
  read(line(pos(1):pos(2)),*) l
  call position_quote(line,pos,NFIL,iErr)
  read(line(pos(1):pos(2)),*) m
  call position_quote(line,pos,NFIL,iErr)
  read(line(pos(1):pos(2)),*) iz
  call find_string('<data>',6,line,NFIL,.false.,iErr)
end subroutine read_orb_signature

subroutine find_orbital(iorb0,iorb,NFIL,iErr)
!..................................................................
! it jumps to the first line of the heading of the orbital iorb,
! being currently at orbital iorb0 
!..................................................................
implicit none
integer iorb0,iorb,NFIL,iErr,i
character line*200
  do i=iorb0+1,iorb
     call find_string('<orbital',8,line,NFIL,.false.,iErr)
  end do
end subroutine find_orbital

subroutine do_pdos(listat,index,pdos,ien,NFIL,jspin,iErr)
!..................................................................
! calculates the PDOS for the given selection of atoms in listat(i)
!..................................................................
use param
implicit none
integer iorb0,iorb1,iorb2,NFIL,iErr,nat,iorb,i,ien,index(NIONS,2),jspin
real*8 proj,pdos(ien),pr
logical listat(NIONS)
character line*200

   rewind(NFIL)
   pdos=0.0d0 ; iorb0=0
   do nat=1,NIONS
      if(listat(nat)) then
         iorb1=index(nat,1) ; iorb2=index(nat,2)
         call find_orbital(iorb0,iorb1,NFIL,iErr)
         if(iErr.eq.1) go to 59
         do iorb=iorb1,iorb2
            call find_string('<data>',6,line,14,.false.,iErr)
            if(iErr.eq.1) go to 59
            do i=1,ien
               if(jspin.eq.1) then
                  read(NFIL,*,err=59) proj
               else
                  read(NFIL,*,err=59) pr,proj
               end if
               pdos(i)=pdos(i)+proj
            end do
         end do
         iorb0=iorb2
      end if
   end do
   return
!......... error
59 write(*,*)'ERROR while reading PDOS file!'
   iErr=1
end subroutine do_pdos
