module siesta_EIG
use param
use code
use atoms
implicit none
real*8,dimension(:,:),allocatable :: energiesf
integer N_of_states
real*8 fermi
logical fermi_found

CONTAINS

subroutine read_siesta_fermi(Err,jspin)
!.................................................................
!  reading SIESTA output information: number of occ states and
!  the Fermi energy
!.................................................................
implicit none
character Line*200,cha1,cha2*2
integer LinEnd(100),LinPos(100),isp,NumLin,iErr,is,tot,i,nkp0,jspin,nk,nb
real*8,dimension(:),allocatable :: Num_el
logical,dimension(:),allocatable :: found_ps
logical appr,en_alloc,Err
real*8 fermi1,e
    
    allocate(Num_el(NSPEC)) ; allocate(found_ps(NSPEC))
    Err=.false. ; found_ps=.false.
    
    write(*,*)'Reading OCC states & Fermi information ...'
    open(19,file=trim(nameout),status='old',form='formatted',err=100)

!_____________________ read number of valence electrons in every species

    do is=1,NSPEC
       call find_string('atom: Called for',16,line,19,.false.,iErr)
       if(iErr.eq.1) go to 100
       call CutStr(Line,NumLin,LinPos,LinEnd,0,0,iErr)
       if(LinEnd(4)-LinPos(4)+1.eq.1) then
          read(Line(LinPos(4):LinEnd(4)),'(a1)') cha1
          cha2=' '//cha1
       else if(LinEnd(4)-LinPos(4)+1.eq.2) then
          read(Line(LinPos(4):LinEnd(4)),'(a2)') cha2
       else
          go to 100
       end if
       do i=1,NSPEC
          if(Species(i).eq.cha2) then
             isp=i ; found_ps(isp)=.true. ; go to 10
          end if
       end do
       write(*,*)'ERROR! Species '//cha2//' is NOT recognised!'
       go to 100
10     write(*,'(3a,i2,a)') &
         '... Reading pseudopotential info for species: ',Species(isp),' [',isp,']'
       call find_string('Total valence charge:',21,line,19,.false.,iErr)
       if(iErr.eq.1) go to 100
       call CutStr(Line,NumLin,LinPos,LinEnd,0,0,iErr)
       read(Line(LinPos(4):LinEnd(4)),*,err=100) Num_el(isp)
       write(*,'(a,f10.5,a)')'... found ',Num_el(isp), ' valence electrons for this species'
    end do
    do i=1,NSPEC
       if(.not.found_ps(i)) then
          write(*,'(3a)')'ERROR: species ',Species(i),' was NOT found!'
          Err=.true.
       end if
    end do
    if(Err) go to 100

!_____________ calculate the number of occupied states
    tot=0
    do i=1,NSPEC
       tot=tot+Num_el(i)*NspN(i)
    end do
    appr=.false.
    if(tot/2*2.ne.tot) then
       write(*,*)'WARNING: you have an ODD number of electrons!'
       if(ispin.eq.1) write(*,*)'ERROR: this is NOT spin-polarised calculation!'
       appr=.true.
    else
       write(*,*)'OK, you have an EVEN number of electrons!'
    end if
    if(ispin.eq.2) then
       write(*,*)' =====>     Spin polarisation is IGNORED'
       appr=.true.
    end if
    if(appr) write(*,*)' =====>     HOMO/LUMO numbers are approximate!'
    tot=tot/2

!_________________ read the Fermi energy

    fermi_found=.false.
    call find_string('siesta: Fermi energy =',22,line,19,.false.,iErr)
    if(iErr.eq.1) then
       write(*,*)'WARNING: Fermi energy is NOT found in '//trim(nameout)
    else
       call CutStr(Line,NumLin,LinPos,LinEnd,0,0,iErr)
       read(Line(LinPos(5):LinEnd(5)),*) fermi
       write(*,'(a,f10.5,a)')'     Fermi energy = ',fermi, ' eV'
       write(*,'(a)')'.......................................................'
       fermi_found=.true.
    end if
    close(19)
100 write(*,*)'WARNING! File '//trim(nameout)//' not found or bad!'

!_____________ read Fermi energy (again!) and all other energies

    en_alloc=.false.
    write(*,'(a)')'................ HOMO/LUMO information ................'
    write(*,*)'Reading from '//trim(seed)//'.EIG ...'
    open(54,file=trim(seed)//'.EIG',form='formatted',status='old',err=20) 
    read(54,*,err=20) fermi1
    write(*,'(a,f10.5,a)')'     Fermi energy = ',fermi1, ' eV'
    if(fermi_found .and. abs(fermi1-fermi).gt.0.0001 ) then
       write(*,*)'WARNING: The two Fermi energies are different!'
       write(*,*)'         Accepting the last one anyway' 
    end if
    fermi_found=.true. ; fermi=fermi1
    read(54,*,err=20) N_of_states,isp,nkp0
    write(*,'(a,i5,a,i1,a)')'... found ',nkp0,' k-points and ',isp,' spins'

    if(isp.lt.jspin) then
       write(*,*)'ERROR! The file only contains the 1st spin' ; Err=.true. ; return
    end if
    allocate(energiesf(N_of_states,nkp0)) ; en_alloc=.true.
    do nk=1,nkp0
       if(jspin.eq.1 .and. ispin.eq.1) then
          read(54,*,err=20) isp,(energiesf(nb,nk),nb=1,N_of_states)
       else if(jspin.eq.1 .and. ispin.eq.2) then
          read(54,*,err=20) isp,(energiesf(nb,nk),nb=1,N_of_states),(e,nb=1,N_of_states)
       else
          read(54,*,err=20) isp,(e,nb=1,N_of_states),(energiesf(nb,nk),nb=1,N_of_states)
       end if
    end do
    close(54) ;  go to 25
20  close(54)
    write(*,*)'WARNING: energies from '//trim(seed)//'.EIG are NOT available ...'
    write(*,'(a,i5)')'    HOMO state is number ',tot 
    write(*,'(a,i5)')'    LUMO state is number ',tot+1 
    go to 30
25  write(*,'(a,i5,a,f10.5,a)') &
         '    HOMO state is number ',tot,' and energy= ',energiesf(tot,1),' eV' 
    write(*,'(a,i5,a,f10.5,a)') &
         '    LUMO state is number ',tot+1,' and energy= ',energiesf(tot+1,1),' eV' 
30  write(*,'(a)')'.......................................................'
       
!_______________ finish/errors

    deallocate(Num_el) ; deallocate(found_ps)
    if(en_alloc) deallocate(energiesf)
    Err=.true. ; if(fermi_found) Err=.false.
end subroutine read_siesta_fermi

subroutine dos_from_EIG(jspin)
!...............................................................
! plot DOS (smeared with Gaussians) from the eigenvalues only
! containied in the [seed].EIG file for spin jspin
!...............................................................
! All available k-points in the EIG files will be used, but with
! identical weights which is not quite correct. Cannot do better
! here as SIESTA uses its own k-points. A better method is to use
! tetr > brill.dat > lev00 route and BandPoints of SIESTA.
!...............................................................
implicit none
logical en_alloc
integer :: isp,np=1000,i,j,lenght(20),jspin,nkp0,nk,nb
real*8 :: sigm=0.2,en1,en2,Up_lim=0.0,Up_lim_E=0.0,step,e,s,en10,en20
character Title*50,title_pl(20)*8,filen(20)*12,item*2
logical DoneDos

!___________ read energiesf() from the EIG file

    en_alloc=.false.
    write(*,*)'Reading from '//trim(seed)//'.EIG ...'
    open(54,file=trim(seed)//'.EIG',form='formatted',status='old',err=20) 
    read(54,*,err=20) fermi
    read(54,*,err=20) N_of_states,isp,nkp0
    write(*,'(a,i5,a,i1,a)')'... found ',nkp0,' k-points and ',isp,' spins'
    if(isp.lt.jspin) then
       write(*,*)'ERROR! The file only contains the 1st spin' ; return
    end if
    allocate(energiesf(N_of_states,nkp0)) ; en_alloc=.true.
    do nk=1,nkp0
       if(jspin.eq.1 .and. ispin.eq.1) then
          read(54,*,err=20) isp,(energiesf(nb,nk),nb=1,N_of_states)
       else if(jspin.eq.1 .and. ispin.eq.2) then
          read(54,*,err=20) isp,(energiesf(nb,nk),nb=1,N_of_states),(e,nb=1,N_of_states)
       else
          read(54,*,err=20) isp,(e,nb=1,N_of_states),(energiesf(nb,nk),nb=1,N_of_states)
       end if
    end do
    go to 25
20  close(54)
    write(*,*)'WARNING: energies from '//trim(seed)//'.EIG are NOT available ...'
    go to 50

25  close(54)

!____________ find the boundaries
    en1=1.0e10 ; en2=-1.0e10
    do nk=1,nkp0
       if(energiesf(1,nk).lt.en1) en1=energiesf(1,nk)
       if(energiesf(N_of_states,nk).gt.en2) en2=energiesf(N_of_states,nk)
    end do
    en10=en1 ; en20=en2

!____________ run the menu
    DoneDos=.false.
21  write(*,*)'.....................................................'
    write(*,'(a,f10.5,a)')'   Fermi energy = ',Fermi,' eV'
    write(*,*)'Choose the necessary options:'
    write(*,*)
    write(*,'(a,f10.5)')'Sg. Sigma for the Gaussian smearing (eV):',sigm
    if(DoneDos) then
       write(*,'(a)')'CD. Calculate the DOS <= Done!'
       write(*,'(a)')'    Data file: dos.dat'
       write(*,'(a)')' P. Plot the DOS'
    else
       write(*,'(a)')'CD. Calculate the DOS'
    end if
    write(*,'(/a)')'=========== Settings for the plot ================'
    write(*,'(a,2(f10.5,a))')' E. Change energy interval [',en1,',',en2,']'
    write(*,'(a,i5)')' N. Number of points in the graph: ',np
    step=(en2-en1)/np
    if(Up_lim.eq.0.0) then
       write(*,'(a)')' Y. Upper limit to the dos-range of the plot: DISABLED'
    else
       write(*,'(a,f10.5)') ' Y. Upper limit to the dos-range of the plot: ',Up_lim
    end if
    if(Up_lim_E.eq.0.0) then
       write(*,'(a)')' X. Upper limit to the energy-range of the plot: DISABLED'
    else
       write(*,'(a,f10.5)') &
            ' X. Upper limit to the energy-range of the plot: ',Up_lim_E
    end if
    write(*,'(a)')' Q. Quit'
    write(*,*)
    write(*,*)'------------>'
    read(*,'(a)') item

![N]............ number of points

    IF(trim(item).eq.'N') then

13     write(*,*)'Specify number of points in the graph:'
       read(*,*,err=13) np
       if(np.lt.1) go to 13
       DoneDos=.false.

![Sg]............ assign sigma for the Gaussian smearing

    ELSE IF(trim(item).eq.'Sg') THEN
14     write(*,*)'Specify the Gaussian smearing (eV):'
       read(*,*,err=14) sigm
       if(sigm.lt.0.01) go to 14
       DoneDos=.false.

![Y]............ assign some upper limit on the y-axis for the plot

    ELSE IF(trim(item).eq.'Y') THEN
78     write(*,*)'Specify the upper limit for y-axis on your plot:'
       read(*,*,err=78) Up_lim

![X]............ assign some upper limit on the x-axis for the plot

    ELSE IF(trim(item).eq.'X') THEN
79     write(*,*)'Specify the upper limit for x-axis on your plot:'
       read(*,*,err=79) Up_lim_E

![E]............. choose energy interval

    ELSE IF(trim(item).eq.'E') then

30     write(*,'(2(a,f10.5),a)') 'Specify energy interval for the DOS between ',&
           en10-3*sigm,' and ',en20+3*sigm,' eV'
       read(*,*,err=30) en1,en2
       if(en1.ge.en2) go to 30
       if(en1.lt.en10-3*sigm) en1=en10-3*sigm ; if(en2.gt.en20+3*sigm) en2=en20+3*sigm
       if(en2.gt.Up_lim_E) Up_lim_E=0.0
       DoneDos=.false.

![CD].............. calculate PDOS

    ELSE IF(trim(item).eq.'CD') then

       open(24,file='dos.dat',form='formatted')
       e=en1-step
       do i=1,np
          s=0.0
          do j=1,N_of_states
             do nk=1,nkp0
                s=s+Gauss(e-energiesf(j,nk),sigm)
             end do
          end do
          write(24,'(2(f10.5,x))') e,s/nkp0
          e=e+step
       end do
       close (24)
       filen(1)='dos.dat' ; lenght(1)=7 ; title_pl(1)='   '
       write(*,*)'File '//filen(1)(1:lenght(1))//' has been written!' 
       DoneDos=.true.

![P].............. plot PDOS

    ELSE IF(trim(item).eq.'P' .and. DoneDos) then
      
       Title='                                             '
       call Plot2(1,filen,lenght,Title,title_pl, &
            'Energy (eV)         ', &
            '     DOS (arb.units)', 'Screen', 33, &
            Up_lim,Up_lim_E,'O')
       
![Q]................ quit

    ELSE IF(trim(item).eq.'Q') then
       go to 50
    ELSE
       go to 21
    END IF
    go to 21

50  close (14)
    if(en_alloc) deallocate(energiesf)
end subroutine dos_from_EIG

real*8 function Gauss(e,s)
implicit none
real*8 e,s,s2
real*8,parameter :: pi=3.141592654
  s2=s*s
  Gauss=1./sqrt(2*pi*s2) * exp(-e*e/(2*s2))
end function Gauss

end module siesta_EIG
