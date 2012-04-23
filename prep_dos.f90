subroutine prep_dos(iFlag)
!.....................................................................
!    It calculates "band group DOS" from PSI2(k-point,band,task)
!
! The "band group DOS" means the sum of individual bands contributions
! to the DOS for bands which enter the group i=[Nbands0(i),Nbands1(i)].
!.....................................................................
! iFlag = 0 - only total DOS needs to be calculated;
! iFlag = 1 - both total and projected DOS need to be calculated.
!.....................................................................
use param
use atoms
use dos
implicit none
!.....................................................................
real*8, parameter :: pi=3.1415927
real*8,dimension (:,:), allocatable :: E
real*8,dimension (:), allocatable :: B_min,B_max,B_step,B0_min,B0_max,B_step1
integer,dimension(:),allocatable :: N_min,N_max, N0_min,N0_max,Nbands0,Nbands1
integer lenght(99),k_pnt_ref(4),Num_B0,i0,notet,i,iTask,ijk,iFlag,n_energy
integer n_en,i1,j,i2,k,len,it
integer NKP,NB,iErr,ngroup,Num_B
real*8 factor,v,E_step,Dispers2,Up_lim,Up_lim_E,Emin,Emax
character cha1*1,cha2*2,answer*1,filen(99)*12,cha,line*16,cha13*13
character Title*50,title_pl(99)*7,item*2,chal*2
character IndivPrev,IndivPstsc,cha16*35,Yes_Job,Smeared
logical SmearedPrev,Missing,Yes_Fast,Smeared_,k_ref_alloc,Err
character task_list*200
data E_step/0.1/,IndivPrev/'N'/,IndivPstsc/'N'/, &
     ijk/1/,Yes_Job/'N'/, &
     Smeared/'N'/,SmearedPrev/.false./, &
     Up_lim/0.0/,Yes_Fast/.false./,Up_lim_E/0.0/
data Title/'                                                  '/
!
!.....................................................................
!...................... PREPARATION ..................................
!.....................................................................
!

!_________ read energies from the file band.out
!
    allocate(E(NKPTS,NBANDS))
    call get_energies(E,Emin,Emax,jspin,Err)
    if(.not.Err) then
       write(*,'(a23,e12.6,a5,e12.6,a3)')'Energies E are between ', &
                                               Emin,' and ',Emax,' eV'
    else
       deallocate(E) ; return
    end if

!.... allocate memory
    allocate (B_step(NBANDS)) ; allocate (B_step1(NBANDS)) ; allocate (B0_min(NBANDS))
    allocate (B0_max(NBANDS)) ; allocate (B_min(NBANDS)) ;   allocate (B_max(NBANDS))
    allocate (N_min(NBANDS)) ;  allocate (N_max(NBANDS)) ;   allocate (N0_min(NBANDS))
    allocate (N0_max(NBANDS)) ; allocate(Nbands0(NBANDS)) ;  allocate(Nbands1(NBANDS))
    if(iFlag.eq.1) then
       allocate(band_yes(NBANDS)) ; allocate(which_task(Ntask)) ; which_task=.false.
    end if
    k_ref_alloc=.false.
    ngroup=0
!
!_________ and  analyze the structure of physical bands:
!      {B0_min(i),B0_max(i)}, i=1,...,Num_B0 <--- energies
!      {N0_min(i),N0_max(i)}, i=1,...,Num_B0 <--- states
!
    call phys_band(E,B0_min,B0_max,Num_B0)
    call phys_band_nmb(E,B0_min,B0_max,Num_B0,N0_min,N0_max)
!
!_________ a loop over all the tetrahedra stored in the file 'brill.dat'
!  k_ref(k-point,edge) - gives the k-point number (in the common list)
!                        for every particular edge of the tetrahedra;
!
    open(44,file='brill.dat',status='old',form='formatted',err=100)

!... find the number of tetrahedra first
    read(44,'(a)',err=120,end=120) line
    Ntetra=0
10  read(44,*,err=120,end=30) notet
    if(notet-Ntetra.ne.1) go to 140
    Ntetra=notet
    read(44,*,err=120) (k_pnt_ref(k),k=1,4)
    do k=1,4
       if(k_pnt_ref(k).gt.NKPTS) then
          write(*,*)'ERROR! <brill.dat> file does NOT correspond to this setup!'
          write(*,*)'       Cannot proceed with the DOS. Leaving ...'
          go to 120
       end if
    end do
    go to 10
30  write(*,'(a,i8)')'Total number of tetragedra is ',Ntetra

!.... allocate memory
    allocate(k_ref(Ntetra,4)) ; k_ref_alloc=.true.

!.... read k-points and tetrahedra information
    rewind(44)
    read(44,'(a)',err=120,end=120) line
    i0=INDEX(line,'yes')
    if(i0.eq.1) write(*,*) 'WARNING! Your data belong to the compressed z-direction!'

!________ fill in k_ref(notetra,edge) using the reference aray
! k_pnt_ref(edge); nkp - is used to count k-vectors;
!
    nkp=0
    do i=1,Ntetra
       read(44,*,err=120,end=30) notet,(k_pnt_ref(k),k=1,4)
       do k=1,4
          k_ref(notet,k)=k_pnt_ref(k)
          if(k_pnt_ref(k).gt.nkp) then
             nkp=nkp+1
             if(nkp.ne.k_pnt_ref(k)) go to 140
          end if
       end do
    end do
    close (44)

!________ calculate a factor = 2* (UC volume) / (2*pi)**3
    factor = 2*VOLC/(2*pi)**3
!
!________ update the factor if calculation is spin-polarised
    if(ispin.eq.2) factor=factor*0.5
!
!______ volume of any small cube (in the reciprocal space)
!
    v = (2*pi)**3/VOLC / (Ntetra/6.0) ; factor = factor*v
    write(*,'(a22,f10.6)')'Small cube volume v = ',v
    write(*,'(a22,f10.6)')'Final value of factor = ',factor
!
!.....................................................................
!............. START MAIN MENU HERE ..................................
!.....................................................................
!  ijk - counts different cycles of calculations (no more than 9);
!        for instance, one can choose a different problem set 
!
!_________  analyze the structure of physical bands:
!      {B_min(i),B_max(i)}, i=1,...,Num_B
!
    iTask=0 
21  if(Smeared.eq.'S') then
       call phys_band_sm(B_min,B_max,Num_B,B0_min,B0_max,Num_B0)
    else
       Num_B=Num_B0
       do i=1,Num_B
          B_min(i)=B0_min(i) ; B_max(i)=B0_max(i)
       end do
    end if
    call phys_band_nmb(E,B_min,B_max,Num_B,N_min,N_max)
    write(*,*)' '
    write(*,*)'_______> Structure of "physical" bands <________'
    do i=1,Num_B
       write(*,'(a6,i5,a25,f10.5,a3,f10.5,a14,i5,a1,i5,a1)') &
            ' Band ',i,' spans energy interval: [',B_min(i),' , ', &
            B_max(i),'] and states {',N_min(i),',',N_max(i),'}'
    end do
    write(*,*)' '
    
    write(*,*)'..............MENU DOS/LDOS .......................'
    write(*,*)'..... Change these parameters if necessary:...'
    write(*,*)

    if(iFlag.eq.0) then
       write(*,*)'  ______________> Only total DOS <_______________'
    else
       write(*,*)'  _____________> LDOS + total DOS <______________'
    end if
    
    write(*,'(a30,i2)')'   0. THE CURRENT PICTURE IS: ',ijk
    if(ijk.lt.10) then
       write(cha,'(i1)') ijk
    else
       write(cha2,'(i2)') ijk
    end if

    IF (iFlag.eq.0) THEN
       iTask=1
    ELSE
       write(*,'(a)')'  ST. Show all tasks [chosen tasks are indicated with Y]'
       if(iTask.eq.0) then
          write(*,'(a)')'  CT. Current PDOS task(s): ... undefined ...'
       else
          write(*,'(a)')'  CT. Add to the current PDOS task(s):'
          write(*,'(a)')'    '//task_list(1:len)
       end if
       write(*,'(a)')'   1. Choose/Add a single PDOS task'
    END IF
    
    if(E_step.eq.0.0) then
       write(*,'(a)') &
            '   2. Maximum energy step for the plotting: ... undefined ...'
    else
       write(*,'(a44,f10.5)') &
            '   2. Maximum energy step for the plotting: ',E_step
       n_energy=0
       do i=1,Num_B
          n_en=(B_max(i)-B_min(i))/E_step+1
          if(n_en.lt.10) then
             n_en=10
             B_step(i)=(B_max(i)-B_min(i))/n_en
          else
             B_step(i)=E_step
          end if
          n_energy=n_energy+n_en
       end do
       write(*,'(a42,i5)') '     > Number of points for the DOS plot: ',n_energy
       write(*,'(a)') '     > Steps for each "physical" band:'
       do i=1,Num_B,6
          i1=i+5 ; if(i1.gt.Num_B) i1=Num_B
          write(*,'(7x,6(f8.5,1x))') (B_step(j),j=i,i1)
       end do
       if(Smeared.eq.'S') then
          Missing=.false.
          call step_ph_band(B_min,B_max,Num_B,E_step,B_step1,B0_min,B0_max,Num_B0)
          do i=1,Num_B
             if(B_step(i) .gt. B_step1(i)) Missing=.true.
          end do
          if(Missing) write(*,'(a54)') &
               '     > WARNING! Some of original bands may be missing!'
       end if
    end if
    
    if(Smeared.eq.'S') then
       write(*,'(a)') '   3. Smeared DOS/LDOS: ENABLED slow method'
       Smeared_=.true.
    else if(Smeared.eq.'F') then
       write(*,'(a)') '   3. Smeared DOS/LDOS: ENABLED fast method'
       Smeared_=.false.
    else
       write(*,'(a)') '   3. Smeared DOS/LDOS: DISABLED'
       Smeared_=.false.
    end if
    
    if(Smeared.ne.'N') then
       write(*,'(a33,f6.3,a12)') &
            '  BB. Broaden band boundaries by ',Broad_Band,' dispersions'
       if(Dispers.le.0.0001) Dispers=0.01
       write(*,'(a35,f10.5)') '   D. Dispersion for the smearing: ',Dispers
       
       Dispers2=Dispers*sqrt(2.0)
       if(SmearedPrev) then
          write(*,'(a)') '  Pr. For previewing: smeared DOS/LDOS '
       else
          write(*,'(a)') '  Pr. For previewing: original DOS/LDOS '
       end if
    end if

    if(ngroup.eq.0) then
       write(*,'(a)') '  Gr. Number of groups of bands: ... undefined ...'
    else
       write(*,'(a,i5)') '  Gr. Number of groups of bands: ',ngroup
       write(*,'(a)')    '     > The following groups of bands are defined:'
       if(ngroup.gt.4) then
          i2=ngroup-4
          write(cha16,2) ngroup-i2
          write(*,cha16)(i,Nbands0(i),Nbands1(i),i=1,4)
          write(cha16,2) i2
          write(*,cha16)(i,Nbands0(i),Nbands1(i),i=5,ngroup)
       else
          write(cha16,2) ngroup
          write(*,cha16)(i,Nbands0(i),Nbands1(i),i=1,ngroup)
       end if
2      format('(10x,',i1,'(''{'',i1,''} '',i5,''_'',i5,''; ''))')
       do i=1,ngroup
          write(cha1,'(i1)') i
          if(ijk.lt.10) then
             filen(i)='dos.dat'//cha1//'_'//cha
             lenght(i)=10
          else
             filen(i)='dos.dat'//cha1//'_'//cha2
             lenght(i)=11
          end if
       end do
       if(ngroup.eq.1) then
          write(*,'(a)') '     > File for the plot: '//filen(1)
          if(Smeared.eq.'F') write(*,'(a)') '     > File with smearing plot: '// &
               filen(1)(:lenght(1))//'_'
          write(*,'(a)') '     > Name of the PostScript file: '// &
               filen(1)(:lenght(1))//'.ps'
       else
          write(*,'(a)') '     > Files for the plot: '// &
               filen(1)//' , ... , '//filen(ngroup)
          if(Smeared.eq.'F') write(*,'(a)') &
               '     > Files with smearing plot: '//filen(1)(:lenght(1))// &
               '_ , ... , '//filen(ngroup)(:lenght(ngroup))//'_'
          write(*,'(a)') '     > Names of the PostScript files: '// &
               filen(1)(:lenght(1))//'.ps' &
               //' , ... , '//filen(ngroup)(:lenght(ngroup))//'.ps'
       end if
    end if
    
    if(Yes_Job.eq.'N') then
       write(*,'(a)') '   C. Adopt present setting and calculate DOS/LDOS'
    else
       write(*,'(a)') '   C. Adopt present setting and calculate DOS/LDOS <-- DONE!'
       if(Smeared.eq.'F'.and.Yes_Fast) then
          write(*,'(a)') '  CS. Adopt present setting and smear the DOS/LDOS <-- DONE!'
       else if(Smeared.eq.'F'.and. (.not.Yes_Fast)) then
          write(*,'(a)') '  CS. Adopt present setting and smear the DOS/LDOS'
       end if
    end if
    
    if(ngroup.gt.1) then
       if(IndivPrev.eq.'Y') then
          write(*,'(a)') '   9. Do you want to preview each individual DOS: YES'
       else
          write(*,'(a)') '   9. Do you want to preview each individual DOS: NO'
       end if
       if(IndivPstsc.eq.'Y') then
          write(*,'(a)') '  10. To make PostScript for each individual DOS: YES'
       else
          write(*,'(a)') '  10. To make PostScript for each individual DOS: NO'
       end if
       write(*,'(a)')    '   P. Preview ALL DOS curves for the current group'
       write(*,'(a)')    '  Ps. Create a PostScript file dos.dat0'// &
            filen(1)(9:lenght(1))//' for ALL DOS curves'
    else
       write(*,'(a)')    '   P. Preview the current DOS/LDOS'
       write(*,'(a)')    '  Ps. PostScript file for the current DOS/LDOS'
    end if
    
    if(Up_lim.eq.0.0) then
       write(*,'(a)')    '  Ud. Upper limit to the dos-range of the plot: DISABLED'
    else
       write(*,'(a46,f10.5)') '  Ud. Upper limit to the dos-range of the plot: ',Up_lim
    end if
    if(Up_lim_E.eq.0.0) then
       write(*,'(a)')   '  Ue. Upper limit to the energy-range of the plot: DISABLED'
    else
       write(*,'(a46,f10.5)') &
            '  Ue. Upper limit to the energy-range of the plot: ',Up_lim_E
    end if
    write(*,'(a)')'   Q. Quit: do not proceed.'
    write(*,*)
    write(*,*)'------> Choose the item and press ENTER:'
    read(*,'(a)',err=101) item
!
!............ start a new job ijk
!
     IF(trim(item).eq.'0') THEN
        if(ijk.eq.99) then
           write(*,*) 'PREP_DOS: You cannot increment the job number anymore!'
           write(*,*) 'PREP_DOS: Press ENTER!'
           read(*,*)
        else
           ijk=ijk+1 ; iTask=0 ; ngroup=0 ; Yes_Job='N' ; Yes_Fast=.false.
        end if
!
!............ show all tasks indicating if any were already chosen
!
     ELSE IF(iFlag.eq.1 .and. trim(item).eq.'ST') THEN
        write(*,*)'...... You have the following tasks: ......'
        do it=1,Ntask
           chal=' N' ; if(which_task(it)) chal=' Y'
           if(CaseDos.eq.'Sphere') then
              if(mod(it+3,4).eq.0) then
                 write(*,'(2(a,i5),a,a,a,f5.2,a,3(f7.3,1x),a)') &
                      'at= ',v_atm(it),' '//Species(v_spec(it))//' ',&
                       it,' task {'//Flag(it)//'}', &
                       chal,' Radius= ',Radius(it), ' center= (', &
                      (Point(i,it),i=1,3),') '//method(it)
              else
                 write(*,'(13x,i5,2a)') it,' task {'//Flag(it)//'}',chal
              end if
           end if
        end do
        write(*,*)'Press ENTER when done ...' ; read(*,*)
!
!............ choose the PDOS problems 
!
     ELSE IF(iFlag.eq.1 .and. trim(item).eq.'CT') THEN
        
        call choose_tasks(iTask,task_list,len)
        Yes_Job='N' ; Yes_Fast=.false.
!
!............ choose a single PDOS task
!
     ELSE IF(iFlag.eq.1 .and. trim(item).eq.'1') THEN

109     write(*,'(a25,i5,a2)')'Choose one PDOS task between 1 and ',Ntask,' :'
        write(*,*) '[ It will be added to the current list of tasks]'
        read(*,*,err=109) it
        if(it.lt.1.or.it.gt.Ntask) go to 109
        which_task(it)=.true.
!
!............ specify maximum allowed energy step E_step
!
     ELSE IF(trim(item).eq.'2') THEN

31      write(*,*)'Enter max(STEP) on the energy scale:'
        read(*,*,err=31) E_step
        Yes_Job='N' ; Yes_Fast=.false.
!
!............ specify if the smearing DOS/LDOS curves should be also
!             calculated and what method is to be used: fast or slow
!
     ELSE IF(trim(item).eq.'3') THEN
        if(Smeared.eq.'N') then
           Smeared='F'
        else if(Smeared.eq.'F') then
           Smeared='S' ; Yes_Job='N'
        else
           Smeared='N' ; Yes_Job='N'
        end if
        Yes_Fast=.false.
!
!............ Broaden band boundaries by a certain number of dispersions
!
     ELSE IF(Smeared.ne.'N' .and. trim(item).eq.'BB') THEN
37      write(*,*)'Enter broadening factor (in terms of dispersion):'
        read(*,*,err=37) Broad_Band
!
!............ specify the dispersion for the Gaussian smearing
!
     ELSE IF(Smeared.ne.'N' .and. trim(item).eq.'D') THEN
 32     write(*,*)'Specify dispersion for the Gaussian smearing:'
        read(*,*,err=32) Dispers
!
!............ specify if preview smeared or original DOS/LDOS curves
!
     ELSE IF(Smeared.ne.'N' .and. trim(item).eq.'Pr') THEN
        if(SmearedPrev) then
           SmearedPrev=.false.
        else
           SmearedPrev=.true.
        end if
!
!................. Choose the bands. All allowed bands are ..........
!  gathered into groups [Nbands0(i),Nbands1(i)] of bands which
!  are treated individually.
!
     ELSE IF(trim(item).eq.'Gr') THEN
108     call group_band(Nbands0,Nbands1,ngroup)
        do i=1,ngroup
           write(*,'(a18,i3,a13)')'..........< Group ',i,' >...........'
           call minmax(E,Emin,Emax,Nbands0(i),Nbands1(i))
           write(*,'(a23,e12.6,a5,e12.6,a3)')'Energies E are between ', &
                Emin,' and ',Emax,' eV'
        end do
!_______ find out whether each band is available or not (band_yes)
        if(iFlag.eq.1) band_yes='y'
        Yes_Job='N'
        Yes_Fast=.false.
!
!............ calculate the DOS
!............ General loop over groups of bands
!
     ELSE if(trim(item).eq.'C') THEN
        if(iTask.eq.0 .or. E_step.eq.0.0 .or. ngroup.eq.0) then
           write(*,*)'ERROR! You still have undefined parameters!'
           go to 21
        end if
        do i=1,ngroup
           write(*,'(a18,i3,a13)')'..........< Group ',i,' >..........'
           call minmax(E,Emin,Emax,Nbands0(i),Nbands1(i))
           write(*,'(a23,e12.6,a5,e12.6,a3)')'Energies E are between ', &
                Emin,' and ',Emax,' eV'
           write(title_pl(i),'(i3,a1,i3)') Nbands0(i),'_',Nbands1(i)
           write(*,'(a)') title_pl(i)
!
!__________ open the file for the output first
!
           open(31,file=filen(i)(:lenght(i)),status='unknown',form='formatted')
           write(*,*) 'The file '//filen(i)//' has been opened for the DOS.'
           if(Smeared.eq.'S') then
              write(31,'(a)') '#---energy---Total DOS--'// &
                   '-Smeared---Projected DOS---Smeared'
           else
              write(31,'(a)') '#---energy---Total DOS---Projected DOS-'
           end if
!
!________ calculate the DOS/LDOS and write it into file UNIT=31 filen(i)
!
           call do_dos(31,iFlag,E,B_step,Num_B,B_min,B_max, &
                 Dispers2,Smeared_,Nbands0(i),Nbands1(i),factor)
           close (31)
           write(*,*)'The file '//filen(i)//' has been created!'
!
!_________ IN the case of the fast method: do smearing here
!
           if(Smeared.eq.'F') then
              call smear_fast(filen(i),lenght(i),filen(i),lenght(i), &
                    Num_B0,B0_min,B0_max,N0_min,N0_max, &
                    Nbands0(i),Nbands1(i),iFlag,E_step,iErr)
              if(iErr.eq.0) Yes_Fast=.true.
           end if
!
!_________ Plot the DOS for the current group
!
           if(IndivPrev.eq.'Y'.and.ngroup.gt.1) &
                call Plot1(filen(i),lenght(i),Title,title_pl(i), &
                 'Energy (eV)         ', &
                 'DOS (arb.units)     ', 'Screen', 33,iFlag, &
                 Smeared,SmearedPrev,Up_lim,Up_lim_E)
           if(IndivPstsc.eq.'Y'.and.ngroup.gt.1) then
              write(*,*)'Give the title:'
              read(*,'(a)') Title
              call Plot1(filen(i),lenght(i),Title,title_pl(i), &
                    'Energy (eV)         ', &
                    'DOS (arb.units)     ', 'Postsc', 33,iFlag, &
                    Smeared,SmearedPrev,Up_lim,Up_lim_E)
           end if
        end do
        Yes_Job='Y'
!
!............ if to preview each time while doing a group?
!
     ELSE IF(ngroup.gt.1.and.trim(item).eq.'9') THEN
        if(IndivPrev.eq.'Y') then
           IndivPrev='N'
        else
           IndivPrev='Y'
        end if
!
!............ if to make a PoastScript file each time while doing a group?
!
      ELSE IF(ngroup.gt.1.and.trim(item).eq.'10') THEN
         if(IndivPstsc.eq.'Y') then
            IndivPstsc='N'
         else
            IndivPstsc='Y'
         end if
!
!............ preview DOS/LDOS for the current group
!
      ELSE IF(trim(item).eq.'P') THEN
         if(Yes_Job.eq.'Y') then
            if(Smeared.eq.'F'.and.(.not.Yes_Fast)) then
               write(*,*)'IGNORED! You should run the item CS first!'
            else
               call Plot(ngroup,filen,lenght,Title,title_pl, &
                    'Energy (eV)         ', &
                    'DOS (arb.units)     ', 'Screen', 33,iFlag, &
                    Smeared,SmearedPrev,Up_lim,Up_lim_E)
            end if
         else
            write(*,*)'IGNORED! You should run the item 8 first!'
         end if
!
!............ make a PostScript file for the current group
!
      ELSE IF(trim(item).eq.'Ps') THEN
        if(Yes_job.eq.'Y') then
          if(Smeared.eq.'F'.and.(.not.Yes_Fast)) then
            write(*,*)'IGNORED! You should run the item CS first!'
          else
            write(*,*)'Give the title:'
            read(*,'(a)') Title
            call Plot(ngroup,filen,lenght,Title,title_pl, &
                'Energy (eV)         ', &
                'DOS (arb.units)     ', 'Postsc', 33,iFlag, &
                 Smeared,SmearedPrev,Up_lim,Up_lim_E)
          end if
        else
          write(*,*)'IGNORED! You should run the item 8 first!'
        end if
!
!............ assign some upper limit on the y-axis for the plot
!
     ELSE IF(trim(item).eq.'Ud') THEN
78      write(*,*)'Specify the upper limit for y-axis on your plot:'
        read(*,*,err=78) Up_lim
!
!............ assign some upper limit on the x-axis for the plot
!
     ELSE IF(trim(item).eq.'Ue') THEN
79      write(*,*)'Specify the upper limit for x-axis on your plot:'
        read(*,*,err=79) Up_lim_E
!
!............ run FAST smearing for the current unsmeared DOS/LDOS
!
     ELSE IF(Smeared.eq.'F'.and.Yes_Job.eq.'Y'.and.trim(item).eq.'CS') THEN
        do i=1,ngroup
           call smear_fast(filen(i),lenght(i),filen(i),lenght(i), &
                   Num_B0,B0_min,B0_max,N0_min,N0_max, &
                   Nbands0(i),Nbands1(i),iFlag,E_step,iErr)
        end do
        if(iErr.eq.0) Yes_Fast=.true.
!
!............ quit: return to the main routine
!
     ELSE IF(trim(item).eq.'Q') THEN
        go to 200
     ELSE
        go to 101
     END IF
     go to 21
!............. error
101  write(*,*)'Incorrect item number! Try again!'
     go to 21
!
!................................................................ 
!.... ........ END OF THE MAIN DOS MENU .........................
!................................................................
!
!............... errors!
100  write(*,*)'FATAL: error while opening <brill.dat> file!'
     go to 200
120  write(*,*)'FATAL: error while reading <brill.dat> file!'
     go to 200
140  write(*,*)'FATAL: internal inconsistency in <brill.dat> file!'
!
!............ finish
200  deallocate (E) ; deallocate (B_step) ; deallocate (B_step1) ; deallocate (B0_min)
     deallocate (B0_max) ; deallocate (B_min) ; deallocate (B_max) ; deallocate (N_min)
     deallocate (N_max) ; deallocate (N0_min) ; deallocate (N0_max) 
     deallocate(Nbands0) ; deallocate(Nbands1) 
     if(k_ref_alloc) deallocate(k_ref)
     if(iFlag.eq.1) then
        deallocate(PSI2) ; deallocate(Point)  ; deallocate(Radius)
        deallocate(Flag) ; deallocate(method) 
        deallocate(v_atm) ; deallocate(v_spec) 
        deallocate(band_yes) ; deallocate(which_task)
     end if
end subroutine prep_dos

subroutine do_dos(NIT,iFlag,E,B_step,Num_B,B_min,B_max, &
          Disp,Smeared,Nbands0,Nbands1,factor)
!.....................................................................
! Calculates DOS and LDOS and writes it into file UNIT=NIT.
!.....................................................................
! iFlag = 0 - only total DOS needs to be calculated;
! iFlag = 1 - both total and projected DOS need to be calculated.
!.....................................................................
! Smeared=.true.  - the smearing is provided as well as the 3rd and 5th
!                   columns in the file (unit=NIT) for the DOS and LDOS,
!                   respectively.
! Smeared=.false. - the smearing is not provided, i.e. we write DOS
!                   and LDOS to the columns 2 and 3, respectively.
!.....................................................................
use param
use dos
implicit none
real*8 E(NKPTS,NBANDS),B_min(NBANDS),B_max(NBANDS),B_step(NBANDS)
integer i_sort(4),iFlag,NIT,Num_B,Nbands0,Nbands1,i,k,iPh,NB,notet,it
integer k0123(4),nk
real*8 E0123(4),A0123(4),Disp,factor,tiny,pos_num
real*8 E0,E1,B_low,B_high,el_num,avr_en,S_tot1,E_step,energy,tot,proj
real*8 S_tot,S_proj,S_tot_sm,S_proj_sm,tot_NB,proj_NB,tot_NB_sm,proj_NB_sm
real*8,dimension(:,:),allocatable :: psi3
logical Smeared
data tiny/0.001/
!
!........ obtain exact energy boundaries [E0,E1] between
!         Nbands0 and Nbands1
!
   E0=1000000.0 ; E1=-E0
   do i=Nbands0,Nbands1
      do k=1,NKPTS
         if(E(k,i).lt.E0) E0=E(k,i)
         if(E(k,i).gt.E1) E1=E(k,i)
      end do
   end do
!
!.......... prepare psi3 as a sum of psi2 over all PDOS tasks
!
   if(iFlag.eq.1) then
      allocate(psi3(NKPTS,NBANDS)) ; psi3=0.0d0
      do nk=1,NKPTS
         do nb=1,NBANDS
            do it=1,Ntask
               if(which_task(it)) psi3(nk,nb)=psi3(nk,nb)+PSI2(nk,nb,it)
            end do
         end do
      end do
   end if
!
!......................................................................
!.............. main loop over "physical" bands (1,...,Num_B) .........
!               with the step E_step=B_step(phys.band)
!......................................................................
!  el_num - number of electrons inside every band (it is calculated
!           by the Simpson's method directly by integrating the total
!           DOS with respect to the energy inside every band)
!  avr_en - average energy for every band found
!  energy - current value of the energy within the physical band
!......................................................................
!
!................. general loop over "physical" bands
!
   DO 100 iPh=1,Num_B
!....... define actual boundaries of the bands
      if(B_min(iPh).lt.E0.and.B_max(iPh).lt.E0) go to 100
      if(B_min(iPh).gt.E1.and.B_max(iPh).gt.E1) go to 100
      if(B_step(iPh).lt.tiny) then
         write(*,*) 'WARNING! Band ',iPh, &
              ' will be missing due to a very small step!'
         go to 100
      end if
      B_low=B_min(iPh) ; B_high=B_max(iPh)
!        if(B_low.lt.E0) B_low=E0
!        if(B_high.gt.E1) B_high=E1
      write(*,'(a14,i3,a14,f10.5)') &
           '.... The band ',iPh,' starts at E= ',B_low-tiny
      if(iFlag.eq.0) then
         if(Smeared) then
            write(NIT,'(f10.5,a)') B_low-tiny,' 0.0 0.0'
         else
            write(NIT,'(f10.5,a)') B_low-tiny,' 0.0'
         end if
      else
         if(Smeared) then
            write(NIT,'(f10.5,a)') B_low-tiny,' 0.0 0.0 0.0 0.0'
         else
            write(NIT,'(f10.5,a)') B_low-tiny,' 0.0 0.0'
         end if
      end if
!....... initialise number of electrons in the band and the average
!        energy
      el_num = 0.0 ; avr_en = 0.0 ; S_tot1 = 0.0
      E_step=B_step(iPh) ; energy=B_low-E_step
      
      write(*,*) energy, E_step,iPh,B_step(iPh)
!
!....... energy loop inside current "physical" band
!
50    energy=energy+E_step
!
!________ get contributions from every tetrahedra to:
!  total DOS:        S_tot  ;  smeared - S_tot_sm
!  projected DOS:    S_proj ;  smeared - S_proj_sm
!
      S_tot=0.0 ; S_proj=0.0 ; S_tot_sm=0.0 ; S_proj_sm=0.0
!
!................. loop over bands inside the group ...................
!
      do 90 NB=Nbands0,Nbands1
!____ asks whether the band NB enters in some island
         if(iFlag.eq.1) then
            if(band_yes(NB).eq.'n') go to 90
         end if
         tot_NB=0.0 ; proj_NB=0.0 ; tot_NB_sm=0.0 ; proj_NB_sm=0.0
!
!............ loop over all the tetrahedra
!
         do 80 notet=1,Ntetra
!
!___________ sort 4 edges of the tetrahedra in increasing order of
!            their energies for the band NB:
! i_sort(new_order_of_edges) - old order of edges of the tetrahedra;
!____ energies will be ordered in such a way that:
!  E(k0,NB)  =< E(k1,NB) =< E(k2,NB) =< E(k3,NB) ):
!
            call sort(E(k_ref(notet,1),NB),E(k_ref(notet,2),NB), &
                     E(k_ref(notet,3),NB),E(k_ref(notet,4),NB),i_sort)
!
!__________ correct order of k-points for corners 0,1,2,3 in k0123;
!           energies on every corner in E01234;
!           the weighting function on the corners in A0123:
!
            do k=1,4
               k0123(k) = k_ref(notet,i_sort(k))
               E0123(k) = E(k0123(k),NB)
!_________________ for PDOS sum up all chosen
               if(iFlag.eq.1) A0123(k) = psi3(k0123(k),NB)
            end do
!
!___________ calculate the contributions to both total (tot) and
!   projected  (proj) DOS:
!
            call do_contr(tot,proj,A0123,E0123,energy,iFlag)
            tot_NB = tot_NB + tot
            proj_NB = proj_NB + proj
!
!___________ calculate smeared DOS and LDOS
!
            if(Smeared) then
               call do_contr_sm(Disp,tot,proj,A0123,E0123,energy,iFlag)
               tot_NB_sm = tot_NB_sm + tot
               proj_NB_sm = proj_NB_sm + proj
            end if
80       end do
!
!_______ add the band values to form the total contributions
         S_tot = S_tot + tot_NB
         S_proj = S_proj + proj_NB
         S_tot_sm = S_tot_sm + tot_NB_sm
         S_proj_sm = S_proj_sm + proj_NB_sm
90    end do
!
!____ total contribution
!
      S_tot = S_tot * factor
      S_tot_sm = S_tot_sm * factor
      S_proj = S_proj * factor
      S_proj_sm = S_proj_sm * factor
!
!________ estimate the number of "electrons" (per cell)
!
      el_num = el_num + 0.5 * ( S_tot1 + S_tot ) * E_step
      avr_en = avr_en + 0.5 * ( S_tot1*(energy-E_step) + &
                                           S_tot*energy ) * E_step
      S_tot1 = S_tot
!
!........ write S_tot, S_proj, S_tot_sm and S_proj_sm to the file
!
      if(iFlag.eq.0) then
         if(Smeared) then
            write(NIT,'(f10.5,2(x,e12.6))') energy, pos_num(S_tot), &
                                                    pos_num(S_tot_sm)
         else
            write(NIT,'(f10.5,2(x,e12.6))') energy, pos_num(S_tot)
         end if
      else
         if(Smeared) then
            write(NIT,'(f10.5,4(x,e16.6))') &
                 energy,pos_num(S_tot),pos_num(S_tot_sm),pos_num(S_proj),&
                                                         pos_num(S_proj_sm)
         else
            write(NIT,'(f10.5,2(x,e12.6))') energy, pos_num(S_tot), &
                                                    pos_num(S_proj)
         end if
      end if
!
!....... repeat the cycle for the next energy untill B_max(iPh)
!
      if(energy+E_step.le.B_high) go to 50
!
!....... finish with this band
!
      if(iFlag.eq.0) then
         if(Smeared) then
            write(NIT,'(f10.5,a)') B_high+tiny,' 0.0 0.0'
         else
            write(NIT,'(f10.5,a)') B_high+tiny,' 0.0'
         end if
      else
         if(Smeared) then
            write(NIT,'(f10.5,a)') B_high+tiny,' 0.0 0.0 0.0 0.0'
         else
            write(NIT,'(f10.5,a)') B_high+tiny,' 0.0 0.0'
         end if
      end if
      write(*,'(a14,i3,a12,f10.5)') &
           '.... The band ',iPh,' ends at E= ',B_high+tiny
      write(*,'(a29,f10.5)') '_____> Number of electrons = ',el_num
      if(nint(el_num).eq.0) then
         write(*,'(a39)') '_____> Average energy = ** undefined **'
      else
         write(*,'(a24,f10.5,a3)') &
              '_____> Average energy = ',avr_en/nint(el_num),' eV'
      end if
100 end do
   if(iFlag.eq.1)  deallocate(psi3)
end subroutine do_dos

!==========================================================================

subroutine smear_fast(filen1,lenght1,filen,lenght, &
                  Num_B0,B0_min,B0_max,N0_min,N0_max, &
                  Nbands0,Nbands1,iFlag,step,iErr)
!.....................................................................
! iFlag = 0 - only total DOS needs to be calculated;
! iFlag = 1 - both total and projected DOS need to be calculated.
!.....................................................................
! filen1(1:lenght1) - file with the input (raw) DOS
! filen(1:lenght)//'_' - file with the output (smeared) DOS
! step - energy step for the final (smeared) DOS
! Dispers - dispersion for the Gaussian used for smearing
! [Emin,Emax] - gives the energy interval
!.....................................................................
! WARNING! Works correctly only for complete physical bands enclosed
! in [Nbands0,Nbands1]!
!.....................................................................
use param
use kpoints
use dos
implicit none
real*8, parameter :: pi=3.1415927,tiny=0.000001
character line*80,filen*12,filen1*12,cha*12
real*8 Emin,Emax,s1,s2,en,tot,proj,s_norm,sl_norm,pos_num
real*8 tot_Ph,proj_Ph,ee,expd,factor
real*8,dimension(:),allocatable :: e,dos0,dosl
real*8,dimension(:),allocatable ::  Wg_tot,Wg_proj
real*8 B0_min(NBANDS),B0_max(NBANDS),step,W,a,s
integer N0_min(NBANDS),N0_max(NBANDS)
integer iErr,Num_B0,lenght,iFlag,iPh,iPh_min,iPh_max,NB0,NB1,NB,NKP,i,j
integer Npoints,n,lenght1,Nbands0,Nbands1,it
logical alloc
   iErr=0
   alloc=.false.
!
!............ multiplication factor if ISPIN=2
   if(ispin.eq.1) then
      factor=2.0
   else
      factor=1.0
   end if
!
!............ calculate prefactors Wg_tot and Wg_proj for every physical
!             band iPh
!
   allocate(Wg_tot(NBANDS))
   allocate(Wg_proj(NBANDS))
   write(*,'(a)')'>>>>>>> Structure of physical bands <<<<<<<'
   do iPh=1,Num_B0
      if(Nbands0.ge.N0_min(iPh) .and. Nbands0.le.N0_max(iPh)) iPh_min=iPh
      if(Nbands1.ge.N0_min(iPh) .and. &
                               Nbands1.le.N0_max(iPh)) iPh_max=iPh
      Wg_tot(iPh)=(N0_max(iPh)-N0_min(iPh)+1)*factor
      if(iFlag.eq.0) then
         write(*,'(a,i2,a,i3,a1,i3,a,e12.6)') &
              'Ph.band= ',iPh,' states= [',N0_min(iPh),',',N0_max(iPh), &
              '] posesses weight(tot)= ',Wg_tot(iPh)
      else
         NB0=Nbands0
         if(N0_min(iPh).gt.NB0) NB0=N0_min(iPh)
         NB1=Nbands1
         if(N0_max(iPh).lt.NB1) NB1=N0_max(iPh)
         W=0.0d0
         do NB=NB0,NB1
            do NKP=1,NKPTS
               s=0.0d0
               do it=1,Ntask
                  if(which_task(it)) s=s+PSI2(NKP,NB,it)
               end do
               W=W+WTKPT(NKP)*s
            end do
         end do
         Wg_proj(iPh)=W*factor
         write(*,'(a,i2,a,i3,a1,i3,a,e12.6,a,e12.6)') &
              'Ph.band= ',iPh,' states= [',N0_min(iPh),',',N0_max(iPh), &
              '], weight(tot)= ',Wg_tot(iPh), &
              ' and weight(proj)=',Wg_proj(iPh)
      end if
   end do
   write(*,'(a,i2,a,i2)') &
           'Active phys. bands are: from ',iPh_min,' to ',iPh_max
!
!...... open data file
!
   open(1,file=filen1(1:lenght1),form='formatted',status='old',err=200)
   read(1,'(a)') line
   i=0
10 i=i+1
   if(iFlag.eq.0) then
      read(1,*,err=201,end=11) a,a
   else
      read(1,*,err=201,end=11) a,a,a
   end if
   go to 10
11 Npoints=i-1
   rewind(1)
   allocate(e(npoints)) ; allocate(dos0(npoints)) ; allocate(dosl(npoints))
   alloc=.true.
!
!.....................................................................
   read(1,'(a)') line
   do i=1,Npoints
      if(iFlag.eq.0) then
         read(1,*,err=201) e(i),dos0(i)
      else
         read(1,*,err=201) e(i),dos0(i),dosl(i)
      end if
   end do
   close(1)
!
!.... set the total number of points and the energy interval
!
   n=i-1
   Emin=e(1)-Dispers*Broad_Band ; Emax=e(n)+Dispers*Broad_Band
   write(*,*)'WARNING! Total # of data points is ',n
   write(*,'(1x,a28,e12.6,a1,e12.6,a1)') &
                 'WARNING! Energy interval = [',Emin,',',Emax,']'
!
!.......... open an output file for the final DOS/LDOS
!
   write(*,'(a)') &
        ' Writing data to the output file '//filen(1:lenght)//'_ ...'
   open(1,file=filen(1:lenght)//'_',form='formatted',status='unknown')
   write(1,'(a)') line
!
!.......... do the smearing
!
   s1=2*Dispers*Dispers ; s2=1./sqrt(s1*pi) ; s1=1./s1
   en=Emin-step
15 en=en+step
   tot=0.0 ; proj=0.0
   do iPh=iPh_min,iPh_max
      s_norm=0.0 ;    sl_norm=0.0 ;    tot_Ph=0.0 ;    proj_Ph=0.0
      do 16 j=1,n
         if(e(j).lt.B0_min(iPh)) go to 16
         if(e(j).gt.B0_max(iPh)) go to 17
         ee=en-e(j) ; expd=dexp(-s1*ee*ee) 
         tot_Ph=tot_Ph+dos0(j)*expd ; s_norm=s_norm+dos0(j)
         if(iFlag.eq.1) then
            proj_Ph=proj_Ph+dosl(j)*expd
            sl_norm=sl_norm+dosl(j)
         end if
16    end do
17    if(s_norm.gt.tiny) tot = tot + Wg_tot(iPh)*tot_Ph/s_norm*s2
      if(iFlag.eq.1.and.sl_norm.gt.tiny) &
           proj = proj + Wg_proj(iPh)*proj_Ph/sl_norm*s2
   end do
   if(iFlag.eq.0) then
      write(1,'(f10.5,2(x,e12.6))') en,pos_num(tot)
   else
      write(1,'(f10.5,2(x,e12.6))') en, pos_num(tot), pos_num(proj)
   end if
   if(en.lt.Emax) go to 15
   close (1)
   write(*,*)' Done!'  ;   go to 400
!
!..... errors
!
 200  write(*,*)'FATAL! Cannot open file '//filen1(1:lenght1)//' !'
      alloc=.false.
      go to 202
 201  write(*,*)'FATAL! The file '//filen1(1:lenght1)//' is wrong!'
 202  iErr=1
!
!............. deallocate
400   deallocate(Wg_tot) ; deallocate(Wg_proj)
      if(alloc) then
         deallocate(e) ; deallocate(dos0) ; deallocate(dosl)
      end if
end subroutine smear_fast

