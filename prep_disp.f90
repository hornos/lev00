subroutine prep_disp()
!.....................................................................
!      Total DOS with different dispersion of smearing
!.....................................................................
use param
use atoms
use dos
implicit none
real*8, parameter :: pi=3.1415927
real*8,dimension (:,:), allocatable :: E
real*8,dimension (:), allocatable :: B_min,B_max,B_step,B0_min,B0_max,B_step1
integer,dimension(:),allocatable :: N_min,N_max, N0_min,N0_max
integer Nbands0,Nbands1,ijk,Num_B0,i0,Num_B,iErr,j0
real*8 Dispers2(9),Disper(9),E_step,Up_lim,Up_lim_E,Emin,Emax,factor,v
integer lenght(99),k_pnt_ref(4),notet,i,k,nkp,ngroup,item,jErr,ii,j,n_energy,n_en,j1
logical Missing,Smeared_,Yes_Fast,Err,k_ref_alloc
character cha1*1,cha2*2,filen(99)*12,cha,line*16
character Title*50,title_pl(99)*7
character Yes_Job,IndivPrev,IndivPstsc,Smeared
data E_step/0.1/,IndivPrev/'N'/,IndivPstsc/'N'/
data ijk/1/,Yes_Job/'N'/,Smeared/'F'/, &
     Yes_Fast/.false./,Up_lim/0.0/,Up_lim_E/0.0/
data Title/'                                                  '/
   iErr=0
!
!.....................................................................
!...................... PREPARATION ..................................
!.....................................................................
!
!__________ read energies from the output file band.out
!
    allocate (E(NKPTS,NBANDS))
    call get_energies(E,Emin,Emax,jspin,Err)
    if(.not.Err) then
       write(*,'(a23,e12.6,a5,e12.6,a3)')'Energies E are between ', &
                                               Emin,' and ',Emax,' eV'
    else
       deallocate(E)
       return
    end if

!.... allocate memory
    allocate (B_step(NBANDS)) ; allocate (B_step1(NBANDS)) ; allocate (B0_min(NBANDS))
    allocate (B0_max(NBANDS)) ; allocate (B_min(NBANDS)) ;   allocate (B_max(NBANDS))
    allocate (N_min(NBANDS)) ;  allocate (N_max(NBANDS)) ;   allocate (N0_min(NBANDS))
    allocate (N0_max(NBANDS))
    k_ref_alloc=.false.
!
!_________ and  analyze the structure of physical bands:
!      {B0_min(i),B0_max(i)}, i=1,...,Num_B0
!      {N0_min(i),N0_max(i)}, i=1,...,Num_B0 <--- states
!
    call phys_band(E,B0_min,B0_max,Num_B0)
    call phys_band_nmb(E,B0_min,B0_max,Num_B0,N0_min,N0_max)
!
!_________ a loop over all the tetrahedra stored in the file 'brill.dat'
!  notetra - is used to count the tetrahedra
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
30  write(*,'(a29,i5)')'Total number of tetragedra is ',Ntetra

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
!
!________ calculate a factor = 2* (UC volume) / (2*pi)**3 * Rotation
    factor = 2*VOLC/(2*pi)**3
!
!________ update the factor if calculation is spin-polarised
    if(ispin.eq.2) factor=factor*0.5
!
!______ volume of any small cube (in the reciprocal space)
!
    v = (2*pi)**3/VOLC / (Ntetra/6.0)
    write(*,'(a22,f10.6)')'Small cube volume v = ',v
    factor = factor*v
    write(*,'(a22,f10.6)')'Final value of factor = ',factor
!
!______ define the complete range of bands
!
    Nbands0=1
    Nbands1=NBANDS
!.....................................................................
!............. START MAIN MENU HERE ..................................
!.....................................................................
!  ijk - counts different cycles of calculations (no more than 9);
!
    ngroup=0
1   write(*,*)' '

    write(*,*)'..............MENU DOS - check SMEARING ...........'
    write(*,*)'..... Change these parameters if necessary:...'
    write(*,*)
    write(*,*)'  ______________> DOS + Smearing <_______________'

    write(*,'(a30,i2)')'   0. THE CURRENT PICTURE IS: ',ijk
    if(ijk.lt.10) then
       write(cha,'(i1)') ijk
    else
       write(cha2,'(i2)') ijk
    end if
    
    if(E_step.eq.0.0) then
       write(*,'(a)') &
            '   1. Maximum energy step for the plotting: ... undefined ...'
    else
       write(*,'(a44,f10.5)') &
            '   1. Maximum energy step for the plotting: ',E_step
    end if

    if(Smeared.eq.'S') then
       write(*,'(a)') '   2. Smeared DOS/LDOS: ENABLED slow method'
       Smeared_=.true.
    else
       write(*,'(a)') '   2. Smeared DOS/LDOS: ENABLED fast method'
       Smeared_=.false.
    end if

    write(*,'(a33,f6.3,a12)')'   3. Broaden band boundaries by ',Broad_Band,' dispersions'

    if(ngroup.eq.0) then
       write(*,'(a)') '   4. Number of dispersions to be examined: ... undefined ...'
    else
       write(*,'(a43,i5)') '   4. Number of dispersions to be examined: ',ngroup
       write(*,'(a)') '     > Dispersions are:'
       if(ngroup.le.5) then
          write(*,'(2(10x,5(f10.5,1x)))') (Disper(i),i=1,ngroup)
       else
          write(*,'(2(10x,5(f10.5,1x)/))') (Disper(i),i=1,ngroup)
       end if
       do i=1,ngroup
          Dispers2(i)=Disper(i)*sqrt(2.0)
          write(cha1,'(i1)') i
          if(ijk.lt.10) then
             filen(i)='dos.dsp'//cha1//'_'//cha
             lenght(i)=10
          else
             filen(i)='dos.dsp'//cha1//'_'//cha2
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
       write(*,'(a)') '   7. Adopt present setting and calculate DOS'
    else
       write(*,'(a)') '   7. Adopt present setting and calculate DOS <-- DONE!'
       if(Smeared.eq.'F'.and.Yes_Fast) then
          write(*,'(a)') '  88. Adopt present setting and smear the DOS/LDOS <-- DONE!'
       else if(Smeared.eq.'F'.and. (.not.Yes_Fast)) then
          write(*,'(a)') '  88. Adopt present setting and smear the DOS/LDOS'
       end if
    end if

    if(ngroup.gt.1) then
       if(IndivPrev.eq.'Y') then
          write(*,'(a)') '   8. Do you want to preview each individual DOS: YES'
       else
          write(*,'(a)') '   8. Do you want to preview each individual DOS: NO'
       end if
       if(IndivPstsc.eq.'Y') then
          write(*,'(a)') '   9. To make PostScript for each individual DOS: YES'
       else
          write(*,'(a)') '   9. To make PostScript for each individual DOS: NO'
       end if
       write(*,'(a)') '  10. Preview ALL DOS curves for the current group'
       write(*,'(a)') '  11. Create a PostScript file dos.dat0'// &
            filen(1)(9:lenght(1))//' for ALL DOS curves'
    else
       write(*,'(a)')'  10. Preview the current DOS'
       write(*,'(a)')'  11. PostScript file for the current DOS'
    end if

    if(Up_lim.eq.0.0) then
       write(*,'(a)') '  12. Upper limit to the dos-range of the plot: DISABLED'
    else
       write(*,'(a46,f10.5)') '  12. Upper limit to the dos-range of the plot: ',Up_lim
    end if
    if(Up_lim_E.eq.0.0) then
       write(*,'(a)') '  13. Upper limit to the energy-range of the plot: DISABLED'
    else
       write(*,'(a46,f10.5)') &
            '  13. Upper limit to the energy-range of the plot: ',Up_lim_E
    end if

    write(*,'(a)')'  14. Quit: do not proceed.'
    write(*,*)
    write(*,*)'------> Choose the item and press ENTER:'
    read(*,*,err=101) item
!
!............ start a new job ijk
!
    IF(item.eq.0) THEN
       if(ijk.eq.99) then
          write(*,*)'PREP_DOS: You cannot increment the job number!'
          write(*,*)'PREP_DOS: Press ENTER!'
          read(*,*)
       else
          ijk=ijk+1
          ngroup=0
          Yes_Job='N'
          Yes_Fast=.false.
       end if
!
!............ specify maximal allowed energy step E_step
!
    ELSE IF(item.eq.1) THEN
31     write(*,*)'Enter max(STEP) on the energy scale:'
       read(*,*,err=31) E_step
       Yes_Job='N'
       Yes_Fast=.false.
!
!............ specify the method for the smearing: fast or slow
!
    ELSE IF(item.eq.2) THEN
       if(Smeared.eq.'F') then
          Smeared='S'
          Yes_Job='N'
          Yes_Fast=.false.
       else
          Smeared='F'
          Yes_Fast=.false.
          if(Yes_Job.eq.'Y') Yes_Fast=.true.
       end if
!
!............ Broaden band boundaries by a certain number of dispersions
!
    ELSE IF(item.eq.3) THEN
37     write(*,*)'Enter broadening factor (in terms of dispersion):'
       read(*,*,err=37) Broad_Band
       if(Smeared.eq.'S') Yes_Job='N'
       Yes_Fast=.false.
!
!............ specify the dispersion for the Gaussian smearing
!
    ELSE IF(item.eq.4) THEN
33     write(*,*) 'Specify the number of dispersions (not more than 9):'
       read(*,*,err=33) ngroup
       if(ngroup.gt.9) go to 33
32     write(*,'(a12,i2,a39)')'Specify all ',ngroup, &
                        ' dispersions for the Gaussian smearing:'
       read(*,*,err=32) (Disper(i),i=1,ngroup)
       if(Smeared.eq.'S') Yes_Job='N'
       Yes_Fast=.false.
!
!............ calculate the DOS
!............ General loop over groups of bands
!
    ELSE if(item.eq.7) THEN
       if(E_step.eq.0.0 .or. ngroup.eq.0) then
          write(*,*)'ERROR! You still have undefined parameters!'
          go to 1
       end if
       jErr=0
       do i=1,ngroup
          write(*,'(a18,i3,a13)')'........< Dispers ',i,' >.........'
          write(title_pl(i),'(a2,f5.3)') 's=',Disper(i)
          write(*,'(a)') title_pl(i)
!
!___________  analyze the structure of physical bands:
!      {B_min(i),B_max(i)}, i=1,...,Num_B
!
          if(Smeared.eq.'F') then
             if(i.gt.1) go to 50
             Num_B=Num_B0
             do ii=1,Num_B
                B_min(ii)=B0_min(ii)
                B_max(ii)=B0_max(ii)
             end do
          else
             Dispers=Disper(i)
             call phys_band_sm(B_min,B_max,Num_B,B0_min,B0_max,Num_B0)
          end if
          call phys_band_nmb(E,B_min,B_max,Num_B,N_min,N_max)
          write(*,*)' '
          write(*,*)'_______> Structure of "physical" bands <________'
          do j=1,Num_B
             write(*,'(a6,i5,a25,f10.5,a3,f10.5,a14,i5,a1,i5,a1)') &
                  ' Band ',j,' spans energy interval: [',B_min(j),' , ', &
                  B_max(j),'] and states {',N_min(j),',',N_max(j),'}'
          end do
!
!__________ set steps for each physical band
!
          n_energy=0
          do j=1,Num_B
             n_en=(B_max(j)-B_min(j))/E_step+1
             if(n_en.lt.10) then
                n_en=10
                B_step(j)=(B_max(j)-B_min(j))/n_en
             else
                B_step(j)=E_step
             end if
             n_energy=n_energy+n_en
          end do
          write(*,'(a42,i5)') &
               '     > Number of points for the DOS plot: ',n_energy
          write(*,'(a)') '     > Steps for each "physical" band:'
          do j=1,Num_B,6
             j1=j+5
             if(j1.gt.Num_B) j1=Num_B
             write(*,'(7x,6(f8.5,1x))') (B_step(j0),j0=j,j1)
          end do
!
!__________ we call a routine to check whether some original bands
!           may be missing
!
          if(Smeared.eq.'S') then
             Missing=.false.
             call step_ph_band(B_min,B_max,Num_B,E_step,B_step1, &
                  B0_min,B0_max,Num_B0)
             do j=1,Num_B
                if(B_step(j) .gt. B_step1(j)) Missing=.true.
             end do
             if(Missing) write(*,'(a54)') &
                  '     > WARNING! Some of original bands may be missing!'
             write(*,*)' '
          end if
!
!__________ open the file for the output first
!
          open(31,file=filen(i)(:lenght(i)),status='unknown',form='formatted')
          write(*,*) 'The file '//filen(i)//' has been opened for the DOS.'
          if(Smeared.eq.'F') then
             write(31,'(a)') '#----energy-----Total DOS-'
          else
             write(31,'(a)') '#----energy-----Total DOS-----Smeared DOS-'
          end if
!
!________ calculate the DOS and write it into file UNIT=31 filen(i)
!
          call do_dos(31,0,E,B_step,Num_B,B_min,B_max, &
               Dispers2(i),Smeared_,Nbands0,Nbands1,factor)
          close (31)
          write(*,*)'The file '//filen(i)//' has been created!'
!
!_________IN the case of the fast method: do smearing here
!
50        if(Smeared.eq.'F') then
             call smear_fast1(filen(1),lenght(1),filen(i),lenght(i), &
                  Num_B0,B0_min,B0_max,N0_min,N0_max, &
                  Nbands0,Nbands1,0,E_step,Disper(i),Broad_Band,iErr)
             if(iErr.eq.1) jErr=1
          end if
!
!_________ Plot the DOS for the current group
!
          if(IndivPrev.eq.'Y'.and.ngroup.gt.1) &
               call Plot_sm(1,filen(1),lenght(1), &
               filen(i),lenght(i),Title,title_pl(i), &
               'Energy (eV)         ','DOS (arb.units)     ',&
               'Screen', 33,Smeared,Up_lim,Up_lim_E)
          if(IndivPstsc.eq.'Y'.and.ngroup.gt.1) then
             write(*,*)'Give the title:'
             read(*,'(a)') Title
             call Plot_sm(1,filen(1),lenght(1), &
                  filen(i),lenght(i),Title,title_pl(i), &
                  'Energy (eV)         ','DOS (arb.units)     ',&
                  'Postsc', 33,Smeared,Up_lim,Up_lim_E)
          end if
       end do
       Yes_Job='Y'
       Yes_Fast=.false.
       if(jErr.eq.0) Yes_Fast=.true.
!
!............ run FAST smearing for the current unsmeared DOS/LDOS
!
    ELSE IF(Smeared.eq.'F'.and.Yes_Job.eq.'Y'.and.item.eq.88) THEN
       jErr=0
       do i=1,ngroup
          write(*,'(a18,i3,a13)')'........< Dispers ',i,' >..........'
          write(title_pl(i),'(a2,f5.3)') 's=',Disper(i)
          write(*,'(a)') title_pl(i)
          call smear_fast1(filen(1),lenght(1),filen(i),lenght(i), &
               Num_B0,B0_min,B0_max,N0_min,N0_max, &
               Nbands0,Nbands1,0,E_step,Disper(i),Broad_Band,iErr)
          if(iErr.eq.1) jErr=1
       end do
       Yes_Fast=.false.
       if(jErr.eq.0) Yes_Fast=.true.
!
!............ if to preview each time while doing a group?
!
    ELSE IF(ngroup.gt.1.and.item.eq.8) THEN
       if(IndivPrev.eq.'Y') then
          IndivPrev='N'
       else
          IndivPrev='Y'
       end if
!
!............ if to make a PoastScript file each time while doing a group?
!
    ELSE IF(ngroup.gt.1.and.item.eq.9) THEN
       if(IndivPstsc.eq.'Y') then
          IndivPstsc='N'
       else
          IndivPstsc='Y'
       end if
!
!............ preview DOS/LDOS for the current group
!
    ELSE IF(item.eq.10) THEN
       if(Yes_Job.eq.'Y') then
          call Plot_sm(ngroup,filen(1),lenght(1), &
               filen,lenght,Title,title_pl, &
               'Energy (eV)         ', &
               'DOS (arb.units)     ', 'Screen', 33,Smeared, &
               Up_lim,Up_lim_E)
       else
          write(*,*)'IGNORED! You should run the item 7 first!'
       end if
!
!............ make a PostScript file for the current group
!
    ELSE IF(item.eq.11) THEN
       if(Yes_job.eq.'Y') then
          write(*,*)'Give the title:'
          read(*,'(a)') Title
          call Plot_sm(ngroup,filen(1),lenght(1), &
               filen,lenght,Title,title_pl, &
               'Energy (eV)         ', &
               'DOS (arb.units)     ', 'Postsc', 33,Smeared, &
               Up_lim,Up_lim_E)
       else
          write(*,*)'IGNORED! You should run the item 7 first!'
       end if
!
!............ assign some upper limit on the y-axis for the plot
!
    ELSE IF(item.eq.12) THEN
78     write(*,*)'Specify the upper limit for y-axis on your plot:'
       read(*,*,err=78) Up_lim
!
!............ assign some upper limit on the x-axis for the plot
!
    ELSE IF(item.eq.13) THEN
79     write(*,*)'Specify the upper limit for x-axis on your plot:'
       read(*,*,err=79) Up_lim_E
!
!............ quit: return to the main routine
!
    ELSE IF(item.eq.14) THEN
       go to 200
    ELSE
       go to 101
    END IF
    go to 1
 !............. error
101 write(*,*)'Incorrect item number! Try again!'
    go to 1
!
!................................................................ 
!.... ........ END OF THE MAIN DOS MENU .........................
!................................................................
!
!............... errors!
 100   write(*,*)'FATAL: error while open <brill.dat> file!'
       go to 200
 120   write(*,*)'FATAL: error while reading <brill.dat> file!'
       go to 200
 140   write(*,*)'FATAL: internal inconsistency in <brill.dat> file!'
       go to 200
!
!................... finish
200    deallocate (E) ;      deallocate (B_step) ; deallocate (B_step1) ; 
       deallocate (B0_min) ; deallocate (B0_max) ; deallocate (B_min)
       deallocate (B_max) ;  deallocate (N_min) ;  deallocate (N_max)
       deallocate (N0_min) ; deallocate (N0_max) 
       if(k_ref_alloc) deallocate(k_ref)
end subroutine prep_disp

subroutine smear_fast1(filen1,lenght1,filen,lenght, &
                  Num_B0,B0_min,B0_max,N0_min,N0_max, &
                  Nbands0,Nbands1,iFlag,step,Dispers,Broad_Band,iErr)
!.....................................................................
! this is a cut down version of smear_fast() (without PSI2)
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
implicit none
real*8, parameter :: pi=3.1415927,tiny=0.000001
character line*80,filen*12,filen1*12,cha*12
real*8 Emin,Emax,s1,s2,en,tot,proj,s_norm,sl_norm,pos_num
real*8 tot_Ph,proj_Ph,ee,expd,factor
real*8,dimension(:),allocatable :: e,dos
real*8,dimension(:),allocatable ::  Wg_tot
real*8 B0_min(NBANDS),B0_max(NBANDS),step,Dispers,Broad_Band,W,a
integer N0_min(NBANDS),N0_max(NBANDS)
integer iErr,Num_B0,lenght,iFlag,iPh,iPh_min,iPh_max,NB0,NB1,NB,NKP,i,j
integer Npoints,n,lenght1,Nbands0,Nbands1
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
   write(*,'(a)')'>>>>>>> Structure of physical bands <<<<<<<'
   do iPh=1,Num_B0
      if(Nbands0.ge.N0_min(iPh) .and. Nbands0.le.N0_max(iPh)) iPh_min=iPh
      if(Nbands1.ge.N0_min(iPh) .and. &
                               Nbands1.le.N0_max(iPh)) iPh_max=iPh
      Wg_tot(iPh)=(N0_max(iPh)-N0_min(iPh)+1)*factor
      write(*,'(a,i2,a,i3,a1,i3,a,e12.6)') &
              'Ph.band= ',iPh,' states= [',N0_min(iPh),',',N0_max(iPh), &
              '] posesses weight(tot)= ',Wg_tot(iPh)
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
   read(1,*,err=201,end=11) a,a ; go to 10
11 Npoints=i-1
   rewind(1)
   allocate(e(npoints)) ; allocate(dos(npoints)) 
   alloc=.true.
!
!.....................................................................
   read(1,'(a)') line
   do i=1,Npoints
      read(1,*,err=201) e(i),dos(i)
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
      s_norm=0.0 ; sl_norm=0.0 ; tot_Ph=0.0 ; proj_Ph=0.0
      do 16 j=1,n
         if(e(j).lt.B0_min(iPh)) go to 16
         if(e(j).gt.B0_max(iPh)) go to 17
         ee=en-e(j) ; expd=dexp(-s1*ee*ee)
         tot_Ph=tot_Ph+dos(j)*expd
         s_norm=s_norm+dos(j)
16    end do
17    if(s_norm.gt.tiny) tot = tot + Wg_tot(iPh)*tot_Ph/s_norm*s2
   end do
   write(1,*) en,pos_num(tot)
   if(en.lt.Emax) go to 15
   close (1)
   write(*,*)' Done!'
   go to 400
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
400   deallocate(Wg_tot)
      if(alloc) then
         deallocate(e) ; deallocate(dos)
      end if
      return
end subroutine smear_fast1
