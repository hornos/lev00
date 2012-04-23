subroutine get_energies(E,Emin,Emax,jspin,Err)
!.....................................................................
!  Gets energies for every k-point from the last iteration (in eV):
!  E(k,i) - the i-th eigenvalue associated with the k-th k-point
!  Emin - the lowest energy;
!  Emax - the higest ebergy;
! jspin =1,2 - spin
!.....................................................................
use param
implicit none
real*8 E(NKPTS,NBANDS),Emin,Emax
integer jspin,NKP,jNB,jNB1,i,nn,NB
character cha
logical Err

   Err=.false.
!
!....... read in energies from the output file file_out.
! WARNING! Lines contained #-characters are ignored!
!
   if(ispin.eq.1) then
      write(*,*)'Reading energies from band.out file ...'
      open (1,file='band.out',status='old',form='formatted',err=15)
   else
      write(cha,'(i1)') jspin
      write(*,*)'Reading energies from band.out.'//cha//' file ...'
      open (1,file='band.out.'//cha,status='old', &
                             form='formatted',err=15)
   end if
   do NKP=1,NKPTS
      do jNB=1,NBANDS,5
         jNB1=jNB+4
         if(jNB1.gt.NBANDS) jNB1=NBANDS
         read(1,*,err=15,end=15) i, (E(i,NB),NB=jNB,jNB1)
         if(i.ne.NKP) go to 15
      end do
   end do
   close (1)

   if(ispin.eq.1) then
      write(*,*)'Energy output file <band.out> was read!'
   else
      write(*,*)'Energy output file <band.out.'//cha//'> was read!'
   end if
!
!_________ find Emin and Emax
   nn=NBANDS
   call minmax(E,Emin,Emax,1,nn)
   return
!
!............. error:
15 write(*,*)"FATAL! Your file <band.out> is bad or absent!"
   Err=.true.
end subroutine get_energies

!===============================================================================

subroutine minmax(E,Emin,Emax,Nbands0,Nbands1)
!.....................................................................
!  Gives energy interval, [Emin,Emax], for the group of bands between
! Nbands0 and Nbands1
!.....................................................................
use param
implicit none
real*8 E(NKPTS,NBANDS),Emin,Emax
integer NB,NKP,Nbands0,Nbands1
!
!_________ find Emin and Emax
Emin=1.0E10
Emax=-Emin
do NB=Nbands0,Nbands1
   do NKP=1,NKPTS
      if(E(NKP,NB).gt.Emax) Emax = E(NKP,NB)
      if(E(NKP,NB).lt.Emin) Emin = E(NKP,NB)
   end do
end do
end subroutine minmax

!===============================================================================

subroutine phys_band(E,B_min,B_max,Num_B)
!.......... analyze the structure of physical bands:.............
!
!                { B_min(i) , B_max(i) }, i=1,...,Num_B
!
!   from the eigenvalues E(k-point,band).
!................................................................
use param
implicit none
real*8 E(NKPTS,NBANDS),B_min(NBANDS),B_max(NBANDS)
real*8,dimension(:),allocatable :: Pmin,Pmax
real*8,parameter :: tiny=0.001
integer NB,N,i,j,NB1,Num_B,NKP
real*8 Tmin,Tmax
!
allocate(Pmin(NBANDS)) ; allocate(Pmax(NBANDS))
!
!____ loop over bands NB first: since all k-points contribute into every
!     band (1,...,NBANDS), we find first boundaries [Pmin,Pmax] of every
!     such a band.
!
      do NB=1,NBANDS
        Pmin(NB)=1.0E10 ; Pmax(NB)=-Pmin(NB)
        do NKP=1,NKPTS
          if(E(NKP,NB).gt.Pmax(NB)) Pmax(NB) = E(NKP,NB)
          if(E(NKP,NB).lt.Pmin(NB)) Pmin(NB) = E(NKP,NB)
        end do
        if(abs(Pmax(NB)-Pmin(NB)).lt.tiny) &
            write(*,'(a,i5)') 'WARNING: zero spread for band ',NB
      end do
!
!....... study overlap of bands: creation of "physical" bands (numbered by N)
!
      N=0
      do 100 NB=1,NBANDS
         Tmin=Pmin(NB) ; Tmax=Pmax(NB)
!___ check if this band has been already included
         if(N.eq.0) go to 10
         do j=1,N
           if(Tmin.ge.B_min(j) .and. Tmax.le.B_max(j)) go to 100
         end do
!___ scan all bands with respect to the current one NB
10       i=0
         do 20 NB1=1,NBANDS
           if(NB.eq.NB1) go to 20
!____________ check if band NB1 has no overlap at all with NB
           if(Pmax(NB1).lt.Tmin) go to 20
           if(Pmin(NB1).gt.Tmax) go to 20
!____________ check if we need to change boundaries of N-th "ph" band
           if(Pmin(NB1).lt.Tmin) then
              Tmin=Pmin(NB1) ; i=i+1
           end if
           if(Pmax(NB1).gt.Tmax) then
              Tmax=Pmax(NB1) ; i=i+1
           end if
20       end do
         if(i.ne.0) go to 10
         N=N+1 ; B_min(N)=Tmin ; B_max(N)=Tmax
100   end do
      Num_B=N
end subroutine phys_band

!===============================================================================

subroutine phys_band_nmb(E,B_min,B_max,Num_B,N_min,N_max)
!.......... ....................................................
!   Find the boundaries (in numbers of states) of physical bands
!
!                { B_min(i) , B_max(i) }, i=1,...,Num_B
!
!  So, for each band, {N_min(i),N_max(i)} gives such a boundary.
!................................................................
use param
implicit none
real*8 E(NKPTS,NBANDS),B_min(NBANDS),B_max(NBANDS)
integer N_min(NBANDS),N_max(NBANDS),i,N1,Num_B,j

!...... it is enough to check just the gamma point
   N1=1
   do 10 i=1,Num_B
      do j=N1,NBANDS
         if(E(1,j).gt.B_max(i)) then
            N_min(i)=N1 ; N_max(i)=j-1 ; N1=j 
            go to 10
         end if
      end do
10 end do
   N_min(Num_B)=N1
   N_max(Num_B)=NBANDS
 end subroutine phys_band_nmb

!===============================================================================

subroutine phys_band_sm(B_min,B_max,Num_B,B0_min,B0_max,Num_B0)
!.......... correct the structure of physical bands:.............
!
!                { B_min(i) , B_max(i) }, i=1,...,Num_B
!
!   taking into account Gaussian smearing (convolution). Each
!   band is widened by the Gaussian dispersion Dispers from both
!   wings.
!................................................................
!  Broad_Band - a number of dispersions to be taken on both sides
!               of any band wings as an artificial broadening
!................................................................
!   { B0_min(i) , B0_max(i) }, i=1,...,Num_B0 -> original bands
!................................................................
use param
use dos
implicit none
real*8  B_min(NBANDS),B_max(NBANDS),B0_min(NBANDS),B0_max(NBANDS)
integer Num_B,Num_B0,i,n_ov,j

!.... shift boundaries of bands
   Num_B=Num_B0
   do i=1,Num_B
      B_min(i)=B0_min(i)-Dispers*Broad_Band
      B_max(i)=B0_max(i)+Dispers*Broad_Band
   end do
!.... correct boundaries and bands in case they overlap
5  if(Num_B.eq.1) return
   n_ov=0
   i=1
10 i=i+1
   if(B_max(i-1).ge.B_min(i)) then
      n_ov=n_ov+1
      
      B_max(i-1)=B_max(i)
      if(B_min(i).lt.B_min(i-1)) B_min(i-1)=B_min(i)
      
      if(i.lt.Num_B) then
         do j=i+1,Num_B
            B_min(j-1)=B_min(j)
            B_max(j-1)=B_max(j)
         end do
      end if
      Num_B=Num_B-1
   end if
   if(i.lt.Num_B) go to 10
   if(n_ov.ne.0) go to 5
end subroutine phys_band_sm

!===============================================================================

subroutine step_ph_band(B_min,B_max,Num_B,E_step,B_step1, &
                                   B0_min,B0_max,Num_B0)
!................................................................
!   It chooses steps, B_step1(physical band), for each phys. band
! in such a way that not to miss any of the actual bands: the step
! is chosen not larger than for the smallest actual band inside the
! interval spanned by the physical band.
!   If the smearing is 'on', then this option allows one not to miss
! any actual narrow peak while showing (previewing) not smeared DOS.
! It is expensive, since the step may become really tiny for some
! physical bands.
!................................................................
!   { B0_min(i) , B0_max(i) }, i=1,...,Num_B0 -> original bands
!   { B_min(i) , B_max(i) }, i=1,...,Num_B    -> physical bands
!   Step - the largest possible step inside any physical band
!................................................................
use param
implicit none
integer n_en,i,Num_B,Num_B0,j
real*8 E_step,step,st
real*8  B_min(NBANDS),B_max(NBANDS),B_step1(NBANDS),B0_min(NBANDS),B0_max(NBANDS)
!
!...... loop over physical bands: we are looking for all actual
! (original) bands matching its interval; then we calculate the step
!
   do i=1,Num_B
!      write(*,*)'----- physical band ',i,'---------'
!      write(*,*)'spans: [',B_min(i),',',B_max(i),'] eV'
!_____ normally, the step would be step:
      n_en=(B_max(i)-B_min(i))/E_step+1
      if(n_en.lt.10) then
         n_en=10
         step=(B_max(i)-B_min(i))/n_en
      else
         step=E_step
      end if
!_____ check however all original bands and choose the smallest step
      do 50 j=1,Num_B0
         if(B0_max(j).lt.B_min(i) .or. B0_min(j).gt.B_max(i)) go to 50
         n_en=(B0_max(j)-B0_min(j))/E_step+1
         if(n_en.lt.10) then
            n_en=10
            st=(B0_max(j)-B0_min(j))/n_en
         else
            st=E_step
         end if
         if(st.lt.step) step=st
50    end do
      B_step1(i)=step
   end do
end subroutine step_ph_band

!===============================================================================

subroutine group_band(Nbands0,Nbands1,ngroup)
!................. Choose the bands. All allowed bands are ..........
!  gathered into groups [Nbands0(i),Nbands1(i)] of bands which
!  are treated individually. 
!....................................................................
use param
implicit none
integer ngroup,i
integer Nbands0(NBANDS),Nbands1(NBANDS)

 16    write(*,*)'How many continuous groups of bands?'
       read(*,*,ERR=17) ngroup
       if(ngroup.gt.8) then
          write(*,*)'GROUP_BAND: Not more than 8 intervals in the DOS!'
          go to 16
       end if
       if(ngroup.ge.0.and.ngroup.le.NBANDS) go to 14
 17    write(*,*)'Error! Try again!'
       go to 16
 14    do 20 i=1,ngroup
 25      write(*,'(a11,i5,a29)') '.......... ',i,'-th group of bands ..........'
         write(*,'(a46,i5,a2)') &
                    'The range of bands in the group between 1 and ',NBANDS,' :'
         read(*,*,ERR=7) Nbands0(i), Nbands1(i)
         if(Nbands0(i).ge.1 .and. Nbands1(i).le.NBANDS .and. &
                                    Nbands0(i).le.Nbands1(i)) goto 20
 7       write(*,*)'Error while giving bands numbers! Try again!'
         go to 25
20    end do

!........... final info
      do i=1,ngroup
         write(*,'(a14,i3,a19,i5,a3,i5)') '... The group ',i, &
                  ' covers bands from ', Nbands0(i),' - ',Nbands1(i)
      end do
end subroutine group_band

!===============================================================================

subroutine sort(E1,E2,E3,E4,i_sort)
!....................................................................
! It sorts out 4 edges of the tetrahedra in increasing order of
!            their energies E0,E1,E2,E3:
! i_sort(new_order_of_edges) - old order of edges of the tetrahedra
! (as in the enrgies E1,E2,E3,E4 in the input).
!....................................................................
!   Energies must be ordered in such a way that:
!  E(new 0)  =< E(new 1) =< E(new 2) =< E(new 3), i.e.
!  E(i_sort(1)) =< E(i_sort(2)) =< E(i_sort(3)) =< E(i_sort(4))
!....................................................................
implicit none
real*8 E(4),E1,E2,E3,E4
integer i_sort(4),k,i,ii
!_________ prepare for sorting procedure
   E(1)=E1
   E(2)=E2
   E(3)=E3
   E(4)=E4
   do i=1,4
      i_sort(i)=i
   end do
!_________ sorting them out; k - counts permutations in every cycle
5  k=0
   do i=1,3
      if( E(i_sort(i)) .gt. E(i_sort(i+1)) ) then
         k=k+1
         ii=i_sort(i+1)
         i_sort(i+1)=i_sort(i)
         i_sort(i)=ii
      end if
   end do
   if(k.ne.0) go to 5
 end subroutine sort

!===============================================================================
subroutine do_contr(tot,proj,A,E,energy,iFlag)
!............calculate the contributions to both total (tot) .........
!   and projected  (proj) DOS for the given 'energy' from the
!   thetrahedra with energies E(0) < E(1) < E(2) < E(3) and with
!   'projected' factors A(0),...,A(3)
!.....................................................................
!   For every case a check is made if there is a pair of equal energies
! in the edges of the tetrahedra. This case happens, for example, when
! z-direction is supressed. In this case some of the cases lead to
! another results (usually, to 0). Note that there can be only one pair
! of equavalent edges because all the tetrahedra are chosen from tiny
! cubes which are oriented in such a way (when z-supression is allowed)
! that their vertical sides are parallel (exactly) to the z-axes.
!   However, if the size of tiny cubes was chosen to be not very small,
! there may be equivalent edges in the cubes by the symmetry. In this
! case more than 2 edges could have the same energy. That is why all
! cases when the enrgies coincide must be considered explicitly.
!.....................................................................
! iFlag = 0 - if only total DOS needs to be calculated;
! iFlag = 1 - if both total and projected DOS need to be calculated.
!.....................................................................
!
implicit none
real*8 E(0:3),A(0:3),tot,proj,energy,E21,E0,E10,E20,E30,A10,A20,A30,E1,E31
real*8 tot0,tot1,A21,A31,proj0,proj1
integer iFlag
real*8,parameter :: tiny=0.00000001

!............ set to 0.0 if the energy is out of range
   tot=0.0 ; proj=0.0
   if(energy.le.E(0).or.energy.ge.E(3)) return
!
!............ the contributions separately for every case .............
!
!_________ E(0) < energy < E(1): these are "f0" and "J0"
!
   IF(energy.le.E(1)) THEN
      
      if( abs(E(0)-E(1)).gt.tiny ) then
         tot = (energy-E(0))**2/((E(1)-E(0))*(E(2)-E(0))*(E(3)-E(0)))
         if(iFlag.eq.1) then
            proj = tot * ( A(0) + (energy-E(0)) * &
                 ( (A(1)-A(0))/(E(1)-E(0)) + (A(2)-A(0))/(E(2)-E(0)) + &
                             (A(3)-A(0))/(E(3)-E(0)) )/3. )
         end if
      end if
      
!
!_________ E(1) < energy < E(2): these are "f0-f1" and "J0-J1"
!
   ELSE IF(energy.le.E(2)) THEN

      E21=E(2)-E(1)
      if( abs(E21).lt.tiny ) return
      E0=energy-E(0)
      E10=E(1)-E(0) ; E20=E(2)-E(0) ; E30=E(3)-E(0) ; A10=A(1)-A(0)
      A20=A(2)-A(0) ; A30=A(3)-A(0)
      if( abs(E10).gt.tiny ) then
         E1=energy-E(1) ; E31=E(3)-E(1)
         tot0 = E0**2/(E10*E20*E30) ; tot1 = E1**2/(E10*E21*E31)
         tot = tot0 - tot1
         if(iFlag.eq.1) then
            A21=A(2)-A(1) ; A31=A(3)-A(1)
            proj0 = tot0*(A(0) + E0*(A10/E10 + A20/E20 + A30/E30)/3.)
            proj1 = tot1*(A(1) + E1*(A10/E10 + A21/E21 + A31/E31)/3.)
            proj = proj0 - proj1
         end if
      else
         tot= E0/(E20*E30) * ( 2.0 - E0*(1./E20 + 1./E30) )
         if(iFlag.eq.1) proj = A(0)*tot + E0*E0/(E20*E30)* &
              (  A20/E20*(1.0-E0/3.*(2./E20 + 1./E30) ) + &
              A30/E30*(1.0-E0/3.*(2./E30 + 1./E20) ) )
      end if
!
!_________ E(2) < energy < E(3): these are "f3" and "J3"
!
   ELSE

      if( abs(E(2)-E(3)).gt.tiny ) then
         tot = (energy-E(3))**2/((E(3)-E(1))*(E(3)-E(2))*(E(3)-E(0)))
         if(iFlag.eq.1) proj = tot * ( A(3) + (energy-E(3)) * &
              ( (A(1)-A(3))/(E(1)-E(3)) + &
              (A(2)-A(3))/(E(2)-E(3)) + &
              (A(0)-A(3))/(E(0)-E(3)) )/3. )
      end if
   END IF
   tot = 0.5*tot
   proj = 0.5*proj
 end subroutine do_contr
!===============================================================================

subroutine do_contr_sm(Disp,tot,proj,A,E,energy,iFlag)
implicit none
!............calculate a smearing contributions to both total (tot)
!   and projected  (proj) DOS for the given 'energy' from the
!   thetrahedra with energies E(0) < E(1) < E(2) < E(3) and with
!   'projected' factors A(0),...,A(3)
!.....................................................................
!   For every case a check is made if there is a pair of equal energies
! in the edges of the tetrahedra. This case happens, for example, when
! z-direction is supressed. In this case some of the cases lead to
! another results (usually, to 0). Note that there can be only one pair
! of equavalent edges because all the tetrahedra are chosen from tiny
! cubes which are oriented in such a way (when z-supression is allowed)
! that their vertical sides are parallel (exactly) to the z-axes.
!   However, if the size of tiny cubes was chosen to be not very small,
! there may be equivalent edges in the cubes by the symmetry. In this
! case more than 2 edges could have the same energy. That is why all
! cases when the enrgies coincide must be considered explicitly.
!.....................................................................
! iFlag = 0 - if only total DOS needs to be calculated;
! iFlag = 1 - if both total and projected DOS need to be calculated.
!.....................................................................
!
real*8 E(0:3),A(0:3),tot,proj,energy,E21,E0,E10,E20,E30,A10,A20,A30,E1,E31
real*8 tot0,tot1,A21,A31,proj0,proj1,B0,Disp,fact,tot2,proj2,fact0,B1,fact1
real*8 El1,El2,tot3,proj3,E23,E13,E03,B3
integer iFlag
real*8,parameter :: tiny=0.00000001
!
!............ do the contributions separately for every case ..........
!
!_________ E(0) < energy < E(1): these are integrals over "f0" and "J0"
!
   tot1=0.0
   proj1=0.0
   E10=E(1)-E(0)
   if( abs(E10).gt.tiny ) then
      E20=E(2)-E(0)
      E30=E(3)-E(0)
      B0=Disp/3.0*((A(1)-A(0))/E10+(A(2)-A(0))/E20+(A(3)-A(0))/E30)
      fact=0.564189583/(E20*E30)*Disp
      call cont_sm(E,A,Disp,energy,iFlag,tot1,proj1,fact,B0,0,0,1)
   end if
!
!_________ E(1) < energy < E(2): these are integrals over "f0-f1" and "J0-J1"
!
   tot2=0.0
   proj2=0.0
   E21=E(2)-E(1)
   if( abs(E21).gt.tiny ) then
      E31=E(3)-E(1)
      if(abs(E(0)-E(1)).gt.tiny ) then
         E10=E(1)-E(0)
         E20=E(2)-E(0)
         E30=E(3)-E(0)
         B0=Disp/3.0*((A(1)-A(0))/E10+(A(2)-A(0))/E20+(A(3)-A(0))/E30)
         fact0=0.564189583*E21/(E10*E20*E30)*Disp
!          call cont_sm(E,A,Disp,energy,iFlag,tota,proja,fact0,B0,0,1,2)
         B1=Disp/3.0*((A(1)-A(0))/E10+(A(2)-A(1))/E21+(A(3)-A(1))/E31)
         fact1=0.564189583/(E10*E31)*Disp
!          call cont_sm(E,A,Disp,energy,iFlag,totb,projb,fact1,B1,1,1,2)
!          tot2=tota-totb
!          proj2=proja-projb
         call cont_sm12b(E,A,Disp,energy,iFlag,tot2,proj2,fact0,B0, &
                                         fact1,B1)
      else
         El1=energy-E(1)
         El2=energy-E(2)
         A21=A(2)-A(1)
         A31=A(3)-A(1)
         call cont_sm12(Disp,iFlag,tot2,proj2,E21,E31,El1,El2, &
                                                          A(1),A21,A31)
      end if
   end if
!
!_________ E(2) < energy < E(3): these are integrals over "f3" and "J3"
!
   tot3=0.0
   proj3=0.0
   E23=E(2)-E(3)
   if( abs(E23).gt.tiny ) then
      E03=E(0)-E(3)
      E13=E(1)-E(3)
      B3=Disp/3.0*((A(0)-A(3))/E03+(A(1)-A(3))/E13+(A(2)-A(3))/E23)
      fact=0.564189583/(E03*E13)*Disp
      call cont_sm(E,A,Disp,energy,iFlag,tot3,proj3,fact,B3,3,2,3)
   end if
!
!_________ sun contributions from these three energy regions
!
   tot = 0.5*(tot1+tot2+tot3)
   proj = 0.5*(proj1+proj2+proj3)
   return
 end subroutine do_contr_sm

!===============================================================================
real*8 function URFC(Y)
!......................................................................
!                 urfc(y)=erfc(y)*exp(y**2),
!  where erfc(y) is the usual error function, erf(y)=1-erfc(y).
!  Thus, if we want erfc(y), we do the following:
!
!         erfc(y) =>  exp(-y*y)*urfc(y)
!
!......................................................................
!  We use here the interpolation given in the Abramovitz and Steagan
!  book for the erfc(y).
!......................................................................
implicit none
real*8 T,Y
  T=1./(1.+0.3275911*Y)
  URFC = T*(0.254829592+T*(-0.284496736+T*(1.421413741+ &
     T*(-1.453152027+T*1.061405429))))
end function URFC

!===============================================================================

subroutine cont_sm12b(E,A,Disp,energy,iFlag,tot,proj,f0,B0,f1,B1)
implicit none
!........................................................................
!  This is:  lambda_0(1,2) - lambda_1(1,2)
!........................................................................
real*8 E(0:3),A(0:3),I0,I1,J0,J1,energy,x0,x1,x2,Disp,tot,proj,f0,B0,f1,B1
real*8 d0,d1,d11,tot1,dd1,d00,dd0
integer iFlag
  x1=(energy-E(1))/Disp
  x2=(energy-E(2))/Disp
  call J0_J1_I0_I1(x1,x2,J0,J1,I0,I1)
  x0=(energy-E(0))/Disp
  tot = I1*(f1-f0) -2*J0*(f1*x1-f0*x0) + (f1*x1*x1-f0*x0*x0)*I0
  if(iFlag.eq.1) then
     d0=f0*A(0)
     d1=f1*A(1)
     tot1 = I1*(d1-d0) -2*J0*(d1*x1-d0*x0) + (d1*x1*x1-d0*x0*x0)*I0
     d0=f0*B0
     d1=f1*B1
     d11=d1*x1
     dd1=d11*x1
     d00=d0*x0
     dd0=d0*x0*x0
     proj = tot1 - J1*(d1-d0) + 3*I1*(d11-d00) - 3*J0*(dd1-dd0) + &
          I0*(dd1*x1 - dd0*x0)
  end if
end subroutine cont_sm12b

!===============================================================================
subroutine cont_sm12(Disp,iFlag,tot,proj,E,E1,El1,El2,AA,A,A1)
implicit none
real*8 I0,I1,J0,J1,M,M1,D,fact,Disp,tot,proj,E,E1,El1,El2,AA,A,A1,x1,x2,ET
integer iFlag

  fact=0.564189583/(Disp*E1)
!
  x1=El1/Disp
  x2=El2/Disp
  call J0_J1_I0_I1(x1,x2,J0,J1,I0,I1)
!
  M=A/E
  M1=A1/E1
  D=M*(2.0/E+1.0/E1)+M1*(1.0/E+2.0/E1)
  M=M+M1
  ET=1.0/E+1.0/E1
!
  tot = -El1*(2.0-El1*ET)*I0+Disp*(2*(1.0-El1*ET)*J0+Disp*ET*I1)
  if(iFlag.eq.1) then
     proj = -El1*El1*I0*(M-El1/3.0*D)+Disp*( El1*J0*(2*M-El1*D) + &
          Disp*(I1*(-M+El1*D)-Disp*J1/3.0*D) )
     proj= (tot*AA + proj)*fact
  end if
  tot = fact*tot
end subroutine cont_sm12
!===============================================================================

subroutine cont_sm(E,A,Disp,energy,iFlag,tot,proj,fact,B,imu,i,j)
implicit none
!........................................................................
!  This is:  lambda_imu(i,j)
!........................................................................
real*8 E(0:3),A(0:3),I0,I1,J0,J1,energy,x1,x2,x3,Disp,tot,proj,B,fact
integer iFlag,i,j,imu

  x1=(energy-E(i))/Disp
  x2=(energy-E(j))/Disp
  call J0_J1_I0_I1(x1,x2,J0,J1,I0,I1)
  x3=(energy-E(imu))/Disp
  tot = -x3*(x3*I0-2*J0)-I1
  if(iFlag.eq.1) then
     proj = (A(imu)*tot - B*(-J1+x3*(3*I1-x3*(3*J0-x3*I0))))*fact
  end if
  tot = fact*tot
end subroutine cont_sm

!===============================================================================
subroutine J0_J1_I0_I1(A,B,J0,J1,I0,I1)
!.......................................................................
!     Functions J0,J1,I0,I1 are calculated for arguments A,B. Special
! care is taken if A is close to B. In this case a series expansion
! is used in every case. - 16.04.97
!.......................................................................
implicit none
real*8 J0,J1,I0,I1,A,B,y,exp1,s,ak,akk,z,Hk,Hkm1,Hk1,exp2,erfc1,urfc,erfc2
real*8 A2,B2,dexp
real*8, parameter :: tiny=0.01,tiny1=0.00000001
integer k

  A2=A*A
  B2=B*B
  exp1=dexp(-A2)
  y=B2-A2
  IF(abs(y).lt.tiny) THEN
!_____________calculating J0
     k=0
     ak=1.0
     s=ak
10   k=k+1
     akk=-y/(k+1)*ak
     s=s+akk
     if(abs(akk-ak).lt.tiny1) go to 11
     ak=akk
     go to 10
11   J0=-0.5*s*exp1*(A+B)
!_____________calculating I0
     z=B-A
     k=1
     s=-1.+A*z
     ak=z*A
     Hk=2*A
     Hkm1=1.
20   Hk1=2*(A*Hk-k*Hkm1)
     k=k+1
     akk=-z/(k+1) * Hk1/Hk * ak
     s=s+akk
     if(abs(akk-ak).lt.tiny1) go to 21
     ak=akk
     Hkm1=Hk
     Hk=Hk1
     go to 20
21   I0=exp1*s
  ELSE
!_____________calculating J0,I0
     z=A-B
     exp2=dexp(-B2)
     J0=0.5*(exp1-exp2)/z
     if(A.ge.0.0) then
        erfc1=exp1*urfc(A)
     else
        erfc1=2.0 - exp1*urfc(-A)
     end if
     if(B.ge.0.0) then
        erfc2=exp2*urfc(B)
     else
        erfc2=2.0 - exp2*urfc(-B)
     end if
     I0=0.886226925*(erfc1-erfc2)/z
  END IF
!_____________calculating J1,I1
  J1=J0*(1.+B2) + 0.5*(A+B)*exp1
  I1=0.5*(I0 +2*B*J0 + exp1)
end subroutine J0_J1_I0_I1

!===============================================================================
