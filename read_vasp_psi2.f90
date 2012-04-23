subroutine get_psi2_from_vasp(Err)
!..................................................................
! VASP: all the information for the projected DOS is read in from
!       PROCAR file or constructed from whatever is available at
!       this stage
!
! Ntask = 4*NIONS
! Point() - positions of spheres from all atomic positions
! Radius() - from RWIGS
! method() = ' from_vasp'
! Flag() - read in from PROCAR
! psi2() - read in from PROCAR
!..................................................................
use param
use atoms
use dos
implicit none
character Line*200,cha(4)
integer LinEnd(100),LinPos(100),isp,it,nat,i,nkp,iErr,NumLin,nkp1,nb,nb1,nat1,k
real*8 s,px,py,pz,dxy,dyz,dz2,dxz,dx2,tot,p,d
logical first,Err
data cha/'t','s','p','d'/
!
   Err=.false.
!........... allocate memory
   Ntask=4*NIONS
   allocate(PSI2(NKPTS,NBANDS,Ntask)) ; allocate(method(Ntask)) 
   allocate(Point(3,Ntask)) ; allocate(Radius(Ntask)) ; allocate(Flag(Ntask))
   allocate(v_atm(Ntask)) ; allocate(v_spec(Ntask)) 

!........... tot,s,p,d jobs for every atom
   CaseDos='Sphere' ; nat=0 ; it=0
   do isp=1,NSPEC
      do k=1,NspN(isp)
         nat=nat+1
         do i=1,4
            it=it+1 ; v_atm(it)=nat ; v_spec(it)=isp
            Point(1:3,it)=TI(1:3,nat) 
            Radius(it)=RWIGS(isp) ; Flag(it)=cha(i)
            method(it) =' from_vasp' ;
         end do
      end do
   end do
!
!........... read in tot,s,p,d-weights of DOS to psi2()
   open(15,file='PROCAR',form='formatted',status='old',err=300)
   first=.true.
15 do nkp=1,NKPTS
      call find_string(' k-point  ',10,line,15,.true.,iErr)
      if(iErr.eq.1) go to 300
      call CutStr(line,NumLin,LinPos,LinEnd,0,0,iErr)
      read(line(LinPos(2):LinEnd(2)),*) nkp1
      if(nkp1.ne.nkp) go to 300
      
      do nb=1,NBANDS
         call find_string('band ',5,line,15,.false.,iErr)
         if(iErr.eq.1) go to 300
         call CutStr(line,NumLin,LinPos,LinEnd,0,0,iErr)
         read(line(LinPos(2):LinEnd(2)),*) nb1
         if(nb1.ne.nb) go to 300
         
         call find_2strings('ion ',4,' s ',3,line,15, &
              .false.,iErr)
         if(iErr.eq.1) go to 300
         call CutStr(line,NumLin,LinPos,LinEnd,0,0,iErr)
         do nat=1,NIONS
            if(NumLin.eq.11) then
               read(15,*) nat1,s,px,py,pz,dxy,dyz,dz2,dxz,dx2,tot
               if(nat.ne.nat1)  go to 300
               psi2(nkp,nb,4*nat-3)=tot
               psi2(nkp,nb,4*nat-2)=s
               psi2(nkp,nb,4*nat-1)=px+py+pz
               psi2(nkp,nb,4*nat  )=dxy+dyz+dz2+dxz+dx2
            else if(NumLin.eq.5) then
               read(15,*) nat1,s,p,d,tot
               if(nat.ne.nat1)  go to 300
               psi2(nkp,nb,4*nat-3)=tot
               psi2(nkp,nb,4*nat-2)=s
               psi2(nkp,nb,4*nat-1)=p
               psi2(nkp,nb,4*nat  )=d
            else
               go to 300
            end if
         end do
      end do
      
   end do
!______________ go one more time if spin down
   if(jspin.eq.2 .and. first) then
      first=.false. ; go to 15
   end if
   close (15)
   write(*,*)'PROCAR read in correctly!'
   return
!
!........... errors
300 write(*,*)'FATAL! The file PROCAR is bad or absent!'
    Err=.true.
    deallocate(PSI2) ; deallocate(Point) ; deallocate(Radius)
    deallocate(Flag) ; deallocate(method)
    deallocate(v_atm) ; deallocate(v_spec) 
end subroutine get_psi2_from_vasp


