PROGRAM LDOS
!===================================================================
! constant current image z(x,y) is calculated
!===================================================================
! It is assumed that the 1st and 2nd lattice vectors lie within the
! surface plane, while the 3rd vector is perpendicular to it
!===================================================================
  implicit none
  integer i,j,NGX,NGY,NGZ,nspins,iz,iy,ix,jX,jY,ix1,iy1,c1,c2
  character*100 :: kappa
  real*16,allocatable,dimension(:,:,:) :: grid
  real*16,allocatable,dimension(:,:)   :: Jmin,Jmax
  real*16 Imin,Imax,Iset,Z0,Z1,Z2,a1,Zset,X,Y,Z,Z00
  type lattice
     real*16,dimension(1:3) :: dir
  end type lattice
  type(lattice) A(1:3)
!
!(1)----- open dialog
!
  write(*,*) 'type the Z-window from the surface to the half vacuum gap: Z0 and Z00'
  read (*,*) Z0,Z00
  write(*,*) 'type the name of file to read' 
  read (*,*) kappa
  open(UNIT=31, Name=trim(kappa))
! 
!(2)----- read the LDOS file
!
!_______________ read the lattice vectors 
!
  do i=1,3
     read(31,*) A(i)%dir(1),A(i)%dir(2),A(i)%dir(3)
     WRITE(*,'(3(f10.5,x))') A(i)%dir(1),A(i)%dir(2),A(i)%dir(3)
  end do
!
!_______________ read the grid
!
  read(31,*) NGX,NGY,NGZ,nspins
  write(*,'(a,3(i5,x),a,i1)') 'grid= ',NGX,NGY,NGZ,' spin=',nspins
!
!_______________ read the LDOS 
!
  allocate (grid(NGX,NGY,NGZ))
  do iz=1,NGZ
     do iy=1,NGY
        read(31,*) (grid(ix,iy,iz),ix=1,NGX)
!        write(*,*)  (grid(ix,iy,iz),ix=1,NGX)
     end do
  end do
  close (31)
  write(*,*) 'type the unit cell replica : c1 and c2'
  read (*,*) c1,c2
! 10 write(*,'(2(a,i5),a)')'plot for which 1 < iX =< ',NGX,' and 1 < iY =< ',NGY,' ?'
!   read(*,*) jX,jY
!   open(1,file='1')
!   do iz=1,NGZ
!      Z=iz*A(3)%dir(3)/NGZ
 !     if(Z > Z0 .and. Z < Z00) then
 !         write(1,'(i5,x,e12.6)')    grid(jx,jy,iz)
 !     endif
 !  end do
 !  close (1)
!    stop
! 
!(3)----- analyse LDOS: min/max currents for each (ix,iy) grid points
!
  allocate (Jmin(NGX,NGY))
  allocate (Jmax(NGX,NGY))
  do ix=1,NGX
     do iy=1,NGY
        Imin=1.e10
        Imax=-1.e10
        do iz=1,NGZ
           Z=iz*A(3)%dir(3)/NGZ
           if(Z > Z0 .and. Z < Z00) then
              if(grid(ix,iy,iz) > Imax ) Imax=grid(ix,iy,iz)
              if(grid(ix,iy,iz) < Imin)  Imin=grid(ix,iy,iz)
!              write(*,*) Imax, Imin
           endif
        end do
        Jmin(ix,iy)=Imin
        Jmax(ix,iy)=Imax
!         write(*,*)       Jmin(ix,iy),Jmax(ix,iy)
     end do
  end do
! 
!(4)----- analyse LDOS: window [Imax:Imin] of currents across 
!                       all (ix,iy) grid points
!
  Imin= 1.e10
  Imax=-1.e10
  do ix=1,NGX
     do iy=1,NGY
        if(Jmin(ix,iy)> Imax) Imax=Jmin(ix,iy)
        if(Jmax(ix,iy)< Imin) Imin=Jmax(ix,iy)
     end do
  end do
! write(*,*) Imax, Imin
  deallocate (Jmin)
  deallocate (Jmax)
!
!(5)----- calculate the image z(ix,iy)
!
1 write(*,'(2(a,f10.5),a)')'Enter current between', Imax, 'and', Imin,' :'
  read(*,*) Iset
!  write(*,*) Iset, Imax, Imin
  if (Iset < Imin .or. Iset > Imax) go to 1

  open(UNIT=32, Name='Zout')

  do ix1=1,c1*NGX
     do iy1=1,c2*NGY
        do iz=NGZ,1,-1
           ix=ix1-(ix1/NGX)*NGX
           iy=iy1-(iy1/NGX)*NGX
           X=ix1*A(1)%dir(1)/NGX+iy1*A(2)%dir(1)/NGY
           Y=ix1*A(1)%dir(2)/NGX+iy1*A(2)%dir(2)/NGY
           Z=iz*A(3)%dir(3)/NGZ
           if(Z > Z0 .and. Z < Z00) then
              if(grid(ix,iy,iz) > Iset) then 
                 if(iz == NGZ) then
                    write(*,*) 'This current is too small!'
                    go to 1
                 endif
                 Z1=Z
                 Z2=Z+A(3)%dir(3)/NGZ
                 a1=(Z1-Z2)/(grid(ix,iy,iz)-grid(ix,iy,iz+1))
                 Zset=a1*(Iset-grid(ix,iy,iz))+Z1
                 exit
              end if
           end if
        end do
        write(32,'(3(f10.5,x))') X,Y,Zset
     end do
     write(32,*)
  end do
  close(32) 
  deallocate(grid)
END PROGRAM LDOS



