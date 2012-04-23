subroutine hat()
  implicit none
  character line*30

      write(*,*)'.....................................................'
      write(*,*)'.....................................................'
      write(*,*)'.......   Postprocesing tool for VASP/SIESTA   ......'
      write(*,*)'.......   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^   ......'
      write(*,*)'.......   It works with:                       ......'
      write(*,*)'.......   ^^^^^^^^^^^^^                        ......'
      write(*,*)'.......  - total or partial charge density     ......'
      write(*,*)'.......  - total and projected DOS             ......'
      write(*,*)'.......  - electrostatic potential             ......'
      write(*,*)'.....................................................'
      write(*,*)'>>>>>>>>  1st release: 25.07.1994'
      include 'build.hat'
      call system('date > tmp.date')
      open(1,file='tmp.date')
      read(1,'(a)') line
      close (1,status='delete')
      write(*,*)'>>>>>>>>  today      : '//line
      write(*,*)'.....................................................'
      write(*,*)'......      inquiries to L. Kantorovich     .........'
      write(*,*)'......      lev.kantorovitch@kcl.ac.uk      .........'
      write(*,*)'.....................................................'

end subroutine hat



