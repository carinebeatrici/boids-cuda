program gerapov
  real,allocatable::x(:),y(:),z(:),vx(:),vy(:)
!  real ::xm,ym,t
  real :: comp
  integer :: nr,nc,n
  character ::a*20,fmtcam*60,fmtesf*60,arquivopov*12
  real::t
  open(10,file="posicoes.dat")
  open(20,file="obs")
  read(20,*)n,nv,nc,comp
  close(20)
nr=nv
!  open(20,file="in")
!  read(20,*)n
  !comp=12.
!  nr=n/4.
!  nc=n-nr
  !read(10,*)nr,nc
  !ladox=12;ladoy=12;ladoz=25
  ladox=comp;ladoy=comp;ladoz=comp
  allocate(x(nr+nc),y(nr+nc),z(nr+nc),vx(nr+nc),vy(nr+nc))
  do lup=1,1000
     !do i=1,7
     read(10,*)t
     !end do
     do i=1,nr
        read(10,*)x(i),y(i),vx(i),vy(i)
!        write(*,*)i,x(i),y(i),vx(i),vy(i)
     end do
!     read(10,*)a
     do i=1+nr,nr+nc
        read(10,*)x(i),y(i),vx(i),vy(i)
!        write(*,*)i,x(i),y(i),vx(i),vy(i)
     end do
     z=2.
!     xm=sum(x)/real(nr+nc)
!     ym=sum(y)/real(nr+nc)
     !write(*,*)xm,ym
!     x=x-xm;y=y-ym
     ! exporta .pov
     write(arquivopov,'(a3,i4.4,a4)')'ray',lup,'.pov'
     !write(*,*)arquivopov
     open(30, file = arquivopov)
     fmtcam = '(a16,f8.3,a1,f8.3,a1,f8.3,a9,f8.3,a1,f8.3,a1,f8.3,a2)'
     write(30,*) '#include "colors.inc"'
     write(30,fmtcam)'camera{location<',ladoy/2.,',',2.*ladoz,',',ladox/2.,'>look_at<',ladoy/2.,',',2.,',',ladox/2.,'>}'
     !write(30,fmtcam)'camera{location<',ladox/2.,',',ladoy/2.,',',2*ladoz,'>look_at<',ladox/2.,',',ladoy/2.,',',ladoz/2.,'>}'
     write(30,'(a14,f8.3,a1,f8.3,a1,f8.3,a13)')'light_source{<',0.,',',4.*ladoz,',',-2.*ladox,'>color White}'
     !write(30,'(a14,f8.3,a1,f8.3,a1,f8.3,a13)')'light_source{<',-2*ladox,',',0.,',',-4.*ladoz,'>color White}'
     write(30,'(a10,f8.3,a35)') 'plane { y,', 0., ' pigment{color rgb<0.2, 0.2, 0.2>}}'
     fmtesf = '(a8,f10.7,a1,f10.7,a1,f10.7,a2,f8.3,a38)'
     !!!  r =  RAIO DE EQUILIBRIO  !!!
     r=0.2
     do i = 1, nr
       !write(30,fmtesf)'sphere{<',x(i),',',y(i),',',z(i),'>,',r,'pigment{color Red}finish{phong 1}}'
       write(30,fmtesf)'sphere{<',x(i),',',z(i),',',y(i),'>,',r,'pigment{color Red}finish{phong 1}}'
     end do
     do i=nr+1,nc+nr
       !write(30,fmtesf)'sphere{<',x(i),',',y(i),',',z(i),'>,',r,'pigment{color Cyan}finish{phong 1}}'
       write(30,fmtesf)'sphere{<',x(i),',',z(i),',',y(i),'>,',r,'pigment{color Cyan}finish{phong 1}}'
     end do
     close(30)
  end do
end program gerapov
