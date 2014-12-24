    Program PPolyn

    implicit none
    real*8 x(129),y1(129),  pk_ref(129)
    real*8 dy1(129)
    integer*4 i,n1, status
    character*1 c

    c = ','

  n1= 1
    open(unit =1,file='bao_wiggles_powerspectrum.txt',status = 'unknown',action = 'read',  iostat=status)
    If(status == 0) then 
    do while(.true.)
       read(1,*,iostat=status) x(n1), y1(n1), pk_ref(n1)
      IF ( status < 0 ) EXIT
      n1 = n1+1
    enddo	
    end if
    close(unit= 1)
    n1 = n1-1
    call derivation_Mikrom(n1,x,y1,dy1)   ! the derivative of the f_bao with respect to k
    !call derivation(n1,x,y1,dy1)
    open(222,file='deriv_mikrom_dfbao_dk.dat',status='unknown')
    do i=1,n1
          write(*, *) x(i),';', y1(i),dy1(i), ';',  pk_ref(i)
       	  write(222,*) x(i), y1(i),dy1(i), pk_ref(i)	  
     enddo
     close(unit=222)  	  

    end

!    --------------------------------------------------------

    Subroutine derivation(n,x,y,dydx)

    implicit none

    real*8 x(100),y(100),dydx(100)
    real*8 r1,r2,r3,r4,s1,s2
    real*8 A,B,C
    integer*4 i,n

    do i=2,n-1
       s1 = (x(i-1)-x(i))*(x(i)**2-x(i+1)**2)
       s2 = (x(i)-x(i+1))*(x(i-1)**2-x(i)**2)
       r1 = (y(i-1)-y(i))*(x(i)-x(i+1))
       r2 = (y(i)-y(i+1))*(x(i-1)-x(i))
       r3 = (y(i-1)-y(i))*(x(i)**2-x(i+1)**2)
       r4 = (y(i)-y(i+1))*(x(i-1)**2-x(i)**2)
       A = (r1-r2)/(s2-s1)
       B = (r3-r4)/(s1-s2)
       C = y(i) - A*x(i)**2 - B*x(i)
       dydx(i) = 2*A*x(i) + B
       if(i.eq.2) dydx(i-1) = 2*A*x(i-1) + B
       if(i.eq.(n-1)) dydx(i+1) = 2*A*x(i+1) + B
    enddo

    return
    end

!    --------------------------------------------------------

    Subroutine derivation_Mikrom(n,x,y,dydx)

    implicit none

    real*8 x(100),y(100),dydx(100)
    real*8 a1,A,B
    integer*4 i,n

    do i=2,n-1
       a1 = (y(i)-y(i-1))/(x(i)-x(i-1))
       A = ((y(i+1)-y(i))/(x(i+1)-x(i)) - a1)/(x(i+1)-x(i-1))
       B = a1 - A*(x(i-1)+x(i))
       dydx(i) = 2*A*x(i) + B
       if(i.eq.2) dydx(i-1) = 2*A*x(i-1) + B
       if(i.eq.(n-1)) dydx(i+1) = 2*A*x(i+1) + B
    enddo

    return
    end

!    --------------------------------------------------------
