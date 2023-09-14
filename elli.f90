  PROGRAM ELLI

    IMPLICIT NONE

    INTEGER, PARAMETER :: n = 10
    REAL v(3,n), x(3,n), dt, t, R2, a(3,n), r, c, phi
    REAL m, w, da, dx(3)
    REAL, PARAMETER :: pi = 3.14159
    COMMON/VAR/ v, x, dt, t, r
    COMMON/COFF/ R2, a, c, phi
    COMMON/STRING/ m, da, dx, w
    INTEGER :: st,stsh
    
    REAL sm_device
    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF(sm_device('X11 -bg BLACK -g 1000x1000') .lt. 0) then
       print *, 'Can''t open output device.'
       stop
    ENDIF



    CALL sm_graphics
    CALL sm_defvar('TeX_strings', '1')
    CALL sm_ptype(1.,1)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    stsh = 10000

    CALL  init


    DO  st = 0,10000000000


       if((st/stsh)*stsh.eq.st) CALL show

       CALL step



    ENDDO

    STOP


  END PROGRAM ELLI


  SUBROUTINE init
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: n = 100
    REAL :: v(3,n), x(3,n), dt, t, R2, a(3,n), r, c, phi
    REAL :: m, w, da, dx(3)
    REAL, PARAMETER :: pi = 3.14159
    COMMON/VAR/ v, x, dt, t, r
    COMMON/COEFF/ R2, a, c, phi
    COMMON/STRING/ m, da, dx, w
    INTEGER :: j, k

    WRITE(*,*) 'dt =?'
    READ(*,*) dt


    t = 0




    m = .1

    DO 1 k = 1,n
       DO 1 j = 1,3
          r = 2*(RAND(0))**(1/3)
          c = (-1 + 2*RAND(0))
          phi = 2*pi*RAND(0)
          x(1,k) = r*SQRT(1-C**2)*COS(phi)
          x(2,k) = r*SQRT(1-C**2)*SIN(phi)
          X(3,K) = r*c



          v(j,k) = ((2*M/r)**(1/2))*(rand(0))**(1/2)
1         CONTINUE



          RETURN
       END DO
  END SUBROUTINE intit
     
     
     
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE step
    IMPLICIT NONE
     INTEGER, PARAMETER :: n = 100
     REAL :: v(3,n), x(3,n), dt, t, R2, a(3,n), r, c, phi
     REAL :: m, w, da, dx(3)
     REAL, PARAMETER :: pi = 3.14159
     COMMON/VAR/ v, x, dt, t, r
     COMMON/COEFF/ R2, a, c, phi
     COMMON/STRING/ m, da, dx, w
     INTEGER :: j, k, k1, k2


     t+dt
     DO 1 k = 1, n
        DO 1 j = 1, 3           x(j,k) = x(j,k) + dt*v(j,k)
1          CONTINUE
           DO k = 1, n
              R2 = (x(1,k))**2+(x(2,k))**2+(x(3,k))**2
              w = 15/(R2*SQRT(R2))
              DO j = 1, 3
                 a(j,k) = -w*x(j,k)
              ENDDO
           ENDDO

           DO 2 k1 = 1, n-1
              DO 2 k2 = k1 + 1, n
                 DO j = 1, 3
                    dx(j) = x(j,k1) - x(j,k2)

                 ENDDO

                 R2 = dx(1)**2 + dx(2)**2 + dx(3)**2
                 w = m/(R2*SQRT(R2))
                 DO j = 1, 3
                    da = - w*dx(j)
                    a(j,k1) = a(j,k1) + da
                    a(j,k2) = a(j,k2) - da
                 ENDDO
2                CONTINUE

                 DO 3 k = 1,n
                    DO 3 j = 1,3
                       v(j,k) = v(j,k) + dt*a(j,k)
3                      CONTINUE



                       RETURN
      
                    
                    
                   
  END SUBROUTINE step
                     
                  
  SUBROUTINE show
     INTEGER, PARAMETER :: n = 100
     REAL :: v(3,n), x(3,n), dt, t, R2, a(3,n), r, c, phi
     REAL :: m, w, da, dx(3)
     REAL, PARAMETER :: pi = 3.14159
     COMMON/VAR/ v, x, dt, t, r
     COMMON/COEFF/ R2, a, c, phi
     COMMON/STRING/ m, da, dx, w
     INTEGER ::  k

     CALL sm_limits(-3.,3.,-3.,3.)
     CALL sm_ctype('white')
     CALL sm_lweight(2.)
     CALL sm_ctype('yellow')
     CALL sm_gflush
     DO k = 1,n
        CALL sm_relocate(x(1,k),x(2,k))
        CALL sm_dot
     ENDDO
     CALL sm_gflush





     WRITE(*,*) 't,r=',t,x(1,2)
   END SUBROUTINE show
 
