! =============================================================================
!  sumsq.f90 — Evaluación de una función de base hexaédrica de Nédélec
!  (componente xi de H2) sobre una nube de puntos del hexaedro de referencia.
!
!  Fortran NO permite usar caracteres griegos como identificadores, por lo
!  que hay que escribir xi, eta, zeta con nombres ASCII.
!
!  Compilar y ejecutar:
!      gfortran -O2 sumsq.f90 -o sumsq && ./sumsq
! =============================================================================
program sumsq
    implicit none
    integer, parameter :: n = 200
    integer :: i, j, k
    real(8) :: xi, eta, zeta, total
    real(8) :: t0, t1

    call cpu_time(t0)
    total = 0.0d0
    ! Recorremos una nube n x n x n de puntos del hexaedro de referencia [0,1]^3
    do i = 1, n
        xi = dble(i) / n
        do j = 1, n
            eta = dble(j) / n
            do k = 1, n
                zeta = dble(k) / n
                total = total + xi * eta**2 * zeta**2
            end do
        end do
    end do
    call cpu_time(t1)

    write(*,'(A,I0,A,F0.6,A,F0.2,A)') &
        'Fortran n=', n, '  resultado=', total, &
        '  tiempo=', (t1 - t0) * 1000.0d0, ' ms'
end program sumsq
