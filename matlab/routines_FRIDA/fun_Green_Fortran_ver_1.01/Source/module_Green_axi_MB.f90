module module_Green_axi_MB

contains



! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine sub_Green_Aphi_serial(npt_source,r_source,z_source,I_source, & 
                                     npt_point,r_point,z_point,vec_Aphi)


    implicit none

    integer npt_source, npt_point, ii, jj

    real(kind=8) r_source_jj, z_source_jj, I_source_jj, r_point_ii
    real(kind=8) z_point_ii, vec_Aphi_ii
    real(kind=8) num, den, kk_square, kk, res_aphi, pi
    real(kind=8) J1, J2

    real(kind=8) r_source(npt_source), z_source(npt_source)
    real(kind=8) I_source(npt_source)
    real(kind=8) r_point(npt_point), z_point(npt_point)

    real(kind=8) vec_Aphi(npt_point)

    !     cccccccccccccccccccccccccccccccc

    pi = 3.14159265358979323d0

    do ii=1, npt_point

        r_point_ii = r_point(ii)
        z_point_ii = z_point(ii)

        vec_Aphi_ii = 0.d0

        do jj=1, npt_source

            r_source_jj = r_source(jj)
            z_source_jj = z_source(jj)
            I_source_jj = I_source(jj)

            if (r_point_ii.eq.0) then

                res_aphi = 0

            else

                num = 4*r_source_jj*r_point_ii
                den = (r_source_jj+r_point_ii)**2 + (z_point_ii-z_source_jj)**2
                kk_square=(num/den)
                kk = sqrt(kk_square)

                call ellipke(kk_square,J1,J2)

                res_aphi =4.d-7/kk*sqrt(r_source_jj/r_point_ii)*((1-kk**2/2)*J1-J2);

            end if

            vec_Aphi_ii = vec_Aphi_ii + res_aphi*I_source_jj

        end do


        vec_Aphi(ii) = vec_Aphi_ii

    end do


    return
    end subroutine sub_Green_Aphi_serial



! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine sub_Green_Aphi_parallel(npt_source,r_source,z_source,I_source, & 
                                     npt_point,r_point,z_point,n_thread,vec_Aphi)


    implicit none

    integer npt_source, npt_point, ii, jj, n_thread

    real(kind=8) r_source_jj, z_source_jj, I_source_jj, r_point_ii
    real(kind=8) z_point_ii, vec_Aphi_ii
    real(kind=8) num, den, kk_square, kk, res_aphi, pi
    real(kind=8) J1, J2

    real(kind=8) r_source(npt_source), z_source(npt_source)
    real(kind=8) I_source(npt_source)
    real(kind=8) r_point(npt_point), z_point(npt_point)

    real(kind=8) vec_Aphi(npt_point)

    !     cccccccccccccccccccccccccccccccc

    pi = 3.14159265358979323d0

    call omp_set_num_threads(n_thread)

    !$OMP PARALLEL SHARED(npt_point,r_point,z_point,r_source,z_source,I_source,vec_Aphi)
    !$OMP DO PRIVATE(ii,r_point_ii,z_point_ii,vec_Aphi_ii, &
        r_source_jj,z_source_jj,I_source_jj,res_aphi,num,den,kk_square, &
        kk,J1,J2)
    do ii=1, npt_point

        r_point_ii = r_point(ii)
        z_point_ii = z_point(ii)

        vec_Aphi_ii = 0.d0

        do jj=1, npt_source

            r_source_jj = r_source(jj)
            z_source_jj = z_source(jj)
            I_source_jj = I_source(jj)

            if (r_point_ii.eq.0) then

                res_aphi = 0

            else

                num = 4*r_source_jj*r_point_ii
                den = (r_source_jj+r_point_ii)**2 + (z_point_ii-z_source_jj)**2
                kk_square=(num/den)
                kk = sqrt(kk_square)

                call ellipke(kk_square,J1,J2)

                res_aphi =4.d-7/kk*sqrt(r_source_jj/r_point_ii)*((1-kk**2/2)*J1-J2);

            end if

            vec_Aphi_ii = vec_Aphi_ii + res_aphi*I_source_jj

        end do


        vec_Aphi(ii) = vec_Aphi_ii

    end do
    !$OMP END DO NOWAIT 
    !$OMP END PARALLEL 

    return
    end subroutine sub_Green_Aphi_parallel



! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine sub_Green_flux_serial(npt_source,r_source,z_source,I_source, & 
                                     npt_point,r_point,z_point,vec_flux)


    implicit none

    integer npt_source, npt_point, ii, jj

    real(kind=8) r_source_jj, z_source_jj, I_source_jj, r_point_ii
    real(kind=8) z_point_ii, vec_flux_ii
    real(kind=8) num, den, kk_square, kk, res_psi, pi
    real(kind=8) J1, J2

    real(kind=8) r_source(npt_source), z_source(npt_source)
    real(kind=8) I_source(npt_source)
    real(kind=8) r_point(npt_point), z_point(npt_point)

    real(kind=8) vec_flux(npt_point)

    !     cccccccccccccccccccccccccccccccc

    pi = 3.14159265358979323d0

    do ii=1, npt_point

        r_point_ii = r_point(ii)
        z_point_ii = z_point(ii)

        vec_flux_ii = 0.d0

        do jj=1, npt_source

            r_source_jj = r_source(jj)
            z_source_jj = z_source(jj)
            I_source_jj = I_source(jj)

            if (r_point_ii.eq.0) then

                res_psi = 0

            else

                num = 4*r_source_jj*r_point_ii
                den = (r_source_jj+r_point_ii)**2 + (z_point_ii-z_source_jj)**2
                kk_square=(num/den)
                kk = sqrt(kk_square)

                call ellipke(kk_square,J1,J2)

                res_psi =4.d-7/kk*sqrt(r_source_jj/r_point_ii)*((1-kk**2/2)*J1-J2)*(2*pi*r_point_ii);

            end if

            vec_flux_ii = vec_flux_ii + res_psi*I_source_jj

        end do


        vec_flux(ii) = vec_flux_ii

    end do


    return
    end subroutine sub_Green_flux_serial



! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine sub_Green_flux_parallel(npt_source,r_source,z_source,I_source, & 
                                     npt_point,r_point,z_point,n_thread,vec_flux)


    implicit none

    integer npt_source, npt_point, ii, jj, n_thread

    real(kind=8) r_source_jj, z_source_jj, I_source_jj, r_point_ii
    real(kind=8) z_point_ii, vec_flux_ii
    real(kind=8) num, den, kk_square, kk, res_psi, pi
    real(kind=8) J1, J2

    real(kind=8) r_source(npt_source), z_source(npt_source)
    real(kind=8) I_source(npt_source)
    real(kind=8) r_point(npt_point), z_point(npt_point)

    real(kind=8) vec_flux(npt_point)

    !     cccccccccccccccccccccccccccccccc

    pi = 3.14159265358979323d0

    call omp_set_num_threads(n_thread)

    !$OMP PARALLEL SHARED(npt_point,r_point,z_point,r_source,z_source,I_source,vec_flux)
    !$OMP DO PRIVATE(ii,r_point_ii,z_point_ii,vec_flux_ii, &
        r_source_jj,z_source_jj,I_source_jj,res_psi,num,den,kk_square, &
        kk,J1,J2)
    do ii=1, npt_point

        r_point_ii = r_point(ii)
        z_point_ii = z_point(ii)

        vec_flux_ii = 0.d0

        do jj=1, npt_source

            r_source_jj = r_source(jj)
            z_source_jj = z_source(jj)
            I_source_jj = I_source(jj)

            if (r_point_ii.eq.0) then

                res_psi = 0

            else

                num = 4*r_source_jj*r_point_ii
                den = (r_source_jj+r_point_ii)**2 + (z_point_ii-z_source_jj)**2
                kk_square=(num/den)
                kk = sqrt(kk_square)

                call ellipke(kk_square,J1,J2)

                res_psi =4.d-7/kk*sqrt(r_source_jj/r_point_ii)*((1-kk**2/2)*J1-J2)*(2*pi*r_point_ii);

            end if

            vec_flux_ii = vec_flux_ii + res_psi*I_source_jj

        end do


        vec_flux(ii) = vec_flux_ii

    end do
    !$OMP END DO NOWAIT 
    !$OMP END PARALLEL 

    return
    end subroutine sub_Green_flux_parallel



! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine sub_Green_BrBz_serial(npt_source,r_source,z_source,I_source, & 
                                     npt_point,r_point,z_point,vec_Br,vec_Bz)


    implicit none

    integer npt_source, npt_point, ii, jj

    real(kind=8) r_source_jj, z_source_jj, I_source_jj, r_point_ii
    real(kind=8) z_point_ii, vec_Br_ii, vec_Bz_ii
    real(kind=8) num, den, kk_square, kk, res_br, res_bz, pi
    real(kind=8) J1, J2

    real(kind=8) r_source(npt_source), z_source(npt_source)
    real(kind=8) I_source(npt_source)
    real(kind=8) r_point(npt_point), z_point(npt_point)

    real(kind=8) vec_Br(npt_point), vec_Bz(npt_point)

    !     cccccccccccccccccccccccccccccccc

    pi = 3.14159265358979323d0

    do ii=1, npt_point

        r_point_ii = r_point(ii)
        z_point_ii = z_point(ii)

        vec_Br_ii = 0.d0
        vec_Bz_ii = 0.d0

        do jj=1, npt_source

            r_source_jj = r_source(jj)
            z_source_jj = z_source(jj)
            I_source_jj = I_source(jj)

            if (r_point_ii.eq.0) then

                res_br = 0
                res_bz = 4.d-7*pi*r_source_jj**2/(2*(r_source_jj**2 + (z_point_ii - & 
                         z_source_jj)**2)**1.5d0)

            else

                num = 4*r_source_jj*r_point_ii
                den = (r_source_jj + r_point_ii)**2 + (z_point_ii - & 
                      z_source_jj)**2
                kk_square=(num/den)
                kk = sqrt(kk_square)

                call ellipke(kk_square,J1,J2)

                res_br= +1.d-7*kk*(z_point_ii-z_source_jj)/(r_point_ii*sqrt(r_source_jj*r_point_ii))*(& 
                        -J1+(r_source_jj**2+r_point_ii**2+(z_point_ii - & 
                        z_source_jj)**2)/((r_source_jj - r_point_ii)**2+(z_point_ii - &
                        z_source_jj)**2)*J2)

                res_bz= +1.d-7*kk/sqrt(r_source_jj*r_point_ii)*(J1+(r_source_jj**2 - &
                        r_point_ii**2-(z_point_ii-z_source_jj)**2)/((r_source_jj-r_point_ii) &
                        **2+(z_point_ii-z_source_jj)**2)*J2)

            end if

            vec_Br_ii = vec_Br_ii + res_br*I_source_jj
            vec_Bz_ii = vec_Bz_ii + res_bz*I_source_jj

        end do


        vec_Br(ii) = vec_Br_ii
        vec_Bz(ii) = vec_Bz_ii

    end do


    return
    end subroutine sub_Green_BrBz_serial



! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine sub_Green_BrBz_parallel(npt_source,r_source,z_source,I_source, & 
                                     npt_point,r_point,z_point,n_thread,vec_Br,vec_Bz)


    implicit none

    integer npt_source, npt_point, ii, jj, n_thread

    real(kind=8) r_source_jj, z_source_jj, I_source_jj, r_point_ii
    real(kind=8) z_point_ii, vec_Br_ii, vec_Bz_ii
    real(kind=8) num, den, kk_square, kk, res_br, res_bz, pi
    real(kind=8) J1, J2

    real(kind=8) r_source(npt_source), z_source(npt_source)
    real(kind=8) I_source(npt_source)
    real(kind=8) r_point(npt_point), z_point(npt_point)

    real(kind=8) vec_Br(npt_point), vec_Bz(npt_point)

    !     cccccccccccccccccccccccccccccccc

    pi = 3.14159265358979323d0

    call omp_set_num_threads(n_thread)

    !$OMP PARALLEL SHARED(npt_point,r_point,z_point,r_source,z_source,I_source,vec_Br,vec_Bz)
    !$OMP DO SCHEDULE(DYNAMIC,1) PRIVATE(ii,r_point_ii,z_point_ii,vec_Br_ii,vec_Bz_ii, &
        r_source_jj,z_source_jj,I_source_jj,res_Br,res_Bz,num,den,kk_square, &
        kk,J1,J2)
    do ii=1, npt_point

        r_point_ii = r_point(ii)
        z_point_ii = z_point(ii)

        vec_Br_ii = 0.d0
        vec_Bz_ii = 0.d0

        do jj=1, npt_source

            r_source_jj = r_source(jj)
            z_source_jj = z_source(jj)
            I_source_jj = I_source(jj)

            if (r_point_ii.eq.0) then

                res_br = 0
                res_bz = 4.d-7*pi*r_source_jj**2/(2*(r_source_jj**2 + (z_point_ii - & 
                         z_source_jj)**2)**1.5d0)

            else

                num = 4*r_source_jj*r_point_ii
                den = (r_source_jj + r_point_ii)**2 + (z_point_ii - & 
                      z_source_jj)**2
                kk_square=(num/den)
                kk = sqrt(kk_square)

                call ellipke(kk_square,J1,J2)

                res_br= +1.d-7*kk*(z_point_ii-z_source_jj)/(r_point_ii*sqrt(r_source_jj*r_point_ii))*(& 
                        -J1+(r_source_jj**2+r_point_ii**2+(z_point_ii - & 
                        z_source_jj)**2)/((r_source_jj - r_point_ii)**2+(z_point_ii - &
                        z_source_jj)**2)*J2)

                res_bz= +1.d-7*kk/sqrt(r_source_jj*r_point_ii)*(J1+(r_source_jj**2 - &
                        r_point_ii**2-(z_point_ii-z_source_jj)**2)/((r_source_jj-r_point_ii) &
                        **2+(z_point_ii-z_source_jj)**2)*J2)

            end if

            vec_Br_ii = vec_Br_ii + res_br*I_source_jj
            vec_Bz_ii = vec_Bz_ii + res_bz*I_source_jj

        end do


        vec_Br(ii) = vec_Br_ii
        vec_Bz(ii) = vec_Bz_ii

    end do
    !$OMP END DO NOWAIT 
    !$OMP END PARALLEL 


    return
    end subroutine sub_Green_BrBz_parallel



! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine ellipke(m,k,e)

    implicit none

    real(kind=8) m, e, k
    real(kind=8) pi, a0, b0, c0, s0, mm, tol, iter, itmax
    real(kind=8) a1, b1, c1, i1, w1




    !      print*, "num =", num
    !      print*, "den =", den
    !      print*, "m =", m

    !      pause

    if (m.eq.1) then
        e = 1
        k = 1d+10

    else

        pi = 3.14159265358979323d0

        a0 = 1;
        b0 = sqrt(1-m);
        c0 = 0;
        s0 = m;
        i1 = 0;
        mm = 1d+10;
        itmax = 100
        iter = 0

        tol = 2.2204d-16

        do while ((mm.gt.tol).and.(iter.le.itmax))

            iter = iter+1
            a1 = (a0+b0)/2;
            b1 = sqrt(a0*b0);
            c1 = (a0-b0)/2;
            i1 = i1 + 1;
            w1 = (2**i1)*(c1**2);
            mm = w1;

            s0 = s0 + w1;
            a0 = a1;
            b0 = b1;
            c0 = c1;

        end do

        k = pi/(2*a1);
        e = k*(1-s0/2);

    end if

    return
    end subroutine ellipke

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

end module



