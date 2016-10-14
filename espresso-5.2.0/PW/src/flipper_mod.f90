MODULE flipper_mod

  USE kinds, ONLY: DP

  IMPLICIT NONE
  SAVE


  REAL(DP), ALLOCATABLE     :: flipper_ewald_force(:,:)         ! LEONID: Forces due to ewald field (classic, flipper_force_ewald.f90)
  REAL(DP), ALLOCATABLE     :: flipper_forcelc(:,:)             ! LEONID: Forces from local potentials, flipper_force_lc.f90
  REAL(DP), ALLOCATABLE     :: flipper_forcenl(:,:)             ! LEONID: Forces from local potentials, flipper_force_lc.f90
  REAL(DP), ALLOCATABLE     :: total_force(:,:)                 ! LEONID: Force as calculated by flipper_force_lc + flipper_ewald_force
  REAL(DP)                  ::  flipper_energy_external, &      ! LEONID: The external energy coming from the charge density interacting the the local pseudopotential
                                flipper_ewld_energy, &          ! LEONID: Real for the ewald energy
                                flipper_energy_kin, &           ! LEONID: Kinetic energy, only of the pinballs
                                flipper_nlenergy, &             ! LEONID: flipper_ewld_energy+flipper_energy_external+flipper_energy_kin
                                flipper_epot, &                 ! LEONID: flipper_energy_external + flipper_ewld_energy
                                flipper_cons_qty                ! LEONID: flipper_ewld_energy+flipper_energy_external+flipper_energy_kin

  REAL(DP)                  :: flipper_temp
  INTEGER                   :: nr_of_pinballs                                  ! LEONID

  COMPLEX(DP), ALLOCATABLE  :: charge_g(:)                ! LEONID
  INTEGER                   :: flipper_ityp                     ! LEONID the integer that corresponds to the type of atom that is the pinball
  
  
  CONTAINS
    SUBROUTINE flipper_forces_potener()
        USE io_global,        ONLY : stdout, ionode
        USE force_mod,        ONLY : force
        USE lsda_mod,         ONLY : nspin
        USE fft_base,         ONLY : dfftp
        USE cell_base,        ONLY : omega
        USE vlocal,           ONLY : strf, vloc
        USE ions_base,        ONLY : nat, tau, ntyp => nsp, zv, ityp
        USE gvect,            ONLY : ngm, gstart, g, eigts1, &    !strf 
                                    eigts2, eigts3, gcutm, &
                                    gg,ngl, nl, igtongl
        USE cell_base,        ONLY : bg, at, alat
        USE control_flags,    ONLY : gamma_only
        USE klist,            only : xk
        USE uspp,             ONLY : vkb
        USE wvfct,            ONLY : npw, igk
        
        IMPLICIT NONE
        
        REAL(DP), EXTERNAL :: ewald
        INTEGER  :: iat, ntyp_save, ipol, na
        
        ntyp_save = ntyp
        ntyp = 1

        CALL struc_fact( nat, tau, ntyp, ityp, ngm, g, bg, &
                   dfftp%nr1, dfftp%nr2, dfftp%nr3, strf, eigts1, eigts2, eigts3 )

        CALL flipper_setlocal()

        ! set ntyp back to what it originally was. Is that necessary?
        ntyp = ntyp_save
        
        flipper_ewld_energy = ewald( alat, nat, ntyp, ityp, zv, at, bg, tau, &
                omega, g, gg, ngm, gcutm, gstart, gamma_only, strf )

        CALL flipper_force_ewald( alat, nat, ntyp, ityp, zv, at, bg, tau, omega, g, &
                 gg, ngm, gstart, gamma_only, gcutm, strf )


        CALL flipper_force_lc( nat, tau, ityp, alat, omega, ngm, ngl, igtongl, &
                 g, nl, nspin, gstart, gamma_only, vloc )

        CALL init_us_2( npw, igk, xk(1,1), vkb )

        CALL flipper_force_energy_us (flipper_forcenl, flipper_nlenergy)

        DO ipol = 1, 3
             DO na = 1, nr_of_pinballs
                total_force(ipol,na) = &
                                flipper_ewald_force(ipol,na)    + &
                                flipper_forcelc(ipol,na)        + &
                                flipper_forcenl(ipol, na)
            END DO
        END DO
        force(:,:)=0.d0
        do iat=1,nr_of_pinballs
            force(1:3,iat)=total_force(1:3,iat)
        end do


! 9036 FORMAT(5X,'atom ',I4,' type ',I2,'   force = ',3F14.8)

    END SUBROUTINE flipper_forces_potener




END MODULE flipper_mod

