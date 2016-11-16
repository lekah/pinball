MODULE pinball

  USE kinds, ONLY: DP

  IMPLICIT NONE
  SAVE


  REAL(DP), ALLOCATABLE     :: flipper_ewald_force(:,:)         ! LEONID: Forces due to ewald field (classic, flipper_force_ewald.f90)
  REAL(DP), ALLOCATABLE     :: flipper_ewald_force_pinball(:,:)         ! LEONID: Forces due to ewald field (classic, flipper_force_ewald.f90)
  REAL(DP), ALLOCATABLE     :: flipper_ewald_force_rigid(:,:)         ! LEONID: Forces due to ewald field (classic, flipper_force_ewald.f90)
  REAL(DP), ALLOCATABLE     :: flipper_forcelc(:,:)             ! LEONID: Forces from local potentials, flipper_force_lc.f90
  REAL(DP), ALLOCATABLE     :: flipper_forcenl(:,:)             ! LEONID: Forces from local potentials, flipper_force_lc.f90
  !   REAL(DP), ALLOCATABLE     :: total_force(:,:)                 ! LEONID: Force as calculated by flipper_force_lc + flipper_ewald_force
  REAL(DP)                  ::  flipper_energy_external, &      ! LEONID: The external energy coming from the charge density interacting the the local pseudopotential
                                flipper_ewld_energy_rigid, &    ! LEONID: Real for the ewald energy
                                flipper_ewld_energy_pinball, &  ! LEONID: Real for the ewald energy
                                flipper_ewld_energy_total, &  ! LEONID: Real for the ewald energy
                                flipper_ewld_energy, &  ! LEONID: Real for the ewald energy
                                flipper_energy_kin, &           ! LEONID: Kinetic energy, only of the pinballs
                                flipper_nlenergy, &             ! LEONID: flipper_ewld_energy+flipper_energy_external+flipper_energy_kin
                                flipper_epot, &                 ! LEONID: flipper_energy_external + flipper_ewld_energy
                                flipper_cons_qty                ! LEONID: flipper_ewld_energy+flipper_energy_external+flipper_energy_kin

  REAL(DP)                  :: flipper_nonlocal_correction, &
                                flipper_local_factor, &
                                flipper_ewald_pinball_factor, &
                                flipper_ewald_rigid_factor
  INTEGER                   :: nr_of_pinballs                                  ! LEONID

  COMPLEX(DP), ALLOCATABLE  :: charge_g(:)                ! LEONID
  COMPLEX(DP), allocatable  :: evc_grad(:,:,:)
  INTEGER                   :: flipper_ityp                     ! LEONID the integer that corresponds to the type of atom that is the pinball
  

  LOGICAL :: lflipper
  LOGICAL :: flipper_do_nonloc = .true.
  LOGICAL, allocatable       :: flipper_ewald_list(:)

  CONTAINS

    SUBROUTINE init_flipper()
        USE scf,                ONLY : rho, charge_density
        USE fft_base,           ONLY : dfftp
        USE ions_base,          ONLY : ityp, nat, nsp
        use wvfct,              only : nbnd, wg
        USE mp_bands,           ONLY : intra_bgrp_comm
        USE io_global,          ONLY : stdout, ionode
        USE fft_base,           ONLY : dffts
        USE gvecs,              ONLY : nls, nlsm
        USE wvfct,              ONLY : npw, igk
        USE fft_interfaces,     ONLY : invfft
        USE cell_base,          ONLY : omega
        use io_files,           ONLY : nwordwfc, diropn, iunwfc
        USE io_files,           ONLY : prefix
        USE input_parameters,   ONLY : prefix_flipper_charge
        USE lsda_mod,           ONLY : nspin
        USE wavefunctions_module,  ONLY: psic, evc
        USE gvect,              ONLY :  ngm, gstart, g, eigts1, &
                                    eigts2, eigts3, gcutm, &
                                    gg,ngl, nl, igtongl
        USE fft_interfaces,     ONLY : fwfft
        USE mp,                 ONLY : mp_sum
        USE io_rho_xml,         ONLY : read_rho

        IMPLICIT NONE
        
        INTEGER  :: igrid,iatom,igm
        INTEGER  :: counter
        INTEGER  :: na
        INTEGER  :: ipol, ibnd, ig
        INTEGER :: iv, i, iun
        INTEGER, EXTERNAL :: find_free_unit
        LOGICAL  :: exst
        
        REAL(DP) :: q_tot, q_tot2, q_tot3
        REAL(DP), allocatable :: charge2(:)




        nr_of_pinballs = 0
        wg(:,:) = 2.0d0
        DO counter = 1, nat
            IF (ityp(counter) == 1) THEN
                nr_of_pinballs = nr_of_pinballs + 1
            END IF
        END DO

        if (ionode) THEN
            print*, '    THIS IS A FLIPPER CALCULATION'
            print*, ' nr_of_pinballs               =  ', nr_of_pinballs
            print*, ' flipper_nonlocal_correation  = ', flipper_nonlocal_correction
            print*, ' flipper_local_factor         = ', flipper_local_factor
            print*, ' flipper_ewald_rigid_factor   = ', flipper_ewald_rigid_factor
            print*, ' flipper_ewald_pinball_factor = ', flipper_ewald_pinball_factor
        end if
        call allocate_flipper()   


        CALL read_rho(charge_density, nspin)

        iun=find_free_unit()


        call diropn(iun, 'wfc', 2*nwordwfc, exst )
        call davcio (evc, 2*nwordwfc, iun, 1,-1) 


! This whole bunch of code below provided a check to see whether the charge sums up
! correctly from the wavefunctions.

!~         allocate(charge2(dffts%nnr))
!~         charge2(:)=0.d0
!~         do iv=1,nbnd !,2
!~            psic=0.d0
!~ 
!~            psic(nls(1:npw))=evc(1:npw,iv)
!~            psic(nlsm(1:npw))=CONJG(evc(1:npw,iv))
!~ 
!~            call invfft ('Wave', psic, dffts)
!~            charge2(1:dffts%nnr)=charge2(1:dffts%nnr)+dble(psic(1:dffts%nnr))**2.0 ! + dimag(psic(1:dffts%nnr))**2.0
!~ 
!~         end do
!~         q_tot=0.
!~         q_tot2=0.
!~         q_tot3=0.
!~         do i=1,dffts%nnr
!~            q_tot=q_tot + ( 2.d0 * charge2(i) - omega * charge_density%of_r(i, 1) )**2
!~            q_tot2=q_tot2 + 2.d0 *charge2(i)
!~            q_tot3=q_tot3 + charge_density%of_r(i, 1) * omega
!~         end do
!~         q_tot = q_tot / (dffts%nr1*dffts%nr2*dffts%nr3)
!~         q_tot2 = q_tot2 / (dffts%nr1*dffts%nr2*dffts%nr3)
!~         q_tot3 = q_tot3 /    (dffts%nr1*dffts%nr2*dffts%nr3)
!~         call mp_sum(q_tot, intra_bgrp_comm)
!~         call mp_sum(q_tot2, intra_bgrp_comm)
!~         call mp_sum(q_tot3, intra_bgrp_comm)
!~         print*,'q_tot',q_tot
!~         print*,'q_tot2',q_tot2
!~         print*,'q_tot3',q_tot3
!~         deallocate(charge2)
!~         ! temp END
!~          print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

         psic(:)=(0.d0,0.d0)

         psic(1:dfftp%nnr)=dcmplx(charge_density%of_r(1:dfftp%nnr,1)) 
         CALL fwfft ('Dense', psic, dfftp)
         DO igm=1,ngm 
            charge_g(igm)=psic(nl(igm))
         END DO


        DO ipol = 1, 3
            DO ibnd = 1, nbnd
                DO ig=1, npw
                    evc_grad(ipol, ig, ibnd) = evc(ig, ibnd) * (0.D0,1.D0) * g(ipol,igk(ig))
                END DO
            END DO
        END DO


    flipper_ewald_list(1:nsp) = .true.
!~     flipper_ewald_list(1) = .false.
    END SUBROUTINE init_flipper
    
    
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
        
        REAL(DP), EXTERNAL :: flipper_energy_ewald
        INTEGER  :: iat, ntyp_save, ipol, na
        
        
        ntyp_save = ntyp
        ntyp = 1
        CALL struc_fact( nat, tau, ntyp, ityp, ngm, g, bg, &
                   dfftp%nr1, dfftp%nr2, dfftp%nr3, strf, eigts1, eigts2, eigts3 )
        ntyp = ntyp_save


    !!!!!!!!!!!!!!!!!!!!!!!!!!! LOCAL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Calculating the local contribution to potential energy:
        CALL flipper_setlocal()
        !
        ! Adding the factor:
        flipper_energy_external = flipper_local_factor*flipper_energy_external
        CALL start_clock( 'pb_force_lc' )
        CALL flipper_force_lc( nat, tau, ityp, alat, omega, ngm, ngl, igtongl, &
                 g, nl, nspin, gstart, gamma_only, vloc )
        ! Applying the factor:
        flipper_forcelc(:,:) = flipper_local_factor*flipper_forcelc(:,:)
        CALL stop_clock( 'pb_force_lc' )




    !!!!!!!!!!!!!!!!!!!!!!!!!!! EWALD !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! ENERGY :
        ! ewald for all atoms PB-PB, PB-R, R-R:
        flipper_ewald_list(1:ntyp) = .true.
        flipper_ewld_energy_total = flipper_energy_ewald( alat, nat, ntyp, ityp, zv, at, bg, tau, &
                omega, g, gg, ngm, gcutm, gstart, gamma_only, strf )

        flipper_ewald_list(2:ntyp) = .false.
        ! This is only pinball-pinball ionic interaction
        flipper_ewld_energy_pinball = flipper_energy_ewald( alat, nat, ntyp, ityp, zv, at, bg, tau, &
                omega, g, gg, ngm, gcutm, gstart, gamma_only, strf )


        ! For phonons??
        ! PB-R, R-R
        ! Put factors
        flipper_ewld_energy = &
            flipper_ewald_rigid_factor*(flipper_ewld_energy_total - flipper_ewld_energy_pinball)  &
            + flipper_ewald_pinball_factor*flipper_ewld_energy_pinball
        ! FORCES:
        CALL start_clock( 'pb_ewald' )
        ! Calculating forces for PB-R
        flipper_ewald_list(2:ntyp) = .true.
        flipper_ewald_list(1) = .false.
        CALL flipper_force_ewald( alat, nat, ntyp, ityp, zv, at, bg, tau, omega, g, &
                 gg, ngm, gstart, gamma_only, gcutm, strf )
        flipper_ewald_force_rigid(:,:) = flipper_ewald_force(:,:)
        ! Calculating forces for PB-PB
        flipper_ewald_list(2:ntyp) = .false.
        flipper_ewald_list(1) = .true.
        CALL flipper_force_ewald( alat, nat, ntyp, ityp, zv, at, bg, tau, omega, g, &
                 gg, ngm, gstart, gamma_only, gcutm, strf )
        flipper_ewald_force_pinball(:,:) = flipper_ewald_force(:,:)

        flipper_ewald_force(:,:) =  &
                flipper_ewald_pinball_factor*flipper_ewald_force_pinball(:,:) &
                +  flipper_ewald_rigid_factor*flipper_ewald_force_rigid(:,:)
        CALL stop_clock( 'pb_ewald' )



    !!!!!!!!!!!!!!!!!!!!!!!! NONLOCAL  !!!!!!!!!!!!!!!!!!!!!!!!!!

    
        CALL start_clock( 'pb_nonloc' )
        IF ( flipper_do_nonloc ) THEN
            CALL init_us_2( npw, igk, xk(1,1), vkb )
            CALL flipper_force_energy_us (flipper_forcenl, flipper_nlenergy)
            ! Applying the correction to the nonlocal term based on a linear factor
            flipper_forcenl(:,:) = flipper_nonlocal_correction * flipper_forcenl(:,:)
            flipper_nlenergy     = flipper_nonlocal_correction*flipper_nlenergy
        ELSE
            flipper_forcenl(:,:) =  0.D0
            flipper_nlenergy     = 0.D0
        ENDIF
        CALL stop_clock( 'pb_nonloc' )


        ! flipper_forcelc(:,:) = 0.99057635*flipper_forcelc(:,:)

    !!!!!!!!!!!!!!!!!!!!!!! SUMMING THE FORCES !!!!!!!!!!!!!!!!!!!!!!!


        DO ipol = 1, 3
             DO na = 1, nr_of_pinballs
                force(ipol,na) = &
                                flipper_ewald_force(ipol,na)    + &
                                flipper_forcelc(ipol,na)        + &
                                flipper_forcenl(ipol, na)
            END DO
            ! In the static picture, all other forces are 0:
            DO na = nr_of_pinballs+1, nat
                force(ipol, na) = 0.D0
            END DO
        END DO




! 9036 FORMAT(5X,'atom ',I4,' type ',I2,'   force = ',3F14.8)

    END SUBROUTINE flipper_forces_potener



    subroutine deallocate_flipper


        DEALLOCATE(flipper_ewald_force)                 ! LEONID
        DEALLOCATE(flipper_ewald_force_pinball)                 ! LEONID
        DEALLOCATE(flipper_ewald_force_rigid)                 ! LEONID
        DEALLOCATE(flipper_forcelc)                 ! LEONID
        ! DEALLOCATE(total_force)                 ! LEONID
        DEALLOCATE(charge_g)
        DEALLOCATE(evc_grad)
        DEALLOCATE(flipper_ewald_list)

    end subroutine deallocate_flipper



    subroutine allocate_flipper
        use gvect,                  ONLY :ngm
        use ions_base,              ONLY : nat, nsp
        USE wvfct,                  ONLY : nbnd, npwx
        USE noncollin_module,       ONLY : npol
        ALLOCATE(flipper_ewald_force( 3, nr_of_pinballs ))
        ALLOCATE(flipper_ewald_force_pinball( 3, nr_of_pinballs ))
        ALLOCATE(flipper_ewald_force_rigid( 3, nr_of_pinballs ))
        ALLOCATE(flipper_forcelc( 3, nr_of_pinballs))
        ! ALLOCATE(total_force( 3, nr_of_pinballs))
        ALLOCATE(charge_g(ngm)) 
        ALLOCATE(flipper_forcenl( 3, nat ))
        ALLOCATE(evc_grad(3, npwx*npol, nbnd ) )
        ALLOCATE(flipper_ewald_list(nsp))

    end subroutine allocate_flipper



    SUBROUTINE flipper_setlocal()
      !----------------------------------------------------------------------
      !
      !    This routine computes the local potential in real space vltot(ir)
      !
      USE kinds,     ONLY : DP
      USE constants, ONLY : eps8
      USE ions_base, ONLY : zv, ntyp => nsp
      USE cell_base, ONLY : omega
      USE extfield,  ONLY : tefield, dipfield, etotefield
      USE gvect,     ONLY : igtongl, gg
      USE scf,       ONLY : rho, v_of_0, vltot
      USE vlocal,    ONLY : strf, vloc
      USE fft_base,  ONLY : dfftp
      USE fft_interfaces,ONLY : invfft
      USE gvect,     ONLY : nl, nlm, ngm, gstart
      USE control_flags, ONLY : gamma_only
      USE mp_bands,  ONLY : intra_bgrp_comm
      USE mp,        ONLY : mp_sum
      USE martyna_tuckerman, ONLY : wg_corr_loc, do_comp_mt
      USE esm,       ONLY : esm_local, esm_bc, do_comp_esm
      USE qmmm,      ONLY : qmmm_add_mm_field

      IMPLICIT NONE

      COMPLEX(DP), ALLOCATABLE :: vltot_g (:)

      INTEGER :: nt, ng, igrid,igm

      ALLOCATE (vltot_g( ngm))

      
      vltot_g(:)=(0.d0,0.d0)

      DO nt = 1, 1   !  ntyp
          DO ng = 1, ngm
              vltot_g (ng)                    = vltot_g(ng) + vloc (igtongl (ng), nt) * strf (ng, nt)
                 ! Useful to get the real part of the VLTOT
                 ! vltot_g_to_get_vltot_r (nl(ng)) = vltot_g_to_get_vltot_r(nl(ng)) + vloc (igtongl (ng), nt) * strf (ng, nt)
          END DO
      END DO

      flipper_energy_external = 0.d0
      IF (gamma_only) THEN
          DO igm = gstart, ngm
            flipper_energy_external = flipper_energy_external &
                     + 2.d0 * DBLE(charge_g(igm))*DBLE(vltot_g(igm)) &
                     + 2.d0 * DIMAG(charge_g(igm))*DIMAG(vltot_g(igm)) 
          END DO
          if (gstart==2) flipper_energy_external = flipper_energy_external + charge_g(1)*conjg(vltot_g(1))
      ELSE
          DO igm = 1, ngm !aris
              flipper_energy_external = flipper_energy_external &
                     + DBLE(charge_g(igm)*CONJG(vltot_g(igm)))
          END DO
      END IF    
      
      flipper_energy_external = flipper_energy_external * omega
      
      call mp_sum(flipper_energy_external,intra_bgrp_comm)

      DEALLOCATE(vltot_g)

    END SUBROUTINE flipper_setlocal



END MODULE pinball

