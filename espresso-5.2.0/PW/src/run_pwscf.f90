!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE run_pwscf ( exit_status ) 
  !----------------------------------------------------------------------------
  !
  ! ... Run an instance of the Plane Wave Self-Consistent Field code 
  ! ... MPI initialization and input data reading is performed in the 
  ! ... calling code - returns in exit_status the exit code for pw.x, 
  ! ... returned in the shell. Values are:
  ! ... * 0: completed successfully
  ! ... * 1: an error has occurred (value returned by the errore() routine)
  ! ... * 2-127: convergence error
  ! ...   * 2: scf convergence error
  ! ...   * 3: ion convergence error
  ! ... * 128-255: code exited due to specific trigger
  !       * 255: exit due to user request, or signal trapped,
  !              or time > max_seconds
  ! ...     (note: in the future, check_stop_now could also return a value
  ! ...     to specify the reason of exiting, and the value could be used
  ! ..      to return a different value for different reasons)
  ! ... Will be eventually merged with NEB
  !
  USE io_global,        ONLY : stdout, ionode, ionode_id
  USE parameters,       ONLY : ntypx, npk, lmaxx
  USE cell_base,        ONLY : fix_volume, fix_area
  USE control_flags,    ONLY : conv_elec, gamma_only, lscf, twfcollect
  USE control_flags,    ONLY : conv_ions, istep, nstep, restart, lmd, lbfgs
  USE force_mod,        ONLY : lforce, lstres, sigma, force
  USE check_stop,       ONLY : check_stop_init, check_stop_now
  USE mp_images,        ONLY : intra_image_comm
  USE qmmm,             ONLY : qmmm_initialization, qmmm_shutdown, &
                               qmmm_update_positions, qmmm_update_forces
  USE pinball,          ONLY : deallocate_flipper, lflipper, flipper_forces_potener
  USE hustler,          ONLY : init_hustler, end_hustler, lhustle

  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: exit_status


  
  IF (ionode .and. lflipper) THEN
    print*, ' ##########################################'
    print*, '      THIS IS A FLIPPER CALCULATION'
    print*, ' ##########################################'
    
    
    
  END IF
  ! END LEONID
  
  !
  exit_status = 0
  IF ( ionode ) WRITE( unit = stdout, FMT = 9010 ) ntypx, npk, lmaxx
  !
  IF (ionode) CALL plugin_arguments()
  CALL plugin_arguments_bcast( ionode_id, intra_image_comm )
  !
  ! ... needs to come before iosys() so some input flags can be
  !     overridden without needing to write PWscf specific code.
  ! 
  CALL qmmm_initialization()
  !
  ! ... convert to internal variables
  !
  CALL iosys()
  !
  IF ( gamma_only ) WRITE( UNIT = stdout, &
     & FMT = '(/,5X,"gamma-point specific algorithms are used")' )
  !
  ! call to void routine for user defined / plugin patches initializations
  !
  CALL plugin_initialization()
  !
  CALL check_stop_init()
  !
  CALL setup ()
  !
  CALL qmmm_update_positions()
  !
  CALL init_run()
  !
  ! ... dry run: code will stop here if called with exit file present
  ! ... useful for a quick and automated way to check input data
  !
  IF ( check_stop_now() ) THEN
     CALL punch( 'config' )
     exit_status = 255
     RETURN
  ENDIF


  
  
  main_loop: DO
     !
     ! ... electronic self-consistency or band structure calculation
     !
     
     IF ( lflipper) THEN
        IF ( check_stop_now() ) THEN
            exit_status = 255
            ! CALL punch( 'config' )
            RETURN
        ENDIF
     ELSE
         IF ( .NOT. lscf) THEN
            CALL non_scf ()
         ELSE
            CALL electrons()
         END IF
         !
         ! ... code stopped by user or not converged
         !
         IF (( check_stop_now() .OR. .NOT. conv_elec ) .AND. ( .NOT. lhustle)) THEN
            IF ( check_stop_now() ) exit_status = 255
            !
            ! I will only exit the calculation if this is not a hustler
            IF ( .NOT. conv_elec ) exit_status =  2
            ! workaround for the case of a single k-point
            twfcollect = .FALSE.
            CALL punch( 'config' )
            RETURN
         ENDIF
     ENDIF
     !
     ! ... ionic section starts here
     !
     
     CALL start_clock( 'ions' )
     conv_ions = .TRUE.
     !
     ! ... recover from a previous run, if appropriate
     !
     !IF ( restart .AND. lscf ) CALL restart_in_ions()
     !
     ! ... file in CASINO format written here if required
     !
     IF ( lmd ) CALL pw2casino()
     !
     ! ... force calculation
     !
     
     ! LEONID
     IF (lflipper) THEN

        CALL flipper_forces_potener()
        CALL start_clock( 'move_ions' )
        CALL move_ions()
        CALL stop_clock( 'move_ions' )

      ELSE 
         IF ( lforce ) CALL forces()          
         !
         ! ... stress calculation
         !
         IF ( lstres ) CALL stress ( sigma )
         !
         ! ... send out forces to MM code in QM/MM run
         !
         CALL qmmm_update_forces(force)
         !
         IF ( lmd .OR. lbfgs ) THEN
            !
            if (fix_volume) CALL impose_deviatoric_stress(sigma)
            !
            if (fix_area)  CALL  impose_deviatoric_stress_2d(sigma)
            !
            ! ... ionic step (for molecular dynamics or optimization)
            !
            CALL move_ions()
            !
            ! ... then we save restart information for the new configuration
            IF ( istep < nstep .AND. .NOT. conv_ions ) &
               CALL punch( 'config' )
            !
         END IF
     END IF
     !
     CALL stop_clock( 'ions' )
     !
     ! ... exit condition (ionic convergence) is checked here
     !
     IF ( conv_ions ) EXIT main_loop
     !
     ! ... receive new positions from MM code in QM/MM run
     !
     IF ( .NOT. lflipper) THEN
        CALL qmmm_update_positions()
         !
         ! ... terms of the hamiltonian depending upon nuclear positions
         ! ... are reinitialized here
         !
         IF ( lmd .OR. lbfgs ) CALL hinit1()
         !
     END IF
  END DO main_loop

  IF (lflipper)  then 
      call deallocate_flipper()
  END IF
  
  if ( lhustle )  CALL end_hustler()
  !
  ! ... save final data file
  !
  IF ( .not. lmd) CALL pw2casino()
  
  IF ( .NOT. lflipper) THEN
    CALL punch('all')
    CALL qmmm_shutdown()
  ENDIF
  
  !
  IF ( .NOT. conv_ions )  exit_status =  3
  RETURN
  !
9010 FORMAT( /,5X,'Current dimensions of program PWSCF are:', &
           & /,5X,'Max number of different atomic species (ntypx) = ',I2,&
           & /,5X,'Max number of k-points (npk) = ',I6,&
           & /,5X,'Max angular momentum in pseudopotentials (lmaxx) = ',i2)
           
9035 FORMAT(5X,'atom ',I4,' type ',I2,'   force = ',3F14.8)  ! LEONID
END SUBROUTINE run_pwscf
