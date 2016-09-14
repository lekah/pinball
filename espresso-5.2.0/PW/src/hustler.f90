
MODULE hustler

    USE basic_algebra_routines  ! norm
    USE input_parameters,       only : hustlerfile
    USE io_files,               only : iunhustle, tmp_dir
    USE io_global,              ONLY : stdout, ionode
    USE kinds,                  ONLY : DP
    USE constants,              ONLY : au_ps, eps8 , ry_to_kelvin ,amu_ry
    USE input_parameters,       ONLY : lflipper
    USE dynamics_module,    ONLY : write_traj, tau_new, tau_old, vel, mass, dt, ndof, vel_defined
    USE flipper_mod,        ONLY : flipper_ewld_energy, &
                          flipper_energy_external, &
                          flipper_energy_kin, &
                          flipper_cons_qty, nr_of_pinballs
    USE control_flags,      ONLY : iprint, istep
    USE ener,               ONLY : etot
    
    USE ions_base,      ONLY : nat, nsp, ityp, tau, if_pos, atm, amass
    USE cell_base,      ONLY : alat, omega
    USE ener,           ONLY : etot
    USE force_mod,      ONLY : force, lstres
    USE control_flags,  ONLY : istep, nstep, conv_ions, lconstrain, tv0rd, &
                         iverbosity  !LEONID
    !
    USE constraints_module, ONLY : nconstr, check_constraint
    USE constraints_module, ONLY : remove_constr_force, remove_constr_vec

    REAL(DP) :: total_mass, temp_new, temp_av, elapsed_time


    INTEGER :: istep0 = 0 
    INTEGER :: mynat
    INTEGER :: i !, j
    CHARACTER*3 :: atom
    CHARACTER*256 :: buffer
    REAL(DP) :: x, y, z

!    CHARACTER(LEN=256) :: hustlerfile



CONTAINS 

    SUBROUTINE init_hustler()
    
    
        IF ( ionode ) THEN
            write(UNIT=stdout, FMT=*) "   OPENING HUSTLER file ", hustlerfile
            OPEN(UNIT=iunhustle,FILE=trim(tmp_dir)//'../'//trim(hustlerfile),FORM="FORMATTED",STATUS="OLD",ACTION="READ")
        
            
            istep0 = 0
            elapsed_time = 0.D0
              IF ( lflipper) THEN
                ndof = 3*nr_of_pinballs
                mynat = nr_of_pinballs
              ELSE
                  IF ( ANY( if_pos(:,:) == 0 ) ) THEN
                     ndof = 3*nat - count( if_pos(:,:) == 0 ) - nconstr
                  ELSE
                     ndof = 3*nat - 3 - nconstr
                  ENDIF
                  mynat = nat
              END IF

            READ(UNIT=iunhustle, FMT=*) buffer  ! Read the line of stuff
            DO i=1, mynat
                READ(UNIT=iunhustle, FMT=101) atom, x, y ,z
                tau(1,i) = x / alat
                tau(2,i) = y / alat
                tau(3,i) = z / alat

            END DO        
            DO na = 1, nat
                mass(na) = amass( ityp(na) ) * amu_ry
            ENDDO
            tau_new(:,:)= tau(:,:)

        END IF
101 format(A,3(1X,f16.10))

    END SUBROUTINE init_hustler





    SUBROUTINE hustle_ions()

        

      

      ! USE flipper_mod
      ! USE input_parameters,     ONLY : nstep_thermo
      !
      IMPLICIT NONE
      !
      REAL(DP) :: ekin, etotold
      
      REAL(DP) :: delta(3), ml(3), mlt
      INTEGER  :: na
      ! istep0 counts MD steps done during this run
      ! (istep counts instead all MD steps, including those of previous runs)
      



      ! REAL(DP) :: kstress(3,3)
      ! LEONID
      REAL(DP) :: walltime_s
      REAL(DP), EXTERNAL :: get_clock
      INTEGER :: io


      print*, "ISTEP0", istep0

      IF ( istep0 > nstep ) conv_ions = .true.

        READ(UNIT=iunhustle, FMT=*, IOSTAT=io) buffer  ! Read the line of stuff
        if ( io .ne. 0 ) THEN
            conv_ions = .true.
            
        ElSE IF (lflipper) THEN
            DO i=1, nr_of_pinballs
                READ(UNIT=iunhustle, FMT=102) atom, x, y ,z
                tau_new(1,i) = x / alat
                tau_new(2,i) = y / alat
                tau_new(3,i) = z / alat

            END DO
        ELSE
            DO i=1, nat
                READ(UNIT=iunhustle, FMT=102) atom, x, y ,z
                tau_new(1,i) = x / alat
                tau_new(2,i) = y / alat
                tau_new(3,i) = z / alat
            END DO
        END IF
        IF ( istep0 > 0 ) THEN
            vel = ( tau_new - tau_old ) / ( 2.D0*dt )
        ELSE
            vel(:,:) = 0.0D0 ! TODO what if velocities are in input file
        END IF

    ! LEONID: Computing kinetic energy and temperature for the flipper
      if (lflipper) then
          flipper_energy_kin= 0.D0  ! LEONID: Here we calculate the kinetic energy of the flipper
          ! LEONID: Only the pinballs contribute:
          DO na = 1, nr_of_pinballs
             !
             ! ml(:) = ml(:) + vel(:,na) * mass(na)
             flipper_energy_kin  = flipper_energy_kin + 0.5D0 * mass(na) * &
                            ( vel(1,na)**2 + vel(2,na)**2 + vel(3,na)**2 )
            ! print*, flipper_energy_kin, mass(na)
          ENDDO
          
          flipper_energy_kin = flipper_energy_kin * alat**2

          ekin = flipper_energy_kin
          etot = flipper_energy_external + flipper_ewld_energy
          flipper_cons_qty  =  ekin + etot
      else
          ml   = 0.D0
          ekin = 0.D0
 !         kstress = 0.d0
          !
          DO na = 1, nat
             !
!             ml(:) = ml(:) + vel(:,na) * mass(na)
             ekin  = ekin + 0.5D0 * mass(na) * &
                            ( vel(1,na)**2 + vel(2,na)**2 + vel(3,na)**2 )
!~              do i = 1, 3
!~                  do j = 1, 3
!~                      kstress(i,j) = kstress(i,j) + mass(na)*vel(i,na)*vel(j,na)
!~                  enddo
!~              enddo
!~              !
          ENDDO
          !
          ekin = ekin*alat**2
      end if    

      temp_new = 2.D0 / dble( ndof ) * ekin * ry_to_kelvin
      temp_av = temp_av + temp_new
      walltime_s = get_clock( 'PWSCF' )

      IF (mod(istep, iprint) .eq. 0) THEN
          IF (lflipper) THEN
              call write_traj(                                              &
                    istep, elapsed_time, tau, vel, force, ekin, etot,       &
                    flipper_cons_qty, temp_new, walltime_s, nr_of_pinballs  &
                ) !aris ! LEONID: added force and walltime
          ELSE
              call write_traj(                                              &
                      istep, elapsed_time, tau, vel, force, ekin, etot,     &
                      ekin+etot, temp_new, walltime_s, nat                  &
                  ) !aris ! LEONID: added force and walltime
          END IF
      END IF

      elapsed_time = elapsed_time + dt*2.D0*au_ps
      !
      istep0= istep0+ 1
      istep = istep + 1


      ! ... here the tau are shifted
      ! LEONID: I moved this to this point, since before, since tau was overwritten with
      ! future positions before being printed
      ! 
      tau_old(:,:) = tau(:,:)
      tau(:,:) = tau_new(:,:)
      
      
102 format(A,3(1X,f16.10))


    END SUBROUTINE hustle_ions

    SUBROUTINE end_hustler()
    
    
        IF ( ionode ) THEN
            write(UNIT=stdout, FMT=*) "   CLOSING HUSTLER file ", hustlerfile
            CLOSE(UNIT=iunhustle, STATUS="KEEP")
        END IF

    END SUBROUTINE end_hustler



END MODULE