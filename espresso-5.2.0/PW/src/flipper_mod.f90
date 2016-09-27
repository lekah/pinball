MODULE flipper_mod

  USE kinds, ONLY: DP

  SAVE


  REAL(DP), ALLOCATABLE     :: flipper_ewald_force(:,:)         ! LEONID: Forces due to ewald field (classic, flipper_force_ewald.f90)
  REAL(DP), ALLOCATABLE     :: flipper_forcelc(:,:)             ! LEONID: Forces from local potentials, flipper_force_lc.f90
  REAL(DP), ALLOCATABLE     :: flipper_forcenl(:,:)             ! LEONID: Forces from local potentials, flipper_force_lc.f90
  REAL(DP), ALLOCATABLE     :: total_force(:,:)                 ! LEONID: Force as calculated by flipper_force_lc + flipper_ewald_force
  REAL(DP)                  ::  flipper_energy_external, &      ! LEONID: The external energy coming from the charge density interacting the the local pseudopotential
                                flipper_ewld_energy, &          ! LEONID: Real for the ewald energy
                                flipper_energy_kin, &           ! LEONID: Kinetic energy, only of the pinballs
                                flipper_epot, &                 ! LEONID: flipper_energy_external + flipper_ewld_energy
                                flipper_cons_qty, &             ! LEONID: flipper_ewld_energy+flipper_energy_external+flipper_energy_kin
                                flipper_nlenergy                ! LEONID: flipper_ewld_energy+flipper_energy_external+flipper_energy_kin
  REAL(DP)                  :: flipper_temp
  INTEGER                   :: nr_of_pinballs                                  ! LEONID
!~   LOGICAL                   :: lflipper                                  ! LEONID
  COMPLEX(DP), ALLOCATABLE  :: charge_g(:)                ! LEONID
  INTEGER                   :: flipper_ityp                     ! LEONID the integer that corresponds to the type of atom that is the pinball
  
  

END MODULE flipper_mod

