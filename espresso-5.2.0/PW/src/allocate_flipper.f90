subroutine allocate_flipper
use flipper_mod
use gvect, only :ngm
!
     ALLOCATE(flipper_ewald_force( 3, nr_of_pinballs ))
     ALLOCATE(flipper_forcelc( 3, nr_of_pinballs))
     ALLOCATE(total_force( 3, nr_of_pinballs))
     ALLOCATE(charge_g(ngm)) 
!
end subroutine allocate_flipper


