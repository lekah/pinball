subroutine deallocate_flipper
use pinball

    DEALLOCATE(flipper_ewald_force)                 ! LEONID
    DEALLOCATE(flipper_forcelc)                 ! LEONID
    DEALLOCATE(total_force)                 ! LEONID
    DEALLOCATE(charge_g)
!
end subroutine deallocate_flipper


