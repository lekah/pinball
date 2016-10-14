!
! Copyright (C) 2001-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
!subroutine flipper_force_lc (nat, nr_of_pinballs, tau, ityp, alat, omega, ngm, ngl, &
!     igtongl, g, aux_rho_of_g, nl, nspin, gstart, gamma_only, vloc, forcelc)
subroutine flipper_force_lc (nat, tau, ityp, alat, omega, ngm, ngl, &
     igtongl, g, nl, nspin, gstart, gamma_only, vloc)

  !----------------------------------------------------------------------
  !
  USE kinds
  USE constants, ONLY : tpi
  USE mp_bands,  ONLY : intra_bgrp_comm
  USE mp,        ONLY : mp_sum
  USE fft_base,  ONLY : dfftp
  
  USE esm,       ONLY : esm_force_lc, esm_bc ! do_comp_esm,
  USE pinball
  implicit none
  !
  !   first the dummy variables
  !
  integer, intent(in) :: nat, ngm, nspin, ngl, gstart, &
                         igtongl (ngm), nl (ngm), ityp (nat)
!, nr_of_pinballs
  ! nat:    number of atoms in the cell
  ! ngm:    number of G vectors
  ! nspin:  number of spin polarizations
  ! ngl:    number of shells
  ! igtongl correspondence G <-> shell of G
  ! nl:     correspondence fft mesh <-> G vec
  ! ityp:   types of atoms

  logical, intent(in) :: gamma_only

  real(DP), intent(in) :: tau (3, nat), g (3, ngm), vloc (ngl, * ), &
         alat, omega
  ! tau:  coordinates of the atoms
  ! g:    coordinates of G vectors
  ! vloc: local potential
  ! rho:  valence charge
  ! alat: lattice parameter
  ! omega: unit cell volume

!  real(DP), intent(out) :: forcelc (3, nr_of_pinballs)
  ! the local-potential contribution to forces on atoms

  integer :: ipol, ig, na
  ! counter on polarizations
  ! counter on G vectors
  ! counter on atoms

  ! complex(DP), intent(in) :: aux_rho_of_g (dfftp%nnr)
  ! auxiliary space for FFT
  real(DP) :: arg, fact
  !
  ! contribution to the force from the local part of the bare potential
  ! F_loc = Omega \Sum_G n*(G) d V_loc(G)/d R_i
  !

  if (gamma_only) then
     fact = 2.d0
  else
     fact = 1.d0
  end if
  do na = 1, nr_of_pinballs
     do ipol = 1, 3
        flipper_forcelc (ipol, na) = 0.d0
     enddo
  enddo
  do na = 1, nr_of_pinballs
     ! contribution from G=0 is zero
     do ig = gstart, ngm
        arg = (g (1, ig) * tau (1, na) + g (2, ig) * tau (2, na) + &
               g (3, ig) * tau (3, na) ) * tpi
        do ipol = 1, 3
!~            flipper_forcelc (ipol, na) = flipper_forcelc (ipol, na) + &
!~                 g (ipol, ig) * vloc (igtongl (ig), ityp (na) ) * &
!~                 (sin(arg)*DBLE(aux_rho_of_g(nl(ig))) + cos(arg)*AIMAG(aux_rho_of_g(nl(ig))) )
           flipper_forcelc (ipol, na) = flipper_forcelc (ipol, na) + &
                g (ipol, ig) * vloc (igtongl (ig), ityp (na) ) * &
                (sin(arg)*DBLE(charge_g(ig)) + cos(arg)*AIMAG(charge_g(ig)) )
        enddo
     enddo
     do ipol = 1, 3
        flipper_forcelc (ipol, na) = fact * flipper_forcelc (ipol, na) * omega * tpi / alat
     enddo
  enddo
!~   IF ( do_comp_esm .and. ( esm_bc .ne. 'pbc' ) ) THEN
!~      !
!~      ! ... Perform corrections for ESM method (add long-range part)
!~      !
!~      CALL esm_force_lc ( aux_rho_of_g, flipper_forcelc )
!~   ENDIF
  !
  call mp_sum(  flipper_forcelc, intra_bgrp_comm )
  !
  ! deallocate (aux)
  return
end subroutine flipper_force_lc
