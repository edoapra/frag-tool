! This file is part of frag.
! SPDX-Identifier: Apache-2.0
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

!> Implementation of the dimer projection method for extended tight binding methods.
module frag
   use mctc_env, only : error_type, fatal_error, get_argument, wp
   use mctc_io, only : structure_type, new, write_structure
   use mctc_io_math, only : matinv_3x3
   use mctc_io_convert, only : autoev
   use frag_bondorder, only : get_wiberg_bondorder
   use frag_fragment, only : get_wiberg_fragment, unfold_fragment
   use frag_output, only : format_list, to_string
   use frag_rad, only : guess_bonds
   use frag_xtb, only : get_calculator
   use tblite_basis_type, only : get_cutoff, basis_type
   use tblite_blas, only : dot, gemv, gemm
   use tblite_context_type, only : context_type
   use tblite_cutoff, only : get_lattice_points
   use tblite_integral_overlap, only : get_overlap
   use tblite_output_format, only : format_string
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction, &
      & get_density_matrix
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   implicit none
   private

   public :: get_jab, jab_input

   !> Configuration data for calculation
   type :: jab_input
      !> Name of the requested tight binding method
      character(len=:), allocatable :: method
      !> List of fragments, generated if not given here
      integer, allocatable :: fragment(:)
      !> Threshold for generating fragments from bond orders
      real(wp) :: thr = 0.1_wp
      !> Calculation accuracy
      real(wp) :: accuracy = 1.0_wp
      !> Output verbosity
      integer :: verbosity = 2
      !> Electronic temperature in Kelvin
      real(wp) :: etemp = 300.0_wp
   end type jab_input

   !> Conversion factor from temperature to energy (Boltzmann's constant in atomic units)
   real(wp), parameter :: ktoau = 3.166808578545117e-06_wp

contains

!> Entry point for calculation of dipole projection related properties
subroutine get_jab(input, mol, error)
   !> Input data for the calculation
   type(jab_input), intent(in) :: input
   !> Molecular structure data
   type(structure_type), intent(inout) :: mol
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: spin, charge, stat, unit, ifr, nfrag, nao
   logical :: exist
   real(wp) :: energy, cutoff, jab, sab, jeff, invlat(3, 3)
   real(wp), allocatable :: gradient(:, :), sigma(:, :)
   type(context_type) :: ctx
   type(xtb_calculator) :: calc, fcalc
   type(structure_type), allocatable :: mfrag(:)
   type(wavefunction_type) :: wfn
   type(wavefunction_type), allocatable :: wfx(:)
   real(wp), allocatable :: overlap(:, :), trans(:, :), wbo(:, :), abc(:, :)
   real(wp), allocatable :: orbital(:, :), scmat(:, :), fdim(:, :), scratch(:), efrag(:)
   integer, allocatable :: fragment(:)

   call get_calculator(calc, mol, input%method, error)
   if (allocated(error)) return

   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, &
      & input%etemp * ktoau)

   call xtb_singlepoint(ctx, mol, calc, wfn, input%accuracy, energy, &
      & verbosity=input%verbosity-1)
   if (ctx%failed()) then
      call ctx%get_error(error)
      return
   end if

   allocate(overlap(calc%bas%nao, calc%bas%nao))
   cutoff = get_cutoff(calc%bas)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans)
   call get_overlap(mol, trans, cutoff, calc%bas, overlap)


   if (allocated(input%fragment)) then
      fragment = input%fragment
   else
      allocate(wbo(mol%nat, mol%nat))
      call get_wiberg_bondorder(calc%bas, overlap, wfn%density, wbo)

      allocate(fragment(mol%nat))
      call get_wiberg_fragment(fragment, wbo, input%thr)
   end if

   nfrag = maxval(fragment)
   select case(nfrag)
   case(:1)
      !call fatal_error(error, "Found no fragments in the input structure")
      !return
   case(2:)
      call ctx%message("Found "//to_string(nfrag)//" fragments:")
      do ifr = 1, nfrag
         call ctx%message("Fragment "//to_string(ifr)//":  "//format_list(fragment == ifr))
      end do
      call ctx%message("")
   end select

   call guess_bonds(mol, wbo)

   call get_wiberg_fragment(fragment, wbo, input%thr)

   nfrag = maxval(fragment)
   select case(nfrag)
   case(:1)
      !call fatal_error(error, "Found no fragments in the input structure")
      !return
   case(2:)
      call ctx%message("Found "//to_string(nfrag)//" fragments:")
      do ifr = 1, nfrag
         call ctx%message("Fragment "//to_string(ifr)//":  "//format_list(fragment == ifr))
      end do
      call ctx%message("")
   end select

   if (any(mol%periodic)) then
      invlat = matinv_3x3(mol%lattice)
      abc = matmul(invlat, mol%xyz)
      call unfold_fragment(abc, wbo, input%thr)
      mol%xyz(:, :) = matmul(mol%lattice, abc)
   end if

   allocate(mfrag(nfrag))
   call get_structure_fragment(mfrag(1), mol, fragment == 1)
   call write_structure(mfrag(1), "final.xyz", error)

end subroutine get_jab


!> Unpack coefficients from a fragment orbital space in the full orbital space
subroutine unpack_coeff(full_bas, frag_bas, full, frag, mask)
   !> Basis set information for full system
   type(basis_type), intent(in) :: full_bas
   !> Basis set information for fragment
   type(basis_type), intent(in) :: frag_bas
   !> Quantity in full orbital space
   real(wp), intent(out) :: full(:)
   !> Quantity in fragment orbital space
   real(wp), intent(in) :: frag(:)
   !> Atom resolved mask for this fragment
   logical, intent(in) :: mask(:)

   integer :: iat, jat, ish, ii, jj, iao, jao, nao

   jat = 0
   do iat = 1, size(mask)
      if (mask(iat)) then
         jat = jat + 1
         ii = full_bas%ish_at(iat)
         jj = frag_bas%ish_at(jat)
         do ish = 1, full_bas%nsh_at(iat)
            iao = full_bas%iao_sh(ii+ish)
            jao = frag_bas%iao_sh(jj+ish)
            nao = full_bas%nao_sh(ii+ish)
            full(iao+1:iao+nao) = frag(jao+1:jao+nao)
         end do
      else
         ii = full_bas%ish_at(iat)
         do ish = 1, full_bas%nsh_at(iat)
            iao = full_bas%iao_sh(ii+ish)
            nao = full_bas%nao_sh(ii+ish)
            full(iao+1:iao+nao) = 0.0_wp
         end do
      end if
   end do
end subroutine unpack_coeff


!> Extract a the fragment structure from the full structure
!>
!> Todo: This routine can currently only create neutral fragments,
!>       charge and spin information is not partition or transferred to fragments
subroutine get_structure_fragment(frag, mol, mask)
   !> Molecular structure data of the fragment
   type(structure_type), intent(out) :: frag
   !> Molecular structure data of the full system
   type(structure_type), intent(in) :: mol
   !> Atom resolved mask for this fragment
   logical, intent(in) :: mask(:)

   integer :: nat
   integer, allocatable :: num(:)
   real(wp), allocatable :: xyz(:, :)

   nat = count(mask)
   num = pack(mol%num(mol%id), mask)
   xyz = reshape(pack(mol%xyz, spread(mask, 1, 3)), [3, nat])
   call new(frag, num, xyz)
end subroutine get_structure_fragment


end module frag
