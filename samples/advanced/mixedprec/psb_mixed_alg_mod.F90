module psb_mixed_alg_mod
  use psb_base_mod
  use psb_prec_mod
  use psb_d_linsolve_conv_mod
  use psb_linsolve_mod
  use psb_mixed_support_mod

  implicit none

contains

  subroutine psb_cg_mixed(a, prec, b, x, eps, desc_a, info, &
       & itmax, iter, err, itrace, istop, cond)
    implicit none
    type(psb_dspmat_type), intent(in)    :: a
    Type(psb_desc_type), Intent(in)      :: desc_a
    class(psb_sprec_type), intent(inout) :: prec
    type(psb_d_vect_type), Intent(inout) :: b
    type(psb_d_vect_type), Intent(inout) :: x
    Real(psb_dpk_), Intent(in)           :: eps
    integer(psb_ipk_), intent(out)                 :: info
    integer(psb_ipk_), Optional, Intent(in)        :: itmax, itrace, istop
    integer(psb_ipk_), Optional, Intent(out)       :: iter
    Real(psb_dpk_), Optional, Intent(out) :: err, cond
		! =   Local data
    real(psb_dpk_), allocatable, target   :: aux(:), td(:), tu(:), eig(:), ewrk(:)
    integer(psb_mpk_), allocatable :: ibl(:), ispl(:), iwrk(:)
    type(psb_d_vect_type), allocatable, target :: wwrk(:)
    type(psb_d_vect_type), pointer  :: q, p, r, z, w
		! =   Single precision workspace for the preconditioner application
    real(psb_spk_), allocatable, target   :: saux(:)
    type(psb_s_vect_type)                 :: r_single, z_single
    real(psb_dpk_)   :: alpha, beta, rho, rho_old, sigma, alpha_old, beta_old
    integer(psb_ipk_) :: itmax_, istop_, naux, it, itx, itrace_, &
         &  n_col, n_row, err_act, ieg, nspl, istebz
    integer(psb_lpk_) :: mglob
    integer(psb_ipk_) :: debug_level, debug_unit
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: np, me
    real(psb_dpk_)     :: derr
    type(psb_itconv_type)       :: stopdat
    logical                     :: do_cond
    character(len=20)           :: name
    character(len=*), parameter :: methdname='CG-MIXED'

    info = psb_success_
    name = 'psb_cg_mixed'
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()

    ctxt = desc_a%get_context()

    call psb_info(ctxt, me, np)
    if (.not.allocated(b%v)) then
      info = psb_err_invalid_vect_state_
      call psb_errpush(info, name)
      goto 9999
    endif
    if (.not.allocated(x%v)) then
      info = psb_err_invalid_vect_state_
      call psb_errpush(info, name)
      goto 9999
    endif


    mglob = desc_a%get_global_rows()
    n_row = desc_a%get_local_rows()
    n_col = desc_a%get_local_cols()


    if (present(istop)) then
      istop_ = istop
    else
      istop_ = psb_get_istop_default()
    endif
    if (.not.psb_is_valid_istop(istop_)) then
      info=psb_err_invalid_istop_
      err=info
      call psb_errpush(info, name, i_err=(/istop_/))
      goto 9999
    end if
    !
    !  istop_ = 1:  normwise backward error, infinity norm
    !  istop_ = 2:  ||r||/||b||   norm 2
    !
    select case(istop_)
    case(psb_istop_ani_, psb_istop_bn2_, &
         & psb_istop_rn2_abs_, psb_istop_rrn2_)
      ! nothing needed
    case default
      ! should never get here
      info=psb_err_internal_error_
      err=info
      call psb_errpush(info, name, a_err="invalid istop_")
      goto 9999
    end select

    call psb_chkvect(mglob, lone, x%get_nrows(), lone, lone, desc_a, info)
    if (info == psb_success_)&
         & call psb_chkvect(mglob, lone, b%get_nrows(), lone, lone, desc_a, info)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      call psb_errpush(info, name, a_err='psb_chkvect on X/B')
      goto 9999
    end if

    naux=4*n_col
    allocate(aux(naux), saux(naux), stat=info)
    if (info == psb_success_) call psb_geall(wwrk, desc_a, info, n=5_psb_ipk_)
    if (info == psb_success_) call psb_geasb(wwrk, desc_a, info, mold=x%v, scratch=.true.)
    !
    !  The descriptor is precision independent: the single precision
    !  vectors handed to the preconditioner are allocated on the very
    !  same descriptor as the double precision ones.
    !
    if (info == psb_success_) call psb_geall(r_single, desc_a, info)
    if (info == psb_success_) call psb_geasb(r_single, desc_a, info, scratch=.true.)
    if (info == psb_success_) call psb_geall(z_single, desc_a, info)
    if (info == psb_success_) call psb_geasb(z_single, desc_a, info, scratch=.true.)
    if (info /= psb_success_) then
      info=psb_err_from_subroutine_non_
      call psb_errpush(info, name)
      goto 9999
    end if

    p  => wwrk(1)
    q  => wwrk(2)
    r  => wwrk(3)
    z  => wwrk(4)
    w  => wwrk(5)
    !
    !  The vectors crossing the precision boundary are converted in
    !  full, halo included; the halo of R is only written by the
    !  preconditioner, hence the scratch storage is cleaned here to
    !  avoid converting uninitialized values.
    !
    call r%zero()
    call z_single%zero()


    if (present(itmax)) then
      itmax_ = itmax
    else
      itmax_ = 1000
    endif

    if (present(itrace)) then
      itrace_ = itrace
    else
      itrace_ = 0
    end if

    do_cond=present(cond)
    if (do_cond) then
      istebz = 0
      allocate(td(itmax_), tu(itmax_), eig(itmax_), &
           & ibl(itmax_), ispl(itmax_), iwrk(3*itmax_), ewrk(4*itmax_), &
           & stat=info)
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_non_
        call psb_errpush(info, name)
        goto 9999
      end if
    end if
    itx=0
    alpha = dzero

    restart: do
			! =
			! =    r0 = b-Ax0
			! =
      if (itx>= itmax_) exit restart

      it = 0
      call psb_geaxpby(done, b, dzero, r, desc_a, info)
      if (info == psb_success_) call psb_spmm(-done, a, x, done, r, desc_a, info, work=aux)
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_non_
        call psb_errpush(info, name)
        goto 9999
      end if

      rho = dzero

      call psb_init_conv(methdname, istop_, itrace_, itmax_, a, x, b, eps, desc_a, stopdat, info)
      if (info /= psb_success_) Then
        call psb_errpush(psb_err_from_subroutine_non_, name)
        goto 9999
      End If

      iteration:  do

        it   = it + 1
        itx = itx + 1

        !
        !  This is the only step performed in single precision:
        !  convert the residual, apply the single precision
        !  preconditioner, then convert the result back to double
        !  and carry on with the double precision recurrence.
        !
        call psb_d2s_vect(r, r_single, info)
        if (info == psb_success_) call prec%apply(r_single, z_single, desc_a, info, work=saux)
        if (info == psb_success_) call psb_s2d_vect(z_single, z, info, mold=x%v)
        if (info /= psb_success_) then
          info=psb_err_from_subroutine_non_
          call psb_errpush(info, name, a_err='single precision preconditioner')
          goto 9999
        end if

        rho_old = rho
        rho     = psb_gedot(r, z, desc_a, info)

        if (it == 1) then
          call psb_geaxpby(done, z, dzero, p, desc_a, info)
        else
          if (rho_old == dzero) then
            if (debug_level >= psb_debug_ext_)&
                 & write(debug_unit, *) me, ' ', trim(name), &
                 & ': CG Iteration breakdown rho'
            exit iteration
          endif
          beta = rho/rho_old
          call psb_geaxpby(done, z, beta, p, desc_a, info)
        end if

        call psb_spmm(done, a, p, dzero, q, desc_a, info, work=aux)
        sigma = psb_gedot(p, q, desc_a, info)
        if (sigma == dzero) then
          if (debug_level >= psb_debug_ext_)&
               & write(debug_unit, *) me, ' ', trim(name), &
               & ': CG Iteration breakdown sigma'
          exit iteration
        endif
        alpha_old = alpha
        alpha = rho/sigma

        if (do_cond) then
          istebz = istebz + 1
          if (istebz == 1) then
            td(istebz) = done/alpha
          else
            td(istebz) = done/alpha + beta/alpha_old
            tu(istebz-1) = sqrt(beta)/alpha_old
          end if
        end if


        call psb_geaxpby(alpha, p, done, x, desc_a, info)
        call psb_geaxpby(-alpha, q, done, r, desc_a, info)

        if (psb_check_conv(methdname, itx, x, r, desc_a, stopdat, info)) exit restart
        if (info /= psb_success_) Then
          call psb_errpush(psb_err_from_subroutine_non_, name)
          goto 9999
        End If

      end do iteration
    end do restart
    if (do_cond) then
      if (me == psb_root_) then
				! #if defined(PSB_HAVE_LAPACK)
        call dstebz('A', 'E', istebz, dzero, dzero, 0, 0, -done, td, tu, &
             & ieg, nspl, eig, ibl, ispl, ewrk, iwrk, info)
        if (info < 0) then
          call psb_errpush(psb_err_from_subroutine_ai_, name, &
               & a_err='dstebz', i_err=(/info/))
          info=psb_err_from_subroutine_ai_
          goto 9999
        end if
        cond = eig(ieg)/eig(1)
				! #else
				! 				cond = dzero
				! #endif
        info=psb_success_
      end if
      call psb_bcast(ctxt, cond)
    end if


    call psb_end_conv(methdname, itx, desc_a, stopdat, info, derr, iter)
    if (present(err)) err = derr

    if (info == psb_success_) call psb_gefree(wwrk, desc_a, info)
    if (info == psb_success_) call psb_gefree(r_single, desc_a, info)
    if (info == psb_success_) call psb_gefree(z_single, desc_a, info)
    if (info == psb_success_) deallocate(aux, saux, stat=info)
    if (info /= psb_success_) then
      call psb_errpush(info, name)
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

  9999 call psb_error_handler(err_act)
    return
  end subroutine psb_cg_mixed

	subroutine psb_mixed_pMPK_split(spmat, prec, vec_in, Z, Q, s, vec_single, desc, info,  & 
                            & base_type, alpha, beta, gamma, mvec_temp, farr_temp, sarr_temp)
    use psb_base_mod
    use psb_prec_mod
    implicit none
    type(psb_dspmat_type), intent(in)           :: spmat
    class(psb_sprec_type), intent(inout)        :: prec
    type(psb_d_vect_type), intent(inout)        :: vec_in
    type(psb_d_multivect_type), intent(inout)   :: Z, Q
    integer(psb_ipk_), intent(in)               :: s
    type(psb_s_vect_type), intent(inout)        :: vec_single
    type(psb_desc_type), intent(in)             :: desc
    integer(psb_ipk_), intent(out)              :: info
    character, optional, intent(in)                             :: base_type
    real(psb_dpk_), optional, intent(in)                        :: alpha, beta, gamma
    type(psb_d_multivect_type), optional, target, intent(inout) :: mvec_temp
    real(psb_dpk_), optional, target, intent(inout)             :: farr_temp(:)
    real(psb_spk_), optional, target, intent(inout)             :: sarr_temp(:)

    integer(psb_ipk_) :: err_act
    real(psb_dpk_)    :: gamma_
    logical           :: save_r
    character         :: base_type_
    character(len=20) :: name = "psb_mixed_pMPK"
    real(psb_dpk_), pointer :: aux(:)
		real(psb_spk_), pointer :: saux(:)

    call psb_erractionsave(err_act)
    info = psb_success_

    !Check s value
    if(s <= 0)  then
        info = psb_err_iarg_invalid_value_
        call psb_errpush(info, name)
        goto 9999
    endif

    !Dimension checks. Z = n*s Q = n*(s+1) or n*s
    if((Z%get_ncols() /= s) .or. ((Q%get_ncols() /= s) .and. (Q%get_ncols() /= s + 1))) then
        info = psb_err_invalid_mvect_size_
        call psb_errpush(info, name)
        goto 9999
    endif

    !Select if save r as the first vector of Q using is dimension
    save_r = (Q%get_ncols() == s + 1);

    if(present(farr_temp)) then
        !TO DO: check dimension
        aux => farr_temp
    else
        allocate(aux(4*desc%get_local_cols()))
    end if

		if(present(sarr_temp)) then
        !TO DO: check dimension
        saux => sarr_temp
    else
        allocate(saux(4*desc%get_local_cols()))
    end if

    !Check type base
    if(present(base_type)) then 
        base_type_ = psb_toupper(base_type)
    else
        base_type_ = "M"
    endif

    select case(base_type_)
        case ("M")
            call psb_mixed_pMPK_split_monomial() 
        
        case ("C")
            ! Check presence of Chebyshev parameters alfa and beta (error or default values?)
            if((.not. present(alpha)) .or. (.not. present(beta))) then
                info = psb_err_invalid_args_combination_
                call psb_errpush(info, name)
                goto 9999
            endif

            !Chebyshev parameters gamma defaulted to 1 if not present
            if(present(gamma)) then
                gamma_ = gamma
            else
                gamma_ = done
            endif

            call psb_mixed_pMPK_split_chebyshev()
        
        case default
            info = psb_err_invalid_input_
            call psb_errpush(info, name)
            goto 9999
    end select

    if(.not. present(farr_temp)) deallocate(aux)

    call psb_erractionrestore(err_act)
    return

	9999 call psb_error_handler(err_act)
			return

	contains
			subroutine psb_mixed_pMPK_split_monomial()
					implicit none
					character(len=20) :: name = "psb_mixed_pMPK_monomial"
					integer(psb_ipk_) :: i, idx_Z, idx_Q

					! Copy r in the first column of Q
					idx_Q = 1
					call psb_geaxpby(done, vec_in, dzero, Q, idx_Q, desc, info)


					! First iteration
					idx_Z = 1
					
					call psb_d2s_mvect_col(Q, idx_Q, vec_single, desc, info)
					call prec%apply(vec_single, desc, info, work = saux)
					call psb_s2d_mvect_col(vec_single, Z, idx_Z, desc, info)

					idx_Q = merge(idx_Q + 1, idx_Q, save_r) !If not save_r overwrite it with the first vector
					call psb_spmm(done, spmat, Z, idx_Z, dzero, Q, idx_Q, desc, info, work = aux)
					
					if(s == 1) return

					do i = 2, s
							idx_Z = idx_Z + 1
							call psb_d2s_mvect_col(Q, idx_Q, vec_single, desc, info)
							call prec%apply(vec_single, desc, info, work = saux)
							call psb_s2d_mvect_col(vec_single, Z, idx_Z, desc, info)
							
							idx_Q = idx_Q + 1
							call psb_spmm(done, spmat, Z, idx_Z, dzero, Q, idx_Q, desc, info, work = aux)
					end do
			end subroutine psb_mixed_pMPK_split_monomial

			subroutine psb_mixed_pMPK_split_chebyshev()
					implicit none
					character(len=20) :: name = "psb_mixed_pMPK_chebyshev"
					integer(psb_ipk_) :: i, idx_Z, idx_Q, ind_tmp
					type(psb_d_multivect_type), pointer :: mvec_tmp

					if(present(mvec_temp)) then
							! TO DO: check dimensions
							mvec_tmp => mvec_temp
					else
							call psb_geall(mvec_tmp, desc, info, n = 3)    ! I need only tre stored temp vectors
							call psb_geasb(mvec_tmp, desc, info)
					endif
					
					! Copy r in the first column of Q
					idx_Q = 1
					call psb_geaxpby(done, vec_in, dzero, Q, idx_Q, desc, info)
					
					! Inizialize first column of temp multivector
					ind_tmp = 1
					call psb_geaxpby(done, vec_in, dzero, mvec_tmp, ind_tmp, desc, info)

					! First iteration
					idx_Z = 1 
					! call prec%apply(mvec_tmp, ind_tmp, Z, idx_Z, desc, info, work = aux)
					call psb_d2s_mvect_col(mvec_tmp, ind_tmp, vec_single, desc, info)
					call prec%apply(vec_single, desc, info, work = saux)
					call psb_s2d_mvect_col(vec_single, Z, idx_Z, desc, info)

					idx_Q = merge(idx_Q + 1, idx_Q, save_r) !If not save_r overwrite it with the first vector
					call psb_spmm(done, spmat, Z, idx_Z, dzero, Q, idx_Q, desc, info, work = aux)

					!Check first early exit
					if(s == 1) goto 9998

					! First update of temp multivector + second iteration
					ind_tmp = ind_tmp + 1
					call psb_geaxpby(alpha, Q, idx_Q, -beta, mvec_tmp, ind_tmp - 1, dzero, mvec_tmp, ind_tmp, desc, info)
					idx_Z = idx_Z + 1

					!call prec%apply(mvec_tmp, ind_tmp, Z, idx_Z, desc, info, work = aux)
					call psb_d2s_mvect_col(mvec_tmp, ind_tmp, vec_single, desc, info)
					call prec%apply(vec_single, desc, info, work = saux)
					call psb_s2d_mvect_col(vec_single, Z, idx_Z, desc, info)
					
					idx_Q = idx_Q + 1
					call psb_spmm(done, spmat, Z, idx_Z, dzero, Q, idx_Q, desc, info, work = aux)

					!Check second early exit
					if(s == 2) goto 9998
					
					! Loop for s > 2
					do i = 3, s
							ind_tmp = modulo(ind_tmp, 3) + 1
							call psb_geaxpby(2*alpha, Q, idx_Q, -2*beta, mvec_tmp, modulo(ind_tmp - 2, 3) + 1, &
																	& -gamma_, mvec_tmp, modulo(ind_tmp - 3, 3) + 1, mvec_tmp, ind_tmp, desc, info)
							idx_Z = idx_Z + 1
							!call prec%apply(mvec_tmp, ind_tmp, Z, idx_Z, desc, info, work = aux)
							call psb_d2s_mvect_col(mvec_tmp, ind_tmp, vec_single, desc, info)
							call prec%apply(vec_single, desc, info, work = saux)
							call psb_s2d_mvect_col(vec_single, Z, idx_Z, desc, info)

							idx_Q = idx_Q + 1
							call psb_spmm(done, spmat, Z, idx_Z, dzero, Q, idx_Q, desc, info, work = aux)
					end do
			
			9998 if(.not. present(mvec_temp)) call psb_gefree(mvec_tmp, desc, info)
					return
			end subroutine psb_mixed_pMPK_split_chebyshev
	end subroutine psb_mixed_pMPK_split

	subroutine psb_dscg_mixed(a, prec, b, x, s, eps, desc_a, info, &
                      & itmax, iter, err, itrace, istop, &
                      & base_type, eigext, Gram_solver, FGS_sweeps)
		use psb_base_mod
		use psb_prec_mod
		use psb_d_linsolve_conv_mod
		use psb_linsolve_mod
		implicit none
		type(psb_dspmat_type), intent(in)     :: a
		class(psb_sprec_type), intent(inout)  :: prec
		type(psb_d_vect_type), intent(inout)  :: b, x
		integer(psb_ipk_), intent(in)         :: s
		real(psb_dpk_), intent(in)            :: eps
		type(psb_desc_type), intent(in)       :: desc_a
		integer(psb_ipk_), intent(out)        :: info
		integer(psb_ipk_), optional, intent(in)   :: itmax, itrace, istop
		integer(psb_ipk_), optional, intent(out)  :: iter
		real(psb_dpk_), optional, intent(out)     :: err
		character, optional, intent(in)           :: base_type
		real(psb_dpk_), optional, intent(in)      :: eigext(2)
		character(len=3), optional, intent(in)    :: Gram_solver
		integer(psb_ipk_), optional, intent(in)   :: FGS_sweeps

		! Local vars
		type(psb_ctxt_type) :: ctxt
		integer(psb_ipk_)   :: istop_, itmax_, itrace_, FGS_sweeps_
		character(len=3)    :: base_type_, Gram_solver_
		integer(psb_ipk_)   :: err_act, np, me, debug_level, debug_unit, n_col, n_row
		integer(psb_lpk_)   :: mglob
		character(len=20)           :: name = 'psb_dscg'
		character(len=*), parameter :: methdbasename = 'sStepCG'
		character(len=20)           :: methdfullname

		real(psb_dpk_), allocatable     :: alpha(:), beta(:, :), W(:, :), temp_fa(:, :)
		integer(psb_ipk_), allocatable  :: pW(:)
		type(psb_d_vect_type)           :: r  
		type(psb_s_vect_type)           :: v_single  
		type(psb_d_multivect_type)      :: Z, Q, P, V, temp_mv
		real(psb_dpk_)                  :: cheb_coeff(3)
		integer(psb_ipk_)               :: itidx
		
		type(psb_itconv_type)         :: stopdat
		real(psb_dpk_)                :: derr

		type(psb_d_multivect_type), target  :: aux_mv
		real(psb_dpk_), allocatable, target :: aux_fa(:)
		real(psb_spk_), allocatable, target :: aux_sa(:)

		character(len=3), parameter   :: forwardGS = "FGS"
		character(len=3), parameter   :: lapackLU = "LLU"
		character(len=3), parameter   :: lapackCC = "LCC"

		info = psb_success_
		call psb_erractionsave(err_act)
		
		debug_unit  = psb_get_debug_unit()
		debug_level = psb_get_debug_level()

		ctxt = desc_a%get_context()
		call psb_info(ctxt, me, np)
		if(debug_level >= psb_debug_ext_) &
				& write(debug_unit, *) me, ' ', trim(name), ': from psb_info ', np
		
		if(s < 1) then
			info = psb_err_iarg_invalid_value_
			call psb_errpush(info, name)
			goto 9999
		endif

		write(methdfullname, '(A, "(", I0, ")")') methdbasename, s

		if((.not. allocated(b%v)) .or. (.not.allocated(x%v))) then 
			info = psb_err_invalid_vect_state_
			call psb_errpush(info, name)
			goto 9999
		endif

		istop_ = 2
		if(present(istop)) istop_ = istop
		
		!  ISTOP_ = 1:  Normwise backward error, infinity norm 
		!  ISTOP_ = 2:  ||r||/||b||, 2-norm 
		if((istop_ < 1 ) .or. (istop_ > 2 )) then
			info = psb_err_invalid_istop_
			err = info
			call psb_errpush(info, name, i_err = (/istop_/))
			goto 9999
		endif

		itmax_ = 1000
		if(present(itmax)) itmax_ = itmax

		itrace_ = 0
		if(present(itrace)) itrace_ = itrace

		base_type_ = "C"
		if(present(base_type)) base_type_ = base_type

		Gram_solver_ = lapackCC
		if(present(Gram_solver)) Gram_solver_ = Gram_solver
		
		FGS_sweeps_ = 30
		if(present(FGS_sweeps)) FGS_sweeps_ = FGS_sweeps

		mglob = desc_a%get_global_rows()
		n_row = desc_a%get_local_rows()
		n_col = desc_a%get_local_cols()

		call psb_chkvect(mglob, lone, x%get_nrows(), lone, lone, desc_a, info)
		if(info /= psb_success_) then
			info = psb_err_from_subroutine_
			call psb_errpush(info, name,  a_err = 'psb_chkvect on x')
			goto 9999
		end if

		call psb_chkvect(mglob, lone, b%get_nrows(), lone, lone, desc_a, info)
		if(info /= psb_success_) then
			info = psb_err_from_subroutine_    
			call psb_errpush(info, name, a_err='psb_chkvect on b')
			goto 9999
		end if

		!Allocate and assembly data structure
		allocate(alpha(s), beta(s, s), W(s, s), pW(s), temp_fa(s, s + 1), aux_fa(4*n_col), aux_sa(4*n_col), stat = info)
		if(info == psb_success_) call psb_geall(r, desc_a, info)
		if(info == psb_success_) call psb_geall(v_single, desc_a, info)
		if(info == psb_success_) call psb_geall(Z, desc_a, info, n = s)
		if(info == psb_success_) call psb_geall(Q, desc_a, info, n = s)
		if(info == psb_success_) call psb_geall(P, desc_a, info, n = s)
		if(info == psb_success_) call psb_geall(V, desc_a, info, n = s)
		if(info == psb_success_) call psb_geall(temp_mv, desc_a, info, n = s)
		if(info == psb_success_) call psb_geall(aux_mv, desc_a, info, n = 3)
		if(info == psb_success_) call psb_geasb(r, desc_a, info)
		if(info == psb_success_) call psb_geasb(v_single, desc_a, info)
		if(info == psb_success_) call psb_geasb(Z, desc_a, info)
		if(info == psb_success_) call psb_geasb(Q, desc_a, info)
		if(info == psb_success_) call psb_geasb(P, desc_a, info)
		if(info == psb_success_) call psb_geasb(V, desc_a, info)
		if(info == psb_success_) call psb_geasb(temp_mv, desc_a, info)
		if(info == psb_success_) call psb_geasb(aux_mv, desc_a, info)

		if(info /= psb_success_) then 
			info = psb_err_from_subroutine_ 
			call psb_errpush(info, name)
			goto 9999
		end if

		! First residual calculation
		call psb_geaxpby(done, b, dzero, r, desc_a, info)
		if(info == psb_success_) call psb_spmm(-done, a, x, done, r, desc_a, info)
		if(info /= psb_success_) then 
			info = psb_err_from_subroutine_ 
			call psb_errpush(info, name)
			goto 9999
		end if

		! Init converence
		call psb_init_conv(methdfullname, istop_, itrace_, itmax_, a, x, b, eps, desc_a, stopdat, info)
		if(info /= psb_success_) then 
			info = psb_err_from_subroutine_ 
			call psb_errpush(info, name)
			goto 9999
		end if

		! check convergence here?

		! Chebyshev coefficient calculation
		select case (base_type_)
			case("M")
				cheb_coeff = dzero
			case("C")
				cheb_coeff = psb_d_chebyshev_coefficients(a, prec, desc_a, info, eigext)
			case default
				info = psb_err_invalid_input_ 
				call psb_errpush(info, name)
				goto 9999
		end select

		if(info /= psb_success_) then 
			info = psb_err_from_subroutine_ 
			call psb_errpush(info, name)
			goto 9999
		end if
		
		! First matrix power kernel
		call psb_mixed_pMPK_split(a, prec, r, P, V, s, v_single, desc_a, info, base_type = base_type_, &
										& alpha = cheb_coeff(1), beta = cheb_coeff(2), gamma = cheb_coeff(3), &
										& mvec_temp = aux_mv, farr_temp = aux_fa, sarr_temp = aux_sa)
		if(info /= psb_success_) then 
			info = psb_err_from_subroutine_ 
			call psb_errpush(info, name)
			goto 9999
		end if

		! Loop until convergence (or maxiter)
		do itidx = 1, itmax_
			! Compute matrix W and rhs for alpha
			call psb_gedots(P, V, temp_fa(:, 1 : s), desc_a, info, global = .false.)
			call psb_gedots(P, r, temp_fa(:, s + 1), desc_a, info, global = .false.)
			call psb_sum(desc_a%get_context(), temp_fa)
			W = temp_fa(:, 1 : s)
			alpha = temp_fa(:, s + 1)

			! Factor matrix W (if soving with LU or Cholesky factorization)
			if(Gram_solver_ == lapackLU) call dgetrf(s, s, W, s, pW, info)
			if(Gram_solver_ == lapackCC) call dpotrf('L', s, W, s, info)

			! Solve for alpha
			select case(Gram_solver_)
				case(forwardGS);  call inner_solver_fgs_1D(W, alpha, FGS_sweeps_)
				case(lapackLU);   call dgetrs('N', s, 1, W, s, pW, alpha, s, info)
				case(lapackCC);   call dpotrs('L', s, 1, W, s, alpha, s, info)
				case default
					info = psb_err_invalid_input_ 
					call psb_errpush(info, name)
					goto 9999
			end select

			! Update solution and residual
			call psb_geaxpby(P, alpha, x, desc_a, info, upd_flag = .true.)
			call psb_geaxpby(V, -alpha, r, desc_a, info, upd_flag = .true.)

			! Check convergence. 
			if(psb_check_conv(methdfullname, itidx, x, r, desc_a, stopdat, info)) exit
			
			! Matrix power kernel
			call psb_mixed_pMPK_split(a, prec, r, Z, Q, s, v_single, desc_a, info, base_type = base_type_, &
										& alpha = cheb_coeff(1), beta = cheb_coeff(2), gamma = cheb_coeff(3), &
										& mvec_temp = aux_mv, farr_temp = aux_fa, sarr_temp = aux_sa)

			! Compute rhs for beta
			call psb_gedots(P, Q, beta, desc_a, info, .true.)
			beta = -beta;

			! Solve for beta
			select case(Gram_solver_)
				case(forwardGS);  call inner_solver_fgs_2D(W, beta, FGS_sweeps_)
				case(lapackLU);   call dgetrs('N', s, s, W, s, pW, beta, s, info)
				case(lapackCC);   call dpotrs('L', s, s, W, s, beta, s, info)
				case default
					info = psb_err_invalid_input_ 
					call psb_errpush(info, name)
					goto 9999
			end select

			! Update P and V. Use of temp_mv in needed because internal dgemm constraint
			call psb_geaxpby(P, beta, temp_mv, desc_a, info, .false.)
			call psb_geaxpby(done, Z, done, temp_mv, P, desc_a, info)
			call psb_geaxpby(V, beta, temp_mv, desc_a, info, .false.)
			call psb_geaxpby(done, Q, done, temp_mv, V, desc_a, info)
		end do

		call psb_end_conv(methdfullname, itidx, desc_a, stopdat, info, derr, iter)
		if(present(err)) err = derr
		if(present(iter)) iter = iter * s

		if(info == psb_success_) call psb_gefree(r, desc_a, info)
		if(info == psb_success_) call psb_gefree(v_single, desc_a, info)
		if(info == psb_success_) call psb_gefree(Z, desc_a, info)
		if(info == psb_success_) call psb_gefree(Q, desc_a, info)
		if(info == psb_success_) call psb_gefree(P, desc_a, info)
		if(info == psb_success_) call psb_gefree(V, desc_a, info)
		if(info == psb_success_) call psb_gefree(temp_mv, desc_a, info)
		if(info == psb_success_) call psb_gefree(aux_mv, desc_a, info)

		if(info == psb_success_) deallocate(alpha, beta, W, pW, temp_fa, aux_fa, aux_sa, stat = info)
		if(info /= psb_success_) then
			call psb_errpush(info,name)
			goto 9999
		end if

		call psb_erractionrestore(err_act)
		return

	9999 call psb_error_handler(err_act)
		return
	contains
		function psb_d_chebyshev_coefficients(a, prec, desc, info, eigext) result(coeff)
			type(psb_dspmat_type), intent(in)     :: a
			class(psb_sprec_type), intent(inout)  :: prec
			type(psb_desc_type), intent(in)       :: desc
			integer(psb_ipk_), intent(out)        :: info
			real(psb_dpk_), optional, intent(in)  :: eigext(2)
			
			real(psb_dpk_) :: lambda_max, lambda_min
			real(psb_dpk_) :: coeff(3)
	
			info = psb_success_

			if(present(eigext)) then
				lambda_min = eigext(1)
				lambda_max = eigext(2)
			else 
				call psb_mixed_powermethod(a, prec, lambda_max, desc, info)
				lambda_min = dzero
			end if

			if(info /= psb_success_) then 
				info = psb_err_from_subroutine_
				return
			end if

			if(lambda_min < 0) then 
				info = psb_err_fatal_ ! TODO: set the correct error
				return
			end if

			if(lambda_max < lambda_min) then 
				info = psb_err_fatal_ ! TODO: set the correct error
				return
			end if

			coeff(1) = 2_psb_dpk_ / (lambda_max - lambda_min)
			coeff(2) = (lambda_max + lambda_min) / (lambda_max - lambda_min)
			coeff(3) = done
		end function psb_d_chebyshev_coefficients

		subroutine inner_solver_fgs_1D(M, rhs, num_iter)
			real(psb_dpk_), intent(in)    :: M(:, :)
			real(psb_dpk_), intent(inout) :: rhs(:)
			integer(psb_ipk_), intent(in) :: num_iter

			integer(psb_ipk_) :: iter_idx, i, j, n
			real(psb_dpk_)    :: sol(size(rhs))

			sol = dzero
			n = size(sol)
			do iter_idx = 1, num_iter
				do i = 1, n
					sol(i) = rhs(i)
					do j = 1, n
						if(j /= i) sol(i) = sol(i) - M(i, j) * sol(j)
					end do
					sol(i) = sol(i) / M(i, i)
				end do
			end do
			rhs = sol
		end subroutine inner_solver_fgs_1D

		subroutine inner_solver_fgs_2D(M, rhs, num_iter)
			real(psb_dpk_), intent(in)    :: M(:, :)
			real(psb_dpk_), intent(inout) :: rhs(:, :)
			integer(psb_ipk_), intent(in) :: num_iter

			integer(psb_ipk_) :: iter_idx, i, j, n
			real(psb_dpk_)    :: sol(size(rhs, 1), size(rhs, 2))

			sol = dzero
			n = size(sol, 1)
			do iter_idx = 1, num_iter
				do i = 1, n
					sol(i, :) = rhs(i, :)
					do j = 1, n
						if(j /= i) sol(i, :) = sol(i, :) - M(i, j) * sol(j, :)
					end do
					sol(i, :) = sol(i, :) / M(i, i)
				end do
			end do
			rhs = sol
		end subroutine inner_solver_fgs_2D
	end subroutine psb_dscg_mixed

	subroutine psb_dscg2_mixed(a, prec, b, x, s, eps, desc_a, info, &
													& itmax, iter, err, itrace, istop, &
													& base_type, eigext, Gram_solver, FGS_sweeps)
		use psb_base_mod
		use psb_prec_mod
		use psb_d_linsolve_conv_mod
		use psb_linsolve_mod

		implicit none
		type(psb_dspmat_type), intent(in)     :: a
		class(psb_sprec_type), intent(inout)  :: prec
		type(psb_d_vect_type), intent(inout)  :: b, x
		integer(psb_ipk_), intent(in)         :: s
		real(psb_dpk_), intent(in)            :: eps
		type(psb_desc_type), intent(in)       :: desc_a
		integer(psb_ipk_), intent(out)        :: info
		integer(psb_ipk_), optional, intent(in)   :: itmax, itrace, istop
		integer(psb_ipk_), optional, intent(out)  :: iter
		real(psb_dpk_), optional, intent(out)     :: err
		character, optional, intent(in)           :: base_type
		real(psb_dpk_), optional, intent(in)      :: eigext(2)
		character(len=3), optional, intent(in)    :: Gram_solver
		integer(psb_ipk_), optional, intent(in)   :: FGS_sweeps

		! Local vars
		type(psb_ctxt_type) :: ctxt
		integer(psb_ipk_)   :: istop_, itmax_, itrace_, FGS_sweeps_
		character(len=3)    :: base_type_, Gram_solver_
		integer(psb_ipk_)   :: err_act, np, me, debug_level, debug_unit, &
														& n_col, n_row
		integer(psb_lpk_)   :: mglob
		character(len=20)           :: name = 'psb_dscg'
		character(len=*), parameter :: methdbasename = 'sStepCGv2'
		character(len=20)           :: methdfullname

		real(psb_dpk_), allocatable     :: alpha(:), beta(:, :), W(:, :), temp_fa(:, :), B2(:, :), c0(:)
		integer(psb_ipk_), allocatable  :: pW(:)
		type(psb_d_vect_type)           :: r  
		type(psb_d_multivect_type)      :: Z, Q, P, V, temp_mv
		real(psb_dpk_)                  :: cheb_coeff(3)
		integer(psb_ipk_)               :: itidx
		
		type(psb_itconv_type)         :: stopdat
		real(psb_dpk_)                :: derr 

		type(psb_d_multivect_type), target  :: aux_mv
		real(psb_dpk_), allocatable, target :: aux_fa(:)

		character(len=3), parameter   :: forwardGS = "FGS"
		character(len=3), parameter   :: lapackLU = "LLU"
		character(len=3), parameter   :: lapackCC = "LCC"

		info = psb_success_
		call psb_erractionsave(err_act)
		
		debug_unit  = psb_get_debug_unit()
		debug_level = psb_get_debug_level()

		ctxt = desc_a%get_context()
		call psb_info(ctxt, me, np)
		if(debug_level >= psb_debug_ext_) &
				& write(debug_unit, *) me, ' ', trim(name), ': from psb_info ', np
		
		if(s < 1) then
			info = psb_err_iarg_invalid_value_
			call psb_errpush(info, name)
			goto 9999
		endif

		write(methdfullname, '(A, "(", I0, ")")') methdbasename, s

		if((.not. allocated(b%v)) .or. (.not.allocated(x%v))) then 
			info = psb_err_invalid_vect_state_
			call psb_errpush(info, name)
			goto 9999
		endif

		istop_ = 2
		if(present(istop)) istop_ = istop
		
		!  ISTOP_ = 1:  Normwise backward error, infinity norm 
		!  ISTOP_ = 2:  ||r||/||b||, 2-norm 
		if((istop_ < 1 ) .or. (istop_ > 2 )) then
			info = psb_err_invalid_istop_
			err = info
			call psb_errpush(info, name, i_err = (/istop_/))
			goto 9999
		endif

		itmax_ = 1000
		if(present(itmax)) itmax_ = itmax

		itrace_ = 0
		if(present(itrace)) itrace_ = itrace

		base_type_ = "C"
		if(present(base_type)) base_type_ = base_type

		Gram_solver_ = lapackCC
		if(present(Gram_solver)) Gram_solver_ = Gram_solver
		
		FGS_sweeps_ = 30
		if(present(FGS_sweeps)) FGS_sweeps_ = FGS_sweeps

		mglob = desc_a%get_global_rows()
		n_row = desc_a%get_local_rows()
		n_col = desc_a%get_local_cols()

		call psb_chkvect(mglob, lone, x%get_nrows(), lone, lone, desc_a, info)
		if(info /= psb_success_) then
			info = psb_err_from_subroutine_
			call psb_errpush(info, name,  a_err = 'psb_chkvect on x')
			goto 9999
		end if

		call psb_chkvect(mglob, lone, b%get_nrows(), lone, lone, desc_a, info)
		if(info /= psb_success_) then
			info = psb_err_from_subroutine_    
			call psb_errpush(info, name, a_err='psb_chkvect on b')
			goto 9999
		end if

		!Allocate and assembly data structure
		allocate(alpha(s), beta(s, s), W(s, s), pW(s), temp_fa(s, 2*s + 1), B2(s, s), c0(s), aux_fa(4*n_col), stat = info)
		if(info == psb_success_) call psb_geall(r, desc_a, info)
		if(info == psb_success_) call psb_geall(Z, desc_a, info, n = s)
		if(info == psb_success_) call psb_geall(Q, desc_a, info, n = s)
		if(info == psb_success_) call psb_geall(P, desc_a, info, n = s)
		if(info == psb_success_) call psb_geall(V, desc_a, info, n = s)
		if(info == psb_success_) call psb_geall(temp_mv, desc_a, info, n = s)
		if(info == psb_success_) call psb_geall(aux_mv, desc_a, info, n = 3)
		if(info == psb_success_) call psb_geasb(r, desc_a, info)
		if(info == psb_success_) call psb_geasb(Z, desc_a, info)
		if(info == psb_success_) call psb_geasb(Q, desc_a, info)
		if(info == psb_success_) call psb_geasb(P, desc_a, info)
		if(info == psb_success_) call psb_geasb(V, desc_a, info)
		if(info == psb_success_) call psb_geasb(temp_mv, desc_a, info)
		if(info == psb_success_) call psb_geasb(aux_mv, desc_a, info)

		if(info /= psb_success_) then 
			info = psb_err_from_subroutine_ 
			call psb_errpush(info, name)
			goto 9999
		end if

		! First residual calculation
		call psb_geaxpby(done, b, dzero, r, desc_a, info)
		if(info == psb_success_) call psb_spmm(-done, a, x, done, r, desc_a, info)
		if(info /= psb_success_) then 
			info = psb_err_from_subroutine_ 
			call psb_errpush(info, name)
			goto 9999
		end if

		! Init converence
		call psb_init_conv(methdfullname, istop_, itrace_, itmax_, a, x, b, eps, desc_a, stopdat, info)
		if(info /= psb_success_) then 
			info = psb_err_from_subroutine_ 
			call psb_errpush(info, name)
			goto 9999
		end if

		! check convergence here?

		! Chebyshev coefficient calculation
		select case (base_type_)
			case("M")
				cheb_coeff = dzero
			case("C")
				cheb_coeff = psb_d_chebyshev_coefficients(a, prec, desc_a, info, eigext)
			case default
				info = psb_err_invalid_input_ 
				call psb_errpush(info, name)
				goto 9999
		end select

		if(info /= psb_success_) then 
			info = psb_err_from_subroutine_ 
			call psb_errpush(info, name)
			goto 9999
		end if
		
		! First matrix power kernel
		call psb_mixed_pMPK_split(a, prec, r, P, V, s, desc_a, info, base_type = base_type_, &
										& alpha = cheb_coeff(1), beta = cheb_coeff(2), gamma = cheb_coeff(3), &
										& mvec_temp = aux_mv, farr_temp = aux_fa)
		if(info /= psb_success_) then 
			info = psb_err_from_subroutine_ 
			call psb_errpush(info, name)
			goto 9999
		end if

		! Compute first Gram system components
		call psb_gedots(P, V, temp_fa(:, 1 : s), desc_a, info, global = .false.)
		call psb_gedots(P, r, temp_fa(:, 2*s + 1), desc_a, info, global = .false.)
		call psb_sum(desc_a%get_context(), temp_fa)

		W = temp_fa(:, 1 : s)
		alpha = temp_fa(:, 2*s + 1)

		! Loop until convergence (or maxiter)
		do itidx = 1, itmax_
			! Factor matrix W (if soving with LU or Cholesky factorization)
			if(Gram_solver_ == lapackLU) call dgetrf(s, s, W, s, pW, info)
			if(Gram_solver_ == lapackCC) call dpotrf('L', s, W, s, info)

			! Solve for alpha
			select case(Gram_solver_)
				case(forwardGS);  call inner_solver_fgs_1D(W, alpha, FGS_sweeps_)
				case(lapackLU);   call dgetrs('N', s, 1, W, s, pW, alpha, s, info)
				case(lapackCC);   call dpotrs('L', s, 1, W, s, alpha, s, info)
				case default
					info = psb_err_invalid_input_ 
					call psb_errpush(info, name)
					goto 9999
			end select

			! Update solution and residual
			call psb_geaxpby(P, alpha, x, desc_a, info, upd_flag = .true.)
			call psb_geaxpby(V, -alpha, r, desc_a, info, upd_flag = .true.)

			! Check convergence
			if(psb_check_conv(methdfullname, itidx, x, r, desc_a, stopdat, info)) exit

			! Matrix power kernel
			call psb_mixed_pMPK_split(a, prec, r, Z, Q, s, desc_a, info, base_type = base_type_, &
										& alpha = cheb_coeff(1), beta = cheb_coeff(2), gamma = cheb_coeff(3), &
										& mvec_temp = aux_mv, farr_temp = aux_fa)

			! Compute dot products
			call psb_gedots(P, Q, temp_fa(:, 1 : s), desc_a, info, global = .false.)
			call psb_gedots(Z, Q, temp_fa(:, s+1 : 2*s), desc_a, info, global = .false.)
			call psb_gedots(Z, r, temp_fa(:, 2*s + 1), desc_a, info, global = .false.)
			call psb_sum(desc_a%get_context(), temp_fa)

			!Compute new rhs for beta
			B2 = -temp_fa(:, 1 : s)

			! Solve for beta
			beta = B2
			select case(Gram_solver_)
				case(forwardGS);  call inner_solver_fgs_2D(W, beta, FGS_sweeps_)
				case(lapackLU);   call dgetrs('N', s, s, W, s, pW, beta, s, info)
				case(lapackCC);   call dpotrs('L', s, s, W, s, beta, s, info)
				case default
					info = psb_err_invalid_input_ 
					call psb_errpush(info, name)
					goto 9999
			end select

			! Update P and V. Use of temp_mv in needed because internal dgemm constraint
			call psb_geaxpby(P, beta, temp_mv, desc_a, info, upd_flag = .false.)
			call psb_geaxpby(done, Z, done, temp_mv, P, desc_a, info)
			call psb_geaxpby(V, beta, temp_mv, desc_a, info, upd_flag = .false.)
			call psb_geaxpby(done, Q, done, temp_mv, V, desc_a, info)

			!Compute new Gram matrix
			W = temp_fa(:, s+1 : 2*s)
			call dgemm('T', 'N', s, s, s, -done, beta, s, B2, s, done, W, s)
			
			!Compute new rhs for alpha
			alpha = temp_fa(:, 2*s + 1)
		end do

		call psb_end_conv(methdfullname, itidx, desc_a, stopdat, info, derr, iter)
		if(present(err)) err = derr
		if(present(iter)) iter = iter * s

		if(info == psb_success_) call psb_gefree(r, desc_a, info)
		if(info == psb_success_) call psb_gefree(Z, desc_a, info)
		if(info == psb_success_) call psb_gefree(Q, desc_a, info)
		if(info == psb_success_) call psb_gefree(P, desc_a, info)
		if(info == psb_success_) call psb_gefree(V, desc_a, info)
		if(info == psb_success_) call psb_gefree(temp_mv, desc_a, info)
		if(info == psb_success_) call psb_gefree(aux_mv, desc_a, info)

		if(info == psb_success_) deallocate(alpha, beta, W, pW, temp_fa, B2, c0, aux_fa, stat = info)
		if(info /= psb_success_) then
			call psb_errpush(info, name)
			goto 9999
		end if

		call psb_erractionrestore(err_act)
		return

	9999 call psb_error_handler(err_act)
		return
	contains 
		function psb_d_chebyshev_coefficients(a, prec, desc, info, eigext) result(coeff)
			type(psb_dspmat_type), intent(in)     :: a
			class(psb_sprec_type), intent(inout)  :: prec
			type(psb_desc_type), intent(in)       :: desc
			integer(psb_ipk_), intent(out)        :: info
			real(psb_dpk_), optional, intent(in)  :: eigext(2)
			
			real(psb_dpk_) :: lambda_max, lambda_min
			real(psb_dpk_) :: coeff(3)
	
			info = psb_success_

			if(present(eigext)) then
				lambda_min = eigext(1)
				lambda_max = eigext(2)
			else 
				call psb_mixed_powermethod(a, prec, lambda_max, desc, info)
				lambda_min = dzero
			end if

			if(info /= psb_success_) then 
				info = psb_err_from_subroutine_
				return
			end if

			if(lambda_min < 0) then 
				info = psb_err_fatal_ ! TODO: set the correct error
				return
			end if

			if(lambda_max < lambda_min) then 
				info = psb_err_fatal_ ! TODO: set the correct error
				return
			end if

			coeff(1) = 2_psb_dpk_ / (lambda_max - lambda_min)
			coeff(2) = (lambda_max + lambda_min) / (lambda_max - lambda_min)
			coeff(3) = done
		end function psb_d_chebyshev_coefficients

		subroutine inner_solver_fgs_1D(M, rhs, num_iter)
			real(psb_dpk_), intent(in)    :: M(:, :)
			real(psb_dpk_), intent(inout) :: rhs(:)
			integer(psb_ipk_), intent(in) :: num_iter

			integer(psb_ipk_) :: iter_idx, i, j, n
			real(psb_dpk_)    :: sol(size(rhs))

			sol = dzero
			n = size(sol)
			do iter_idx = 1, num_iter
				do i = 1, n
					sol(i) = rhs(i)
					do j = 1, n
						if(j /= i) sol(i) = sol(i) - M(i, j) * sol(j)
					end do
					sol(i) = sol(i) / M(i, i)
				end do
			end do
			rhs = sol
		end subroutine inner_solver_fgs_1D

		subroutine inner_solver_fgs_2D(M, rhs, num_iter)
			real(psb_dpk_), intent(in)    :: M(:, :)
			real(psb_dpk_), intent(inout) :: rhs(:, :)
			integer(psb_ipk_), intent(in) :: num_iter

			integer(psb_ipk_) :: iter_idx, i, j, n
			real(psb_dpk_)    :: sol(size(rhs, 1), size(rhs, 2))

			sol = dzero
			n = size(sol, 1)
			do iter_idx = 1, num_iter
				do i = 1, n
					sol(i, :) = rhs(i, :)
					do j = 1, n
						if(j /= i) sol(i, :) = sol(i, :) - M(i, j) * sol(j, :)
					end do
					sol(i, :) = sol(i, :) / M(i, i)
				end do
			end do
			rhs = sol
		end subroutine inner_solver_fgs_2D
	end subroutine psb_dscg2_mixed


	subroutine psb_mixed_powermethod(a, prec, lambda, desc, info, x, flag, itmax, iter, tol)
    use psb_base_mod
    use psb_prec_mod
    implicit none
    type(psb_dspmat_type), intent(in)     :: a
    class(psb_sprec_type), intent(inout)  :: prec 
    real(psb_dpk_), intent(out)           :: lambda
    type(psb_desc_type), intent(in)       :: desc
    integer(psb_ipk_), intent(out)        :: info
    type(psb_d_vect_type), intent(inout), optional  :: x
    logical, intent(in), optional                   :: flag
    integer(psb_ipk_), intent(in), optional         :: itmax
    integer(psb_ipk_), intent(out), optional        :: iter
    real(psb_dpk_), intent(in), optional            :: tol ! def 10^-3

    type(psb_d_vect_type)   :: z, q 
    type(psb_s_vect_type)   :: q_single 
    integer(psb_ipk_)       :: i, itmax_
    real(psb_dpk_)          :: tol_
    real(psb_dpk_)          :: lambda_old, norm_factor
    
    if(present(itmax)) then
        itmax_ = itmax
    else
        itmax_ = 20_psb_ipk_
    end if

    if(present(tol)) then
        tol_ = tol
    else
        tol_ = real(1.0e-3, psb_dpk_)
    end if

    call psb_geall(z, desc, info)
    call psb_geall(q, desc, info)
    call psb_geall(q_single, desc, info)
    call psb_geasb(z, desc, info)
    call psb_geasb(q, desc, info)
    call psb_geasb(q_single, desc, info)
    call q%set(done)

    if(present(x) .and. present(flag)) then     !TO DO: can we avoid the allocation of one vector in this case?
        if(flag) call psb_geaxpby(done, x, dzero, q, desc, info)
    end if

    lambda_old = dzero
    do i = 1, itmax_
        norm_factor = done / psb_genrm2(q, desc, info)
        call psb_geaxpby(norm_factor, q, dzero, z, desc, info)      ! z_k = q_k / |q_k|_2
        call psb_spmm(done, a, z, dzero, q, desc, info)             ! q_k = A z_k 

				call psb_d2s_vect(q, q_single, info)
        call prec%apply(q_single, desc, info)												! q_k = B q_k
        call psb_s2d_vect(q_single, q, info, mold=x%v)

        lambda = psb_gedot(z, q, desc, info)                        ! lambda = <z_k, q_k> = z_k^T BA z_k
        if(abs(lambda - lambda_old) < tol_ * abs(lambda)) exit
        lambda_old = lambda
    end do

    if(present(iter)) iter = i
    if(present(x)) call psb_geaxpby(done, z, dzero, x, desc, info)

    call psb_gefree(z, desc, info)
    call psb_gefree(q, desc, info)  
    call psb_gefree(q_single, desc, info)    
	end subroutine psb_mixed_powermethod

end module psb_mixed_alg_mod