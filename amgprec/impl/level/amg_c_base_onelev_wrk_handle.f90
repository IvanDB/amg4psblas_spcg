!  
!   
!                             AMG4PSBLAS version 1.0
!    Algebraic Multigrid Package
!               based on PSBLAS (Parallel Sparse BLAS version 3.7)
!    
!    (C) Copyright 2021 
!  
!        Salvatore Filippone  
!        Pasqua D'Ambra   
!        Fabio Durastante        
!   
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the AMG4PSBLAS group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific prior written permission.
!   
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AMG4PSBLAS GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!   
!  
submodule (amg_c_onelev_mod)  amg_c_base_onelev_wrk_handle_impl
  use psb_base_mod
  
contains
  
  module subroutine c_base_onelev_move_alloc(lv, b,info)
    implicit none
    class(amg_c_onelev_type), target, intent(inout) :: lv, b
    integer(psb_ipk_), intent(out) :: info

    call b%free(info)
    b%parms  = lv%parms
    b%szratio = lv%szratio
    if (associated(lv%sm2,lv%sm2a)) then
      call move_alloc(lv%sm,b%sm)
      call move_alloc(lv%sm2a,b%sm2a)
      b%sm2 =>b%sm2a
    else
      call move_alloc(lv%sm,b%sm)
      call move_alloc(lv%sm2a,b%sm2a)
      b%sm2 =>b%sm
    end if

    call move_alloc(lv%aggr,b%aggr)
    if (info == psb_success_) call psb_move_alloc(lv%ac,b%ac,info)
    if (info == psb_success_) call psb_move_alloc(lv%tprol,b%tprol,info)
    if (info == psb_success_) call psb_move_alloc(lv%desc_ac,b%desc_ac,info)
    if (info == psb_success_) call psb_move_alloc(lv%linmap,b%linmap,info)
    if (info == psb_success_) call lv%remap_data%move_alloc(b%remap_data,info)
    b%base_a    => lv%base_a
    b%base_desc => lv%base_desc

  end subroutine c_base_onelev_move_alloc

  module subroutine c_base_onelev_allocate_wrk(lv,info,vmold)
    implicit none
    class(amg_c_onelev_type), target, intent(inout) :: lv
    integer(psb_ipk_), intent(out) :: info
    class(psb_c_base_vect_type), intent(in), optional  :: vmold
    !
    integer(psb_ipk_) :: nwv, i
    info = psb_success_
    nwv = lv%get_wrksz()
    if (.not.allocated(lv%wrk)) allocate(lv%wrk,stat=info)
!!$    write(0,*) 'From allocate_wrk :',lv%remap_data%desc_ac_pre_remap%is_asb()
    if (info == 0) then
      if (lv%remap_data%desc_ac_pre_remap%is_asb()) then
        !
        !  Need to fix this, we need two different allocations
        !
        call lv%wrk%alloc(nwv,lv%base_desc,info,vmold=vmold,&
             & desc2=lv%remap_data%desc_ac_pre_remap)
      else
        call lv%wrk%alloc(nwv,lv%base_desc,info,vmold=vmold)
      end if
    end if

  end subroutine c_base_onelev_allocate_wrk

  module subroutine c_base_onelev_free_wrk(lv,info)
    implicit none
    class(amg_c_onelev_type), target, intent(inout) :: lv
    integer(psb_ipk_), intent(out) :: info
    !
    integer(psb_ipk_) :: nwv,i
    info = psb_success_

    if (allocated(lv%wrk)) then
      call lv%wrk%free(info)
      if (info == 0) deallocate(lv%wrk,stat=info)
    end if
  end subroutine c_base_onelev_free_wrk

  module subroutine c_wrk_alloc(wk,nwv,desc,info,vmold, desc2)
    Implicit None
    ! Arguments
    class(amg_cmlprec_wrk_type), target, intent(inout) :: wk
    integer(psb_ipk_), intent(in)                     :: nwv
    type(psb_desc_type), intent(in)                   :: desc
    integer(psb_ipk_), intent(out)                    :: info
    class(psb_c_base_vect_type), intent(in), optional  :: vmold
    type(psb_desc_type), intent(in), optional          :: desc2
    !
    integer(psb_ipk_) :: i

    info = psb_success_
    call wk%free(info)
!!$    write(0,*) 'wrk_alloc D: "',trim(desc%get_fmt()),'"',&
!!$         & present(desc2),desc%is_valid()

    allocate(wk%wv(nwv),stat=info)    
    if  (present(desc2).and.(desc%is_valid())) then
!!$      write(0,*) 'wrk_alloc D2:',desc2%get_fmt(),desc2%is_asb()
      if (desc2%get_local_cols()>desc%get_local_cols()) then
        call c_inner_do_wrk_alloc(wk,nwv,desc2,vmold=vmold)
      else
        call c_inner_do_wrk_alloc(wk,nwv,desc,vmold=vmold)
      end if
    else if (present(desc2)) then
      call c_inner_do_wrk_alloc(wk,nwv,desc2,vmold=vmold)
    else if (desc%is_valid()) then
      call c_inner_do_wrk_alloc(wk,nwv,desc,vmold=vmold)      
    end if

  contains
  end subroutine c_wrk_alloc

  module subroutine c_inner_do_wrk_alloc(wk,nwv,desc,vmold)
    class(amg_cmlprec_wrk_type), target, intent(inout) :: wk
    integer(psb_ipk_), intent(in)                     :: nwv
    type(psb_desc_type), intent(in)                   :: desc
    class(psb_c_base_vect_type), intent(in), optional  :: vmold
    
    integer(psb_ipk_) :: i
    
    call psb_geasb(wk%vx2l,desc,info,&
         & scratch=.true.,mold=vmold)
    call psb_geasb(wk%vy2l,desc,info,&
         & scratch=.true.,mold=vmold)
    call psb_geasb(wk%vtx,desc,info,&
         & scratch=.true.,mold=vmold)
    call psb_geasb(wk%vty,desc,info,&
         & scratch=.true.,mold=vmold)
    do i=1,nwv
      call psb_geasb(wk%wv(i),desc,info,&
           & scratch=.true.,mold=vmold)
    end do
  end subroutine c_inner_do_wrk_alloc
  
  
  module subroutine c_wrk_free(wk,info)

    Implicit None

    ! Arguments
    class(amg_cmlprec_wrk_type), target, intent(inout) :: wk
    integer(psb_ipk_), intent(out)                    :: info
    !
    integer(psb_ipk_) :: i
    info = psb_success_

    if (allocated(wk%tx)) deallocate(wk%tx, stat=info)
    if (allocated(wk%ty)) deallocate(wk%ty, stat=info)
    if (allocated(wk%x2l)) deallocate(wk%x2l, stat=info)
    if (allocated(wk%y2l)) deallocate(wk%y2l, stat=info)
    call wk%vtx%free(info)
    call wk%vty%free(info)
    call wk%vx2l%free(info)
    call wk%vy2l%free(info)
    if (allocated(wk%wv)) then
      do i=1,size(wk%wv)
        call wk%wv(i)%free(info)
      end do
      deallocate(wk%wv, stat=info)
    end if

  end subroutine c_wrk_free

  module subroutine c_wrk_clone(wk,wkout,info)
    Implicit None

    ! Arguments
    class(amg_cmlprec_wrk_type), target, intent(inout) :: wk
    class(amg_cmlprec_wrk_type), target, intent(inout) :: wkout
    integer(psb_ipk_), intent(out)                    :: info
    !
    integer(psb_ipk_) :: i
    info = psb_success_

    call psb_safe_ab_cpy(wk%tx,wkout%tx,info)
    call psb_safe_ab_cpy(wk%ty,wkout%ty,info)
    call psb_safe_ab_cpy(wk%x2l,wkout%x2l,info)
    call psb_safe_ab_cpy(wk%y2l,wkout%y2l,info)
    call wk%vtx%clone(wkout%vtx,info)
    call wk%vty%clone(wkout%vty,info)
    call wk%vx2l%clone(wkout%vx2l,info)
    call wk%vy2l%clone(wkout%vy2l,info)
    if (allocated(wkout%wv)) then
      do i=1,size(wkout%wv)
        call wkout%wv(i)%free(info)
      end do
      deallocate( wkout%wv)
    end if
    allocate(wkout%wv(size(wk%wv)),stat=info)
    do i=1,size(wk%wv)
      call wk%wv(i)%clone(wkout%wv(i),info)
    end do
    return

  end subroutine c_wrk_clone

  module subroutine c_wrk_move_alloc(wk, b,info)
    implicit none
    class(amg_cmlprec_wrk_type), target, intent(inout) :: wk, b
    integer(psb_ipk_), intent(out) :: info

    call b%free(info)
    call move_alloc(wk%tx,b%tx)
    call move_alloc(wk%ty,b%ty)
    call move_alloc(wk%x2l,b%x2l)
    call move_alloc(wk%y2l,b%y2l)
    !
    ! Should define V%move_alloc....
    call move_alloc(wk%vtx%v,b%vtx%v)
    call move_alloc(wk%vty%v,b%vty%v)
    call move_alloc(wk%vx2l%v,b%vx2l%v)
    call move_alloc(wk%vy2l%v,b%vy2l%v)
    call move_alloc(wk%wv,b%wv)

  end subroutine c_wrk_move_alloc

  module subroutine c_wrk_cnv(wk,info,vmold)
    Implicit None

    ! Arguments
    class(amg_cmlprec_wrk_type), target, intent(inout) :: wk
    integer(psb_ipk_), intent(out)                    :: info
    class(psb_c_base_vect_type), intent(in), optional  :: vmold
    !
    integer(psb_ipk_) :: i

    info = psb_success_
    if (present(vmold)) then
      call wk%vtx%cnv(vmold)
      call wk%vty%cnv(vmold)
      call wk%vx2l%cnv(vmold)
      call wk%vy2l%cnv(vmold)
      if (allocated(wk%wv)) then
        do i=1,size(wk%wv)
          call wk%wv(i)%cnv(vmold)
        end do
      end if
    end if
  end subroutine c_wrk_cnv

  module function c_wrk_sizeof(wk) result(val)
    implicit none
    class(amg_cmlprec_wrk_type), intent(in) :: wk
    integer(psb_epk_) :: val
    integer :: i
    val = 0
    val = val + (1_psb_epk_ * (2*psb_sizeof_sp)) * psb_size(wk%tx)
    val = val + (1_psb_epk_ * (2*psb_sizeof_sp)) * psb_size(wk%ty)
    val = val + (1_psb_epk_ * (2*psb_sizeof_sp)) * psb_size(wk%x2l)
    val = val + (1_psb_epk_ * (2*psb_sizeof_sp)) * psb_size(wk%y2l)
    val = val + wk%vtx%sizeof()
    val = val + wk%vty%sizeof()
    val = val + wk%vx2l%sizeof()
    val = val + wk%vy2l%sizeof()
    if (allocated(wk%wv)) then
      do i=1, size(wk%wv)
        val = val + wk%wv(i)%sizeof()
      end do
    end if
  end function c_wrk_sizeof

  module subroutine c_remap_data_clone(rmp, remap_out, info)
    implicit none
    ! Arguments
    class(amg_c_remap_data_type), target, intent(inout) :: rmp
    class(amg_c_remap_data_type), target, intent(inout) :: remap_out
    integer(psb_ipk_), intent(out)                    :: info
    !
    integer(psb_ipk_) :: i

    info = psb_success_

    call rmp%ac_pre_remap%clone(remap_out%ac_pre_remap,info)
    if (info == psb_success_) &
         & call rmp%desc_ac_pre_remap%clone(remap_out%desc_ac_pre_remap,info)
    remap_out%idest = rmp%idest
    call psb_safe_ab_cpy(rmp%isrc,remap_out%isrc,info)
    call psb_safe_ab_cpy(rmp%nrsrc,remap_out%nrsrc,info)
  end subroutine c_remap_data_clone

  module subroutine c_remap_move_alloc(rmp, remap_out, info)
    implicit none
    ! Arguments
    class(amg_c_remap_data_type), target, intent(inout) :: rmp
    class(amg_c_remap_data_type), target, intent(inout) :: remap_out
    integer(psb_ipk_), intent(out)                    :: info
    !
    integer(psb_ipk_) :: i

    info = psb_success_

    call psb_move_alloc(rmp%ac_pre_remap,remap_out%ac_pre_remap,info)
    if (info == psb_success_) &
         & call psb_move_alloc(rmp%desc_ac_pre_remap,remap_out%desc_ac_pre_remap,info)
    remap_out%idest = rmp%idest
    call move_alloc(rmp%isrc,remap_out%isrc)
    call move_alloc(rmp%nrsrc,remap_out%nrsrc)
    call move_alloc(rmp%naggr,remap_out%naggr)
  end subroutine c_remap_move_alloc
  
end submodule amg_c_base_onelev_wrk_handle_impl
