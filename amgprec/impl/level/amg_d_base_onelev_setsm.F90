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
submodule (amg_d_onelev_mod)  amg_d_base_onelev_setsm_impl
  use psb_base_mod
  
contains
  module subroutine amg_d_base_onelev_setsm(lv,val,info,pos)
    implicit none

    ! Arguments
    class(amg_d_onelev_type), target, intent(inout) :: lv
    class(amg_d_base_smoother_type), intent(in)  :: val
    integer(psb_ipk_), intent(out)               :: info
    character(len=*), optional, intent(in)       :: pos

    ! Local variables
    integer(psb_ipk_)                      :: ipos_
    character(len=*), parameter            :: name='amg_base_onelev_setsm'

    info = psb_success_

    if (present(pos)) then
      select case(psb_toupper(trim(pos)))
      case('PRE')
        ipos_ = amg_smooth_pre_
      case('POST')
        ipos_ = amg_smooth_post_
      case default
        ipos_ = amg_smooth_both_
      end select
    else
      ipos_ = amg_smooth_both_
    end if

    if (ipos_ == amg_smooth_both_) then
      if (allocated(lv%sm2a)) then
        call lv%sm2a%free(info)
        deallocate(lv%sm2a, stat=info)
        lv%sm2 => null()
      end if
    end if

    select case(ipos_)
    case(amg_smooth_pre_, amg_smooth_both_) 
      if (allocated(lv%sm)) then
        if (.not.same_type_as(lv%sm,val)) then
          call lv%sm%free(info)
          deallocate(lv%sm, stat=info)
        end if
      endif
      if (.not.allocated(lv%sm)) then
        allocate(lv%sm,mold=val) 
      end if
      call lv%sm%default()        
      if (ipos_ ==  amg_smooth_both_) lv%sm2 => lv%sm
    case(amg_smooth_post_) 
      if (allocated(lv%sm2a)) then
        if (.not.same_type_as(lv%sm2a,val)) then
          call lv%sm2a%free(info)
          deallocate(lv%sm2a, stat=info)
        endif
      end if
      if (.not.allocated(lv%sm2a)) then
        allocate(lv%sm2a,mold=val) 
      end if
      call lv%sm2a%default()
      lv%sm2 => lv%sm2a
    end select

  end subroutine amg_d_base_onelev_setsm

end submodule amg_d_base_onelev_setsm_impl
