#include "rundeck_opts.h"

C****
C**** Specify horizontal resolution of ocean and
C**** import vertical resolution from olayers module
C****

      module oceanres
      use olayers, only : lmo,lmo_min,dZO
      implicit none

      integer, parameter :: ! number of gridcells in
     &      imo = 288,  ! longitude
     &      jmo = 180   ! latitude

      ! to-do: move these res-dependent ODIFF settings to its init
      real*8, parameter ::
     &    AKHMIN = 1.d5     ! minimum horizontal viscosity (m2/s)
     &  , AKHFAC = 1d0      ! factor to multiply built-in scaling for k_diff

      ! to-do: make the following into rundeck parameters
      integer, parameter ::
     &    NORDER=4      !  order of Alternating Binomial Filter (must be even)
      real*8, parameter ::
     &    OABFUX=.15d0/4**NORDER   ! coef. for divergence filter in X dir. 
     &  , OABFVX=.15d0/4**NORDER   ! coef. for vorticity filter in X dir.
     &  , by4tonv=.15/(4.**norder) ! coef. for divergence filter in Y dir.
     &  , by4tonu=.15/(4.**norder) ! coef. for vorticity filter in Y dir.


      end module oceanres
