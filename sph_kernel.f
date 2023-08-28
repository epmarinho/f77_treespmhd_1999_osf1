#ifndef _NBODY
      function sphKernel(r,hp)
                                                            /*
      spline kernel (monaghan & lattanzio 1985)
                                                            */
      INCLUDE 'SPH_common.h'
      q=r/hp
      if(q.lt.1.)then
          w = 1.5 * (q*q) * (.5*q - 1) + 1
      else
          if(q.lt.2.)then
                w = 2. - q
                w = .25 * (w * w * w)
          else
            w = 0
          endif
      endif
#if defined(_2D)
      w = xnorm/(hp*hp) * w
#elif defined(_TUBE)
      w = xnorm/hp * w
#else
      w = xnorm/(hp*hp*hp) * w
#endif
      sphKernel = w
      end



      function D_kernel(r,hp)
                                                            /*
      this is to be used in the folowing form:

      (grad_i w_ij) = + r_ij^/r_ij * D_kernel(r_ij,h)
                                                            */
      INCLUDE 'SPH_common.h'
      q=r/hp
      if(q.lt.1.)then
          dw =  q * (.75*q-1.)
      else
          if(q.lt.2.)then
                dw = 2. - q
                dw = - .25 * (dw * dw)
          else
            dw = 0
          endif
      endif
      overh2 = 1./(hp*hp)
#ifdef _TUBE
      dw = xnorm * overh2 * dw
#elif defined(_2D)
      overh3 = overh2 / hp
      dw = xnorm * overh3 * dw
#else
      overh4 = overh2 * overh2
      dw = xnorm * overh4 * dw
#endif
      D_kernel = 3. * dw
      end



      function spline(q)
      ! spline kernel (monaghan & lattanzio 1985)
      INCLUDE 'SPH_common.h'
      if(q.lt.1.)then
          w = 1.5 * (q*q) * (.5*q - 1) + 1.
      else
          if(q.lt.2.)then
                w = 2. - q
                w = .25 * (w * w * w)
          else
            w = 0
          endif
      endif
      spline = w
      end
#endif

