/*
*  Copyright (C) 2010 Craig Robinson
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/


/**
 * \author Craig Robinson
 *
 * \brief Function to compute the factorized flux as uses in the new EOBNR_PP
 * model. Flux function given by Phys.Rev.D79:064004,2009.
 */

#include <complex.h>

#ifndef _LALSIMIMRFACTORIZEDFLUX_C
#define _LALSIMIMRFACTORIZEDFLUX_C

/**
 * This function calculates the factorized flux in the EOB dynamics for
 * the EOBNR (and potentially subsequent) models. The flux function
 * is found in Phys.Rev.D79:064004,2009.
 */
static REAL8 XLALSimIMREOBFactorizedFlux(
                      REAL8Vector  *values, /**<< Dynamics r, phi, pr, pphi */
                      const REAL8  omega,   /**<< Angular frequency omega */
                      EOBParams    *ak,     /**<< Structure containing pre-computed parameters */
                      const INT4   lMax     /**<< Maximum l to include when calculating flux (between 2 and 8) */
                     )

{

  REAL8 flux = 0.0;
  REAL8 v;
  REAL8 omegaSq;
  COMPLEX16 hLM;
  INT4 l, m;

  EOBNonQCCoeffs *nqcCoeffs;
  COMPLEX16       hNQC;

  if ( !values || !ak )
  {
    XLAL_ERROR_REAL8( XLAL_EFAULT );
  }

  if ( lMax < 2 )
  {
    XLAL_ERROR_REAL8( XLAL_EINVAL );
  }

  nqcCoeffs = ak->nqcCoeffs;

  /* Omegs is the derivative of phi */
  omegaSq = omega*omega;

  v = cbrt( omega );

  /* We need to apply the NQC for the (2,2) mode */
  /* To avoid having an if statement in the loop we will */
  /* deal with (2,2) and (2,1) separately */
  /* (2,2) */
  l = 2;
  m = 2;

  if ( XLALSimIMREOBNonQCCorrection( &hNQC, values, omega, nqcCoeffs ) == XLAL_FAILURE )
  {
    XLAL_ERROR_REAL8( XLAL_EFUNC );
  }

  if ( XLALSimIMREOBGetFactorizedWaveform( &hLM, values, v, l, m, ak )
           == XLAL_FAILURE )
  {
    XLAL_ERROR_REAL8( XLAL_EFUNC );
  }
  /* For the 2,2 mode, we apply NQC correction to the flux */
  hLM *= hNQC;

  flux = (REAL8)(m * m) * omegaSq * (creal(hLM)*creal(hLM) + cimag(hLM)*cimag(hLM));

  /* (2,1) */
  l = 2;
  m = 1;

  if ( XLALSimIMREOBGetFactorizedWaveform( &hLM, values, v, l, m, ak )
           == XLAL_FAILURE )
  {
    XLAL_ERROR_REAL8( XLAL_EFUNC );
  }

  flux += (REAL8)(m * m) * omegaSq * (creal(hLM)*creal(hLM) + cimag(hLM)*cimag(hLM));

  /* All other modes */
  for ( l = 3; l <= lMax; l++ )
  {
    /*INT4 minM = l-3;
    if ( minM < 1 )
      minM = 1;*/

    for ( m = 1; m <= l; m++ )
    {

      if ( XLALSimIMREOBGetFactorizedWaveform( &hLM, values, v, l, m, ak )
             == XLAL_FAILURE )
      {
        XLAL_ERROR_REAL8( XLAL_EFUNC );
      }

      flux += (REAL8)(m * m) * omegaSq * (creal(hLM)*creal(hLM) + cimag(hLM)*cimag(hLM));
    }
  }

  return flux * LAL_1_PI / 8.0;
}

#endif /*_LALSIMIMRFACTORIZEDFLUX_C*/
