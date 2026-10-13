/*
 *  Copyright (C) 2020 Hector Estelles
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

/* Standard LAL */
#include <lal/Sequence.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/XLALError.h>

/* LAL datatypes and constants */
#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>

/* Time series, frequency series and spherical harmonics */
#include <lal/TimeSeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/SphericalHarmonics.h>
#include <lal/LALSimSphHarmMode.h>
#include <lal/FrequencySeries.h>
#include "LALSimInspiralPrecess.h"

/* LALSimulation */
#include <lal/LALSimIMR.h>
#include <lal/LALSimInspiral.h>

/* Standard C */
#include <math.h>
#include <complex.h>
#include <stdbool.h>

/* Link IMRPhenomTHM routines */
#include "LALSimIMRPhenomTHM_internals.h"
#include "LALSimIMRPhenomTHM.h"
#include "LALSimIMRPhenomX.h"
#include "LALSimIMRPhenomTPHM_FrameRotations.c"
#include "LALSimIMRPhenomTPHM_EulerAngles.c"

#define MAX(a,b) (((a)>(b))?(a):(b))

/**
 * @addtogroup LALSimIMRPhenomT_c
 * @{
 *
 * @name Routines for IMRPhenomTP
 * @{
 *
 * @author Héctor Estellés
 *
 * @brief C code for the IMRPhenomTP phenomenological waveform model.
 *
 * This is a precessing time domain model covering the L=2 sector of Spin-Weighted Spherical Harmonic modes through applying the 'twisting-up' procedure
 * [Phys. Rev. D 91, 024043 (2015), Phys.Rev.D85,084003 (2012), Phys. Rev. D 86, 104063 (2012)] to the dominant non-precessing mode from the IMRPhenomT model.
 *
 * Different versions for the Euler angles employed in the 'twisting-up' of IMRPhenomT can be selected through
 * the WaveformParams option 'PhenomXPrecVersion'.
 *
 * A specific version can be passed with a 5-digit number specification: #####, where tree first digits have to correspond to a valid
 * precessing version of the models IMRPhenomXP/IMRPhenomXPHM, fourth digit selects EulerRDVersion internal option, and 5 digit selects
 * EulerGammaVersion internal option (example, 22310 will select PhenomXPrecVersion=223 for the MSA/NNLO angles calculation, EulerRDVersion=1
 * and EulerGammaVersion=0.)
 *
 * If only three digits are passed, these have to correspond to a valid precessing version of the models IMRPhenomXP/IMRPhenomXPHM, or to be 300
 * for computing Euler angles from a numerical PN evolution of the precessing equations. In the case of numerical PN evolution, other options don't apply
 * (for example, it is invalid to call 30010; version 300 always incorporate an approximation of the Euler angles during the ringdown and a numerical
 *  determination of gamma through the minimal rotation condition.)
 *
 * Final spin options from IMRPhenomXP/XPHM can also be employed, with the same LAL parameter 'PhenomXFinalSpinMod'. Additional final spin option can be employed when PrecVersion=300, calling PhenomXFinalSpinMod=4.
 * This new final spin option computes the total in-plane spin at the coalescence (t=0 in model convention) for constructing the geometrical precessing final spin, and employs the parallel components to L of
 * the individual spins to evaluate the non-precessing final spin fit.
 *
 * Summary of options for analytical PN angles:
 * - EulerRDVersion: 0 (no RD angles attached, MSA/NNLO computed for the whole waveform)
 *                   1 (MSA/NNLO computed up to the peak time of the 22 mode, then RD angles are attached.)
 *                   2 (MSA/NNLO computed up to MECO time, linear continuation performed up to the peak time, then RD angles are attached.)
 *
 * - GammaVersion:   0 (Gamma is computed from the MSA/NNLO PN expressions.)
 *                   1 (Gamma computed numerically from minimal rotation condition.)
 *
 * List of model versions (values of PhenomXPrecVersion LAL parameters):
 * - 22300: MSA during the whole waveform duration, including ringdown. Closest to 223 PhenomXP/HM.
 * - 223/22310: MSA up to the t=0 (peak time of non-precessing 22 mode) and then ringdown angles addition.
 * - 22311: MSA up to the t=0 (peak time of non-precessing 22 mode), but gamma angle computed numerically from minimal rotation condition, and then ringdown angles addition.
 * - 22320: MSA up to t=tMECO, then linear continuation of alpha and beta up to t=0, then ringdown angles addition (gamma computed from minimal rotation condition from tMECO).
 * - 22321: MSA up to t=tMECO, then linear continuation of alpha and beta up to t=0, then ringdown angles addition (gamma computed from minimal rotation during the whole waveform).
 * - 10200: NNLO during the whole waveform duration, including ringdown. Closest to 102 PhenomXP/HM.
 * - 102/10210: NNLO up to the t=0 (peak time of non-precessing 22 mode) and then ringdown angles addition.
 * - 10211: NNLO up to the t=0 (peak time of non-precessing 22 mode), but gamma angle computed numerically from minimal rotation condition, and then ringdown angles addition.
 * - 10220: NNLO up to t=tMECO, then linear continuation of alpha and beta up to t=0, then ringdown angles addition (gamma computed from minimal rotation condition from tMECO).
 * - 10221: NNLo up to t=tMECO, then linear continuation of alpha and beta up to t=0, then ringdown angles addition (gamma computed from minimal rotation during the whole waveform).
 * - 300: Numerical evolution of the precessing spin equations, employing PhenomT 22 mode frequency for approximating orbital frequency. Angles determined by the direction evolution of L(t).
 *  This version accepts the same FinalSpinMod options as 102* versions (i.e, FS=3 corresponds to FS=0) but also admits a computation of final spin reading the spin components at merger, this is specified by FinalSpinMod = 4.
 *
 * All analytical versions (MSA/NNLO +  different options) can be called with non-default implementations of MSA/NNLO in XP/XPHM (i.e, 224, 104, 220, etc ...)
 */

/**
 * Routine to compute time-domain polarizations for the IMRPhenomTP model.
 *
 * This is a precessing time domain model covering the L=2 sector of Spin-Weighted Spherical Harmonic modes through
 * applying the 'twisting-up' procedure [Phys. Rev. D 91, 024043 (2015), Phys.Rev.D85,084003 (2012), Phys. Rev. D 86, 104063 (2012)] to the dominant non-precessing mode from the IMRPhenomT/HM models.
 *
 */
int XLALSimIMRPhenomTP(
  REAL8TimeSeries **hp,       /**< [out] TD waveform for plus polarisation */
  REAL8TimeSeries **hc,       /**< [out] TD waveform for cross polarisation */
  REAL8 m1_SI,                /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,                /**< Mass of companion 2 (kg) */
  REAL8 chi1x,                /**< x component of primary spin*/
  REAL8 chi1y,                /**< y component of primary spin*/
  REAL8 chi1z,                /**< z component of primary spin */
  REAL8 chi2x,                /**< x component of secondary spin*/
  REAL8 chi2y,                /**< y component of secondary spin*/
  REAL8 chi2z,                /**< z component of secondary spin */
  REAL8 distance,             /**< Luminosity distance (m) */
  REAL8 inclination,          /**< inclination (in rad) */
  REAL8 deltaT,               /**< sampling interval (s) */
  REAL8 fmin,               /**< starting GW frequency (Hz) */
  REAL8 fRef,               /**< reference GW frequency (Hz) */
  REAL8 phiRef,               /**< reference orbital phase (rad) */
  LALDict *lalParams            /**< LAL dictionary containing accessory parameters */
  )
{

  UINT4 status;

  /* Check if aligned spin limit is reached for version 300.
  If in plane spin is below 10^-8, then fallback to IMRPhenomT model. */
  REAL8 chi_in_plane = sqrt(pow(chi1x + chi2x,2) + pow(chi1y + chi2y,2));
  if(chi_in_plane<1e-8 && XLALSimInspiralWaveformParamsLookupPhenomXPrecVersion(lalParams)==300)
  {
  		status = XLALSimIMRPhenomT(hp, hc, m1_SI, m2_SI, chi1z, chi2z, distance, inclination, deltaT, fmin, fRef, phiRef, lalParams);
  		return status;
  }


  /* Compute modes in the L0 frame */

  SphHarmTimeSeries *hlm = NULL;

  status = XLALSimIMRPhenomTPHM_L0Modes(&hlm, m1_SI, m2_SI, chi1x,chi1y,chi1z,chi2x,chi2y,chi2z,distance,inclination,deltaT,fmin,fRef,phiRef,lalParams, 1);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: function XLALSimIMRPhenomTPHM_L0Modes has failed.");

  LIGOTimeGPS epoch = (hlm)->mode->epoch;
  size_t length = (hlm)->mode->data->length;

  /* Compute polarisations */

  /* Auxiliary time series for setting the polarisations to zero */
  REAL8TimeSeries *hplus = XLALCreateREAL8TimeSeries ("hplus", &epoch, 0.0, deltaT, &lalStrainUnit,length);
  REAL8TimeSeries *hcross = XLALCreateREAL8TimeSeries ("hcross", &epoch, 0.0, deltaT, &lalStrainUnit,length);
  memset( hplus->data->data, 0, hplus->data->length * sizeof(REAL8) );
  memset( hcross->data->data, 0, hcross->data->length * sizeof(REAL8) );

  /* Add the modes to the polarisations using lalsim routines.
  Negative modes are explicitely passed, instead of using the symmetry flag */

  SphHarmTimeSeries *hlms_temp = hlm;
  while ( hlms_temp )
  {
      XLALSimAddMode(hplus, hcross, hlms_temp->mode, inclination, LAL_PI/2. - phiRef, hlms_temp->l, hlms_temp->m, 0);
      hlms_temp = hlms_temp->next;
  }

  /* Point the output pointers to the relevant time series */
  (*hp) = hplus;
  (*hc) = hcross;

  /* Destroy intermediate time series */
  XLALDestroySphHarmTimeSeries(hlm);
  XLALDestroySphHarmTimeSeries(hlms_temp);

  return status;
}

/** @} */
/**
* @name Routines for IMRPhenomTPHM
* @{
*
* @author Héctor Estellés
*
* @brief C Code for the IMRPhenomTPHM phenomenological waveform model.
*
* This is a precessing time domain model covering up to the L=5 sector of Spin-Weighted Spherical Harmonic modes through
* applying the 'twisting-up' procedure [Phys. Rev. D 91, 024043 (2015), Phys.Rev.D85,084003 (2012), Phys. Rev. D 86, 104063 (2012)] to the non-precessing modes from the IMRPhenomTHM model.
*
* @note The same precessing versions specified for the IMRPhenomTP model are available here.
*/


/**
 * Routine to compute time-domain polarizations for the IMRPhenomTPHM model.
 *
 */
int XLALSimIMRPhenomTPHM(
  REAL8TimeSeries **hp,       /**< [out] TD waveform for plus polarisation */
  REAL8TimeSeries **hc,       /**< [out] TD waveform for cross polarisation */
  REAL8 m1_SI,                /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,                /**< Mass of companion 2 (kg) */
  REAL8 chi1x,                /**< x component of primary spin*/
  REAL8 chi1y,                /**< y component of primary spin*/
  REAL8 chi1z,                /**< z component of primary spin */
  REAL8 chi2x,                /**< x component of secondary spin*/
  REAL8 chi2y,                /**< y component of secondary spin*/
  REAL8 chi2z,                /**< z component of secondary spin */
  REAL8 distance,             /**< Luminosity distance (m) */
  REAL8 inclination,          /**< inclination (in rad) */
  REAL8 deltaT,               /**< sampling interval (s) */
  REAL8 fmin,                 /**< starting GW frequency (Hz) */
  REAL8 fRef,                 /**< reference GW frequency (Hz) */
  REAL8 phiRef,               /**< reference orbital phase (rad) */
  LALDict *lalParams          /**< LAL dictionary containing accessory parameters */
  )
{

  UINT4 status;

  /* Check if aligned spin limit is reached for version 300.
  If in plane spin is below 10^-8, then fallback to IMRPhenomTHM model. */
  REAL8 chi_in_plane = sqrt(pow(chi1x + chi2x,2) + pow(chi1y + chi2y,2));
  if(chi_in_plane<1e-8 && XLALSimInspiralWaveformParamsLookupPhenomXPrecVersion(lalParams)==300)
  {
  		status = XLALSimIMRPhenomTHM(hp, hc, m1_SI, m2_SI, chi1z, chi2z, distance, inclination, deltaT, fmin, fRef, phiRef, lalParams);
  		return status;
  }

  /* Compute modes in the L0 frame */

  SphHarmTimeSeries *hlm = NULL;

  status = XLALSimIMRPhenomTPHM_L0Modes(&hlm, m1_SI, m2_SI, chi1x,chi1y,chi1z,chi2x,chi2y,chi2z,distance,inclination,deltaT,fmin,fRef,phiRef,lalParams,0);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: function XLALSimIMRPhenomTPHM_L0Modes has failed.");

  LIGOTimeGPS epoch = (hlm)->mode->epoch;
  size_t length = (hlm)->mode->data->length;

  /* Compute polarisations */

  /* Auxiliary time series for setting the polarisations to zero */
  REAL8TimeSeries *hplus = XLALCreateREAL8TimeSeries ("hplus", &epoch, 0.0, deltaT, &lalStrainUnit,length);
  REAL8TimeSeries *hcross = XLALCreateREAL8TimeSeries ("hcross", &epoch, 0.0, deltaT, &lalStrainUnit,length);
  memset( hplus->data->data, 0, hplus->data->length * sizeof(REAL8) );
  memset( hcross->data->data, 0, hcross->data->length * sizeof(REAL8) );

  /* Add the modes to the polarisations using lalsim routines.
  Negative modes are explicitely passed, instead of using the symmetry flag */

  SphHarmTimeSeries *hlms_temp = hlm;
  while ( hlms_temp )
  {
      XLALSimAddMode(hplus, hcross, hlms_temp->mode, inclination, LAL_PI/2. - phiRef, hlms_temp->l, hlms_temp->m, 0);
      hlms_temp = hlms_temp->next;
  }

  /* Point the output pointers to the relevant time series */
  (*hp) = hplus;
  (*hc) = hcross;

  /* Destroy intermediate time series */
  XLALDestroySphHarmTimeSeries(hlm);
  XLALDestroySphHarmTimeSeries(hlms_temp);

  return status;
}

/**
 * Routine to be used by ChooseTDModes, it returns a list of the time-domain modes of the IMRPhenomTPHM model in the inertial L0-frame.
 This is a wrapper for ChooseTDModes with the following justification: it is desirable to mantain a formal dependence on reference phase
 and inclination for producing the L0 modes employed in the polarisations, in the function XLALSimIMRPhenomTPHM_L0Modes. Although by construction
 L0 frame modes are independent of these parameters, this will allow to check if this actually satisfied if a new angle description is included.
 Since for ChooseTDModes one cannot pass these parameters, they are set to 0 in this wrapper.
 */
SphHarmTimeSeries* XLALSimIMRPhenomTPHM_ChooseTDModes(
  REAL8 m1_SI,                /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,                /**< Mass of companion 2 (kg) */
  REAL8 chi1x,                /**< x component of primary spin*/
  REAL8 chi1y,                /**< y component of primary spin*/
  REAL8 chi1z,                /**< z component of primary spin */
  REAL8 chi2x,                /**< x component of secondary spin*/
  REAL8 chi2y,                /**< y component of secondary spin*/
  REAL8 chi2z,                /**< z component of secondary spin */
  REAL8 distance,             /**< Luminosity distance (m) */
  REAL8 deltaT,               /**< sampling interval (s) */
  REAL8 fmin,                 /**< starting GW frequency (Hz) */
  REAL8 fRef,                 /**< reference GW frequency (Hz) */
  LALDict *lalParams          /**< LAL dictionary containing accessory parameters */
  )
{

  INT4 status;

  /* L0 frame modes are independent of reference phase and inclination */
  REAL8 phiRef = 0.0;
  REAL8 inclination = 0.0;

  /* Compute modes (if a subset of modes is desired, it should be passed through lalParams */
  SphHarmTimeSeries *hlms = NULL;

  status = XLALSimIMRPhenomTPHM_L0Modes(&hlms, m1_SI, m2_SI, chi1x,chi1y,chi1z,chi2x,chi2y,chi2z,distance,inclination,deltaT,fmin,fRef,phiRef,lalParams,0);
  XLAL_CHECK_NULL(status != XLAL_FAILURE, XLAL_EFUNC, "Error: Internal function LALSimIMRPhenomTPHM_L0Modes has failed producing the modes.");

  return hlms;
}

/**
 * Routine to compute a list of the time-domain modes of the IMRPhenomTPHM model in the inertial L0-frame.
 * This function calls XLALSimIMRPhenomTPHM_JModes for producing the modes in the J-frame and the Euler angles employed
 * in the rotation between the co-precessing frame and the J-frame. Then it reads the angles value corresponding to the
 * specified reference frequency and performs and inverse global rotation with these angles from the J-frame to the L0-frame.
 */
int XLALSimIMRPhenomTPHM_L0Modes(
  SphHarmTimeSeries **hlmI,   /**< [out] Modes in the intertial L0=z frame*/
  REAL8 m1_SI,                /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,                /**< Mass of companion 2 (kg) */
  REAL8 chi1x,                /**< x component of primary spin*/
  REAL8 chi1y,                /**< y component of primary spin*/
  REAL8 chi1z,                /**< z component of primary spin */
  REAL8 chi2x,                /**< x component of secondary spin*/
  REAL8 chi2y,                /**< y component of secondary spin*/
  REAL8 chi2z,                /**< z component of secondary spin */
  REAL8 distance,             /**< Luminosity distance (m) */
  REAL8 inclination,          /**< inclination (in rad) */
  REAL8 deltaT,               /**< sampling interval (s) */
  REAL8 fmin,                 /**< starting GW frequency (Hz) */
  REAL8 fRef,                 /**< reference GW frequency (Hz) */
  REAL8 phiRef,               /**< reference orbital phase (rad) */
  LALDict *lalParams,         /**< LAL dictionary containing accessory parameters */
  UINT4 only22                /**< Flag for calling only IMRPhenomTP (dominant 22 coprec mode only) */
  )
{

  UINT4 status;


  *hlmI = NULL;
  REAL8TimeSeries *alphaTS = NULL;
  REAL8TimeSeries *cosbetaTS= NULL;
  REAL8TimeSeries *gammaTS = NULL;
  REAL8 af;

  /* Compute J frame modes */
  status = XLALSimIMRPhenomTPHM_JModes(hlmI, &alphaTS, &cosbetaTS, &gammaTS, &af, m1_SI, m2_SI, chi1x,chi1y,chi1z,chi2x,chi2y,chi2z,distance,inclination,deltaT,fmin,fRef,phiRef,lalParams,only22);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: function XLALSimIMRPhenomTPHM_JModes has failed.");

  /* Set up structs for computing tref index */
  IMRPhenomTWaveformStruct *pWF;
  pWF    = XLALMalloc(sizeof(IMRPhenomTWaveformStruct));
  status = IMRPhenomTSetWaveformVariables(pWF, m1_SI, m2_SI, chi1z, chi2z, distance, deltaT, fmin, fRef, phiRef, lalParams);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: Internal function IMRPhenomTSetWaveformVariables has failed.");

  IMRPhenomTPhase22Struct *pPhase;
  pPhase = XLALMalloc(sizeof(IMRPhenomTPhase22Struct));
  status   = IMRPhenomTSetPhase22Coefficients(pPhase, pWF);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: Internal function IMRPhenomTSetPhase22Coefficients has failed.");

  // Compute tref index

  if(pPhase->tmin<0.0 && pPhase->tRef>=0.0) // Sanity check for reference frequency, to be consistent with _JModes procedure. Proper warning message raised there.
  {
    pPhase->tRef = MAX(-8*pWF->dtM, pPhase->tmin);
  }


  REAL8 tend = -pPhase->tmin;
  REAL8 tint1 = tend + pPhase->tRef;
  size_t length1 = floor(fabs(tint1)/pWF->dtM);

  // Angles for the global rotation from J-frame to L0-frame: Euler angles evaluated at tref.
  REAL8 alphaJtoI, cosbetaJtoI, gammaJtoI;
  alphaJtoI = (alphaTS)->data->data[length1];
  cosbetaJtoI = (cosbetaTS)->data->data[length1];
  gammaJtoI = (gammaTS)->data->data[length1];

  // Perform global rotation from J frame to L0 frame
  PhenomT_precomputed_sqrt *SQRT;
  SQRT  = XLALMalloc(sizeof(PhenomT_precomputed_sqrt));
  IMRPhenomTPHM_SetPrecomputedSqrt(SQRT);

  size_t length = (*hlmI)->mode->data->length;
  status = PhenomTPHM_RotateModes_Global(*hlmI, -alphaJtoI, cosbetaJtoI, -gammaJtoI, length, SQRT);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: Internal function PhenomTPHM_RotateModes_Global has failed.");

  // Destroy time series
  XLALDestroyREAL8TimeSeries(alphaTS);
  XLALDestroyREAL8TimeSeries(cosbetaTS);
  XLALDestroyREAL8TimeSeries(gammaTS);

  // Free structs
  LALFree(pWF);
  LALFree(pPhase);
  LALFree(SQRT);

  return status;
}



/**
 * Routine to compute a list of the time-domain modes of the IMRPhenomTPHM model in the inertial J-frame.
 * It also outputs the time series for the precessing Euler angles in the J-frame.
 */
int XLALSimIMRPhenomTPHM_JModes(
  SphHarmTimeSeries **hlmJ,   /**< [out] Modes in the intertial J0=z frame*/
  REAL8TimeSeries **alphaTS,  /**< [out] Precessing Euler angle alpha */
  REAL8TimeSeries **cosbetaTS,   /**< [out] Precessing Euler angle beta */
  REAL8TimeSeries **gammaTS,  /**< [out] Precessing Euler angle gamma */
  REAL8 *af,                  /**< [out] Final spin */
  REAL8 m1_SI,                /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,                /**< Mass of companion 2 (kg) */
  REAL8 chi1x,                /**< x component of primary spin*/
  REAL8 chi1y,                /**< y component of primary spin*/
  REAL8 chi1z,                /**< z component of primary spin */
  REAL8 chi2x,                /**< x component of secondary spin*/
  REAL8 chi2y,                /**< y component of secondary spin*/
  REAL8 chi2z,                /**< z component of secondary spin */
  REAL8 distance,             /**< Luminosity distance (m) */
  REAL8 inclination,          /**< inclination (in rad) */
  REAL8 deltaT,               /**< sampling interval (s) */
  REAL8 fmin,               /**< starting GW frequency (Hz) */
  REAL8 fRef,               /**< reference GW frequency (Hz) */
  REAL8 phiRef,               /**< reference orbital phase (rad) */
  LALDict *lalParams,       /**< LAL dictionary containing accessory parameters */
  UINT4 only22              /**< Flag for calling only IMRPhenomTP (dominant 22 coprec mode only) */
  )
{

  UINT4 status;

  /* Compute co-precessing modes and Euler angles */
  *hlmJ = NULL;
  status = XLALSimIMRPhenomTPHM_CoprecModes(hlmJ, alphaTS, cosbetaTS, gammaTS, af, m1_SI, m2_SI, chi1x,chi1y,chi1z,chi2x,chi2y,chi2z,distance,inclination,deltaT,fmin,fRef,phiRef,lalParams,only22);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: function XLALSimIMRPhenomTPHM_JModes has failed.");

  /* Rotate coprecessing modes to J frame */

  PhenomT_precomputed_sqrt *SQRT;
  SQRT  = XLALMalloc(sizeof(PhenomT_precomputed_sqrt));
  IMRPhenomTPHM_SetPrecomputedSqrt(SQRT);

  status = PhenomTPHM_RotateModes(*hlmJ, *gammaTS, *cosbetaTS, *alphaTS, SQRT);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: Internal function PhenomTPHM_RotateModes has failed.");

  LALFree(SQRT);

  return status;
}

/**
 * Routine to compute a list of the time-domain modes of the IMRPhenomTPHM model in the co-precessing frame.
 * It also outputs the time series for the precessing Euler angles in the J-frame.
 */
int XLALSimIMRPhenomTPHM_CoprecModes(
  SphHarmTimeSeries **hlmJ,   /**< [out] Modes in the intertial J0=z frame*/
  REAL8TimeSeries **alphaTS,  /**< [out] Precessing Euler angle alpha */
  REAL8TimeSeries **cosbetaTS,   /**< [out] Precessing Euler angle beta */
  REAL8TimeSeries **gammaTS,  /**< [out] Precessing Euler angle gamma */
  REAL8 *af,                  /**< [out] Final spin */
  REAL8 m1_SI,                /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,                /**< Mass of companion 2 (kg) */
  REAL8 chi1x,                /**< x component of primary spin*/
  REAL8 chi1y,                /**< y component of primary spin*/
  REAL8 chi1z,                /**< z component of primary spin */
  REAL8 chi2x,                /**< x component of secondary spin*/
  REAL8 chi2y,                /**< y component of secondary spin*/
  REAL8 chi2z,                /**< z component of secondary spin */
  REAL8 distance,             /**< Luminosity distance (m) */
  REAL8 inclination,          /**< inclination (in rad) */
  REAL8 deltaT,               /**< sampling interval (s) */
  REAL8 fmin,               /**< starting GW frequency (Hz) */
  REAL8 fRef,               /**< reference GW frequency (Hz) */
  REAL8 phiRef,               /**< reference orbital phase (rad) */
  LALDict *lalParams,       /**< LAL dictionary containing accessory parameters */
  UINT4 only22              /**< Flag for calling only IMRPhenomTP (dominant 22 coprec mode only) */
  )
{

  /* Sanity checks */
  if(fRef  <  0.0) { XLAL_ERROR(XLAL_EDOM, "fRef_In must be positive or set to 0 to ignore.\n");  }
  if(deltaT   <= 0.0) { XLAL_ERROR(XLAL_EDOM, "deltaT must be positive.\n");                         }
  if(m1_SI    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "m1 must be positive.\n");                             }
  if(m2_SI    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "m2 must be positive.\n");                             }
  if(fmin    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "f_min must be positive.\n");                          }
  if(distance <= 0.0) { XLAL_ERROR(XLAL_EDOM, "Distance must be positive and greater than 0.\n");    }
  if(fRef > 0.0 && fRef < fmin){ XLAL_ERROR(XLAL_EDOM, "fRef must be >= f_min or =0 to use f_min.\n"); }

  /* Swap components if m2>m1 */
  if(m1_SI < m2_SI)
  {
    REAL8 auxswap;

    auxswap = m2_SI;
    m2_SI = m1_SI;
    m1_SI = auxswap;

    auxswap = chi2x;
    chi2x = chi1x;
    chi1x = auxswap;

    auxswap = chi2y;
    chi2y = chi1y;
    chi1y = auxswap;

    auxswap = chi2z;
    chi2z = chi1z;
    chi1z = auxswap;
  }

  REAL8 mass_ratio = m1_SI / m2_SI;
  REAL8 chi1L = chi1z;
  REAL8 chi2L = chi2z;

  /*
    Perform a basic sanity check on the region of the parameter space in which model is evaluated. Behaviour is as follows:
      - For mass ratios <= 20.0 and spins <= 0.99: no warning messages.
      - For 200 > mass ratio > 20 and spins <= 0.9: print a warning message that model can be pathological.
      - For mass ratios > 200: throw a hard error that model is not valid.
      - For spins > 0.99: throw a warning that we are extrapolating the model to extremal
  */

  if(mass_ratio > 20.0 && chi1L < 0.9 && m2_SI/LAL_MSUN_SI >= 0.5  ) { XLAL_PRINT_INFO("Warning: Extrapolating outside of Numerical Relativity calibration domain."); }
  if(mass_ratio > 20.0 && (chi1L >= 0.9 || m2_SI/LAL_MSUN_SI < 0.5) ) { XLAL_PRINT_INFO("Warning: Model can be pathological at these parameters"); }
  if(mass_ratio > 200. && fabs(mass_ratio - 200) > 1e-12) { XLAL_ERROR(XLAL_EDOM, "ERROR: Model not valid at mass ratios beyond 200."); } // The 1e-12 is to avoid rounding errors
  if(fabs(chi1L) > 0.99 || fabs(chi2L) > 0.99) { XLAL_PRINT_INFO("Warning: Extrapolating to extremal spins, model is not trusted."); }

  /* Set up model options specified in precessing version */
  UNUSED INT4 precVer, precVerX, EulerRDVersion, EulerGammaVersion, EulerInspVersion, FSVer;
  UNUSED INT4 evolvedFS = 0;
  precVer = XLALSimInspiralWaveformParamsLookupPhenomXPrecVersion(lalParams);
  FSVer = XLALSimInspiralWaveformParamsLookupPhenomXPFinalSpinMod(lalParams);
  UINT4 reconsMerger = XLALSimInspiralWaveformParamsLookupPhenomTPHMMergerVersion(lalParams); // FIXME:Internal option, will be deleted for release
  char digits[6];
  INT4 lendigits = snprintf(digits, 6, "%d", precVer);

  if(lendigits == 5)
  {
    sprintf(digits, "%d", precVer);
    sscanf(digits, "%3d%1d%1d", &precVerX, &EulerRDVersion, &EulerGammaVersion);
    EulerInspVersion = 0;
    if(EulerRDVersion>2 || EulerRDVersion<0){XLAL_ERROR(XLAL_EDOM, "Invalid EulerRDVersion. Choose between 0, 1, and 2.\n");}
    if(EulerGammaVersion>1 || EulerGammaVersion<0){XLAL_ERROR(XLAL_EDOM, "Invalid EulerGammaVersion. Choose between 0, 1.\n");}

  }
  else if(lendigits==3 && precVer==300)
  {
    EulerInspVersion = 1;
    precVerX = 102; // Default version for IMRPhenomX struct
    if(FSVer==4)
    {
      evolvedFS = 1;
      FSVer = 2;
    }
  }
  else if(lendigits==3 && precVer!=300)
  {
    EulerInspVersion = 0;
    precVerX = precVer;
    EulerRDVersion = 1; // Default ringdown version for analytical angles
    EulerGammaVersion = 0;
  }
  else
  {
    XLAL_ERROR(XLAL_EDOM, "Incorrect format/values for PhenomXPrecVersion option.\n");
  }


  INT4 status;

  /* Use an auxiliar laldict to not overwrite the input argument */
  LALDict *lalParams_aux;
  if (lalParams == NULL)
  {
      lalParams_aux = XLALCreateDict();
  }
  else{
      lalParams_aux = XLALDictDuplicate(lalParams);
  }
  /* Setup mode array */
  UINT4 LMAX;
  LALValue *ModeArray = NULL;

  if(only22==0)
  {
    XLAL_CHECK(check_input_mode_array_THM(lalParams_aux) == XLAL_SUCCESS, XLAL_EFAULT, "Not available mode chosen.\n");
    lalParams_aux = IMRPhenomTHM_setup_mode_array(lalParams_aux);
    ModeArray = XLALSimInspiralWaveformParamsLookupModeArray(lalParams_aux);
  }
  else if(only22==1)
  {
    XLAL_CHECK(check_input_mode_array_22_THM(lalParams_aux) == XLAL_SUCCESS, XLAL_EFAULT, "Not available mode chosen.\n");
    ModeArray = XLALSimInspiralWaveformParamsLookupModeArray(lalParams_aux);
    if(ModeArray==NULL)
    {
      ModeArray = XLALSimInspiralCreateModeArray();
      /* Activate (2,2) and (2,-2) modes*/
      XLALSimInspiralModeArrayActivateMode(ModeArray, 2, 2);
      XLALSimInspiralModeArrayActivateMode(ModeArray, 2, -2);
      /* Insert ModeArray into lalParams */
      XLALSimInspiralWaveformParamsInsertModeArray(lalParams_aux, ModeArray);
    }
    LMAX = 2;
  }

  /* Select maximum L for which modes have to be activated */

  if(XLALSimInspiralModeArrayIsModeActive(ModeArray, 5, 5)==1||XLALSimInspiralModeArrayIsModeActive(ModeArray, 5, -5)==1)
  {
    LMAX = 5;
  }
  else if(XLALSimInspiralModeArrayIsModeActive(ModeArray, 4, 4)==1||XLALSimInspiralModeArrayIsModeActive(ModeArray, 4, -4)==1)
  {
    LMAX = 4;
  }
  else if(XLALSimInspiralModeArrayIsModeActive(ModeArray, 3, 3)==1||XLALSimInspiralModeArrayIsModeActive(ModeArray, 3, -3)==1)
  {
    LMAX = 3;
  }
  else{
    LMAX=2;}

  /* Add PrecVersion and FinalSpin options */
  status = XLALSimInspiralWaveformParamsInsertPhenomXPrecVersion(lalParams_aux, precVerX);
  status = XLALSimInspiralWaveformParamsInsertPhenomXPFinalSpinMod(lalParams_aux, FSVer);

  /* Initialise waveform struct and phase/frequency struct for the 22, needed for any mode */
  IMRPhenomTWaveformStruct *pWF;
  pWF    = XLALMalloc(sizeof(IMRPhenomTWaveformStruct));
  status = IMRPhenomTSetWaveformVariables(pWF, m1_SI, m2_SI, chi1L, chi2L, distance, deltaT, fmin, fRef, phiRef, lalParams_aux);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: Internal function IMRPhenomTSetWaveformVariables has failed.");

  /* Initialize IMR PhenomX Waveform struct and check that it initialized correctly */
  IMRPhenomXWaveformStruct *pWFX;
  pWFX    = XLALMalloc(sizeof(IMRPhenomXWaveformStruct));
  status = IMRPhenomXSetWaveformVariables(pWFX, m1_SI, m2_SI, chi1L, chi2L, 0.1, fRef, phiRef, fmin, 1./deltaT, distance, inclination, lalParams_aux, 0);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXSetWaveformVariables failed.\n");

  /* Initialize IMR PhenomX Precession struct and check that it generated successfully */
  IMRPhenomXPrecessionStruct *pPrec;
  pPrec  = XLALMalloc(sizeof(IMRPhenomXPrecessionStruct));
  status = IMRPhenomXGetAndSetPrecessionVariables(pWFX, pPrec, m1_SI, m2_SI, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z, lalParams_aux, 0);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXSetPrecessionVariables failed.\n");

  pWF->afinal_prec = pWFX->afinal; // Set final spin to its precessing value
  if(reconsMerger!=1)
  {
        pWF->afinal = pWFX->afinal; // FIXME: Internal option, will be removed before merging to master.
  }

  IMRPhenomTPhase22Struct *pPhase;
  pPhase = XLALMalloc(sizeof(IMRPhenomTPhase22Struct));
  status   = IMRPhenomTSetPhase22Coefficients(pPhase, pWF);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: Internal function IMRPhenomTSetPhase22Coefficients has failed.");

  /* Length of the required sequence and length of the different inspiral regions of the 22 phase. */
  size_t length = pPhase->wflength;
  size_t length_insp_early = pPhase->wflength_insp_early;
  size_t length_insp_late = pPhase->wflength_insp_late;
  size_t length_insp = length_insp_early + length_insp_late;

  /*Initialize time */
  REAL8 t, thetabar, theta, w22, ph22;
  REAL8 factheta = pow(5.0,1./8);

  /* Initialize REAL8 sequences for storing the phase and the frequency of the 22 */
  UNUSED REAL8Sequence *phi22 = NULL;
  REAL8Sequence *xorb = NULL;
  xorb = XLALCreateREAL8Sequence(length);
  REAL8Sequence *timesVec = NULL;
  timesVec = XLALCreateREAL8Sequence(length);

  phi22 = XLALCreateREAL8Sequence(length);

  LIGOTimeGPS ligotimegps_zero = LIGOTIMEGPSZERO; // = {0,0}
  XLALGPSAdd(&ligotimegps_zero, pPhase->tminSec);

  /* Compute 22 phase and x=(0.5*omega_22)^(2/3), needed for all modes.
  Depending on the inspiral region, theta is defined with a fitted t0 parameter or with t0=0. */

  if(pWF->inspVersion!=0) // If reconstruction is non-default, this will compute frequency, phase and PN expansion parameter x from TaylorT3 with fitted t0 for the early inspiral region-
  {
    for(UINT4 jdx = 0; jdx < length_insp_early; jdx++)
    {
      t = pPhase->tmin + jdx*pWF->dtM;
      thetabar = pow(pWF->eta*(pPhase->tt0-t),-1./8);
      theta = factheta*thetabar;
      w22 = IMRPhenomTomega22(t, theta, pWF, pPhase);
      xorb->data[jdx] = pow(0.5*w22,2./3);
      ph22 = IMRPhenomTPhase22(t, thetabar, pWF, pPhase);
      phi22->data[jdx] = ph22;
    }

    for(UINT4 jdx = length_insp_early; jdx < length_insp; jdx++) // For times later than the early-late inspiral boundary, it computes phase, frequency and x with the extended TaylorT3 (with tt0=0)
    {
      t = pPhase->tmin + jdx*pWF->dtM;
      thetabar = pow(-pWF->eta*t,-1./8);
      theta = factheta*thetabar;
      w22 = IMRPhenomTomega22(t, theta, pWF, pPhase);
      xorb->data[jdx] = pow(0.5*w22,2./3);
      ph22 = IMRPhenomTPhase22(t, thetabar, pWF, pPhase);
      phi22->data[jdx] = ph22;
    }
  }

  else // If default reconstruction, only one inspiral region with extended TaylorT3 (with tt0=0) is employed up to the inspiral-merger boundary time
  {
      for(UINT4 jdx = 0; jdx < length_insp; jdx++)
      {
        t = pPhase->tmin + jdx*pWF->dtM;
        thetabar = pow(-pWF->eta*t,-1./8);
        theta = factheta*thetabar;
        w22 = IMRPhenomTomega22(t, theta, pWF, pPhase);
        xorb->data[jdx] = pow(0.5*w22,2./3);
        ph22 = IMRPhenomTPhase22(t, thetabar, pWF, pPhase);
        phi22->data[jdx] = ph22;
      }

  }

  /* During the 22 phase merger region, phase and x are also computed because the higher modes end the inspiral at a fixed time,
  which can occur later than the end of the 22 inspiral. */
  for(UINT4 jdx = length_insp; jdx < length; jdx++)
  {
    t = pPhase->tmin + jdx*pWF->dtM;
    w22 = IMRPhenomTomega22(t, 0.0, pWF, pPhase);
    xorb->data[jdx] = pow(0.5*w22,2./3);
    ph22 = IMRPhenomTPhase22(t, 0.0, pWF, pPhase);
    phi22->data[jdx] = ph22;
  }

  /****** Compute Euler angles *******/

  *alphaTS = XLALCreateREAL8TimeSeries( "alphaPtoJ", &ligotimegps_zero, 0.0, deltaT, &lalStrainUnit,length);
  *cosbetaTS = XLALCreateREAL8TimeSeries( "cosbetaPtoJ", &ligotimegps_zero, 0.0, deltaT, &lalStrainUnit,length);
  *gammaTS = XLALCreateREAL8TimeSeries( "epsilonPtoJ", &ligotimegps_zero, 0.0, deltaT, &lalStrainUnit,length);

  /* Check validity of fRef and fmin. If fmin>f_peak, proceed with analytical angles with EulerRDVersion=1, to just have ringdown waveforms.
  If f_min<f_peak but fRef>fpeak, and numerical angles are selected, then default fRef to some maximal allowed value, close to the peak, warning the user. */

  if(pPhase->tmin>=-10*pWF->dtM)
  {
    EulerInspVersion = 0;
    EulerRDVersion = 1;
    EulerGammaVersion = 0;

    XLAL_PRINT_WARNING("Waveform only contains ringdown, version 300 not meaningful. Defaulting to analytical ringdown angles (EulerInspVersion=0, EulerRDVersion=1.)\n");
  }
  else if(pPhase->tRef>=-8*pWF->dtM)
  {
    pPhase->tRef = MAX(-8*pWF->dtM, pPhase->tmin);
    XLAL_PRINT_WARNING("Waveform contains pre-peak cycles but reference frequency is specified post-peak, which may be not meaningful. Reference frequency will be set close before the peak. \n");
  }

  if(EulerInspVersion==1)
  {
    REAL8 af_evolved;
    status = IMRPhenomTPHM_NumericalEulerAngles(alphaTS,cosbetaTS,gammaTS, &af_evolved, xorb, pPhase->dtM,m1_SI, m2_SI, chi1x,chi1y,chi1z,chi2x,chi2y,chi2z,pPhase->tmin, pWF, pPhase, pPrec);
    XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: Internal function IMRPhenomTPHM_NumericalEulerAngles has failed.");

    if(evolvedFS==1)
    {
      pWF->afinal_prec = af_evolved;
      if(!(reconsMerger==1))
      {
        pWF->afinal = af_evolved;
      }

      IMRPhenomTPhase22Struct *pPhase2;
      pPhase2 = XLALMalloc(sizeof(IMRPhenomTPhase22Struct));
      status   = IMRPhenomTSetPhase22Coefficients(pPhase2, pWF);
      XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: Internal function IMRPhenomTSetPhase22Coefficients has failed.");

      for(UINT4 jdx = length_insp; jdx < length; jdx++)
      {
        t = pPhase->tmin + jdx*pWF->dtM;

        w22 = IMRPhenomTomega22(t, 0.0, pWF, pPhase2);
        xorb->data[jdx] = pow(0.5*w22,2./3);
        ph22 = IMRPhenomTPhase22(t, 0.0, pWF, pPhase2);

        phi22->data[jdx] = ph22;
      }

      LALFree(pPhase2);
      }
  }
  else
  {
    status = PNAnalyticalInspiralEulerAngles(alphaTS, cosbetaTS, gammaTS, xorb, pWF, pPhase, pWFX, pPrec, EulerRDVersion, EulerGammaVersion);
    XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: Internal function PNAnalyticalInspiralEulerAngles has failed.");
  }

  *af = pWF->afinal;

  /***** Loop over modes ******/

  *hlmJ = NULL;

  INT4 posMode, negMode;

  for (UINT4 ell = 2; ell <= LMAX; ell++)
  {
    for (UINT4 emm = 0; emm <= ell; emm++)
    { /* Loop over only positive m is intentional.*/

      posMode = XLALSimInspiralModeArrayIsModeActive(ModeArray, ell, emm);
      negMode = XLALSimInspiralModeArrayIsModeActive(ModeArray, ell, -emm);
      if ( posMode != 1 && negMode != 1)
      { /* skip mode */
        COMPLEX16TimeSeries *tmp_mode = XLALCreateCOMPLEX16TimeSeries( "hlm", &ligotimegps_zero, 0.0, deltaT, &lalStrainUnit,length);
        memset( tmp_mode->data->data, 0, tmp_mode->data->length * sizeof(COMPLEX16) );
        *hlmJ = XLALSphHarmTimeSeriesAddMode(*hlmJ, tmp_mode, ell, emm);
        *hlmJ = XLALSphHarmTimeSeriesAddMode(*hlmJ, tmp_mode, ell, -emm);
        XLALDestroyCOMPLEX16TimeSeries(tmp_mode);
        continue;
      } /* else: generate mode */

      COMPLEX16TimeSeries *tmp_mode = NULL;

      status = LALSimIMRPhenomTHM_OneMode(&tmp_mode, pWF, pPhase, phi22, xorb, ell, emm);

      if(posMode==1) /* Modes are generated only for positive m. If positive m is requested, simply add to the SphHarmTimeSeries structure */
      {
        *hlmJ = XLALSphHarmTimeSeriesAddMode(*hlmJ, tmp_mode, ell, emm);
      }

      if(negMode==1) /* If negative m is requested, transform from positive m mode using symmetry property for aligned spin systems */
      {
        for(UINT4 jdx = 0; jdx < length; jdx++)
        {
          ((*tmp_mode).data)->data[jdx] = pow(-1,ell)*conj(((*tmp_mode).data)->data[jdx]);
        }
        *hlmJ = XLALSphHarmTimeSeriesAddMode(*hlmJ, tmp_mode, ell, -emm);
      }

      XLALDestroyCOMPLEX16TimeSeries(tmp_mode); /* Destroy temporary structure once the desired mode is added */
    }

  }

  /*Free structures and destroy sequences and dict */
  XLALDestroyValue(ModeArray);
  LALFree(pPhase);
  LALFree(pWF);
  LALFree(pWFX);
  LALFree(pPrec);

  XLALDestroyREAL8Sequence(phi22);

  XLALDestroyREAL8Sequence(xorb);
  XLALDestroyREAL8Sequence(timesVec);

  XLALDestroyDict(lalParams_aux);

  return status;
}

/** @} */
/** @} */

/* ==========================================================================
 * Model-time functions: IMRPhenomTP / IMRPhenomTPHM modes and polarisations at
 * the times of a user-supplied TimeArray (units of M, t = 0 at the peak of
 * the non-precessing 22 amplitude, any spacing).
 *
 * Same pipeline as XLALSimIMRPhenomTPHM_CoprecModes -> _JModes -> _L0Modes,
 * with the co-precessing modes evaluated at the times of TimeArray (as in
 * XLALSimIMRPhenomTHMModelTimes) and the Euler angles interpolated from the
 * standard time grid (see IMRPhenomTPHM_EulerAnglesModelTimes).
 *
 * - *ModelTimes:          every time of TimeArray must lie within the time span
 *                         of the standard waveform (the Euler angles are only
 *                         defined there), otherwise an error is raised.
 * - *ModelTimesTruncated: samples outside the time span of the standard
 *                         waveform are set to 0.
 *
 * IMRPhenomTP: dominant co-precessing (2,+-2) modes. IMRPhenomTPHM: modes
 * selected in the ModeArray of lalParams, among (2,+-2), (2,+-1), (3,+-3),
 * (4,+-4), (5,+-5) (all of them by default).
 * ========================================================================== */

/* Co-precessing modes at the times of TimeArray, plus the Euler angles at the
times of TimeArray and on the standard grid. As XLALSimIMRPhenomTPHM_CoprecModes. */
static int IMRPhenomTPHM_CoprecModesModelTimes_engine(
  SphHarmTimeSeries **hlmJ,
  REAL8TimeSeries **alphaTS, REAL8TimeSeries **cosbetaTS, REAL8TimeSeries **gammaTS,
  REAL8TimeSeries **alphaNatTS, REAL8TimeSeries **cosbetaNatTS, REAL8TimeSeries **gammaNatTS,
  REAL8 *af,
  REAL8 m1_SI, REAL8 m2_SI,
  REAL8 chi1x, REAL8 chi1y, REAL8 chi1z,
  REAL8 chi2x, REAL8 chi2y, REAL8 chi2z,
  REAL8 distance, REAL8 inclination, REAL8 deltaT,
  REAL8 fmin, REAL8 fRef, REAL8 phiRef,
  REAL8Sequence *TimeArray, LALDict *lalParams,
  UINT4 only22, UINT4 truncate)
{
  /* Sanity checks */
  if(fRef  <  0.0) { XLAL_ERROR(XLAL_EDOM, "fRef_In must be positive or set to 0 to ignore.\n");  }
  if(deltaT   <= 0.0) { XLAL_ERROR(XLAL_EDOM, "deltaT must be positive.\n");                         }
  if(m1_SI    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "m1 must be positive.\n");                             }
  if(m2_SI    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "m2 must be positive.\n");                             }
  if(fmin    <= 0.0) { XLAL_ERROR(XLAL_EDOM, "f_min must be positive.\n");                          }
  if(distance <= 0.0) { XLAL_ERROR(XLAL_EDOM, "Distance must be positive and greater than 0.\n");    }
  if(fRef > 0.0 && fRef < fmin){ XLAL_ERROR(XLAL_EDOM, "fRef must be >= f_min or =0 to use f_min.\n"); }
  XLAL_CHECK(NULL != TimeArray, XLAL_EFAULT, "Error: TimeArray must not be NULL.\n");
  XLAL_CHECK(TimeArray->length > 0, XLAL_EDOM, "Error: TimeArray is empty.\n");

  /* Swap components if m2>m1 */
  if(m1_SI < m2_SI)
  {
    REAL8 auxswap;
    auxswap = m2_SI; m2_SI = m1_SI; m1_SI = auxswap;
    auxswap = chi2x; chi2x = chi1x; chi1x = auxswap;
    auxswap = chi2y; chi2y = chi1y; chi1y = auxswap;
    auxswap = chi2z; chi2z = chi1z; chi1z = auxswap;
  }

  REAL8 mass_ratio = m1_SI / m2_SI;
  REAL8 chi1L = chi1z;
  REAL8 chi2L = chi2z;

  if(mass_ratio > 20.0 && chi1L < 0.9 && m2_SI/LAL_MSUN_SI >= 0.5  ) { XLAL_PRINT_INFO("Warning: Extrapolating outside of Numerical Relativity calibration domain."); }
  if(mass_ratio > 20.0 && (chi1L >= 0.9 || m2_SI/LAL_MSUN_SI < 0.5) ) { XLAL_PRINT_INFO("Warning: Model can be pathological at these parameters"); }
  if(mass_ratio > 200. && fabs(mass_ratio - 200) > 1e-12) { XLAL_ERROR(XLAL_EDOM, "ERROR: Model not valid at mass ratios beyond 200."); }
  if(fabs(chi1L) > 0.99 || fabs(chi2L) > 0.99) { XLAL_PRINT_INFO("Warning: Extrapolating to extremal spins, model is not trusted."); }

  /* Model options specified in the precessing version */
  UNUSED INT4 precVer, precVerX, EulerRDVersion, EulerGammaVersion, EulerInspVersion, FSVer;
  UNUSED INT4 evolvedFS = 0;
  precVer = XLALSimInspiralWaveformParamsLookupPhenomXPrecVersion(lalParams);
  FSVer = XLALSimInspiralWaveformParamsLookupPhenomXPFinalSpinMod(lalParams);
  UINT4 reconsMerger = XLALSimInspiralWaveformParamsLookupPhenomTPHMMergerVersion(lalParams);

  char digits[6];
  INT4 lendigits = snprintf(digits, 6, "%d", precVer);
  if(lendigits == 5)
  {
    sprintf(digits, "%d", precVer);
    sscanf(digits, "%3d%1d%1d", &precVerX, &EulerRDVersion, &EulerGammaVersion);
    EulerInspVersion = 0;
    if(EulerRDVersion>2 || EulerRDVersion<0){XLAL_ERROR(XLAL_EDOM, "Invalid EulerRDVersion. Choose between 0, 1, and 2.\n");}
    if(EulerGammaVersion>1 || EulerGammaVersion<0){XLAL_ERROR(XLAL_EDOM, "Invalid EulerGammaVersion. Choose between 0, 1.\n");}
  }
  else if(lendigits==3 && precVer==300)
  {
    EulerInspVersion = 1;
    precVerX = 102;
    if(FSVer==4)
    {
      evolvedFS = 1;
      FSVer = 2;
    }
  }
  else if(lendigits==3 && precVer!=300)
  {
    EulerInspVersion = 0;
    precVerX = precVer;
    EulerRDVersion = 1;
    EulerGammaVersion = 0;
  }
  else
  {
    XLAL_ERROR(XLAL_EDOM, "Incorrect format/values for PhenomXPrecVersion option.\n");
  }

  INT4 status;

  /* Use an auxiliar laldict to not overwrite the input argument */
  LALDict *lalParams_aux = (lalParams == NULL) ? XLALCreateDict() : XLALDictDuplicate(lalParams);

  /* Mode array */
  UINT4 LMAX;
  LALValue *ModeArray = NULL;
  if(only22==0)
  {
    XLAL_CHECK(check_input_mode_array_THM(lalParams_aux) == XLAL_SUCCESS, XLAL_EFAULT, "Not available mode chosen.\n");
    lalParams_aux = IMRPhenomTHM_setup_mode_array(lalParams_aux);
    ModeArray = XLALSimInspiralWaveformParamsLookupModeArray(lalParams_aux);
  }
  else
  {
    XLAL_CHECK(check_input_mode_array_22_THM(lalParams_aux) == XLAL_SUCCESS, XLAL_EFAULT, "Not available mode chosen.\n");
    ModeArray = XLALSimInspiralWaveformParamsLookupModeArray(lalParams_aux);
    if(ModeArray==NULL)
    {
      ModeArray = XLALSimInspiralCreateModeArray();
      XLALSimInspiralModeArrayActivateMode(ModeArray, 2, 2);
      XLALSimInspiralModeArrayActivateMode(ModeArray, 2, -2);
      XLALSimInspiralWaveformParamsInsertModeArray(lalParams_aux, ModeArray);
    }
  }

  if(XLALSimInspiralModeArrayIsModeActive(ModeArray, 5, 5)==1||XLALSimInspiralModeArrayIsModeActive(ModeArray, 5, -5)==1)
  {
    LMAX = 5;
  }
  else if(XLALSimInspiralModeArrayIsModeActive(ModeArray, 4, 4)==1||XLALSimInspiralModeArrayIsModeActive(ModeArray, 4, -4)==1)
  {
    LMAX = 4;
  }
  else if(XLALSimInspiralModeArrayIsModeActive(ModeArray, 3, 3)==1||XLALSimInspiralModeArrayIsModeActive(ModeArray, 3, -3)==1)
  {
    LMAX = 3;
  }
  else
  {
    LMAX = 2;
  }

  status = XLALSimInspiralWaveformParamsInsertPhenomXPrecVersion(lalParams_aux, precVerX);
  status = XLALSimInspiralWaveformParamsInsertPhenomXPFinalSpinMod(lalParams_aux, FSVer);

  /* Waveform structs */
  IMRPhenomTWaveformStruct *pWF;
  pWF    = XLALMalloc(sizeof(IMRPhenomTWaveformStruct));
  status = IMRPhenomTSetWaveformVariables(pWF, m1_SI, m2_SI, chi1L, chi2L, distance, deltaT, fmin, fRef, phiRef, lalParams_aux);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: Internal function IMRPhenomTSetWaveformVariables has failed.");

  IMRPhenomXWaveformStruct *pWFX;
  pWFX    = XLALMalloc(sizeof(IMRPhenomXWaveformStruct));
  status = IMRPhenomXSetWaveformVariables(pWFX, m1_SI, m2_SI, chi1L, chi2L, 0.1, fRef, phiRef, fmin, 1./deltaT, distance, inclination, lalParams_aux, 0);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXSetWaveformVariables failed.\n");

  IMRPhenomXPrecessionStruct *pPrec;
  pPrec  = XLALMalloc(sizeof(IMRPhenomXPrecessionStruct));
  status = IMRPhenomXGetAndSetPrecessionVariables(pWFX, pPrec, m1_SI, m2_SI, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z, lalParams_aux, 0);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXSetPrecessionVariables failed.\n");

  pWF->afinal_prec = pWFX->afinal;
  if(reconsMerger!=1)
  {
        pWF->afinal = pWFX->afinal;
  }

  IMRPhenomTPhase22Struct *pPhase;
  pPhase = XLALMalloc(sizeof(IMRPhenomTPhase22Struct));
  status   = IMRPhenomTSetPhase22Coefficients(pPhase, pWF);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: Internal function IMRPhenomTSetPhase22Coefficients has failed.");

  /* 22 phase and x at the times of TimeArray */
  size_t length = TimeArray->length;
  REAL8Sequence *phi22 = XLALCreateREAL8Sequence(length);
  REAL8Sequence *xorb  = XLALCreateREAL8Sequence(length);
  status = LALSimIMRPhenomT_ModelTimes_Phase22(phi22, xorb, TimeArray, pWF, pPhase);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: LALSimIMRPhenomT_ModelTimes_Phase22 has failed.");

  /* Euler angles (pPhase must hold the standard tmin and wflength here) */
  status = IMRPhenomTPHM_EulerAnglesModelTimes(alphaTS, cosbetaTS, gammaTS, alphaNatTS, cosbetaNatTS, gammaNatTS, af, TimeArray,
             m1_SI, m2_SI, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z, deltaT,
             pWF, pPhase, pWFX, pPrec, EulerInspVersion, EulerRDVersion, EulerGammaVersion,
             evolvedFS, reconsMerger, truncate);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: Internal function IMRPhenomTPHM_EulerAnglesModelTimes has failed.");

  /* Final spin from the evolved spins: merger-ringdown 22 phase and x with the updated coefficients */
  if(EulerInspVersion==1 && evolvedFS==1 && pPhase->tmin < -10*pWF->dtM)
  {
    IMRPhenomTPhase22Struct *pPhase2;
    pPhase2 = XLALMalloc(sizeof(IMRPhenomTPhase22Struct));
    status   = IMRPhenomTSetPhase22Coefficients(pPhase2, pWF);
    XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: Internal function IMRPhenomTSetPhase22Coefficients has failed.");
    for(UINT4 jdx = 0; jdx < length; jdx++)
    {
      REAL8 t = TimeArray->data[jdx];
      if(t > pPhase->tCut22)
      {
        xorb->data[jdx] = pow(0.5*IMRPhenomTomega22(t, 0.0, pWF, pPhase2),2./3);
        phi22->data[jdx] = IMRPhenomTPhase22(t, 0.0, pWF, pPhase2);
      }
    }
    LALFree(pPhase2);
  }

  /* Epoch of the output: first requested time. pPhase->tmin keeps the start of the standard waveform. */
  pPhase->tminSec = pWF->M_sec*TimeArray->data[0];
  LIGOTimeGPS epoch = LIGOTIMEGPSZERO;
  XLALGPSAdd(&epoch, pPhase->tminSec);

  /***** Loop over modes ******/
  *hlmJ = NULL;
  INT4 posMode, negMode;
  for (UINT4 ell = 2; ell <= LMAX; ell++)
  {
    for (UINT4 emm = 0; emm <= ell; emm++)
    { /* Loop over only positive m is intentional.*/
      posMode = XLALSimInspiralModeArrayIsModeActive(ModeArray, ell, emm);
      negMode = XLALSimInspiralModeArrayIsModeActive(ModeArray, ell, -emm);
      if ( posMode != 1 && negMode != 1)
      { /* zero mode: the rotation mixes the m's of a given l */
        COMPLEX16TimeSeries *tmp_mode = XLALCreateCOMPLEX16TimeSeries("hlm", &epoch, 0.0, deltaT, &lalStrainUnit, length);
        memset(tmp_mode->data->data, 0, tmp_mode->data->length * sizeof(COMPLEX16));
        *hlmJ = XLALSphHarmTimeSeriesAddMode(*hlmJ, tmp_mode, ell, emm);
        *hlmJ = XLALSphHarmTimeSeriesAddMode(*hlmJ, tmp_mode, ell, -emm);
        XLALDestroyCOMPLEX16TimeSeries(tmp_mode);
        continue;
      }

      COMPLEX16TimeSeries *tmp_mode = NULL;
      if(truncate==1)
      {
        status = LALSimIMRPhenomTHM_OneModeModelTimesTruncated(&tmp_mode, pWF, pPhase, phi22, xorb, TimeArray, ell, emm);
      }
      else
      {
        status = LALSimIMRPhenomTHM_OneModeModelTimes(&tmp_mode, pWF, pPhase, phi22, xorb, TimeArray, ell, emm);
      }
      XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: generating mode (%d,%d) has failed.\n", ell, emm);

      if(posMode==1)
      {
        *hlmJ = XLALSphHarmTimeSeriesAddMode(*hlmJ, tmp_mode, ell, emm);
      }
      if(negMode==1) /* negative m from the aligned spin symmetry (co-precessing frame) */
      {
        for(UINT4 jdx = 0; jdx < length; jdx++)
        {
          tmp_mode->data->data[jdx] = pow(-1,ell)*conj(tmp_mode->data->data[jdx]);
        }
        *hlmJ = XLALSphHarmTimeSeriesAddMode(*hlmJ, tmp_mode, ell, -emm);
      }
      XLALDestroyCOMPLEX16TimeSeries(tmp_mode);
    }
  }

  XLALDestroyValue(ModeArray);
  LALFree(pPhase);
  LALFree(pWF);
  LALFree(pWFX);
  LALFree(pPrec);
  XLALDestroyREAL8Sequence(phi22);
  XLALDestroyREAL8Sequence(xorb);
  XLALDestroyDict(lalParams_aux);

  return status;
}

/* Euler angles (co-precessing -> J frame) at the times of TimeArray, as used by the functions below. */
static int IMRPhenomTPHM_EulerAnglesModelTimes_engine(
  REAL8TimeSeries **alphaTS, REAL8TimeSeries **cosbetaTS, REAL8TimeSeries **gammaTS,
  REAL8 m1_SI, REAL8 m2_SI,
  REAL8 chi1x, REAL8 chi1y, REAL8 chi1z,
  REAL8 chi2x, REAL8 chi2y, REAL8 chi2z,
  REAL8 distance, REAL8 deltaT, REAL8 fmin, REAL8 fRef, REAL8 phiRef,
  REAL8Sequence *TimeArray, LALDict *lalParams, UINT4 truncate)
{
  XLAL_CHECK(NULL != alphaTS && NULL != cosbetaTS && NULL != gammaTS, XLAL_EFAULT);
  XLAL_CHECK(*alphaTS == NULL && *cosbetaTS == NULL && *gammaTS == NULL, XLAL_EFAULT);

  SphHarmTimeSeries *hlm = NULL;
  REAL8 af;
  UINT4 status = IMRPhenomTPHM_CoprecModesModelTimes_engine(&hlm, alphaTS, cosbetaTS, gammaTS, NULL, NULL, NULL, &af,
             m1_SI, m2_SI, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z,
             distance, 0.0, deltaT, fmin, fRef, phiRef, TimeArray, lalParams, 0, truncate);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: function IMRPhenomTPHM_CoprecModesModelTimes_engine has failed.");
  XLALDestroySphHarmTimeSeries(hlm);

  return status;
}

/** Euler angles alpha, cos(beta), gamma of IMRPhenomTP/TPHM (co-precessing -> J frame) at the times of TimeArray,
interpolated from the standard time grid as in the TP/TPHM ModelTimes functions (see IMRPhenomTPHM_EulerAnglesModelTimes).
Same angles as returned by XLALSimIMRPhenomTPHM_JModes on the standard grid. */
int XLALSimIMRPhenomTPHMEulerAnglesModelTimes(
  REAL8TimeSeries **alphaTS,   /**< [out] alpha at the times of TimeArray */
  REAL8TimeSeries **cosbetaTS, /**< [out] cos(beta) at the times of TimeArray */
  REAL8TimeSeries **gammaTS,   /**< [out] gamma at the times of TimeArray */
  REAL8 m1_SI,      /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,      /**< Mass of companion 2 (kg) */
  REAL8 chi1x,      /**< x component of primary spin*/
  REAL8 chi1y,      /**< y component of primary spin*/
  REAL8 chi1z,      /**< z component of primary spin */
  REAL8 chi2x,      /**< x component of secondary spin*/
  REAL8 chi2y,      /**< y component of secondary spin*/
  REAL8 chi2z,      /**< z component of secondary spin */
  REAL8 distance,   /**< Luminosity distance (m) */
  REAL8 deltaT,     /**< Sampling interval (s) */
  REAL8 fmin,       /**< Starting GW frequency (Hz) */
  REAL8 fRef,       /**< Reference GW frequency (Hz) */
  REAL8 phiRef,     /**< Reference orbital phase (rad) */
  REAL8Sequence *TimeArray, /**< Times (units of M, t=0 at the 22 amplitude peak) */
  LALDict *lalParams  /**< LAL dictionary containing accessory parameters */
  )
{
  return IMRPhenomTPHM_EulerAnglesModelTimes_engine(alphaTS, cosbetaTS, gammaTS, m1_SI, m2_SI, chi1x, chi1y, chi1z,
           chi2x, chi2y, chi2z, distance, deltaT, fmin, fRef, phiRef, TimeArray, lalParams, 0);
}

/** As XLALSimIMRPhenomTPHMEulerAnglesModelTimes; outside the time span of the standard waveform the angles take their boundary values. */
int XLALSimIMRPhenomTPHMEulerAnglesModelTimesTruncated(
  REAL8TimeSeries **alphaTS,   /**< [out] alpha at the times of TimeArray */
  REAL8TimeSeries **cosbetaTS, /**< [out] cos(beta) at the times of TimeArray */
  REAL8TimeSeries **gammaTS,   /**< [out] gamma at the times of TimeArray */
  REAL8 m1_SI,      /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,      /**< Mass of companion 2 (kg) */
  REAL8 chi1x,      /**< x component of primary spin*/
  REAL8 chi1y,      /**< y component of primary spin*/
  REAL8 chi1z,      /**< z component of primary spin */
  REAL8 chi2x,      /**< x component of secondary spin*/
  REAL8 chi2y,      /**< y component of secondary spin*/
  REAL8 chi2z,      /**< z component of secondary spin */
  REAL8 distance,   /**< Luminosity distance (m) */
  REAL8 deltaT,     /**< Sampling interval (s) */
  REAL8 fmin,       /**< Starting GW frequency (Hz) */
  REAL8 fRef,       /**< Reference GW frequency (Hz) */
  REAL8 phiRef,     /**< Reference orbital phase (rad) */
  REAL8Sequence *TimeArray, /**< Times (units of M, t=0 at the 22 amplitude peak) */
  LALDict *lalParams  /**< LAL dictionary containing accessory parameters */
  )
{
  return IMRPhenomTPHM_EulerAnglesModelTimes_engine(alphaTS, cosbetaTS, gammaTS, m1_SI, m2_SI, chi1x, chi1y, chi1z,
           chi2x, chi2y, chi2z, distance, deltaT, fmin, fRef, phiRef, TimeArray, lalParams, 1);
}

/* Euler angles of the global rotation from the J frame to the L0 frame: the angles of the standard grid at the
index of tref (as in XLALSimIMRPhenomTPHM_L0Modes). */
static int IMRPhenomTPHM_JtoL0AnglesModelTimes(
  REAL8 *alphaJtoI, REAL8 *cosbetaJtoI, REAL8 *gammaJtoI,
  REAL8TimeSeries *alphaNat, REAL8TimeSeries *cosbetaNat, REAL8TimeSeries *gammaNat,
  REAL8 m1_SI, REAL8 m2_SI, REAL8 chi1z, REAL8 chi2z,
  REAL8 distance, REAL8 deltaT, REAL8 fmin, REAL8 fRef, REAL8 phiRef, LALDict *lalParams)
{
  UINT4 status;
  IMRPhenomTWaveformStruct *pWF;
  pWF    = XLALMalloc(sizeof(IMRPhenomTWaveformStruct));
  status = IMRPhenomTSetWaveformVariables(pWF, m1_SI, m2_SI, chi1z, chi2z, distance, deltaT, fmin, fRef, phiRef, lalParams);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: Internal function IMRPhenomTSetWaveformVariables has failed.");

  IMRPhenomTPhase22Struct *pPhase;
  pPhase = XLALMalloc(sizeof(IMRPhenomTPhase22Struct));
  status   = IMRPhenomTSetPhase22Coefficients(pPhase, pWF);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: Internal function IMRPhenomTSetPhase22Coefficients has failed.");

  if(pPhase->tmin<0.0 && pPhase->tRef>=0.0)
  {
    pPhase->tRef = MAX(-8*pWF->dtM, pPhase->tmin);
  }

  REAL8 tend = -pPhase->tmin;
  REAL8 tint1 = tend + pPhase->tRef;
  size_t length1 = floor(fabs(tint1)/pWF->dtM);

  *alphaJtoI = alphaNat->data->data[length1];
  *cosbetaJtoI = cosbetaNat->data->data[length1];
  *gammaJtoI = gammaNat->data->data[length1];

  LALFree(pWF);
  LALFree(pPhase);

  return status;
}

/* Modes in the inertial L0 frame at the times of TimeArray. As XLALSimIMRPhenomTPHM_L0Modes
(co-precessing -> J frame with the time dependent Euler angles -> L0 frame with the
Euler angles of the standard grid at the reference time). */
static int IMRPhenomTPHM_L0ModesModelTimes_engine(
  SphHarmTimeSeries **hlmI,
  REAL8 m1_SI, REAL8 m2_SI,
  REAL8 chi1x, REAL8 chi1y, REAL8 chi1z,
  REAL8 chi2x, REAL8 chi2y, REAL8 chi2z,
  REAL8 distance, REAL8 inclination, REAL8 deltaT,
  REAL8 fmin, REAL8 fRef, REAL8 phiRef,
  REAL8Sequence *TimeArray, LALDict *lalParams,
  UINT4 only22, UINT4 truncate)
{
  XLAL_CHECK(NULL != hlmI, XLAL_EFAULT);
  XLAL_CHECK(*hlmI == NULL, XLAL_EFAULT);

  UINT4 status;
  REAL8TimeSeries *alphaTS = NULL, *cosbetaTS = NULL, *gammaTS = NULL;
  REAL8TimeSeries *alphaNat = NULL, *cosbetaNat = NULL, *gammaNat = NULL;
  REAL8 af;

  /* Co-precessing modes, rotated to the J frame (as XLALSimIMRPhenomTPHM_JModes) */
  status = IMRPhenomTPHM_CoprecModesModelTimes_engine(hlmI, &alphaTS, &cosbetaTS, &gammaTS, &alphaNat, &cosbetaNat, &gammaNat, &af,
             m1_SI, m2_SI, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z,
             distance, inclination, deltaT, fmin, fRef, phiRef, TimeArray, lalParams, only22, truncate);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: function IMRPhenomTPHM_CoprecModesModelTimes_engine has failed.");

  PhenomT_precomputed_sqrt *SQRT;
  SQRT  = XLALMalloc(sizeof(PhenomT_precomputed_sqrt));
  IMRPhenomTPHM_SetPrecomputedSqrt(SQRT);

  status = PhenomTPHM_RotateModes(*hlmI, gammaTS, cosbetaTS, alphaTS, SQRT);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: Internal function PhenomTPHM_RotateModes has failed.");

  /* Global rotation from the J frame to the L0 frame: Euler angles at tref */
  REAL8 alphaJtoI, cosbetaJtoI, gammaJtoI;
  status = IMRPhenomTPHM_JtoL0AnglesModelTimes(&alphaJtoI, &cosbetaJtoI, &gammaJtoI, alphaNat, cosbetaNat, gammaNat,
             m1_SI, m2_SI, chi1z, chi2z, distance, deltaT, fmin, fRef, phiRef, lalParams);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: function IMRPhenomTPHM_JtoL0AnglesModelTimes has failed.");

  size_t length = (*hlmI)->mode->data->length;
  status = PhenomTPHM_RotateModes_Global(*hlmI, -alphaJtoI, cosbetaJtoI, -gammaJtoI, length, SQRT);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: Internal function PhenomTPHM_RotateModes_Global has failed.");

  XLALDestroyREAL8TimeSeries(alphaTS);
  XLALDestroyREAL8TimeSeries(cosbetaTS);
  XLALDestroyREAL8TimeSeries(gammaTS);
  XLALDestroyREAL8TimeSeries(alphaNat);
  XLALDestroyREAL8TimeSeries(cosbetaNat);
  XLALDestroyREAL8TimeSeries(gammaNat);
  LALFree(SQRT);

  return status;
}

/* Polarisations of IMRPhenomTP (only22=1) or IMRPhenomTPHM (only22=0) at the times of TimeArray.
As XLALSimIMRPhenomTP/TPHM, including the fallback to IMRPhenomT/THM for aligned spins with version 300. */
static int IMRPhenomTPHM_ModelTimes_polarizations(
  REAL8TimeSeries **hp, REAL8TimeSeries **hc,
  REAL8 m1_SI, REAL8 m2_SI,
  REAL8 chi1x, REAL8 chi1y, REAL8 chi1z,
  REAL8 chi2x, REAL8 chi2y, REAL8 chi2z,
  REAL8 distance, REAL8 inclination, REAL8 deltaT,
  REAL8 fmin, REAL8 fRef, REAL8 phiRef,
  REAL8Sequence *TimeArray, LALDict *lalParams,
  UINT4 only22, UINT4 truncate)
{
  XLAL_CHECK(NULL != hp, XLAL_EFAULT);
  XLAL_CHECK(NULL != hc, XLAL_EFAULT);
  XLAL_CHECK(*hp == NULL, XLAL_EFAULT);
  XLAL_CHECK(*hc == NULL, XLAL_EFAULT);

  UINT4 status;

  REAL8 chi_in_plane = sqrt(pow(chi1x + chi2x,2) + pow(chi1y + chi2y,2));
  if(chi_in_plane<1e-8 && XLALSimInspiralWaveformParamsLookupPhenomXPrecVersion(lalParams)==300)
  {
    if(only22==1)
    {
      return (truncate==1) ?
        XLALSimIMRPhenomTModelTimesTruncated(hp, hc, m1_SI, m2_SI, chi1z, chi2z, distance, inclination, deltaT, fmin, fRef, phiRef, TimeArray, lalParams) :
        XLALSimIMRPhenomTModelTimes(hp, hc, m1_SI, m2_SI, chi1z, chi2z, distance, inclination, deltaT, fmin, fRef, phiRef, TimeArray, lalParams);
    }
    return (truncate==1) ?
      XLALSimIMRPhenomTHMModelTimesTruncated(hp, hc, m1_SI, m2_SI, chi1z, chi2z, distance, inclination, deltaT, fmin, fRef, phiRef, TimeArray, lalParams) :
      XLALSimIMRPhenomTHMModelTimes(hp, hc, m1_SI, m2_SI, chi1z, chi2z, distance, inclination, deltaT, fmin, fRef, phiRef, TimeArray, lalParams);
  }

  SphHarmTimeSeries *hlm = NULL;
  status = IMRPhenomTPHM_L0ModesModelTimes_engine(&hlm, m1_SI, m2_SI, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z,
             distance, inclination, deltaT, fmin, fRef, phiRef, TimeArray, lalParams, only22, truncate);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: function IMRPhenomTPHM_L0ModesModelTimes_engine has failed.");

  status = LALSimIMRPhenomT_ModelTimes_AddModes(hp, hc, hlm, inclination, phiRef, deltaT);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: LALSimIMRPhenomT_ModelTimes_AddModes has failed.");

  XLALDestroySphHarmTimeSeries(hlm);

  return status;
}

/** IMRPhenomTP modes in the inertial L0 frame at the times of TimeArray (see XLALSimIMRPhenomTPHM_L0Modes). */
int XLALSimIMRPhenomTPModesModelTimes(
  SphHarmTimeSeries **hlms, /**< [out] Modes (L0 frame) */
  REAL8 m1_SI,      /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,      /**< Mass of companion 2 (kg) */
  REAL8 chi1x,      /**< x component of primary spin*/
  REAL8 chi1y,      /**< y component of primary spin*/
  REAL8 chi1z,      /**< z component of primary spin */
  REAL8 chi2x,      /**< x component of secondary spin*/
  REAL8 chi2y,      /**< y component of secondary spin*/
  REAL8 chi2z,      /**< z component of secondary spin */
  REAL8 distance,   /**< Luminosity distance (m) */
  REAL8 deltaT,     /**< Sampling interval (s) */
  REAL8 fmin,       /**< Starting GW frequency (Hz) */
  REAL8 fRef,       /**< Reference GW frequency (Hz) */
  REAL8 phiRef,     /**< Reference orbital phase (rad) */
  REAL8Sequence *TimeArray, /**< Times (units of M, t=0 at the 22 amplitude peak) */
  LALDict *lalParams  /**< LAL dictionary containing accessory parameters */
  )
{
  return IMRPhenomTPHM_L0ModesModelTimes_engine(hlms, m1_SI, m2_SI, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z,
           distance, 0.0, deltaT, fmin, fRef, phiRef, TimeArray, lalParams, 1, 0);
}

/** As XLALSimIMRPhenomTPModesModelTimes, zero outside the time span of the standard waveform. */
int XLALSimIMRPhenomTPModesModelTimesTruncated(
  SphHarmTimeSeries **hlms, /**< [out] Modes (L0 frame) */
  REAL8 m1_SI,      /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,      /**< Mass of companion 2 (kg) */
  REAL8 chi1x,      /**< x component of primary spin*/
  REAL8 chi1y,      /**< y component of primary spin*/
  REAL8 chi1z,      /**< z component of primary spin */
  REAL8 chi2x,      /**< x component of secondary spin*/
  REAL8 chi2y,      /**< y component of secondary spin*/
  REAL8 chi2z,      /**< z component of secondary spin */
  REAL8 distance,   /**< Luminosity distance (m) */
  REAL8 deltaT,     /**< Sampling interval (s) */
  REAL8 fmin,       /**< Starting GW frequency (Hz) */
  REAL8 fRef,       /**< Reference GW frequency (Hz) */
  REAL8 phiRef,     /**< Reference orbital phase (rad) */
  REAL8Sequence *TimeArray, /**< Times (units of M, t=0 at the 22 amplitude peak) */
  LALDict *lalParams  /**< LAL dictionary containing accessory parameters */
  )
{
  return IMRPhenomTPHM_L0ModesModelTimes_engine(hlms, m1_SI, m2_SI, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z,
           distance, 0.0, deltaT, fmin, fRef, phiRef, TimeArray, lalParams, 1, 1);
}

/** IMRPhenomTP polarisations at the times of TimeArray (see XLALSimIMRPhenomTP). */
int XLALSimIMRPhenomTPModelTimes(
  REAL8TimeSeries **hp, /**< [out] TD waveform for plus polarisation */
  REAL8TimeSeries **hc, /**< [out] TD waveform for cross polarisation */
  REAL8 m1_SI,      /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,      /**< Mass of companion 2 (kg) */
  REAL8 chi1x,      /**< x component of primary spin*/
  REAL8 chi1y,      /**< y component of primary spin*/
  REAL8 chi1z,      /**< z component of primary spin */
  REAL8 chi2x,      /**< x component of secondary spin*/
  REAL8 chi2y,      /**< y component of secondary spin*/
  REAL8 chi2z,      /**< z component of secondary spin */
  REAL8 distance,   /**< Luminosity distance (m) */
  REAL8 inclination,  /**< inclination of source (rad) */
  REAL8 deltaT,     /**< sampling interval (s) */
  REAL8 fmin,     /**< starting GW frequency (Hz) */
  REAL8 fRef,     /**< reference GW frequency (Hz) */
  REAL8 phiRef,     /**< reference orbital phase (rad) */
  REAL8Sequence *TimeArray, /**< Times (units of M, t=0 at the 22 amplitude peak) */
  LALDict *lalParams  /**< LAL dictionary containing accessory parameters */
  )
{
  return IMRPhenomTPHM_ModelTimes_polarizations(hp, hc, m1_SI, m2_SI, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z,
           distance, inclination, deltaT, fmin, fRef, phiRef, TimeArray, lalParams, 1, 0);
}

/** As XLALSimIMRPhenomTPModelTimes, zero outside the time span of the standard waveform. */
int XLALSimIMRPhenomTPModelTimesTruncated(
  REAL8TimeSeries **hp, /**< [out] TD waveform for plus polarisation */
  REAL8TimeSeries **hc, /**< [out] TD waveform for cross polarisation */
  REAL8 m1_SI,      /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,      /**< Mass of companion 2 (kg) */
  REAL8 chi1x,      /**< x component of primary spin*/
  REAL8 chi1y,      /**< y component of primary spin*/
  REAL8 chi1z,      /**< z component of primary spin */
  REAL8 chi2x,      /**< x component of secondary spin*/
  REAL8 chi2y,      /**< y component of secondary spin*/
  REAL8 chi2z,      /**< z component of secondary spin */
  REAL8 distance,   /**< Luminosity distance (m) */
  REAL8 inclination,  /**< inclination of source (rad) */
  REAL8 deltaT,     /**< sampling interval (s) */
  REAL8 fmin,     /**< starting GW frequency (Hz) */
  REAL8 fRef,     /**< reference GW frequency (Hz) */
  REAL8 phiRef,     /**< reference orbital phase (rad) */
  REAL8Sequence *TimeArray, /**< Times (units of M, t=0 at the 22 amplitude peak) */
  LALDict *lalParams  /**< LAL dictionary containing accessory parameters */
  )
{
  return IMRPhenomTPHM_ModelTimes_polarizations(hp, hc, m1_SI, m2_SI, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z,
           distance, inclination, deltaT, fmin, fRef, phiRef, TimeArray, lalParams, 1, 1);
}

/** IMRPhenomTPHM modes in the inertial L0 frame at the times of TimeArray (see XLALSimIMRPhenomTPHM_L0Modes). */
int XLALSimIMRPhenomTPHMModesModelTimes(
  SphHarmTimeSeries **hlms, /**< [out] Modes (L0 frame) */
  REAL8 m1_SI,      /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,      /**< Mass of companion 2 (kg) */
  REAL8 chi1x,      /**< x component of primary spin*/
  REAL8 chi1y,      /**< y component of primary spin*/
  REAL8 chi1z,      /**< z component of primary spin */
  REAL8 chi2x,      /**< x component of secondary spin*/
  REAL8 chi2y,      /**< y component of secondary spin*/
  REAL8 chi2z,      /**< z component of secondary spin */
  REAL8 distance,   /**< Luminosity distance (m) */
  REAL8 deltaT,     /**< Sampling interval (s) */
  REAL8 fmin,       /**< Starting GW frequency (Hz) */
  REAL8 fRef,       /**< Reference GW frequency (Hz) */
  REAL8 phiRef,     /**< Reference orbital phase (rad) */
  REAL8Sequence *TimeArray, /**< Times (units of M, t=0 at the 22 amplitude peak) */
  LALDict *lalParams  /**< LAL dictionary containing accessory parameters */
  )
{
  return IMRPhenomTPHM_L0ModesModelTimes_engine(hlms, m1_SI, m2_SI, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z,
           distance, 0.0, deltaT, fmin, fRef, phiRef, TimeArray, lalParams, 0, 0);
}

/** As XLALSimIMRPhenomTPHMModesModelTimes, zero outside the time span of the standard waveform. */
int XLALSimIMRPhenomTPHMModesModelTimesTruncated(
  SphHarmTimeSeries **hlms, /**< [out] Modes (L0 frame) */
  REAL8 m1_SI,      /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,      /**< Mass of companion 2 (kg) */
  REAL8 chi1x,      /**< x component of primary spin*/
  REAL8 chi1y,      /**< y component of primary spin*/
  REAL8 chi1z,      /**< z component of primary spin */
  REAL8 chi2x,      /**< x component of secondary spin*/
  REAL8 chi2y,      /**< y component of secondary spin*/
  REAL8 chi2z,      /**< z component of secondary spin */
  REAL8 distance,   /**< Luminosity distance (m) */
  REAL8 deltaT,     /**< Sampling interval (s) */
  REAL8 fmin,       /**< Starting GW frequency (Hz) */
  REAL8 fRef,       /**< Reference GW frequency (Hz) */
  REAL8 phiRef,     /**< Reference orbital phase (rad) */
  REAL8Sequence *TimeArray, /**< Times (units of M, t=0 at the 22 amplitude peak) */
  LALDict *lalParams  /**< LAL dictionary containing accessory parameters */
  )
{
  return IMRPhenomTPHM_L0ModesModelTimes_engine(hlms, m1_SI, m2_SI, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z,
           distance, 0.0, deltaT, fmin, fRef, phiRef, TimeArray, lalParams, 0, 1);
}

/** IMRPhenomTPHM polarisations at the times of TimeArray (see XLALSimIMRPhenomTPHM). */
int XLALSimIMRPhenomTPHMModelTimes(
  REAL8TimeSeries **hp, /**< [out] TD waveform for plus polarisation */
  REAL8TimeSeries **hc, /**< [out] TD waveform for cross polarisation */
  REAL8 m1_SI,      /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,      /**< Mass of companion 2 (kg) */
  REAL8 chi1x,      /**< x component of primary spin*/
  REAL8 chi1y,      /**< y component of primary spin*/
  REAL8 chi1z,      /**< z component of primary spin */
  REAL8 chi2x,      /**< x component of secondary spin*/
  REAL8 chi2y,      /**< y component of secondary spin*/
  REAL8 chi2z,      /**< z component of secondary spin */
  REAL8 distance,   /**< Luminosity distance (m) */
  REAL8 inclination,  /**< inclination of source (rad) */
  REAL8 deltaT,     /**< sampling interval (s) */
  REAL8 fmin,     /**< starting GW frequency (Hz) */
  REAL8 fRef,     /**< reference GW frequency (Hz) */
  REAL8 phiRef,     /**< reference orbital phase (rad) */
  REAL8Sequence *TimeArray, /**< Times (units of M, t=0 at the 22 amplitude peak) */
  LALDict *lalParams  /**< LAL dictionary containing accessory parameters */
  )
{
  return IMRPhenomTPHM_ModelTimes_polarizations(hp, hc, m1_SI, m2_SI, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z,
           distance, inclination, deltaT, fmin, fRef, phiRef, TimeArray, lalParams, 0, 0);
}

/** As XLALSimIMRPhenomTPHMModelTimes, zero outside the time span of the standard waveform. */
int XLALSimIMRPhenomTPHMModelTimesTruncated(
  REAL8TimeSeries **hp, /**< [out] TD waveform for plus polarisation */
  REAL8TimeSeries **hc, /**< [out] TD waveform for cross polarisation */
  REAL8 m1_SI,      /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,      /**< Mass of companion 2 (kg) */
  REAL8 chi1x,      /**< x component of primary spin*/
  REAL8 chi1y,      /**< y component of primary spin*/
  REAL8 chi1z,      /**< z component of primary spin */
  REAL8 chi2x,      /**< x component of secondary spin*/
  REAL8 chi2y,      /**< y component of secondary spin*/
  REAL8 chi2z,      /**< z component of secondary spin */
  REAL8 distance,   /**< Luminosity distance (m) */
  REAL8 inclination,  /**< inclination of source (rad) */
  REAL8 deltaT,     /**< sampling interval (s) */
  REAL8 fmin,     /**< starting GW frequency (Hz) */
  REAL8 fRef,     /**< reference GW frequency (Hz) */
  REAL8 phiRef,     /**< reference orbital phase (rad) */
  REAL8Sequence *TimeArray, /**< Times (units of M, t=0 at the 22 amplitude peak) */
  LALDict *lalParams  /**< LAL dictionary containing accessory parameters */
  )
{
  return IMRPhenomTPHM_ModelTimes_polarizations(hp, hc, m1_SI, m2_SI, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z,
           distance, inclination, deltaT, fmin, fRef, phiRef, TimeArray, lalParams, 0, 1);
}

/* Co-precessing modes and the Euler angles of the rotation to the L0 frame at the times of TimeArray.
For each l the L0-frame modes of XLALSimIMRPhenomTP(HM)ModesModelTimes are the co-precessing modes rotated with
PhenomTPHM_RotateModes(h, gamma, cosbeta, alpha) (co-precessing -> J frame, time dependent) and then with
PhenomTPHM_RotateModes_Global(h, -alphaJtoL0, cosbetaJtoL0, -gammaJtoL0) (J -> L0 frame, constant). */
static int IMRPhenomTPHM_CoprecModesAnglesModelTimes_engine(
  SphHarmTimeSeries **hlms,
  REAL8TimeSeries **alphaTS, REAL8TimeSeries **cosbetaTS, REAL8TimeSeries **gammaTS,
  REAL8 *alphaJtoL0, REAL8 *cosbetaJtoL0, REAL8 *gammaJtoL0,
  REAL8 m1_SI, REAL8 m2_SI,
  REAL8 chi1x, REAL8 chi1y, REAL8 chi1z,
  REAL8 chi2x, REAL8 chi2y, REAL8 chi2z,
  REAL8 distance, REAL8 deltaT, REAL8 fmin, REAL8 fRef, REAL8 phiRef,
  REAL8Sequence *TimeArray, LALDict *lalParams, UINT4 only22, UINT4 truncate)
{
  XLAL_CHECK(NULL != hlms && NULL != alphaTS && NULL != cosbetaTS && NULL != gammaTS, XLAL_EFAULT);
  XLAL_CHECK(*hlms == NULL && *alphaTS == NULL && *cosbetaTS == NULL && *gammaTS == NULL, XLAL_EFAULT);
  XLAL_CHECK(NULL != alphaJtoL0 && NULL != cosbetaJtoL0 && NULL != gammaJtoL0, XLAL_EFAULT);

  UINT4 status;
  REAL8TimeSeries *alphaNat = NULL, *cosbetaNat = NULL, *gammaNat = NULL;
  REAL8 af;

  status = IMRPhenomTPHM_CoprecModesModelTimes_engine(hlms, alphaTS, cosbetaTS, gammaTS, &alphaNat, &cosbetaNat, &gammaNat, &af,
             m1_SI, m2_SI, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z,
             distance, 0.0, deltaT, fmin, fRef, phiRef, TimeArray, lalParams, only22, truncate);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: function IMRPhenomTPHM_CoprecModesModelTimes_engine has failed.");

  status = IMRPhenomTPHM_JtoL0AnglesModelTimes(alphaJtoL0, cosbetaJtoL0, gammaJtoL0, alphaNat, cosbetaNat, gammaNat,
             m1_SI, m2_SI, chi1z, chi2z, distance, deltaT, fmin, fRef, phiRef, lalParams);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: function IMRPhenomTPHM_JtoL0AnglesModelTimes has failed.");

  XLALDestroyREAL8TimeSeries(alphaNat);
  XLALDestroyREAL8TimeSeries(cosbetaNat);
  XLALDestroyREAL8TimeSeries(gammaNat);

  return status;
}

/** Co-precessing (2,+-2) modes of IMRPhenomTP at the times of TimeArray (units of M, t=0 at the 22 amplitude peak),
with the Euler angles alpha, cos(beta), gamma (co-precessing -> J frame) at the same times and the angles of the
constant J -> L0 rotation (the Euler angles at the reference time). */
int XLALSimIMRPhenomTPCoprecModesModelTimes(
  SphHarmTimeSeries **hlms,                /**< [out] Co-precessing modes */
  REAL8TimeSeries **alphaTS,               /**< [out] Euler angle alpha (co-precessing -> J frame) */
  REAL8TimeSeries **cosbetaTS,             /**< [out] cos(beta) */
  REAL8TimeSeries **gammaTS,               /**< [out] Euler angle gamma */
  REAL8 *alphaJtoL0,                       /**< [out] alpha of the J -> L0 rotation */
  REAL8 *cosbetaJtoL0,                     /**< [out] cos(beta) of the J -> L0 rotation */
  REAL8 *gammaJtoL0,                       /**< [out] gamma of the J -> L0 rotation */
  REAL8 m1_SI,                             /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,                             /**< Mass of companion 2 (kg) */
  REAL8 chi1x,                             /**< x component of primary spin */
  REAL8 chi1y,                             /**< y component of primary spin */
  REAL8 chi1z,                             /**< z component of primary spin */
  REAL8 chi2x,                             /**< x component of secondary spin */
  REAL8 chi2y,                             /**< y component of secondary spin */
  REAL8 chi2z,                             /**< z component of secondary spin */
  REAL8 distance,                          /**< Luminosity distance (m) */
  REAL8 deltaT,                            /**< Sampling interval (s) of the standard waveform */
  REAL8 fmin,                              /**< Starting GW frequency (Hz) */
  REAL8 fRef,                              /**< Reference GW frequency (Hz) */
  REAL8 phiRef,                            /**< Reference orbital phase (rad) */
  REAL8Sequence *TimeArray,                /**< Times (units of M, t=0 at the 22 amplitude peak) */
  LALDict *lalParams                       /**< LAL dictionary containing accessory parameters */
  )
{
  return IMRPhenomTPHM_CoprecModesAnglesModelTimes_engine(hlms, alphaTS, cosbetaTS, gammaTS, alphaJtoL0, cosbetaJtoL0, gammaJtoL0,
           m1_SI, m2_SI, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z, distance, deltaT, fmin, fRef, phiRef, TimeArray, lalParams, 1, 0);
}

/** As XLALSimIMRPhenomTPCoprecModesModelTimes; zero modes and boundary angles outside the time span of the standard waveform. */
int XLALSimIMRPhenomTPCoprecModesModelTimesTruncated(
  SphHarmTimeSeries **hlms, REAL8TimeSeries **alphaTS, REAL8TimeSeries **cosbetaTS, REAL8TimeSeries **gammaTS,
  REAL8 *alphaJtoL0, REAL8 *cosbetaJtoL0, REAL8 *gammaJtoL0,
  REAL8 m1_SI, REAL8 m2_SI, REAL8 chi1x, REAL8 chi1y, REAL8 chi1z, REAL8 chi2x, REAL8 chi2y, REAL8 chi2z,
  REAL8 distance, REAL8 deltaT, REAL8 fmin, REAL8 fRef, REAL8 phiRef, REAL8Sequence *TimeArray, LALDict *lalParams)
{
  return IMRPhenomTPHM_CoprecModesAnglesModelTimes_engine(hlms, alphaTS, cosbetaTS, gammaTS, alphaJtoL0, cosbetaJtoL0, gammaJtoL0,
           m1_SI, m2_SI, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z, distance, deltaT, fmin, fRef, phiRef, TimeArray, lalParams, 1, 1);
}

/** Co-precessing modes of IMRPhenomTPHM (modes selected in the ModeArray of lalParams) at the times of TimeArray, with the
Euler angles at the same times and the angles of the constant J -> L0 rotation. See XLALSimIMRPhenomTPCoprecModesModelTimes. */
int XLALSimIMRPhenomTPHMCoprecModesModelTimes(
  SphHarmTimeSeries **hlms, REAL8TimeSeries **alphaTS, REAL8TimeSeries **cosbetaTS, REAL8TimeSeries **gammaTS,
  REAL8 *alphaJtoL0, REAL8 *cosbetaJtoL0, REAL8 *gammaJtoL0,
  REAL8 m1_SI, REAL8 m2_SI, REAL8 chi1x, REAL8 chi1y, REAL8 chi1z, REAL8 chi2x, REAL8 chi2y, REAL8 chi2z,
  REAL8 distance, REAL8 deltaT, REAL8 fmin, REAL8 fRef, REAL8 phiRef, REAL8Sequence *TimeArray, LALDict *lalParams)
{
  return IMRPhenomTPHM_CoprecModesAnglesModelTimes_engine(hlms, alphaTS, cosbetaTS, gammaTS, alphaJtoL0, cosbetaJtoL0, gammaJtoL0,
           m1_SI, m2_SI, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z, distance, deltaT, fmin, fRef, phiRef, TimeArray, lalParams, 0, 0);
}

/** As XLALSimIMRPhenomTPHMCoprecModesModelTimes; zero modes and boundary angles outside the time span of the standard waveform. */
int XLALSimIMRPhenomTPHMCoprecModesModelTimesTruncated(
  SphHarmTimeSeries **hlms, REAL8TimeSeries **alphaTS, REAL8TimeSeries **cosbetaTS, REAL8TimeSeries **gammaTS,
  REAL8 *alphaJtoL0, REAL8 *cosbetaJtoL0, REAL8 *gammaJtoL0,
  REAL8 m1_SI, REAL8 m2_SI, REAL8 chi1x, REAL8 chi1y, REAL8 chi1z, REAL8 chi2x, REAL8 chi2y, REAL8 chi2z,
  REAL8 distance, REAL8 deltaT, REAL8 fmin, REAL8 fRef, REAL8 phiRef, REAL8Sequence *TimeArray, LALDict *lalParams)
{
  return IMRPhenomTPHM_CoprecModesAnglesModelTimes_engine(hlms, alphaTS, cosbetaTS, gammaTS, alphaJtoL0, cosbetaJtoL0, gammaJtoL0,
           m1_SI, m2_SI, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z, distance, deltaT, fmin, fRef, phiRef, TimeArray, lalParams, 0, 1);
}


/* Harmonic factors Z_lm(t) of the co-precessing modes h_lm(t) at the times of TimeArray, such that
h+ - i hx = sum over the modes of Z_lm(t) h_lm(t). With D^J(t) and D^G the Wigner-D matrices of
PhenomTPHM_RotateModes(h, gamma, cosbeta, alpha) and PhenomTPHM_RotateModes_Global(h, -alphaJtoL0, cosbetaJtoL0,
-gammaJtoL0), and Y_lm = -2Y_lm(inclination, pi/2 - phiRef) as in LALSimIMRPhenomT_ModelTimes_AddModes:
Z_lm(t) = sum_{m''} D^J_{m m''}(t) sum_{m'} D^G_{m'' m'} Y_lm'. One series per mode of hlms. */
static int IMRPhenomTPHM_CoprecHarmonicsModelTimes(
  SphHarmTimeSeries **harmonics, SphHarmTimeSeries *hlms,
  REAL8TimeSeries *alphaTS, REAL8TimeSeries *cosbetaTS, REAL8TimeSeries *gammaTS,
  REAL8 alphaJtoL0, REAL8 cosbetaJtoL0, REAL8 gammaJtoL0, REAL8 inclination, REAL8 phiRef)
{
  XLAL_CHECK(NULL != harmonics && NULL == *harmonics && NULL != hlms, XLAL_EFAULT);
  INT4 lmax = XLALSphHarmTimeSeriesGetMaxL(hlms);
  XLAL_CHECK(lmax >= 2 && lmax <= 5, XLAL_EINVAL, "Modes with l from 2 to 5 are supported.");
  size_t length = hlms->mode->data->length;

  PhenomT_precomputed_sqrt *SQRT = XLALMalloc(sizeof(PhenomT_precomputed_sqrt));
  IMRPhenomTPHM_SetPrecomputedSqrt(SQRT);
  PhenomTPWignerStruct *wS = XLALCalloc(1, sizeof(PhenomTPWignerStruct));

  /* Y rotated with the constant J -> L0 rotation: Yg[l-2][m''+l] = sum_{m'} D^G_{m'' m'} Y_lm' */
  COMPLEX16 Yg[4][11];
  IMRPhenomTPHM_SetWignerDStruct(wS, SQRT, cosbetaJtoL0, -alphaJtoL0, -gammaJtoL0, lmax, -1, 1);
  for(INT4 l = 2; l <= lmax; l++)
  {
    for(INT4 mpp = -l; mpp <= l; mpp++)
    {
      Yg[l-2][mpp+l] = 0.0;
      for(INT4 m = -l; m <= l; m++)
      {
        Yg[l-2][mpp+l] += PhenomTWignerDMatrix(l, mpp, m, wS) *
          XLALSpinWeightedSphericalHarmonic(inclination, LAL_PI/2. - phiRef, -2, l, m);
      }
    }
  }

  for(SphHarmTimeSeries *h = hlms; h; h = h->next)
  {
    COMPLEX16TimeSeries *Z = XLALCreateCOMPLEX16TimeSeries("Z", &h->mode->epoch, 0.0, h->mode->deltaT, &lalDimensionlessUnit, length);
    memset(Z->data->data, 0, length * sizeof(COMPLEX16));
    *harmonics = XLALSphHarmTimeSeriesAddMode(*harmonics, Z, h->l, h->m);
    XLALDestroyCOMPLEX16TimeSeries(Z);
  }

  /* time dependent co-precessing -> J rotation */
  for(size_t i = 0; i < length; i++)
  {
    IMRPhenomTPHM_SetWignerDStruct(wS, SQRT, cosbetaTS->data->data[i], gammaTS->data->data[i], alphaTS->data->data[i], lmax, 1, 0);
    for(SphHarmTimeSeries *Z = *harmonics; Z; Z = Z->next)
    {
      COMPLEX16 z = 0.0;
      for(INT4 mpp = -(INT4)Z->l; mpp <= (INT4)Z->l; mpp++)
      {
        z += PhenomTWignerDMatrix(Z->l, Z->m, mpp, wS) * Yg[Z->l-2][mpp+Z->l];
      }
      Z->mode->data->data[i] = z;
    }
  }

  LALFree(wS);
  LALFree(SQRT);
  return XLAL_SUCCESS;
}

static int IMRPhenomTPHM_CoprecModesHarmonicsModelTimes_engine(
  SphHarmTimeSeries **hlms, SphHarmTimeSeries **harmonics,
  REAL8 m1_SI, REAL8 m2_SI,
  REAL8 chi1x, REAL8 chi1y, REAL8 chi1z,
  REAL8 chi2x, REAL8 chi2y, REAL8 chi2z,
  REAL8 distance, REAL8 inclination, REAL8 deltaT, REAL8 fmin, REAL8 fRef, REAL8 phiRef,
  REAL8Sequence *TimeArray, LALDict *lalParams, UINT4 only22, UINT4 truncate)
{
  REAL8TimeSeries *alphaTS = NULL, *cosbetaTS = NULL, *gammaTS = NULL;
  REAL8 alphaJtoL0, cosbetaJtoL0, gammaJtoL0;
  int status = IMRPhenomTPHM_CoprecModesAnglesModelTimes_engine(hlms, &alphaTS, &cosbetaTS, &gammaTS,
             &alphaJtoL0, &cosbetaJtoL0, &gammaJtoL0, m1_SI, m2_SI, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z,
             distance, deltaT, fmin, fRef, phiRef, TimeArray, lalParams, only22, truncate);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: function IMRPhenomTPHM_CoprecModesAnglesModelTimes_engine has failed.");
  status = IMRPhenomTPHM_CoprecHarmonicsModelTimes(harmonics, *hlms, alphaTS, cosbetaTS, gammaTS,
             alphaJtoL0, cosbetaJtoL0, gammaJtoL0, inclination, phiRef);
  XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: function IMRPhenomTPHM_CoprecHarmonicsModelTimes has failed.");
  XLALDestroyREAL8TimeSeries(alphaTS);
  XLALDestroyREAL8TimeSeries(cosbetaTS);
  XLALDestroyREAL8TimeSeries(gammaTS);
  return status;
}

/** Co-precessing (2,+-2) modes of IMRPhenomTP at the times of TimeArray (units of M, t=0 at the 22 amplitude peak)
and their harmonic factors Z_lm(t), such that h+ - i hx = sum over the modes of Z_lm(t) h_lm(t). Z includes the
rotation from the co-precessing frame to the L0 frame and the spherical harmonics at (inclination, pi/2 - phiRef). */
int XLALSimIMRPhenomTPCoprecModesHarmonicsModelTimes(
  SphHarmTimeSeries **hlms,                /**< [out] Co-precessing modes */
  SphHarmTimeSeries **harmonics,           /**< [out] Harmonic factors Z_lm of the same modes */
  REAL8 m1_SI,                             /**< Mass of companion 1 (kg) */
  REAL8 m2_SI,                             /**< Mass of companion 2 (kg) */
  REAL8 chi1x,                             /**< x component of primary spin */
  REAL8 chi1y,                             /**< y component of primary spin */
  REAL8 chi1z,                             /**< z component of primary spin */
  REAL8 chi2x,                             /**< x component of secondary spin */
  REAL8 chi2y,                             /**< y component of secondary spin */
  REAL8 chi2z,                             /**< z component of secondary spin */
  REAL8 distance,                          /**< Luminosity distance (m) */
  REAL8 inclination,                       /**< Inclination of source (rad) */
  REAL8 deltaT,                            /**< Sampling interval (s) of the standard waveform */
  REAL8 fmin,                              /**< Starting GW frequency (Hz) */
  REAL8 fRef,                              /**< Reference GW frequency (Hz) */
  REAL8 phiRef,                            /**< Reference orbital phase (rad) */
  REAL8Sequence *TimeArray,                /**< Times (units of M, t=0 at the 22 amplitude peak) */
  LALDict *lalParams                       /**< LAL dictionary containing accessory parameters */
  )
{
  return IMRPhenomTPHM_CoprecModesHarmonicsModelTimes_engine(hlms, harmonics, m1_SI, m2_SI, chi1x, chi1y, chi1z,
           chi2x, chi2y, chi2z, distance, inclination, deltaT, fmin, fRef, phiRef, TimeArray, lalParams, 1, 0);
}

/** As XLALSimIMRPhenomTPCoprecModesHarmonicsModelTimes; zero modes outside the time span of the standard waveform. */
int XLALSimIMRPhenomTPCoprecModesHarmonicsModelTimesTruncated(
  SphHarmTimeSeries **hlms, SphHarmTimeSeries **harmonics,
  REAL8 m1_SI, REAL8 m2_SI, REAL8 chi1x, REAL8 chi1y, REAL8 chi1z, REAL8 chi2x, REAL8 chi2y, REAL8 chi2z,
  REAL8 distance, REAL8 inclination, REAL8 deltaT, REAL8 fmin, REAL8 fRef, REAL8 phiRef,
  REAL8Sequence *TimeArray, LALDict *lalParams)
{
  return IMRPhenomTPHM_CoprecModesHarmonicsModelTimes_engine(hlms, harmonics, m1_SI, m2_SI, chi1x, chi1y, chi1z,
           chi2x, chi2y, chi2z, distance, inclination, deltaT, fmin, fRef, phiRef, TimeArray, lalParams, 1, 1);
}

/** Co-precessing modes of IMRPhenomTPHM (modes selected in the ModeArray of lalParams) at the times of TimeArray and
their harmonic factors. See XLALSimIMRPhenomTPCoprecModesHarmonicsModelTimes. */
int XLALSimIMRPhenomTPHMCoprecModesHarmonicsModelTimes(
  SphHarmTimeSeries **hlms, SphHarmTimeSeries **harmonics,
  REAL8 m1_SI, REAL8 m2_SI, REAL8 chi1x, REAL8 chi1y, REAL8 chi1z, REAL8 chi2x, REAL8 chi2y, REAL8 chi2z,
  REAL8 distance, REAL8 inclination, REAL8 deltaT, REAL8 fmin, REAL8 fRef, REAL8 phiRef,
  REAL8Sequence *TimeArray, LALDict *lalParams)
{
  return IMRPhenomTPHM_CoprecModesHarmonicsModelTimes_engine(hlms, harmonics, m1_SI, m2_SI, chi1x, chi1y, chi1z,
           chi2x, chi2y, chi2z, distance, inclination, deltaT, fmin, fRef, phiRef, TimeArray, lalParams, 0, 0);
}

/** As XLALSimIMRPhenomTPHMCoprecModesHarmonicsModelTimes; zero modes outside the time span of the standard waveform. */
int XLALSimIMRPhenomTPHMCoprecModesHarmonicsModelTimesTruncated(
  SphHarmTimeSeries **hlms, SphHarmTimeSeries **harmonics,
  REAL8 m1_SI, REAL8 m2_SI, REAL8 chi1x, REAL8 chi1y, REAL8 chi1z, REAL8 chi2x, REAL8 chi2y, REAL8 chi2z,
  REAL8 distance, REAL8 inclination, REAL8 deltaT, REAL8 fmin, REAL8 fRef, REAL8 phiRef,
  REAL8Sequence *TimeArray, LALDict *lalParams)
{
  return IMRPhenomTPHM_CoprecModesHarmonicsModelTimes_engine(hlms, harmonics, m1_SI, m2_SI, chi1x, chi1y, chi1z,
           chi2x, chi2y, chi2z, distance, inclination, deltaT, fmin, fRef, phiRef, TimeArray, lalParams, 0, 1);
}
