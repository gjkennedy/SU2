/*!
 * \file nemo_turb_diffusion.hpp
 * \brief Declarations of numerics classes for discretization of
 *        viscous fluxes in turbulence problems.
 * \author F. Palacios, T. Economon
 * \version 7.2.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "../scalar/scalar_diffusion.hpp"


/*!
 * \class CNEMOAvgGrad_TurbSA
 * \brief Class for computing viscous term using average of gradients (Spalart-Allmaras Turbulence model).
 * \ingroup ViscDiscr
 * \author A. Bueno.
 */
class CNEMOAvgGrad_TurbSA final : public CNEMOAvgGrad_Scalar {
private:
  const su2double sigma = 2.0/3.0;

  /*!
   * \brief Adds any extra variables to AD
   */
  void ExtraADPreaccIn(void) override;

  /*!
   * \brief SA specific steps in the ComputeResidual method
   * \param[in] config - Definition of the particular problem.
   */
  void FinishResidualCalc(const CConfig* config) override;

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] correct_grad - Whether to correct gradient for skewness.
   * \param[in] config - Definition of the particular problem.
   */
  CNEMOAvgGrad_TurbSA(unsigned short val_nDim, unsigned short val_nVar,
                      unsigned short val_nVar_NEMO,
                      unsigned short val_nPrimVar,
                      unsigned short val_nPrimVarGrad,
                      bool correct_grad, const CConfig* config);

};

/*!
 * \class CNEMOAvgGrad_TurbSA_Neg
 * \brief Class for computing viscous term using average of gradients (Spalart-Allmaras Turbulence model).
 * \ingroup ViscDiscr
 * \author F. Palacios
 */
class CNEMOAvgGrad_TurbSA_Neg final : public CNEMOAvgGrad_Scalar {
private:
  const su2double sigma = 2.0/3.0;
  const su2double cn1 = 16.0;

  /*!
   * \brief Adds any extra variables to AD
   */
  void ExtraADPreaccIn(void) override;

  /*!
   * \brief SA specific steps in the ComputeResidual method
   * \param[in] config - Definition of the particular problem.
   */
  void FinishResidualCalc(const CConfig* config) override;

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] correct_grad - Whether to correct gradient for skewness.
   * \param[in] config - Definition of the particular problem.
   */
  CNEMOAvgGrad_TurbSA_Neg(unsigned short val_nDim, unsigned short val_nVar,
                          unsigned short val_nVar_NEMO,
                          unsigned short val_nPrimVar,
                          unsigned short val_nPrimVarGrad,
                          bool correct_gradient, const CConfig* config);

};

/*!
 * \class CNEMOAvgGrad_TurbSST
 * \brief Class for computing viscous term using average of gradient with correction (Menter SST turbulence model).
 * \ingroup ViscDiscr
 * \author A. Bueno.
 */
class CNEMOAvgGrad_TurbSST final : public CNEMOAvgGrad_Scalar {
private:
  const su2double
  sigma_k1 = 0.0, /*!< \brief Constants for the viscous terms, k-w (1), k-eps (2)*/
  sigma_k2 = 0.0,
  sigma_om1 = 0.0,
  sigma_om2 = 0.0;

  su2double F1_i, F1_j; /*!< \brief Menter's first blending function */

  /*!
   * \brief Adds any extra variables to AD
   */
  void ExtraADPreaccIn(void) override;

  /*!
   * \brief SST specific steps in the ComputeResidual method
   * \param[in] config - Definition of the particular problem.
   */
  void FinishResidualCalc(const CConfig* config) override;

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] constants - Constants of the model.
   * \param[in] correct_grad - Whether to correct gradient for skewness.
   * \param[in] config - Definition of the particular problem.
   */
  CNEMOAvgGrad_TurbSST(unsigned short val_nDim, unsigned short val_nVar,
                       unsigned short val_nVar_NEMO,
                       unsigned short val_nPrimVar,
                       unsigned short val_nPrimVarGrad,
                       const su2double* constants, bool correct_grad,
                       const CConfig* config);

  /*!
   * \brief Sets value of first blending function.
   */
  void SetF1blending(su2double val_F1_i, su2double val_F1_j) override {
    F1_i = val_F1_i; F1_j = val_F1_j;
  }

};
