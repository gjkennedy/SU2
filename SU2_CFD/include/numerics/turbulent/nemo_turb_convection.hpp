﻿/*!
 * \file nemo_turb_convection.hpp
 * \brief Delarations of numerics classes for discretization of
 *        convective fluxes in turbulence problems.
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

#include "../scalar/scalar_convection.hpp"

/*!
 * \class CUpwSca_TurbSA
 * \brief Class for doing a scalar upwind solver for the Spalar-Allmaras turbulence model equations.
 * \ingroup ConvDiscr
 * \author A. Bueno.
 */
class CNEMOUpwSca_TurbSA final : public CNEMOUpwScalar {
private:
  /*!
   * \brief Adds any extra variables to AD
   */
  void ExtraADPreaccIn() override;

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
   * \param[in] config - Definition of the particular problem.
   */
  CNEMOUpwSca_TurbSA(unsigned short val_nDim, unsigned short val_nVar,
                     unsigned short val_nVar_NEMO,
                     unsigned short val_nPrimVar,
                     unsigned short val_nPrimVarGrad,
                     const CConfig* config);

};

/*!
 * \class CUpwSca_TurbSST
 * \brief Class for doing a scalar upwind solver for the Menter SST turbulence model equations.
 * \ingroup ConvDiscr
 * \author A. Campos.
 */
class CNEMOUpwSca_TurbSST final : public CNEMOUpwScalar {
private:
  /*!
   * \brief Adds any extra variables to AD
   */
  void ExtraADPreaccIn() override;

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
   * \param[in] config - Definition of the particular problem.
   */
  CNEMOUpwSca_TurbSST(unsigned short val_nDim, unsigned short val_nVar,
                      unsigned short val_nVar_NEMO,
                      unsigned short val_nPrimVar,
                      unsigned short val_nPrimVarGrad,
                      const CConfig* config);

};
