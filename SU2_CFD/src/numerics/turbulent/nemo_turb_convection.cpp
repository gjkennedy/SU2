/*!
 * \file nemo_turb_convection.cpp
 * \brief Implementation of numerics classes to compute convective
 *        fluxes in turbulence problems.
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

#include "../../../include/numerics/turbulent/nemo_turb_convection.hpp"

CNEMOUpwScalar_TurbSA::CNEMOUpwScalar_TurbSA(unsigned short val_nDim,
                               unsigned short val_nVar,
                               unsigned short val_nVar_NEMO,
                               unsigned short val_nPrimVar,
                               unsigned short val_nPrimVarGrad,
                               const CConfig* config) :
                CNEMOUpwScalar(val_nDim, val_nVar, val_nVar_NEMO,
                                           val_nPrimVar, val_nPrimVarGrad config) { }



void CNEMOUpwSca_TurbSA::ExtraADPreaccIn() {

  AD::SetPreaccIn(V_i, nDim+1);
  AD::SetPreaccIn(V_j, nDim+1);
}

void CNEMOUpwSca_TurbSA::FinishResidualCalc(const CConfig* config) {

  Flux[0] = a0*ScalarVar_i[0]+a1*ScalarVar_j[0];

  if (implicit) {
    Jacobian_i[0][0] = a0;
    Jacobian_j[0][0] = a1;
  }
}

CNEMOUpwSca_TurbSST::CNEMOUpwSca_TurbSST(unsigned short val_nDim,
                                 unsigned short val_nVar,
                                 unsigned short val_nVar_NEMO,
                                 unsigned short val_nPrimVar, unsigned short val_nPrimVarGrad,
                                         const CConfig* config) :
                                         CNEMOUpwScalar(val_nDim, val_nVar, val_nVar_NEMO,
                                                        val_nPrimVar, val_nPrimVarGrad, config) {}

void CNEMOUpwSca_TurbSST::ExtraADPreaccIn() {

  AD::SetPreaccIn(V_i, nDim+3);
  AD::SetPreaccIn(V_j, nDim+3);
}

void CNEMOUpwSca_TurbSST::FinishResidualCalc(const CConfig* config) {

  Flux[0] = a0*Density_i*ScalarVar_i[0]+a1*Density_j*ScalarVar_j[0];
  Flux[1] = a0*Density_i*ScalarVar_i[1]+a1*Density_j*ScalarVar_j[1];

  if (implicit) {
    Jacobian_i[0][0] = a0;    Jacobian_i[0][1] = 0.0;
    Jacobian_i[1][0] = 0.0;   Jacobian_i[1][1] = a0;

    Jacobian_j[0][0] = a1;    Jacobian_j[0][1] = 0.0;
    Jacobian_j[1][0] = 0.0;   Jacobian_j[1][1] = a1;
  }
}
