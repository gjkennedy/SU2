/*!
 * \file variable_adjoint_levelset.cpp
 * \brief Definition of the solution fields.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 3.2.4 "eagle"
 *
 * SU2, Copyright (C) 2012-2014 Aerospace Design Laboratory (ADL).
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

#include "../include/variable_structure.hpp"


CAdjLevelSetVariable::CAdjLevelSetVariable(void) : CVariable() {}

CAdjLevelSetVariable::CAdjLevelSetVariable(unsigned short val_nDim, unsigned short val_nvar, CConfig *config)
: CVariable(val_nDim, val_nvar, config) {
	
	/*--- Allocate residual structures ---*/
  
	Residual_Sum = new double [nVar]; Residual_Old = new double [nVar];
	
	/*--- Allocate limiter (upwind)---*/
  
	if (config->GetSpatialOrder() == SECOND_ORDER_LIMITER) Limiter = new double [nVar];
	
}

CAdjLevelSetVariable::CAdjLevelSetVariable(double val_levelset, unsigned short val_nDim, unsigned short val_nvar, CConfig *config)
: CVariable(val_nDim, val_nvar,config) {
	
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  
	/*--- Allocate residual structures ---*/
  
	Residual_Sum = new double [nVar]; Residual_Old = new double [nVar];
	
	/*--- Allocate limiter (upwind)---*/
  
	if (config->GetSpatialOrder() == SECOND_ORDER_LIMITER) Limiter = new double [nVar];
	
	/*--- Solution and old solution initialization ---*/
  
	Solution[0] = val_levelset;		Solution_Old[0] = val_levelset;
	
	if (dual_time) {
		Solution_time_n[0] = val_levelset;
		Solution_time_n1[0] = val_levelset;
	}
	
}

CAdjLevelSetVariable::~CAdjLevelSetVariable(void) {
	
	delete [] Residual_Sum; delete [] Residual_Old;
	if (Limiter != NULL) delete [] Limiter;
	
}
