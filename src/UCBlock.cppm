/*--------------------------------------------------------------------------*/
/*--------------------------- File UCBlock.cppm ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Interface of the C++20 module smspp.ucblock, which exports every name that
 * the headers of UCBlock declare at namespace scope, together with the
 * modules it depends on; a translation unit can then write
 *
 *     import smspp.ucblock;
 *
 * in place of the #include of the headers. The macros (e.g., those of the
 * factory in SMSFactory.h) are not exported by a module, so a file that uses
 * them still includes the header that defines them.
 *
 * This file is written by SMS++/cmake/smspp_module_exports.py, which reads
 * the names from the AST of the headers: it is not to be edited by hand.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Ali Ghezelsoflu \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Rafael Durbano Lobato \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Kostas Tavlaridis-Gyparakis \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Ali Ghezelsoflu,
 *                      Rafael Durbano Lobato, Kostas Tavlaridis-Gyparakis,
 *                      Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*------------------------- GLOBAL MODULE FRAGMENT -------------------------*/
/*--------------------------------------------------------------------------*/

module;

#include "ACNetworkBlock.h"
#include "BatteryUnitBlock.h"
#include "ConversionUnitBlock.h"
#include "DesignNetworkBlock.h"
#include "ECNetworkBlock.h"
#include "HydroSystemUnitBlock.h"
#include "IntermittentUnitBlock.h"
#include "NuclearUnitExtDPSolver.h"
#include "OTSNetworkBlock.h"
#include "SlackUnitBlock.h"
#include "ThermalUnitDPSolver.h"
#include "UCBlock.h"
#include "UCBlockModel.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- MODULE PURVIEW -----------------------------*/
/*--------------------------------------------------------------------------*/

export module smspp.ucblock;

export import smspp;

export namespace SMSpp_di_unipi_it {

 // ACNetworkBlock.h
 using SMSpp_di_unipi_it::ACNetworkBlock;
 using SMSpp_di_unipi_it::ACNetworkBlockSolution;

 // BatteryUnitBlock.h
 using SMSpp_di_unipi_it::BatteryUnitBlock;
 using SMSpp_di_unipi_it::BatteryUnitBlockMod;
 using SMSpp_di_unipi_it::BatteryUnitBlockRngdMod;
 using SMSpp_di_unipi_it::BatteryUnitBlockSbstMod;
 using SMSpp_di_unipi_it::BatteryUnitBlockSolution;

 // ConversionUnitBlock.h
 using SMSpp_di_unipi_it::ConversionUnitBlock;
 using SMSpp_di_unipi_it::ConversionUnitBlockMod;
 using SMSpp_di_unipi_it::ConversionUnitBlockSolution;

 // DCNetworkBlock.h
 using SMSpp_di_unipi_it::DCNetworkBlock;
 using SMSpp_di_unipi_it::DCNetworkBlockMod;
 using SMSpp_di_unipi_it::DCNetworkBlockRngdMod;
 using SMSpp_di_unipi_it::DCNetworkBlockSbstMod;
 using SMSpp_di_unipi_it::DCNetworkBlockSolution;
 using SMSpp_di_unipi_it::SpMat;

 // DesignNetworkBlock.h
 using SMSpp_di_unipi_it::DesignNetworkBlock;
 using SMSpp_di_unipi_it::DesignNetworkBlockSolution;

 // ECNetworkBlock.h
 using SMSpp_di_unipi_it::ECNetworkBlock;
 using SMSpp_di_unipi_it::ECNetworkBlockMod;
 using SMSpp_di_unipi_it::ECNetworkBlockRngdMod;
 using SMSpp_di_unipi_it::ECNetworkBlockSbstMod;
 using SMSpp_di_unipi_it::ECNetworkBlockSolution;

 // HydroSystemUnitBlock.h
 using SMSpp_di_unipi_it::HydroSystemUnitBlock;
 using SMSpp_di_unipi_it::HydroSystemUnitBlockSolution;

 // HydroUnitBlock.h
 using SMSpp_di_unipi_it::HydroUnitBlock;
 using SMSpp_di_unipi_it::HydroUnitBlockMod;
 using SMSpp_di_unipi_it::HydroUnitBlockRngdMod;
 using SMSpp_di_unipi_it::HydroUnitBlockSbstMod;
 using SMSpp_di_unipi_it::HydroUnitBlockSolution;

 // IntermittentUnitBlock.h
 using SMSpp_di_unipi_it::IntermittentUnitBlock;
 using SMSpp_di_unipi_it::IntermittentUnitBlockMod;
 using SMSpp_di_unipi_it::IntermittentUnitBlockRngdMod;
 using SMSpp_di_unipi_it::IntermittentUnitBlockSbstMod;
 using SMSpp_di_unipi_it::IntermittentUnitBlockSolution;

 // NetworkBlock.h
 using SMSpp_di_unipi_it::NetworkBlock;
 using SMSpp_di_unipi_it::NetworkBlockMod;
 using SMSpp_di_unipi_it::NetworkBlockRngdMod;
 using SMSpp_di_unipi_it::NetworkBlockSbstMod;
 using SMSpp_di_unipi_it::NetworkBlockSolution;
 using SMSpp_di_unipi_it::deserialize_local_dim;

 // NuclearUnitBlock.h
 using SMSpp_di_unipi_it::NuclearUnitBlock;
 using SMSpp_di_unipi_it::NuclearUnitBlockMod;
 using SMSpp_di_unipi_it::NuclearUnitBlockRngdMod;
 using SMSpp_di_unipi_it::NuclearUnitBlockSbstMod;
 using SMSpp_di_unipi_it::NuclearUnitBlockSolution;

 // NuclearUnitExtDPSolver.h
 using SMSpp_di_unipi_it::NuclearUnitExtDPSolver;

 // OTSNetworkBlock.h
 using SMSpp_di_unipi_it::OTSNetworkBlock;

 // SlackUnitBlock.h
 using SMSpp_di_unipi_it::SlackUnitBlock;
 using SMSpp_di_unipi_it::SlackUnitBlockMod;
 using SMSpp_di_unipi_it::SlackUnitBlockRngdMod;
 using SMSpp_di_unipi_it::SlackUnitBlockSbstMod;

 // ThermalUnitBlock.h
 using SMSpp_di_unipi_it::ThermalUnitBlock;
 using SMSpp_di_unipi_it::ThermalUnitBlockMod;
 using SMSpp_di_unipi_it::ThermalUnitBlockRngdMod;
 using SMSpp_di_unipi_it::ThermalUnitBlockSbstMod;
 using SMSpp_di_unipi_it::ThermalUnitBlockSolution;

 // ThermalUnitDPSolver.h
 using SMSpp_di_unipi_it::ThermalUnitDPSolver;

 // ThermalUnitDPSolverBase.h
 using SMSpp_di_unipi_it::ThermalUnitDPSolverBase;

 // ThermalUnitExtDPSolver.h
 using SMSpp_di_unipi_it::ThermalUnitExtDPSolver;

 // UCBlock.h
 using SMSpp_di_unipi_it::UCBlock;
 using SMSpp_di_unipi_it::UCBlockMod;
 using SMSpp_di_unipi_it::UCBlockRngdMod;
 using SMSpp_di_unipi_it::UCBlockSbstMod;
 using SMSpp_di_unipi_it::UCBlockSolution;

 // UnitBlock.h
 using SMSpp_di_unipi_it::UnitBlock;
 using SMSpp_di_unipi_it::UnitBlockMod;
 using SMSpp_di_unipi_it::UnitBlockSolution;
}

/*--------------------------------------------------------------------------*/
/*------------------------- End File UCBlock.cppm --------------------------*/
/*--------------------------------------------------------------------------*/
