/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2003 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Common Public License. Please see     */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#ifndef _RETURN_VALUES_H
#define _RETURN_VALUES_H

/*****************************************************************************
 *****************************************************************************
 *************                                                      **********
 *************                  Return Values                       **********
 *************                                                      **********
 *****************************************************************************
 *****************************************************************************/

/*----------------------- Global return codes -------------------------------*/
#define FUNCTION_TERMINATED_NORMALLY      0
#define ERROR__USER                      -100

/*-------------- Return codes for sym_parse_comand_line() -------------------*/
#define ERROR__OPENING_PARAM_FILE        -110
#define ERROR__PARSING_PARAM_FILE        -111

/*----------------- Return codes for sym_load_problem() ---------------------*/
#define ERROR__READING_GMPL_FILE         -120
#define ERROR__READING_WARM_START_FILE   -121

/*-------------------- Return codes for sym_solve() -------------------------*/
#define TM_NO_PROBLEM                     225
#define TM_NO_SOLUTION                    226
#define TM_OPTIMAL_SOLUTION_FOUND         227
#define TM_TIME_LIMIT_EXCEEDED            228
#define TM_NODE_LIMIT_EXCEEDED            229
#define TM_TARGET_GAP_ACHIEVED            230
#define TM_ERROR__NO_BRANCHING_CANDIDATE -250
#define TM_ERROR__ILLEGAL_RETURN_CODE    -251
#define TM_ERROR__NUMERICAL_INSTABILITY  -252
#define TM_ERROR__COMM_ERROR             -253
#define TM_ERROR__USER                   -275

#endif
