/****************************************************************************
*
*   Copyright (c) 2017 Windhover Labs, L.L.C. All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions
* are met:
*
* 1. Redistributions of source code must retain the above copyright
*    notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright
*    notice, this list of conditions and the following disclaimer in
*    the documentation and/or other materials provided with the
*    distribution.
* 3. Neither the name Windhover Labs nor the names of its
*    contributors may be used to endorse or promote products derived
*    from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
* FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
* COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
* INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
* BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
* OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
* AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
* ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
* POSSIBILITY OF SUCH DAMAGE.
*
*****************************************************************************/

#include "cfe.h"
#include "pe_app.h"
#include "pe_test_utils.h"
#include "pe_app_test.h"
#include "uttest.h"
#include "ut_osapi_stubs.h"
#include "ut_cfe_sb_stubs.h"
#include "ut_cfe_sb_hooks.h"
#include "ut_cfe_es_stubs.h"
#include "ut_cfe_es_hooks.h"
#include "ut_cfe_evs_stubs.h"
#include "ut_cfe_evs_hooks.h"
#include "ut_cfe_time_stubs.h"
#include "ut_cfe_psp_memutils_stubs.h"
#include "ut_cfe_tbl_stubs.h"
#include "ut_cfe_fs_stubs.h"
#include "ut_cfe_time_stubs.h"
#include <float.h>

int32 hookCalledCount = 0;

/**************************************************************************
 * Tests for PE_InitEvent()
 **************************************************************************/

/**
 * Test PE_InitEvent() with failed CFE_EVS_Register
 */
void Test_PE_InitEvent_Fail_Register(void)
{
    PE oPE;

    /* Set a fail result for EVS */
    int32 result = (CFE_SEVERITY_BITMASK & CFE_SEVERITY_ERROR)
                   | CFE_EVENTS_SERVICE | CFE_EVS_NOT_IMPLEMENTED;
    int32 expected = CFE_EVS_APP_NOT_REGISTERED;

    Ut_CFE_EVS_SetReturnCode(UT_CFE_EVS_REGISTER_INDEX, expected, 1);

    /* Execute the function being tested */
    result = oPE.InitEvent();

    /* Verify results */
    UtAssert_True (result == expected, "InitEvent, failed EVS Register");
}


/**************************************************************************
 * Tests for PE_InitPipe()
 **************************************************************************/
/**
 * Test PE_InitPipe(), fail SCH CFE_SB_CreatePipe
 */
void Test_PE_InitPipe_Fail_CreateSCHPipe(void)
{
    PE oPE;

    /* Set a fail result for SB */
    int32 result = (CFE_SEVERITY_BITMASK & CFE_SEVERITY_ERROR)
                   | CFE_SOFTWARE_BUS_SERVICE | CFE_SB_NOT_IMPLEMENTED;
    int32 expected = CFE_SB_BAD_ARGUMENT;

    Ut_CFE_SB_SetReturnCode(UT_CFE_SB_CREATEPIPE_INDEX, expected, 1);

    /* Execute the function being tested */
    result = oPE.InitPipe();

    /* Verify results */
    UtAssert_True (result == expected, "InitPipe, fail SB create SCH pipe");
}


/**
 * Test PE_InitPipe(), fail CFE_SB_SubscribeEx for wakeup
 */
void Test_PE_InitPipe_Fail_SubscribeWakeup(void)
{
    PE oPE;

    /* Set a fail result for SB */
    int32 result = (CFE_SEVERITY_BITMASK & CFE_SEVERITY_ERROR)
                   | CFE_SOFTWARE_BUS_SERVICE | CFE_SB_NOT_IMPLEMENTED;
    int32 expected = CFE_SB_BAD_ARGUMENT;

    Ut_CFE_SB_SetReturnCode(UT_CFE_SB_SUBSCRIBEEX_INDEX, expected, 1);

    /* Execute the function being tested */
    result = oPE.InitPipe();

    /* Verify results */
    UtAssert_True (result == expected, "InitPipe, fail CFE_SB_SubscribeEx for wakeup");
}


/**
 * Test PE_InitPipe(), fail CFE_SB_SubscribeEx for sendhk
 */
void Test_PE_InitPipe_Fail_SubscribeSendHK(void)
{
    PE oPE;

    /* Set a fail result for SB */
    int32 result = (CFE_SEVERITY_BITMASK & CFE_SEVERITY_ERROR)
                   | CFE_SOFTWARE_BUS_SERVICE | CFE_SB_NOT_IMPLEMENTED;
    int32 expected = CFE_SB_BAD_ARGUMENT;

    Ut_CFE_SB_SetReturnCode(UT_CFE_SB_SUBSCRIBEEX_INDEX, expected, 2);

    /* Execute the function being tested */
    result = oPE.InitPipe();

    /* Verify results */
    UtAssert_True (result == expected, "InitPipe, fail CFE_SB_SubscribeEx for sendhk");
}


/**
 * Test PE_InitPipe(), fail CMD CFE_SB_CreatePipe
 */
void Test_PE_InitPipe_Fail_CreateCMDPipe(void)
{
    PE oPE;

    /* Set a fail result for SB */
    int32 result = (CFE_SEVERITY_BITMASK & CFE_SEVERITY_ERROR)
                   | CFE_SOFTWARE_BUS_SERVICE | CFE_SB_NOT_IMPLEMENTED;
    int32 expected = CFE_SB_BAD_ARGUMENT;

    Ut_CFE_SB_SetReturnCode(UT_CFE_SB_CREATEPIPE_INDEX, expected, 2);

    /* Execute the function being tested */
    result = oPE.InitPipe();

    /* Verify results */
    UtAssert_True (result == expected, "InitPipe, fail SB create CMD pipe");
}


/**
 * Test PE_InitPipe(), fail CFE_SB_Subscribe for CMD msg
 */
void Test_PE_InitPipe_Fail_SubscribeCMD(void)
{
    PE oPE;

    /* Set a fail result for SB */
    int32 result = (CFE_SEVERITY_BITMASK & CFE_SEVERITY_ERROR)
                   | CFE_SOFTWARE_BUS_SERVICE | CFE_SB_NOT_IMPLEMENTED;
    int32 expected = CFE_SB_BAD_ARGUMENT;

    Ut_CFE_SB_SetReturnCode(UT_CFE_SB_SUBSCRIBE_INDEX, expected, 1);

    /* Execute the function being tested */
    result = oPE.InitPipe();

    /* Verify results */
    UtAssert_True (result == expected, "InitPipe, fail CFE_SB_Subscribe for CMD");
}


/**
 * Test PE_InitPipe(), fail DATA CFE_SB_CreatePipe
 */
void Test_PE_InitPipe_Fail_CreateDATAPipe(void)
{
    PE oPE;

    /* Set a fail result for SB */
    int32 result = (CFE_SEVERITY_BITMASK & CFE_SEVERITY_ERROR)
                   | CFE_SOFTWARE_BUS_SERVICE | CFE_SB_NOT_IMPLEMENTED;
    int32 expected = CFE_SB_BAD_ARGUMENT;

    Ut_CFE_SB_SetReturnCode(UT_CFE_SB_CREATEPIPE_INDEX, expected, 2);

    /* Execute the function being tested */
    result = oPE.InitPipe();

    /* Verify results */
    UtAssert_True (result == expected, "InitPipe, fail SB create DATA pipe");
}


/**************************************************************************
 * Tests for PE_InitData()
 **************************************************************************/
/**
 * Test PE_InitData()
 */
void Test_PE_InitData(void)
{
    PE oPE;

    /* Set a fail result */
    int32 expected = CFE_SUCCESS;

    /* Execute the function being tested */
    oPE.InitData();

    /* Verify results */
    //UtAssert_True (result == expected, "InitData");
}


/**************************************************************************
 * Tests for PE_InitApp()
 **************************************************************************/
/**
 * Test PE_InitApp(), fail init event
 */
void Test_PE_InitApp_Fail_InitEvent(void)
{
    PE oPE;

    int32 result = CFE_SUCCESS;
    int32 expected = CFE_EVS_APP_NOT_REGISTERED;

    Ut_CFE_EVS_SetReturnCode(UT_CFE_EVS_REGISTER_INDEX, expected, 1);

    /* Execute the function being tested */
    result = oPE.InitApp();

    /* Verify results */
    UtAssert_True (result == expected, "InitApp, fail init event");
}


/**
 * Test PE_InitApp(), fail init pipe
 */
void Test_PE_InitApp_Fail_InitPipe(void)
{
    PE oPE;

    int32 result = CFE_SUCCESS;
    int32 expected = CFE_SB_BAD_ARGUMENT;

    Ut_CFE_SB_SetReturnCode(UT_CFE_SB_CREATEPIPE_INDEX, expected, 1);

    /* Execute the function being tested */
    result = oPE.InitApp();

    /* Verify results */
    UtAssert_True (result == expected, "InitApp, fail init pipe");
}


/**
 * Test PE_InitApp(), fail init data.
 * NOTE: no current way to fail PE_InitData() in default
 */
void Test_PE_InitApp_Fail_InitData(void)
{
    PE oPE;

    int32 result = CFE_SUCCESS;
    int32 expected = CFE_SUCCESS;

    /* Execute the function being tested */
    result = oPE.InitApp();

    /* Verify results */
    UtAssert_True (result == expected, "InitApp, fail init data");
}


/**
 * Test PE_InitApp(), fail init config table
 */
void Test_PE_InitApp_Fail_InitConfigTbl(void)
{
    PE oPE;

    int32 result = CFE_SUCCESS;
    int32 expected = CFE_TBL_ERR_INVALID_NAME;

    Ut_CFE_TBL_SetReturnCode(UT_CFE_TBL_REGISTER_INDEX, expected, 1);

    /* Execute the function being tested */
    result = oPE.InitApp();

    /* Verify results */
    UtAssert_True (result == expected, "InitApp, fail init config table");
}


/**
 * Test PE_InitApp(), Nominal
 */
void Test_PE_InitApp_Nominal(void)
{
    PE oPE;

    /* Set a fail result for SB */
    int32 result = (CFE_SEVERITY_BITMASK & CFE_SEVERITY_ERROR)
                   | CFE_EXECUTIVE_SERVICE | CFE_ES_ERR_APP_REGISTER;
    int32 expected = CFE_SUCCESS;

    /* Execute the function being tested */
    result = oPE.InitApp();

    /* Verify results */
    UtAssert_True (result == expected, "InitApp, nominal");
}

/**************************************************************************
 * Tests for PE_AppMain()
 **************************************************************************/
/**
 * Test PE_AppMain(), Fail RegisterApp
 */
void Test_PE_AppMain_Fail_RegisterApp(void)
{
    PE oPE;

    /* fail the register app */
    Ut_CFE_ES_SetReturnCode(UT_CFE_ES_REGISTERAPP_INDEX, CFE_ES_ERR_APP_REGISTER, 1);

    /* Execute the function being tested */
    oPE.AppMain();
}


/**
 * Test PE_AppMain(), Fail InitApp
 */
void Test_PE_AppMain_Fail_InitApp(void)
{
    PE oPE;

    /* fail the register app */
    Ut_CFE_EVS_SetReturnCode(UT_CFE_EVS_REGISTER_INDEX, CFE_EVS_APP_NOT_REGISTERED, 1);

    /* Execute the function being tested */
    oPE.AppMain();
}


/**
 * Test PE_AppMain(), Fail AcquireConfigPtrs
 */
void Test_PE_AppMain_Fail_AcquireConfigPtrs(void)
{
    PE oPE;

    /* fail the register app */
    Ut_CFE_TBL_SetReturnCode(UT_CFE_TBL_GETADDRESS_INDEX, CFE_TBL_ERR_INVALID_HANDLE, 2);

    /* Execute the function being tested */
    oPE.AppMain();
}


/**
 * Test PE_AppMain(), Invalid Schedule Message
 */
void Test_PE_AppMain_InvalidSchMessage(void)
{
    PE oPE;

    /* The following will emulate behavior of receiving a SCH message to send HK */
    Ut_CFE_SB_SetReturnCode(UT_CFE_SB_RCVMSG_INDEX, CFE_SUCCESS, 1);
    Ut_CFE_SB_SetReturnCode(UT_CFE_SB_GETMSGID_INDEX, 0, 1);

    Ut_CFE_ES_SetReturnCode(UT_CFE_ES_RUNLOOP_INDEX, FALSE, 2);

    /* Execute the function being tested */
    oPE.AppMain();

}


/**
 * Hook to support: PE_AppMain(), Nominal - SendHK
 */
int32 Test_PE_AppMain_Nominal_SendHK_SendMsgHook(CFE_SB_Msg_t *MsgPtr)
{
    /* TODO:  Test the contents of your HK message here. */

    hookCalledCount++;

    return CFE_SUCCESS;
}

/**
 * Test PE_AppMain(), Nominal - SendHK
 */
void Test_PE_AppMain_Nominal_SendHK(void)
{
    PE oPE;

    /* The following will emulate behavior of receiving a SCH message to WAKEUP */
    Ut_CFE_SB_SetReturnCode(UT_CFE_SB_RCVMSG_INDEX, CFE_SUCCESS, 1);
    Ut_CFE_SB_SetReturnCode(UT_CFE_SB_GETMSGID_INDEX, PE_SEND_HK_MID, 1);

    Ut_CFE_ES_SetReturnCode(UT_CFE_ES_RUNLOOP_INDEX, FALSE, 2);

    /* Used to verify HK was transmitted correctly. */
    hookCalledCount = 0;
    Ut_CFE_ES_SetFunctionHook(UT_CFE_SB_SENDMSG_INDEX, (void*)&Test_PE_AppMain_Nominal_SendHK_SendMsgHook);

    /* Execute the function being tested */
    oPE.AppMain();

    /* Verify results */
    UtAssert_True (hookCalledCount == 1, "AppMain_Nominal_SendHK");

}


/**
 * Test PE_AppMain(), Nominal - Wakeup
 */
void Test_PE_AppMain_Nominal_Wakeup(void)
{
    PE oPE;

    /* The following will emulate behavior of receiving a SCH message to WAKEUP */
    Ut_CFE_SB_SetReturnCode(UT_CFE_SB_RCVMSG_INDEX, CFE_SUCCESS, 1);
    Ut_CFE_SB_SetReturnCode(UT_CFE_SB_GETMSGID_INDEX, PE_WAKEUP_MID, 1);

    Ut_CFE_ES_SetReturnCode(UT_CFE_ES_RUNLOOP_INDEX, FALSE, 2);

    /* Execute the function being tested */
    //oPE.AppMain();

}

/**
 * Test PE_AppMain(), Nominal - Land
 */
void Test_PE_AppMain_Land_Nominal_1(void)
{
    PE oPE;

    math::Matrix10F10 input = {
		{8.014909, -0.000000, -0.000000, 0.000392, -0.000000, -0.000000, -0.000000, 0.000000, 0.000000, 0.000000},
		{-0.000000, 8.014909, -0.000000, 0.000000, 0.000392, -0.000000, -0.000000, -0.000000, -0.000000, 0.000000},
		{-0.000000, -0.000000, 4.002175, -0.000000, -0.000000, 0.002733, 0.000000, 0.000000, -0.000000, 3.999990},
		{0.000392, 0.000000, -0.000000, 0.002367, -0.000000, -0.000000, -0.000000, 0.000001, 0.000000, 0.000000},
		{-0.000000, 0.000392, -0.000000, -0.000000, 0.002367, -0.000000, -0.000001, -0.000000, -0.000000, 0.000000},
		{-0.000000, -0.000000, 0.002733, -0.000000, -0.000000, 0.012086, 0.000000, 0.000000, -0.000002, 0.000835},
		{-0.000000, -0.000000, 0.000000, -0.000000, -0.000001, 0.000000, 0.000003, 0.000000, -0.000000, -0.000000},
		{0.000000, -0.000000, 0.000000, 0.000001, -0.000000, 0.000000, 0.000000, 0.000003, -0.000000, -0.000000},
		{0.000000, -0.000000, -0.000000, 0.000000, -0.000000, -0.000002, -0.000000, -0.000000, 0.000003, 0.000000},
		{0.000000, 0.000000, 3.999990, 0.000000, 0.000000, 0.000835, -0.000000, -0.000000, 0.000000, 3.999844}
    };

	math::Matrix10F10 output = {
		{8.014877, 0.000000, -0.000000, 0.000201, -0.000000, 0.000000, 0.000000, -0.000000, -0.000000, -0.000000},
		{0.000000, 8.014877, -0.000000, 0.000000, 0.000201, 0.000000, 0.000000, 0.000000, 0.000000, -0.000000},
		{-0.000000, -0.000000, 4.000550, 0.000000, 0.000000, 0.001322, 0.000000, 0.000000, -0.000000, 3.999881},
		{0.000201, 0.000000, 0.000000, 0.001216, -0.000000, -0.000000, -0.000000, 0.000000, 0.000000, 0.000000},
		{-0.000000, 0.000201, 0.000000, -0.000000, 0.001216, -0.000000, -0.000000, -0.000000, -0.000000, 0.000000},
		{0.000000, 0.000000, 0.001322, -0.000000, -0.000000, 0.010861, 0.000000, 0.000000, -0.000002, 0.000740},
		{0.000000, 0.000000, 0.000000, -0.000000, -0.000000, 0.000000, 0.000003, 0.000000, -0.000000, -0.000000},
		{-0.000000, 0.000000, 0.000000, 0.000000, -0.000000, 0.000000, 0.000000, 0.000003, -0.000000, -0.000000},
		{-0.000000, 0.000000, -0.000000, 0.000000, -0.000000, -0.000002, -0.000000, -0.000000, 0.000003, 0.000000},
		{-0.000000, -0.000000, 3.999881, 0.000000, 0.000000, 0.000740, -0.000000, -0.000000, 0.000000, 3.999836}
	};

	oPE.m_StateCov = input;

    /* Execute the function being tested */
	oPE.landCorrect();

	for(int i = 0; i < 10; i++)
	{
		for(int j = 0; j < 10; j++)
		{
			UtAssert_DoubleCmpAbs(oPE.m_StateCov[i][j], output[i][j], FLT_EPSILON, "");
		}
	}

}

/**
 * Test PE_AppMain(), Nominal - Land
 */
void Test_PE_AppMain_Land_Nominal_2(void)
{
    PE oPE;


    /* Execute the function being tested */


}

/**
 * Test PE_AppMain(), Nominal - Land
 */
void Test_PE_AppMain_Land_Nominal_3(void)
{
    PE oPE;


    /* Execute the function being tested */


}


/**************************************************************************
 * Rollup Test Cases
 **************************************************************************/
void PE_App_Test_AddTestCases(void)
{
    UtTest_Add(Test_PE_InitEvent_Fail_Register, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_InitEvent_Fail_Register");

    UtTest_Add(Test_PE_InitPipe_Fail_CreateSCHPipe, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_InitPipe_Fail_CreateSCHPipe");
    UtTest_Add(Test_PE_InitPipe_Fail_SubscribeWakeup, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_InitPipe_Fail_SubscribeWakeup");
    UtTest_Add(Test_PE_InitPipe_Fail_SubscribeSendHK, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_InitPipe_Fail_SubscribeSendHK");
    UtTest_Add(Test_PE_InitPipe_Fail_CreateCMDPipe, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_InitPipe_Fail_CreateCMDPipe");
    UtTest_Add(Test_PE_InitPipe_Fail_SubscribeCMD, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_InitPipe_Fail_SubscribeCMD");
    UtTest_Add(Test_PE_InitPipe_Fail_CreateDATAPipe, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_InitPipe_Fail_CreateDATAPipe");

    UtTest_Add(Test_PE_InitData, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_InitData");

    UtTest_Add(Test_PE_InitApp_Fail_InitEvent, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_InitApp_Fail_InitEvent");
    UtTest_Add(Test_PE_InitApp_Fail_InitPipe, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_InitApp_Fail_InitPipe");
    UtTest_Add(Test_PE_InitApp_Fail_InitData, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_InitApp_Fail_InitData");
    UtTest_Add(Test_PE_InitApp_Fail_InitConfigTbl, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_InitApp_Fail_InitConfigTbl");
    UtTest_Add(Test_PE_InitApp_Nominal, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_InitApp_Nominal");

    UtTest_Add(Test_PE_AppMain_Fail_RegisterApp, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_AppMain_Fail_RegisterApp");
    UtTest_Add(Test_PE_AppMain_Fail_InitApp, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_AppMain_Fail_InitApp");
    UtTest_Add(Test_PE_AppMain_Fail_AcquireConfigPtrs, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_AppMain_Fail_AcquireConfigPtrs");
    UtTest_Add(Test_PE_AppMain_InvalidSchMessage, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_AppMain_InvalidSchMessage");
    UtTest_Add(Test_PE_AppMain_Nominal_SendHK, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_AppMain_Nominal_SendHK");
    UtTest_Add(Test_PE_AppMain_Nominal_Wakeup, PE_Test_Setup, PE_Test_TearDown,
               "Test_PE_AppMain_Nominal_Wakeup");

}


